# -*- coding: utf-8 -*-
"""
方法B：材料坐标 + 分域谱元（graded spectral element, nodal Lagrange/Gauss-Lobatto）
弱形式 Galerkin，Robin 为自然边界条件，质量/能量守恒律在离散层面精确成立。

单元划分刻意把表面边界层（厚度 ~ D/(R*h_m)）单独放入最外层单元，
从而避免全局多项式基在薄边界层上的 Gibbs 振铃导致的非线性失稳。

未知量 = 各单元节点上的物理值（温度 / 干基含水率），界面节点共享（C0 连续），
通量连续由弱形式自然给出。
"""
import numpy as np
from numpy.polynomial import polynomial as P
from scipy.special import roots_legendre
from scipy.integrate import solve_ivp

from hb_core import (R0, H_CONV, HM_MASS, T0_C, C0, Ambient, Radius,
                     props_q1, props_q23, props_q4)


def cheb_lobatto(n):
    """[-1,1] 上 n+1 个 Chebyshev-Lobatto 节点（升序）"""
    k = np.arange(n + 1)
    return -np.cos(np.pi * k / n)


class SEM1D:
    """一维（轴对称）分域谱元：内积权重 sigma dsigma，域 sigma in [0,1]。"""

    def __init__(self, breaks=(0.0, 0.6, 0.9, 1.0), n=20, nq=None):
        self.breaks = list(breaks)
        self.n = int(n)
        self.ne = len(self.breaks) - 1
        self.nq = int(nq) if nq else self.n + 3
        self.elems = []
        gid = 0
        for e in range(self.ne):
            a, b = self.breaks[e], self.breaks[e + 1]
            he = b - a
            xi = cheb_lobatto(self.n)                       # n+1 节点
            sig_node = 0.5 * (a + b) + 0.5 * he * xi
            # Lagrange 多项式（升幂系数，变量 xi）
            Lc = []
            for k in range(self.n + 1):
                idx = [m for m in range(self.n + 1) if m != k]
                roots_m = xi[idx]
                denom = np.prod(xi[k] - roots_m)
                Lc.append(P.polyfromroots(roots_m) / denom)
            # Gauss-Legendre 积分点（映射到 [a,b]）
            xq, wq = roots_legendre(self.nq)
            sig_q = 0.5 * (a + b) + 0.5 * he * xq
            w_sig = 0.5 * he * wq * sig_q                   # 含 sigma dsigma
            L = np.array([P.polyval(xq, c) for c in Lc])              # (n+1,nq)
            DL = np.array([P.polyval(xq, P.polyder(c)) * (2.0 / he) for c in Lc])  # d/dsigma
            # 未加权测度 ds 下的权重；CUM[q,p]=int_{-1}^{xi_q} l_p dxi ；TOT[p]=int_{-1}^{1} l_p dxi
            wu = 0.5 * he * wq
            CUM = np.zeros((self.nq, self.n + 1))
            TOT = np.zeros(self.n + 1)
            for p in range(self.n + 1):
                ant = P.polyint(Lc[p])
                CUM[:, p] = (P.polyval(xq, ant) - P.polyval(-1.0, ant)) * 0.5 * he
                TOT[p] = (P.polyval(1.0, ant) - P.polyval(-1.0, ant)) * 0.5 * he
            # 全局编号
            if e == 0:
                g = np.arange(self.n + 1)
                gid = self.n + 1
            else:
                g = np.concatenate([[gid - 1], np.arange(gid, gid + self.n)])
                gid = gid + self.n
            self.elems.append(dict(a=a, b=b, he=he, xi=xi, sig_node=sig_node,
                                   sig_q=sig_q, w=w_sig, wu=wu, CUM=CUM, TOT=TOT,
                                   L=L, DL=DL, g=g))
        self.ndof = gid
        # 输出用 barycentric 权重
        for el in self.elems:
            xi = el['xi']
            n = self.n
            w = np.array([(-1.0) ** k * (0.5 if (k == 0 or k == n) else 1.0)
                          for k in range(n + 1)])
            el['bw'] = w
        self.is_surf = None      # 由外部指定（自由表面 sigma=1 的 DOF）
        # 节点 sigma 坐标（全局）
        self.sig_nodes = np.zeros(self.ndof)
        for el in self.elems:
            self.sig_nodes[el['g']] = el['sig_node']
        order = np.argsort(self.sig_nodes)
        self.surf_dof = int(order[-1])          # sigma=1 对应节点
        assert abs(self.sig_nodes[self.surf_dof] - 1.0) < 1e-14

    # ---------------- 工具 ----------------
    def gather(self, u):
        """把全局节点值 u 插值到所有单元的积分点，返回列表"""
        return [el['L'].T @ u[el['g']] for el in self.elems]

    def assemble(self, coefs, kind='L', wlist=None):
        """组装全局矩阵： sum_q w_q * coef_q * (L 或 DL)[:,q] (.)^T"""
        n = self.ndof
        out = np.zeros((n, n))
        for e, (el, cf) in enumerate(zip(self.elems, coefs)):
            B = el['L'] if kind == 'L' else el['DL']
            ww = el['w'] if wlist is None else wlist[e]
            out[np.ix_(el['g'], el['g'])] += (B * (ww * cf)) @ B.T
        return out

    # ---- 干物质(Lagrangian)坐标下的度量：v(C)=(1+C)/rho(C), I, J ----
    def metrics(self, Cnodes, rho):
        """Cnodes: 各单元节点上的 C 值; rho: 密度函数。返回 (Jq 各积分点, I)"""
        if self.mapping == 'affine':
            return None, 1.0
        vn = [(1.0 + c) / rho(c) for c in Cnodes]
        I = float(sum(el['TOT'] @ v for el, v in zip(self.elems, vn)))
        Jq, acc = [], 0.0
        for el, v in zip(self.elems, vn):
            Jq.append((acc + el['CUM'] @ v) / I)
            acc += float(el['TOT'] @ v)
        return Jq, I

    def interp(self, u, sigma):
        """在任意 sigma（数组）处取值（barycentric, 向量化）"""
        sigma = np.atleast_1d(np.asarray(sigma, float))
        out = np.empty_like(sigma)
        done = np.zeros(sigma.shape, bool)
        for el in self.elems:
            m = (sigma >= el['a'] - 1e-13) & (sigma <= el['b'] + 1e-13) & (~done)
            if not m.any():
                continue
            xi = 2.0 * (sigma[m] - el['a']) / el['he'] - 1.0
            xn = el['xi']
            w = el['bw']
            d = xi[:, None] - xn[None, :]
            ad = np.abs(d)
            t = w[None, :] / np.where(ad < 1e-300, 1e-300, d)
            val = (t @ u[el['g']]) / t.sum(axis=1)
            ex = ad.min(axis=1) < 1e-13
            if ex.any():
                val[ex] = u[el['g'][np.argmin(ad[ex], axis=1)]]
            out[m] = val
            done[m] = True
        if not done.all():
            raise ValueError('sigma out of range: %s' % sigma[~done][:5])
        return out

    def global_quad(self):
        """全局积分点（用于守恒量/诊断）"""
        s = np.concatenate([el['sig_q'] for el in self.elems])
        w = np.concatenate([el['w'] for el in self.elems])
        return s, w


class Model:
    """状态 y = [T 节点值 (ndof); C 节点值 (ndof)]

    mapping='affine' : 仿射(均匀)收缩材料坐标 sigma=r/R(t)（问题1-4基线, 等价于工作区 M0）
    mapping='local'  : 干物质坐标 s + 局部比容 v(C)=(1+C)/rho(C) 的收缩映射
                       r=R(t)sqrt(J(s)),  J=int_0^s v dxi / int_0^1 v dxi
                       （等价于工作区 M1 非均匀收缩自洽模型）
    """

    def __init__(self, sem=None, props=None, amb=None, radius=None, mapping='affine',
                 rho_grad=False):
        self.S = sem or SEM1D()
        self.S.mapping = mapping
        self.mapping = mapping
        self.rho_grad = bool(rho_grad)   # True: 保留 (ln rho_d)_sigma * C_sigma 项
        self.pr = props or props_q23()
        self.amb = amb
        self.rad = radius or Radius(fixed=True)
        self.ndof = self.S.ndof
        self.history = []
        # 常数质量矩阵（仿射模式的质量方程用）及其分解
        from scipy.linalg import lu_factor
        M0 = self.S.assemble([np.ones_like(el['sig_q']) for el in self.S.elems], 'L')
        self.M0 = M0
        self.M0_lu = lu_factor(M0)
        # 常物性快速通道（问题1）：MT、AT 不随时间变化
        self._const = bool(self.pr.get('const', False))
        if self._const:
            cq = [np.full(el['sig_q'].shape, C0) for el in self.S.elems]
            tq = [np.full(el['sig_q'].shape, T0_C) for el in self.S.elems]
            rcp = [self.pr['rho'](c) * self.pr['cp'](c) for c in cq]
            kq = [self.pr['k'](c) for c in cq]
            self.MT_const = self.S.assemble(rcp, 'L')
            self.AT_const = self.S.assemble(kq, 'DL')
            self.MT_lu = lu_factor(self.MT_const)

    # ------------- 物理场 -------------
    def fields(self, y):
        S = self.S
        T = y[:S.ndof]
        C = y[S.ndof:]
        Tq = S.gather(T)
        Cq = S.gather(C)
        return T, C, Tq, Cq

    def rhs(self, t, y):
        """统一弱形式（两种映射）：
        mapping='affine' : 测度 sigma dsigma,  刚度系数 k / D,        边界 h/R, hm/R
        mapping='local'  : 测度 (v/I)ds,       刚度系数 (4I/R2)Jk,    边界 2h/R, 2hm/R
                           水分刚度 (4I/R2)JD/v"""
        S = self.S
        T, C, Tq, Cq = self.fields(y)
        R = float(self.rad(t))
        Tinf, Cinf = self.amb(t)
        if self.mapping == 'affine':
            Jq = [el['sig_q'] for el in S.elems]
            I = 1.0
            wl = [el['w'] for el in S.elems]
            scale = 1.0 / (R * R)
            Dcoef = [self.pr['D'](c, tt) for c, tt in zip(Cq, Tq)]
            kcoef = [self.pr['k'](c) for c in Cq]
            hb = H_CONV / R
            mb = HM_MASS / R
            MT = self.MT_const if self._const else S.assemble(
                [self.pr['rho'](c) * self.pr['cp'](c) for c in Cq], 'L')
            AT = self.AT_const if self._const else S.assemble(kcoef, 'DL')
        else:
            # 干物质坐标 s: 测度 ds;  v=(1+C)/rho;  J=int_0^s v/I;  r=R sqrt(J)
            # 质量: C_t = (4I^2/R^2) d_s((J D/v^2) C_s)
            # 焓  : (1+C)cp T_t = (4I^2/R^2) d_s((J k/v) T_s)
            # 边界: -(2 I hm/(R v1))(C-Cinf),  -(2 I h/R)(T-Tinf)
            Jq, I = S.metrics([C[el['g']] for el in S.elems], self.pr['rho'])
            wl = [el['wu'] for el in S.elems]                      # 测度 ds
            wc = [np.full_like(el['sig_q'], 4.0 * I * I / (R * R)) for el in S.elems]
            vq = [(1.0 + c) / self.pr['rho'](c) for c in Cq]
            v1 = float(vq[-1][-1]) if not isinstance(vq[-1], float) else vq[-1]
            v1 = float(np.asarray(vq[-1]).ravel()[-1])
            Dcoef = [J * self.pr['D'](c, tt) / (v * v)
                     for J, c, tt, v in zip(Jq, Cq, Tq, vq)]
            kcoef = [J * self.pr['k'](c) / v for J, c, v in zip(Jq, Cq, vq)]
            hb = 2.0 * I * H_CONV / R
            mb = 2.0 * I * HM_MASS / (R * v1)
            MT = S.assemble([(1.0 + c) * self.pr['cp'](c) for c in Cq], 'L', wlist=wl)
            AT = S.assemble(kcoef, 'DL', wlist=wc)
            scale = 1.0
        AC = S.assemble(Dcoef, 'DL', wlist=(None if self.mapping == 'affine' else wc))
        nc = S.ndof
        e_s = np.zeros(nc)
        e_s[S.surf_dof] = 1.0
        Ts = T[S.surf_dof]
        Cs = C[S.surf_dof]
        if self._const and self.mapping == 'affine':
            from scipy.linalg import lu_solve
            ra = -scale * (AT @ T) - hb * (Ts - Tinf) * e_s
            dT = lu_solve(self.MT_lu, ra)
        else:
            ra = -scale * (AT @ T) - hb * (Ts - Tinf) * e_s
            dT = np.linalg.solve(MT, ra)
        mc = -scale * (AC @ C) - mb * (Cs - Cinf) * e_s
        if self.rho_grad and self.mapping == 'affine':
            # 严格守恒形式:  d(rho_d C)/dt = (1/R^2 sigma) d_sigma( sigma rho_d D C_sigma )
            #   rho_d = rho(C)/(1+C)；质量矩阵随解变化，边界带表面密度
            #   守恒量: d/dt int rho_d C sigma dsigma = -(rho_d,s hm/R)(C_s - C_inf)
            rhod = [self.pr['rho'](c) / (1.0 + c) for c in Cq]
            Md = S.assemble(rhod, 'L')
            ACd = S.assemble([r * d for r, d in zip(rhod, Dcoef)], 'DL')
            rhod_s = float(np.asarray(rhod[-1]).ravel()[-1])
            mbd = rhod_s * HM_MASS / R
            mcd = -(ACd @ C) / (R * R) - mbd * (Cs - Cinf) * e_s
            dC = np.linalg.solve(Md, mcd)
            return np.concatenate([dT, dC])
        if self.mapping == 'affine':
            dC = self._solve_M0(mc)
        else:
            M1m = S.assemble([np.ones_like(el['sig_q']) for el in S.elems], 'L', wlist=wl)
            dC = np.linalg.solve(M1m, mc)
        return np.concatenate([dT, dC])

    def _solve_M0(self, b):
        from scipy.linalg import lu_solve
        return lu_solve(self.M0_lu, b)

    # ------------- 守恒诊断 -------------
    def totals(self, y):
        S = self.S
        T, C, Tq, Cq = self.fields(y)
        w = [el['w'] for el in S.elems]
        rcp = [self.pr['rho'](c) * self.pr['cp'](c) for c in Cq]
        E = float(sum(np.sum(wi * ri * ti) for wi, ri, ti in zip(w, rcp, Tq)))
        W = float(sum(np.sum(wi * ci) for wi, ci in zip(w, Cq)))
        return E, W, float(T[S.surf_dof]), float(C[S.surf_dof])

    # ------------- 求解 -------------
    def solve(self, t_end, t_eval=None, y0=None, rtol=1e-10, atol=1e-12,
              events=None, dense=False, max_step=np.inf, method='Radau'):
        S = self.S
        if y0 is None:
            y0 = np.concatenate([np.full(S.ndof, T0_C), np.full(S.ndof, C0)])
        if t_eval is not None:
            t_eval = np.asarray(t_eval, float)
            t_eval = t_eval[(t_eval >= 0) & (t_eval <= t_end)]
        sol = solve_ivp(self.rhs, (0.0, t_end), y0, method=method, t_eval=t_eval,
                        rtol=rtol, atol=atol, events=events, dense_output=dense,
                        max_step=max_step)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol

    # ------------- 输出 -------------
    def sample(self, y, t, r_cm=None, frame='material'):
        """frame='material': r_cm 视为初始(材料)坐标标签, sigma=r_cm/R0_cm
           frame='current' : r_cm 视为当前构形距离, sigma=r_cm/R(t)"""
        S = self.S
        if r_cm is None:
            r_cm = np.arange(0, 20.1, 1.0) / 10.0
        r_cm = np.asarray(r_cm, float)
        R = float(self.rad(t))
        if frame == 'material':
            sig = r_cm / 100.0 / R0
        else:
            sig = r_cm / 100.0 / R
        inside = sig <= 1.0 + 1e-12
        sig = np.clip(sig, 0.0, 1.0)
        T = S.interp(y[:S.ndof], sig)
        C = S.interp(y[S.ndof:], sig)
        T = np.where(inside, T, np.nan)
        C = np.where(inside, C, np.nan)
        return T, C


C_TARGET = 0.15


def drying_time(model, t_max=600000.0, rtol=1e-10, atol=1e-12, method='BDF'):
    """事件：中心节点 C(0,t)=0.15（需另行确认 C 沿 sigma 单调递减）"""
    S = model.S
    cdof = int(np.argmin(S.sig_nodes))      # sigma=0 节点
    assert S.sig_nodes[cdof] < 1e-14

    def ev(t, y):
        return y[S.ndof + cdof] - C_TARGET
    ev.terminal = True
    ev.direction = -1
    sol = model.solve(t_max, t_eval=[0.0], rtol=rtol, atol=atol, events=ev,
                      dense=True, method=method)
    if sol.t_events[0].size == 0:
        raise RuntimeError('未达到 C<0.15')
    return float(sol.t_events[0][0]), sol
