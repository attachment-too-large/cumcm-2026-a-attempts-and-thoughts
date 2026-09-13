# -*- coding: utf-8 -*-
"""
方法B（独立实现）：材料坐标 + 偶数Legendre谱Galerkin（弱形式 / 自然Robin边界）

物理模型（圆柱，径向导热 + 干基水分扩散）
  Eulerian:   rho*cp*(dT/dt + v*dT/dr) = (1/r) d/dr ( k r dT/dr ),   v = (Rdot/R) r
              dC/dt + v*dC/dr           = (1/r) d/dr ( D r dC/dr )
  令 sigma = r/R(t) （均匀收缩的材料/Lagrangian坐标），两个对流项对均匀仿射收缩**精确抵消**：
              rho*cp*dT/dt = (k/R^2) (1/sigma) d/dsigma ( sigma dT/dsigma )
              dC/dt        = (D/R^2) (1/sigma) d/dsigma ( sigma dC/dsigma )
  边界(sigma=1)： -k*(1/R)*dT/dsigma = h*(Ts - Tinf) ; -D*(1/R)*dC/dsigma = hm*(Cs - Cinf)
  柱心 sigma=0 自然正则（弱形式无需特殊处理）

空间离散： T(sigma,t)=sum_j a_j P_{2j}(2sigma-1)  （偶数Legendre，tau=2sigma-1）
  质量矩阵 M_ij = int_0^1 phi_i phi_j sigma dsigma = delta_ij / (2(4j+1))  （解析对角！）
  检验函数取弱形式 -> Robin自然边界；phi_j(1)=P_{2j}(1)=1 => 表面值 = sum_j a_j
  守恒律：取 phi_0=1 得 d/dt int rho cp T sigma dsigma = -(h/R)(Ts-Tinf) 精确成立
"""
import os
import numpy as np
import openpyxl
from scipy.special import eval_legendre, roots_legendre
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------------
# 常数
# ----------------------------------------------------------------------------
R0 = 0.02          # m，初始半径
LZ = 0.25          # m，长度
H_CONV = 25.0      # W/(m^2 K)
HM_MASS = 8.0e-7   # m/s
T0_C = 28.0        # 初始温度 ℃
C0 = 2.55          # 初始干基含水率 kg/kg
T_HOLD = 14400.0   # 预热平衡段结束时间（附件1覆盖范围）
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '附件')


# ----------------------------------------------------------------------------
# 数据读取
# ----------------------------------------------------------------------------
def _read_sheet(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb.active
    rows = [r for r in ws.iter_rows(values_only=True)]
    hdr, rows = rows[0], rows[1:]
    rows = [r for r in rows if r[0] is not None]
    return hdr, rows


def load_att1(path=None):
    path = path or os.path.join(DATA_DIR, '附件1.xlsx')
    _, rows = _read_sheet(path)
    t = np.array([float(r[0]) for r in rows])
    T = np.array([float(r[1]) for r in rows])
    C = np.array([float(r[2]) for r in rows])
    return t, T, C


def load_att2(path=None):
    path = path or os.path.join(DATA_DIR, '附件2.xlsx')
    _, rows = _read_sheet(path)
    t = np.array([float(r[0]) for r in rows])
    R = np.array([float(r[1]) for r in rows]) / 100.0   # cm -> m
    return t, R


class Ambient:
    """烘房环境 T_inf(t)[℃], C_inf(t)[kg/kg]。

    0 <= t <= t_hold : 附件1 实测曲线插值（'linear' 或 'pchip'）
    t >  t_hold      : 保持末端值（即恒温干燥段的设定值）
    """

    def __init__(self, t, T, C, mode='linear', t_hold=T_HOLD, smooth=None):
        t = np.asarray(t, float)
        T = np.asarray(T, float)
        C = np.asarray(C, float)
        if smooth:                      # 可选 Savitzky-Golay 平滑（敏感性分析）
            from scipy.signal import savgol_filter
            T = savgol_filter(T, smooth, 2)
            C = savgol_filter(C, smooth, 2)
        self.t_hold = float(t_hold)
        self.mode = mode
        m = t <= self.t_hold
        self.t = t[m]
        self.T = T[m]
        self.C = C[m]
        self.T_hold = float(self.T[-1])
        self.C_hold = float(self.C[-1])
        if mode == 'pchip':
            self._fT = PchipInterpolator(self.t, self.T, extrapolate=False)
            self._fC = PchipInterpolator(self.t, self.C, extrapolate=False)
        elif mode == 'linear':
            self._fT = self._fC = None
        else:
            raise ValueError(mode)

    def __call__(self, t):
        scalar = np.isscalar(t) or np.ndim(t) == 0
        tt = np.atleast_1d(np.asarray(t, float))
        if self.mode == 'pchip':
            T = np.where(tt <= self.t_hold, self._fT(np.minimum(tt, self.t_hold)), self.T_hold)
            C = np.where(tt <= self.t_hold, self._fC(np.minimum(tt, self.t_hold)), self.C_hold)
        else:
            T = np.interp(tt, self.t, self.T, right=self.T_hold)
            C = np.interp(tt, self.t, self.C, right=self.C_hold)
        T = np.where(tt > self.t_hold, self.T_hold, T)
        C = np.where(tt > self.t_hold, self.C_hold, C)
        if scalar:
            return float(T[0]), float(C[0])
        return T, C


class Radius:
    """R(t) [m]。fixed=True 时半径恒定 R0（问题1-3）；否则用附件2插值（问题4）。"""

    def __init__(self, fixed=True, t=None, R=None, mode='pchip'):
        self.fixed = fixed
        if fixed:
            return
        self.t = np.asarray(t, float)
        self.R = np.asarray(R, float)
        self.mode = mode
        if mode == 'pchip':
            self._f = PchipInterpolator(self.t, self.R)

    def __call__(self, t):
        if self.fixed:
            return R0 if np.isscalar(t) or np.ndim(t) == 0 else np.full(np.shape(t), R0)
        t = np.asarray(t, float)
        if self.mode == 'pchip':
            R = self._f(np.clip(t, self.t[0], self.t[-1]))
        else:
            R = np.interp(t, self.t, self.R)
        if np.isscalar(t) or np.ndim(t) == 0:
            return float(R)
        return R


# ----------------------------------------------------------------------------
# 谱基
# ----------------------------------------------------------------------------
class EvenLegendre:
    """phi_j(sigma) = P_{2j}(2 sigma - 1),  sigma in [0,1],  j=0..N-1"""

    def __init__(self, N, nq=None, check=True):
        self.N = int(N)
        self.nq = int(nq) if nq else 2 * self.N + 14
        xi, w = roots_legendre(self.nq)
        self.xi = xi
        self.sig = 0.5 * (1.0 + xi)                 # sigma 节点
        self.wq = 0.25 * w * (1.0 + xi)             # int_0^1 g sigma dsigma = sum wq*g
        P = np.empty((self.N, self.nq))
        G = np.empty((self.N, self.nq))
        for j in range(self.N):
            n = 2 * j
            P[j] = eval_legendre(n, xi)
            c = np.zeros(n + 1)
            c[n] = 1.0
            d = np.polynomial.legendre.legder(c)
            G[j] = 2.0 * np.polynomial.legendre.legval(xi, d)   # d/dsigma
        self.P = P
        self.G = G
        j = np.arange(self.N)
        self.Mana = 1.0 / (2.0 * (4.0 * j + 1.0))   # 解析质量矩阵（对角）
        self.M = np.diag(self.Mana)
        self.A0 = (G * self.wq) @ G.T               # int sigma phi_i' phi_j' dsigma
        self.b = np.ones(self.N)                    # phi_j(1)=1
        self.ones = np.ones(self.N)
        if check:
            Mq = (P * self.wq) @ P.T
            e1 = np.max(np.abs(Mq - self.M))
            # 解析恒等式: int_0^1 sigma phi_i' phi_j' dsigma 对 i=j=0 应为0 等
            assert e1 < 1e-12, ('mass matrix mismatch', e1)
            self.selfcheck = e1

    # 取值： sigma -> 值（1D 系数或多列）
    def eval(self, a, sigma):
        xi = 2.0 * np.asarray(sigma, float) - 1.0
        V = np.polynomial.legendre.legvander(xi, 2 * self.N - 2)[:, 0::2]
        a = np.asarray(a, float)
        if a.ndim == 1:
            return V @ a
        return V @ a

    def nodes_eval(self, a):
        a = np.asarray(a, float)
        return self.P.T @ a


# ----------------------------------------------------------------------------
# 物性
# ----------------------------------------------------------------------------
def props_q1():
    """附录2。D = 7e-9 * exp(-0.89/C)  —— 指数为分式 -0.89 除以 C（已由 A题.pdf 原版面渲染核对）"""
    return dict(
        name='Q1(附录2)',
        rho=lambda C: np.full_like(C, 820.0),
        cp=lambda C: np.full_like(C, 2600.0),
        k=lambda C: np.full_like(C, 0.36),
        drho=lambda C: np.zeros_like(C),
        D=lambda C, T: 7e-9 * np.exp(-0.89 / np.maximum(C, 1e-8)),
    )


def props_q23():
    """附录3。D = 2.4e-3 * exp(-0.45/C) * exp(-3850/T)"""
    return dict(
        name='Q2/Q3(附录3)',
        rho=lambda C: 650.0 + 128.0 * C,
        cp=lambda C: 1450.0 + 2736.0 * C / (C + 1.0),
        k=lambda C: 0.21 + 0.38 * C / (C + 1.0),
        drho=lambda C: np.full_like(C, 128.0),
        D=lambda C, T: 2.4e-3 * np.exp(-0.45 / np.maximum(C, 1e-8))
        * np.exp(-3850.0 / (T + 273.15)),
    )


def props_q4():
    """附录4。D = 4.2e-4 * exp(-0.30/C) * exp(-3850/T)"""
    return dict(
        name='Q4(附录4)',
        rho=lambda C: 760.0 + 90.0 * C,
        cp=lambda C: 1850.0 + 2150.0 * C / (C + 1.0),
        k=lambda C: 0.12 + 0.20 * C / (C + 1.0),
        drho=lambda C: np.full_like(C, 90.0),
        D=lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 1e-8))
        * np.exp(-3850.0 / (T + 273.15)),
    )


# ----------------------------------------------------------------------------
# 模型
# ----------------------------------------------------------------------------
class Model:
    """状态 y = [a(N); c(N)]，a/c 为 T[℃] / C[kg/kg] 的偶数Legendre系数。"""

    def __init__(self, N=28, props=None, amb=None, radius=None, nq=None):
        self.B = EvenLegendre(N, nq)
        self.pr = props or props_q23()
        self.amb = amb
        self.rad = radius or Radius(fixed=True)
        self.N = self.B.N
        self.balance = []          # 记录守恒量

    # --- 由系数得到节点物理量 ---
    def fields(self, y):
        B = self.B
        a = y[:B.N]
        c = y[B.N:]
        Tq = B.P.T @ a
        Cq = B.P.T @ c
        return a, c, Tq, Cq

    def rhs(self, t, y):
        B = self.B
        a, c, Tq, Cq = self.fields(y)
        R = self.rad(t)
        Tinf, Cinf = self.amb(t)
        rho_cp = self.pr['rho'](Cq) * self.pr['cp'](Cq)
        kq = self.pr['k'](Cq)
        Dq = self.pr['D'](Cq, Tq)
        w = B.wq
        MT = (B.P * (w * rho_cp)) @ B.P.T
        AT = (B.G * (w * kq)) @ B.G.T
        AC = (B.G * (w * Dq)) @ B.G.T
        Ts = a.sum()
        Cs = c.sum()
        ra = -(AT @ a) / (R * R) - (H_CONV / R) * (Ts - Tinf) * B.ones
        rc = -(AC @ c) / (R * R) - (HM_MASS / R) * (Cs - Cinf) * B.ones
        da = np.linalg.solve(MT, ra)
        dc = rc * (2.0 * (4.0 * np.arange(B.N) + 1.0))   # M_C = diag(1/(2(4j+1)))
        return np.concatenate([da, dc])

    # --- 守恒量（用于独立校验） ---
    def integrals(self, y):
        """返回 (总焓指标 int rho cp T sigma dsigma, 总水分 int C sigma dsigma,
                 表面温度, 表面浓度)"""
        B = self.B
        a, c, Tq, Cq = self.fields(y)
        w = B.wq
        rho_cp = self.pr['rho'](Cq) * self.pr['cp'](Cq)
        E = float(np.sum(w * rho_cp * Tq))
        W = float(np.sum(w * Cq))
        return E, W, float(a.sum()), float(c.sum())

    # --- 求解 ---
    def solve(self, t_end, t_eval=None, y0=None, rtol=1e-10, atol=1e-12,
              method='Radau', events=None, dense=False, max_step=np.inf):
        B = self.B
        if y0 is None:
            y0 = np.concatenate([np.concatenate([[T0_C], np.zeros(B.N - 1)]),
                                 np.concatenate([[C0], np.zeros(B.N - 1)])])
        if t_eval is None:
            t_eval = np.linspace(0.0, t_end, 201)
        t_eval = np.asarray(t_eval, float)
        t_eval = t_eval[(t_eval >= 0) & (t_eval <= t_end)]
        sol = solve_ivp(self.rhs, (0.0, t_end), y0, method=method, t_eval=t_eval,
                        rtol=rtol, atol=atol, events=events, dense_output=dense,
                        max_step=max_step)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol

    # --- 输出到物理网格 ---
    def sample(self, y, t, r_cm=None):
        """在 r_cm（当前构形，cm）处取值；默认 0..2.0 step 0.1（材料/Lagrangian标签）"""
        B = self.B
        if r_cm is None:
            r_cm = np.arange(0, 20.1, 1.0) / 10.0
        r_cm = np.asarray(r_cm, float)
        R = self.rad(t)
        T = B.eval(y[:B.N], r_cm / 100.0 / R)
        C = B.eval(y[B.N:], r_cm / 100.0 / R)
        return T, C


# ----------------------------------------------------------------------------
# 问题3/4：求烘干结束时间（各处 C < 0.15）
# ----------------------------------------------------------------------------
C_TARGET = 0.15


def drying_time(model, t_max=400000.0, rtol=1e-10, atol=1e-12, verbose=True):
    """事件：中心处 C(0,t) = 0.15（圆柱干燥中心最慢，先确认单调性）。"""
    B = model.B

    def ev(t, y):
        c = y[B.N:]
        return c.sum() - C_TARGET      # 中心浓度 = sum a_j (phi_j(0)=(-1)^j 不是1!)
    # 注意：phi_j(0)=P_{2j}(-1)=1，因此 sum 也是 sigma=0 处的值。核对：
    #   P_n(-1)=(-1)^n, n=2j 偶数 => 1。所以 c.sum() = C(sigma=0) 正确。
    ev.terminal = True
    ev.direction = -1
    sol = model.solve(t_max, t_eval=[0.0], rtol=rtol, atol=atol, events=ev, dense=True)
    if sol.t_events[0].size == 0:
        raise RuntimeError('未在 %.0f s 内达到 C<0.15' % t_max)
    return float(sol.t_events[0][0])
