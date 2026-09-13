# -*- coding: utf-8 -*-
"""M1 非均匀收缩自洽模型的独立复现（干物质坐标 + 三对角有限体积）

干物质坐标 s∈[0,1]（每单位 s 含相同干物质质量），局部比容 v(C)=(1+C)/rho(C)
    I = ∫_0^1 v ds,   J(s) = ∫_0^s v dξ / I,   物理半径 r = R(t) sqrt(J)
控制方程（守恒形式）:
    dC/dt = (4I^2/R^2) d_s( (J D / v^2) d_s C )
    (1+C)cp dT/dt = (4I^2/R^2) d_s( (J k / v) d_s T )
离散（单元中心有限体积，dsm_i = Δs 为干物质质量分数）:
    Δs dC_i/dt = G^C_{i-1/2}(C_{i-1}-C_i) - G^C_{i+1/2}(C_i-C_{i+1})
    G^C_{i+1/2} = 4I^2 J_f D_f /(R^2 v_f^2 Δs)
    G^C_out     = 2 I h_m /(R v_N)      （表面）
    Δs (1+C_i)cp_i dT_i/dt = ... G^T_{i+1/2} = 4I^2 J_f k_f /(R^2 v_f Δs),  G^T_out = 2 I h / R
"""
import numpy as np
from scipy.linalg import solve_banded

R0, LZ = 0.02, 0.25
T0C, C0 = 28.0, 2.55
H_CONV, HM_MASS = 25.0, 8.0e-7


class M1Solver:
    def __init__(self, props, amb, radius, N=400, theta=0.5):
        self.pr, self.amb, self.rad = props, amb, radius
        self.N, self.theta = int(N), theta
        self.ds = 1.0 / int(N)

    def metrics(self, C):
        """由当前 C 场计算 v(单元中心)、I、J(面)"""
        v = (1.0 + C) / self.pr['rho'](C)
        I = self.ds * v.sum()
        Jf = np.concatenate([[0.0], np.cumsum(v) * self.ds / I])   # 面 0..N，J_N=1
        vf = 0.5 * (v[:-1] + v[1:])
        return v, I, Jf, vf

    def _tridiag(self, Yn, G, Gout, Yinf, cap, dt, R, th, drive=0.0):
        """解 Δs*cap*(Y^{n+1}-Y^n)/dt = th*F^{n+1} + (1-th)*F^n + drive"""
        N, ds = self.N, self.ds
        # 内部面 i=0..N-2 对应 G[i]；表面 Gout
        gl = np.concatenate([G, [Gout]])                 # 长度 N：面 1/2..N-1/2 及表面
        diag = ds * cap / dt
        up = np.zeros(N - 1)
        lo = np.zeros(N - 1)
        rhs = ds * cap / dt * Yn
        for i in range(N):
            if i > 0:
                a = gl[i - 1]                            # 左面
                diag[i] += th * a
                lo[i - 1] = -th * a
            if i < N - 1:
                b = gl[i]                                # 右面
                diag[i] += th * b
                up[i] = -th * b
            else:
                diag[i] += th * Gout
                rhs[i] += th * Gout * Yinf
        # 显式部分
        Fn = np.zeros(N)
        Fn[0] = -gl[0] * (Yn[0] - Yn[1]) if N > 1 else 0.0
        for i in range(N):
            f = 0.0
            if i > 0:
                f += gl[i - 1] * (Yn[i - 1] - Yn[i])
            if i < N - 1:
                f -= gl[i] * (Yn[i] - Yn[i + 1])
            else:
                f -= Gout * (Yn[i] - Yinf)
            Fn[i] = f
        rhs += (1 - th) * Fn + th * 0.0
        ab = np.zeros((3, N))
        ab[0, 1:] = up
        ab[1, :] = diag
        ab[2, :-1] = lo
        return solve_banded((1, 1), ab, rhs)

    def step(self, C, T, t, dt):
        R = float(self.rad(t + 0.5 * dt)) if self.rad is not None else R0
        Tinf, Cinf = self.amb(t + 0.5 * dt)
        v, I, Jf, vf = self.metrics(C)
        D = self.pr['D'](C, T)
        k = self.pr['k'](C)
        cp = self.pr['cp'](C)
        Df = 0.5 * (D[:-1] + D[1:])
        kf = 0.5 * (k[:-1] + k[1:])
        Jfc = Jf[1:self.N]                               # 内部面 1..N-1 处的 J
        coef = 4.0 * I * I / (R * R)
        GC = coef * Jfc * Df / (vf ** 2 * self.ds)
        GT = coef * Jfc * kf / (vf * self.ds)
        GCout = 2.0 * I * HM_MASS / (R * v[-1])
        GTout = 2.0 * I * H_CONV / R
        Cn = self._tridiag(C, GC, GCout, Cinf, np.ones(self.N), dt, R, self.theta)
        capT = (1.0 + C) * cp
        Tn = self._tridiag(T, GT, GTout, Tinf, capT, dt, R, self.theta)
        return np.clip(Cn, 0.0, 1e3), np.clip(Tn, -50, 300)

    def run(self, t_end, dt, t_eval=(), target=None):
        N = self.N
        C = np.full(N, C0)
        T = np.full(N, T0C)
        nt = int(np.ceil(t_end / dt))
        want = sorted(set(t_eval))
        store, wi = {}, 0
        cprev = None
        for step in range(nt):
            t = step * dt
            Cn, Tn = self.step(C, T, t, dt)
            cprev = C
            C, T = Cn, Tn
            if target is not None:
                if C[0] <= target:
                    # 线性插值定位
                    if cprev is not None and cprev[0] > target:
                        frac = (cprev[0] - target) / max(cprev[0] - C[0], 1e-30)
                        tf = (step + frac) * dt
                        return tf, C, T
            while wi < len(want) and want[wi] <= (step + 1) * dt + 1e-9:
                store[want[wi]] = (float(C[0]), float(T[0]), float(C[-1]),
                                   float(self.ds * C.sum()))
                wi += 1
            if not np.all(np.isfinite(C)):
                raise RuntimeError('M1 发散于 t=%.1f' % t)
        return store, C, T
