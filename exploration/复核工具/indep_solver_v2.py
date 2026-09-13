# -*- coding: utf-8 -*-
"""自研二阶求解器（用于独立复核的精度提升）

设计要点（均为独立思考后针对"一阶"瓶颈的对症改进）：

1. **时间二阶**：θ 法取 θ=1/2（Crank–Nicolson），配 **Rannacher 启动**
   （前 n_ran 步取 θ=1）。纯 CN 在 t=0 的表面温度/浓度跃变下，最高频空间模态的
   放大因子趋于 −1、几乎不衰减，会产生非物理振荡；Rannacher 用两步全隐式先把
   高频模态压掉再切回 CN，既消振荡又保住二阶。
2. **界面系数可选对数平均** log-mean(D_i,D_{i+1})：D ∝ exp(−a/C) 在 C→0 时跨十余个
   数量级，算术平均高估表层通量；对数平均更贴合指数型剖面（D 的几何平均性质），
   且与调和平均同样稳定。
3. **非线性半步中点线性化**：系数取 (t_n) 与 (t_{n+1} 迭代值) 的中点，使非线性项
   的局部截断误差与 CN 同阶。

**正确性保证**：隐式算子与显式残差共用同一个面通量算子 `_flux_op`，
因此 θ=1、face='arith'、n_ran=∞ 时本类与 `indep_solver.FV` **逐位相同**（见自检）。
"""
import os
import sys
import numpy as np
from scipy.linalg import solve_banded

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_solver import props, Rof, Tinf, Cinf, T0, C0


def harm(a, b):
    return 2.0 * a * b / (a + b + 1e-300)


def logmean(a, b):
    """对数平均 (a−b)/ln(a/b)；a≈b 时退化为算术平均。"""
    a = np.maximum(np.asarray(a, float), 1e-300)
    b = np.maximum(np.asarray(b, float), 1e-300)
    out = 0.5 * (a + b)
    r = np.abs(np.log(a / b))
    m = r > 1e-8
    if np.any(m):
        out[m] = (a[m] - b[m]) / np.log(a[m] / b[m])
    return out


class FV2:
    """节点中心轴对称 FV；θ 法（θ=1/2 为 CN）+ Rannacher 启动。"""

    def __init__(self, mode, shrink, N=400, theta=0.5, n_ran=2, face='log',
                 tol=1e-11, maxit=50, h_conv=25.0, km=8e-7):
        self.mode, self.shrink, self.N = mode, shrink, N
        self.theta, self.n_ran, self.face = theta, n_ran, face
        self.tol, self.maxit = tol, maxit
        self.h_conv, self.km = h_conv, km
        xi = np.linspace(0.0, 1.0, N + 1)
        self.xi = 0.5 * (xi[:-1] + xi[1:])
        self.xf = xi[1:-1]
        lo = np.concatenate([[0.0], self.xf])
        hi = np.concatenate([self.xf, [1.0]])
        self.wV = 0.5 * (hi ** 2 - lo ** 2)
        self.Af = self.xf
        self.dface = np.diff(self.xi)
        self.nstep = 0

    def _faceavg(self, y):
        if self.face == 'log':
            return logmean(y[:-1], y[1:])
        if self.face == 'harm':
            return harm(y[:-1], y[1:])
        return 0.5 * (y[:-1] + y[1:])

    @staticmethod
    def _flux_op(Y, g, Gout, Yinf):
        """与 indep_solver.FV._solve 完全一致的离散通量算子。

        g    : 长度 N-1 的内部面系数
        q[i] = 从节点 i 经其上侧面流出的通量；q[N-1] 为外表面
        返回 F[i] = q[i-1] - q[i]，且 F[0] = -q[0]
        """
        N = len(Y)
        q = np.empty(N)
        q[:N - 1] = g * (Y[:N - 1] - Y[1:])
        q[N - 1] = Gout * (Y[N - 1] - Yinf)
        F = np.empty(N)
        F[0] = -q[0]
        F[1:] = q[:-1] - q[1:]
        return F

    def _solve_theta(self, Yn, cap, g, Gout, Yinf, dt, th):
        N = self.N
        gg = np.concatenate([g, [Gout]])              # 长度 N
        lo = np.concatenate([[0.0], g])               # 长度 N
        diag = cap / dt + th * (lo + gg)
        ab = np.zeros((3, N))
        ab[0, 1:] = -th * gg[:-1]
        ab[1, :] = diag
        ab[2, :-1] = -th * lo[1:]
        Fn = self._flux_op(Yn, g, Gout, Yinf)
        rhs = cap / dt * Yn + (1.0 - th) * Fn
        rhs[N - 1] += th * Gout * Yinf
        return solve_banded((1, 1), ab, rhs)

    def step(self, T, C, t, dt, theta=None):
        R = Rof(t + 0.5 * dt) if self.shrink else 0.02
        V = R ** 2 * self.wV
        Tinf_o, Cinf_o = Tinf(t), Cinf(t)
        th = (1.0 if self.nstep < self.n_ran else self.theta) if theta is None else theta
        Tg, Cg = T.copy(), C.copy()
        for _ in range(self.maxit):
            Tm = 0.5 * (T + Tg)
            Cm = 0.5 * (C + Cg)
            rho, cp, k, D = props(self.mode, Cm, Tm)
            kf = self._faceavg(k)
            Df = self._faceavg(D)
            GT = self.Af * kf / self.dface
            GC = self.Af * Df / self.dface
            Tn = self._solve_theta(T, rho * cp * V, GT, R * self.h_conv, Tinf_o, dt, th)
            Cn2 = self._solve_theta(C, V.copy(), GC, R * self.km, Cinf_o, dt, th)
            d = max(float(np.max(np.abs(Tn - Tg))), float(np.max(np.abs(Cn2 - Cg))))
            Tg, Cg = Tn, Cn2
            if d < self.tol:
                break
        self.nstep += 1
        return Tg, Cg

    def till(self, thr, t_end, dt):
        N = self.N
        T = np.full(N, T0)
        C = np.full(N, C0)
        t = 0.0
        while t < t_end - 1e-12:
            h = min(dt, t_end - t)
            T, C = self.step(T, C, t, h)
            t += h
            if C.max() < thr:
                return t, C, T
        return None, C, T
