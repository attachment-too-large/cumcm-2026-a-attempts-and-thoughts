# -*- coding: utf-8 -*-
"""新解法核心：坐标 u=(r/R)^2 下的 Chebyshev 谱配置 + 自适应隐式时间积分。

与 paper 的差别：
  * 空间：paper 用节点中心有限体积（二阶、待定系数三对角）；本文用 Chebyshev
    谱配置（谱精度、全局耦合）。
  * 坐标：paper 直接在 xi=r/R 上离散；本文取 u=(r/R)^2，使圆柱 Laplace 算子
        (1/r) d/dr ( r D d/dr ) = (4/R^2) d/du ( u D d/du )
    在 u=0 处完全正则（无 1/u 奇异），u=0 与 u=1 无需任何特殊处理。
  * 时间：paper 用固定步长 theta 法（CN + Rannacher 启动）；本文用自适应变阶
    隐式 BDF（SciPy solve_ivp），无固定步长假设。
  * 边界：以"通量形式"配置 —— 内部节点通量 g_j=u_j D_j (dC/du)_j 由谱微分给出，
    表面通量 g_N 直接取 Robin 条件的物理值，因此边界条件以物理通量形式精确进入。
"""
import numpy as np
from scipy.integrate import solve_ivp
import common as cm


def cheb_matrices(N):
    """CGL 节点 x_j = cos(j*pi/N) (从 1 到 -1) 上的 Chebyshev 微分矩阵。"""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * ((-1.0) ** np.arange(N + 1))
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, D @ D, x


class HerbSolver:
    """problems 1-4 统一求解器。

    N      : 谱节点数（节点总数 N+1）
    prob   : 1/2/3/4，决定物性经验式与环境口径
    shrink : 是否使用附件2 的 R(t)（仅问题4）
    """

    def __init__(self, N, prob, shrink=False):
        self.N = N
        self.Np = N + 1
        self.prob = prob
        self.shrink = shrink
        Dx, D2x, x = cheb_matrices(N)
        self.u = (1.0 - x) / 2.0          # u_0 = 0, u_N = 1
        self.Du = -2.0 * Dx               # d/du = -2 d/dx
        self.Np = N + 1

    # ---------------- 几何与环境 ----------------
    def R_of(self, t):
        return float(cm.R_of(t)) if self.shrink else cm.R0

    def T_air(self, t):
        return float(cm.T_air(t))

    def C_air(self, t):
        return float(cm.C_air(t))

    # ---------------- 右端项 ----------------
    def rhs(self, t, y):
        N, Np = self.N, self.Np
        C = y[:Np]
        T = y[Np:]
        R = self.R_of(t)
        rho, cp, k, D = cm.props(self.prob, C, T)

        # 水分：g = u*D*dC/du
        gC = self.u * D * (self.Du @ C)
        gC[N] = 0.5 * cm.hm * R * (self.C_air(t) - C[N])        # Robin 物理通量
        dC = (4.0 / R ** 2) * (self.Du @ gC)

        # 温度：g = u*k*dT/du
        gT = self.u * k * (self.Du @ T)
        gT[N] = 0.5 * cm.h * R * (self.T_air(t) - T[N])
        dT = (4.0 / (R ** 2 * rho * cp)) * (self.Du @ gT)

        out = np.empty_like(y)
        out[:Np] = dC
        out[Np:] = dT
        return out

    # ---------------- 初值 ----------------
    def y0(self):
        return np.concatenate([np.full(self.Np, cm.C0), np.full(self.Np, cm.T0)])

    # ---------------- 求解 ----------------
    def solve(self, t_end, t_eval=None, rtol=1e-10, atol=1e-13, method="BDF",
              events=None, max_step=np.inf, first_step=None, dense=False):
        return solve_ivp(self.rhs, (0.0, t_end), self.y0(), method=method,
                         rtol=rtol, atol=atol, t_eval=t_eval, events=events,
                         max_step=max_step, first_step=first_step, dense_output=dense)

    # ---------------- 采样 ----------------
    def r_nodes(self, t):
        """当前物理半径节点 r_j = R(t)*sqrt(u_j)。"""
        return np.sqrt(self.u) * self.R_of(t)

    def _bary_weights(self):
        N = self.N
        w = (-1.0) ** np.arange(N + 1)
        w[0] *= 0.5
        w[N] *= 0.5
        return w

    def interp_u(self, y, u_t):
        """把状态按 u 变量重心插值到 u_t（CGL 节点，谱精度）。"""
        C = y[:self.Np]
        T = y[self.Np:]
        return (self._bary(C, u_t), self._bary(T, u_t))

    def _bary(self, f, u_t):
        u_t = np.atleast_1d(np.asarray(u_t, dtype=float))
        w = self._bary_weights()
        u = self.u
        out = np.empty(u_t.shape, dtype=float)
        for k, ut in enumerate(u_t):
            diff = ut - u
            hit = np.where(np.abs(diff) < 1e-15)[0]
            if hit.size:
                out[k] = f[hit[0]]
                continue
            out[k] = np.sum(w * f / diff) / np.sum(w / diff)
        return out
