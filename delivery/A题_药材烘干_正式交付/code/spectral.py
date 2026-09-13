# -*- coding: utf-8 -*-
# ==============================================================================
# spectral.py —— 主求解器：u=(r/R)^2 坐标下的 Chebyshev 谱配置 + 自适应隐式时间积分
# ==============================================================================
# 做什么：在"径向坐标平方" u=(r/R)^2 上同时求解水分浓度 C 与温度 T 两个耦合场，四问共用。
#         换元后圆柱 Laplace 算子写成 (4/R^2) d/du ( u D d/du )，在轴心 u=0 处完全正则，
#         故轴心与表面都不需要特殊控制体处理，这是本方法区别于常规 r 坐标离散之处；
#         解在 u 上解析，谱收敛成立。时间方向交给 SciPy 的自适应变阶隐式 BDF：物性在
#         每个右端项求值处按当前 (C, T) 重算，非线性被隐式格式直接吸收，无需外层迭代。
# 输入：不读文件。几何、环境与物性全部取自 common.py：R0、h、hm、T_air(t)、C_air(t)，
#       其中环境式由附件 1 的 241 组实测数据标定，props() 按 prob 选题目附录 2/3/4 的
#       物性系数；问题 4 的收缩半径 R(t) 取 common.R_of（附件 2 实测值，cm 已换算为 m）。
# 输出：不写文件。对外提供 cheb_matrices() 与 HerbSolver 类，被各结果、绘图与校验脚本
#       （make_results.py、make_figures.py、run_all.py、verify_data.py 等）调用取解。
# 关键变量：
#   N     谱节点数（节点总数 N+1），调用方常取 64；空间误差随 N 指数下降
#   u_j   CGL 节点映射后的坐标 u=(r/R)^2，自 0（轴心）升到 1（表面），无量纲
#   prob  1/2/3/4，决定物性经验式；shrink=True 时启用 R(t)（仅问题 4）
#   rtol  自适应时间积分的相对容差，默认 1e-10，调用方多用 1e-11 ~ 1e-12
#   atol  绝对容差，默认 1e-13；C [kg/kg] 与 T [degC] 量纲不同，调用方按状态分段传数组
# ==============================================================================
import numpy as np
from scipy.integrate import solve_ivp
import common as cm


def cheb_matrices(N):
    """构造 CGL 节点 x_j = cos(j*pi/N)（自 1 降到 -1）上的 Chebyshev 微分矩阵。

    N 为多项式阶数（节点共 N+1 个）；返回 (D, D@D, x)：D 为一阶微分矩阵，
    二阶直接取 D@D——在 CGL 节点上用一阶矩阵平方仍给出精确的谱微分矩阵，
    不必另推二阶系数（二阶矩阵随结果一并给出，当前只用一阶）。
    """
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * ((-1.0) ** np.arange(N + 1))
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    # 上面的公式在 i=j 处是 0/0，对角元改用"本行其余元素之和取负"：
    # 这样 D 作用在常向量上严格得到 0（常数的导数为零），不吃浮点残差。
    D -= np.diag(D.sum(axis=1))
    return D, D @ D, x


class HerbSolver:
    """四个问题统一的水分—温度耦合求解器。

    N      谱节点数（节点总数 N+1），空间误差随 N 指数下降
    prob   1/2/3/4，决定 cm.props 用哪套物性经验式与环境口径
    shrink 是否启用附件 2 的实测半径 R(t)，仅问题 4；其余问题 R 恒为 R0
    状态向量 y = [C_0 ... C_N, T_0 ... T_N]，前段单位 [kg/kg]（干基含水率），
    后段单位 [degC]；两段共用同一套 u 节点。
    """

    def __init__(self, N, prob, shrink=False):
        self.N = N
        self.Np = N + 1
        self.prob = prob
        self.shrink = shrink
        # CGL 节点 x_j = cos(j*pi/N) 自 1 降到 -1，与 u 的换元关系是 x = 1 - 2u，
        # 故 u 自 u_0 = 0（轴心）升到 u_N = 1（表面），且 d/du = -2 d/dx。
        Dx, D2x, x = cheb_matrices(N)
        self.u = (1.0 - x) / 2.0
        self.Du = -2.0 * Dx
        self.Np = N + 1

    # ---------------- 几何与环境 ----------------
    # 半径 R(t) [m] 与烘房空气状态 T_air(t) [degC]、C_air(t) [kg/kg]；
    # solve_ivp 需要右端项返回标量系数，故这里统一 float() 化。
    def R_of(self, t):
        """当前半径 R(t) [m]：shrink=True 时取附件 2 实测半径插值，否则恒为 cm.R0。"""
        return float(cm.R_of(t)) if self.shrink else cm.R0

    def T_air(self, t):
        """烘房空气温度 T_air(t) [degC]：t <= 14400 s 用一阶惯性式，其后恒为 Tset。"""
        return float(cm.T_air(t))

    def C_air(self, t):
        """烘房空气含湿量 C_air(t) [kg/kg]：分段口径与 T_air 相同（同一 T_END）。"""
        return float(cm.C_air(t))

    # ---------------- 右端项 ----------------
    # 状态向量 y 的前 N+1 个分量是水分 C，后 N+1 个是温度 T，两者在同一组 u 节点上。
    def rhs(self, t, y):
        """右端项 dy/dt，供 solve_ivp 调用；t [s]，y 为状态向量（含义见类说明）。

        返回与 y 同形状的导数：dC/dt [kg/(kg s)] 与 dT/dt [degC/s]。
        """
        N, Np = self.N, self.Np
        C = y[:Np]
        T = y[Np:]
        R = self.R_of(t)
        rho, cp, k, D = cm.props(self.prob, C, T)

        # 水分：通量形式 g = u*D*dC/du，方程即 dC/dt = (4/R^2) dg/du
        gC = self.u * D * (self.Du @ C)
        # 表面节点不取谱微分算出的通量，而由 Robin 条件 -D dC/dr = hm(C_s - C_air)
        # 除以 dr/du = R/(2 sqrt(u)) 反解：g_N = (R/2)*hm*(C_air - C_N)，
        # 于是边界条件以物理通量的形式精确进入，不引入谱微分的边界误差。
        gC[N] = 0.5 * cm.hm * R * (self.C_air(t) - C[N])        # Robin 物理通量
        dC = (4.0 / R ** 2) * (self.Du @ gC)

        # 温度：同上，只是扩散系数换成 k，且被体积热容 rho*cp 除
        gT = self.u * k * (self.Du @ T)
        gT[N] = 0.5 * cm.h * R * (self.T_air(t) - T[N])
        dT = (4.0 / (R ** 2 * rho * cp)) * (self.Du @ gT)

        out = np.empty_like(y)
        out[:Np] = dC
        out[Np:] = dT
        return out

    # ---------------- 初值 ----------------
    # 状态向量按 [C 段, T 段] 顺序拼接，两场初始都是均匀场。
    def y0(self):
        """初始状态 [C0 ... C0, T0 ... T0]：C0 = 2.55 [kg/kg]、T0 = 28.0 [degC]。"""
        return np.concatenate([np.full(self.Np, cm.C0), np.full(self.Np, cm.T0)])

    # ---------------- 求解 ----------------
    # 这里只做参数透传：method 可选 BDF（默认，变阶自适应）、Radau、LSODA，
    # 换积分器重算同一问题的目的是交叉验证时间离散误差。
    def solve(self, t_end, t_eval=None, rtol=1e-10, atol=1e-13, method="BDF",
              events=None, max_step=np.inf, first_step=None, dense=False):
        """积分到 t_end [s]，返回 scipy 的 OdeResult（sol.t [s]、sol.y = 状态矩阵）。

        t_eval 指定输出时刻 [s]；rtol/atol 控制时间积分误差；events 用于定位烘干
        时刻（判据 C_max 下穿 0.15 kg/kg）；method 换积分器做上述交叉验证。
        """
        return solve_ivp(self.rhs, (0.0, t_end), self.y0(), method=method,
                         rtol=rtol, atol=atol, t_eval=t_eval, events=events,
                         max_step=max_step, first_step=first_step, dense_output=dense)

    # ---------------- 采样 ----------------
    # 谱解活在 u 坐标上，要还原成物理半径剖面或截面平均值都在这一节完成。
    def r_nodes(self, t):
        """当前时刻的物理半径节点 r_j = R(t)*sqrt(u_j) [m]，用于把谱解画到 r 坐标上。"""
        return np.sqrt(self.u) * self.R_of(t)

    def _bary_weights(self):
        """CGL 节点的重心插值权重 w_j（两端点折半），供 _bary 与面积权重共用。"""
        N = self.N
        w = (-1.0) ** np.arange(N + 1)
        w[0] *= 0.5
        w[N] *= 0.5
        return w

    def area_weights(self):
        """截面面积平均 Cbar = int_0^1 C(u,t) du 在 CGL 节点上的权重（缓存）。

        论文 S11.6 给出的截面平均即此定义（u=(r/R)^2 下 du 已含 2r dr 的测度）。
        C 在 u 上是 N 次多项式，故用 Gauss--Legendre 求积对 Lagrange 基函数精确
        积分一次即可得到权重向量；此后任意时刻的平均值都是 weights @ C_nodes，
        与逐时刻插值求积完全等价但快得多。

        注意：CGL 节点的**等权求和**（或若干物理半径采样点的算术平均）都不是
        面积平均——端点附近节点最密，等权重会把近表面薄层放大，结果偏低。
        """
        if getattr(self, "_area_w", None) is None:
            nq = 2 * self.N + 2                  # 被积函数为 N 次多项式，留足余量
            xg, wg = np.polynomial.legendre.leggauss(nq)
            uq = 0.5 * (xg + 1.0)                # 映射到 [0,1]，节点严格在内部
            wq = 0.5 * wg
            w = self._bary_weights()
            diff = uq[:, None] - self.u[None, :]
            basis = (w[None, :] / diff) / np.sum(w[None, :] / diff, axis=1, keepdims=True)
            self._area_w = wq @ basis
        return self._area_w

    def interp_u(self, y, u_t):
        """把状态向量按 u 变量重心插值到 u_t（谱精度），返回 (C, T) 两个数组。

        y 为某一时刻的状态向量；u_t 为目标 u 坐标，与物理半径的换算见 r_nodes。
        """
        C = y[:self.Np]
        T = y[self.Np:]
        return (self._bary(C, u_t), self._bary(T, u_t))

    def _bary(self, f, u_t):
        """单场重心插值：f 为 CGL 节点上的场值，u_t 为目标 u 坐标，返回同形状数组。"""
        u_t = np.atleast_1d(np.asarray(u_t, dtype=float))
        w = self._bary_weights()
        u = self.u
        out = np.empty(u_t.shape, dtype=float)
        for k, ut in enumerate(u_t):
            diff = ut - u
            # 命中判据用相对容差：绝对容差在 u_t 极接近某节点却不严格相等时会漏判，
            # 该分支的 1/diff 会把插值结果污染成无意义的巨值。
            hit = np.where(np.abs(diff) <= 1e-14 * max(1.0, abs(ut)))[0]
            if hit.size:
                out[k] = f[hit[0]]
                continue
            num = np.sum(w * f / diff)
            den = np.sum(w / diff)
            if den == 0.0:
                # 数学上 den = 1/prod(u_t-u_j) 恒不为零，但 N 较大且 u_t 落在
                # [0,1] 之外时，交替符号的项会在浮点下抵消到恰好 0，此时 num/den
                # 变成 inf。兜底取最近节点的值，保证输出有限。
                out[k] = f[int(np.argmin(np.abs(diff)))]
            else:
                out[k] = num / den
        return out
