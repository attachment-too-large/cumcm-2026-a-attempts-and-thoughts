# -*- coding: utf-8 -*-
# ==============================================================================
# fv_uniform.py —— 独立参照解：均匀 r 网格 + 守恒型有限体积 + 后向 Euler + Picard
# ==============================================================================
# 做什么：用一套与主求解器完全不同的离散再解同一个初边值问题，两法一致即说明解的是
#         同一个 PDE、且各自的离散误差都已在目标精度以下。物性在每个时间步内用不动点
#         (Picard) 迭代刷新，故后向 Euler 的非线性不需要牛顿法。
# 输入：不读文件。几何、环境与物性取自 common.py（R0、hm、h、T_air、C_air、按题目
#       附录 2/3/4 取系数的 props、附件 2 的 R(t)）；网格数 N、步长 dt [s] 与报告时刻
#       由调用方 run_all.py、make_figures.py 给出。
# 输出：不写文件。solve_fv() 返回 (xi, C_tab, T_tab)，sample_xi() 做线性插值采样。
# 关键变量：
#   N        径向控制体数（节点共 N+1 个），均匀网格 xi_i = i/N
#   dt       时间步长 [s]；后向 Euler 等步长推进，时间方向一阶
#   n_picard 每步的物性迭代次数，默认 3（系数冻结在上一轮迭代值上）
#   R(t)     不收缩时恒为 R0 = 2.000e-2 [m]；hm 8.0e-7 [m/s]；h 25.0 [W/(m2 K)]
#   体积测度 vol_0 = h^2/8、vol_i = xi_i*h、vol_N = h/2 - h^2/8，与 fv_theta.py 同
# 与主求解器的三处刻意差别：坐标用均匀 r 而非 u = xi^2；离散用守恒型有限体积（面通量
#   加体积测度）而非谱微分；时间用等步长后向 Euler（一阶）而非自适应变阶 BDF。
# ==============================================================================
import numpy as np
import common as cm


def solve_fv(N, dt, t_end, prob, t_report, shrink=False, n_picard=3):
    """积分到 t_end [s] 并采样，返回 (xi, C_tab, T_tab)。

    N 为网格数，dt [s] 为等步长，prob、shrink 决定物性与是否启用 R(t)（见文件头），
    n_picard 为每步的物性迭代次数；t_report 为报告时刻数组 [s]。
      xi    —— 无量纲节点坐标 i/N，长度 N+1；
      C_tab —— 各报告时刻的含水率 [kg/kg]，形状 (len(t_report), N+1)；
      T_tab —— 各报告时刻的温度 [degC]，形状同上；未落在整步上的时刻保持 nan。
    """
    h = 1.0 / N
    xi = np.arange(N + 1) * h
    xif = (np.arange(N) + 0.5) * h           # 内部面 xi_{i+1/2}, i=0..N-1

    # 节点中心格式的体积测度（已约去 2*pi）：节点 i 代表环形控制体
    # int_{xi_{i-1/2}}^{xi_{i+1/2}} xi dxi，故紧贴轴心与表面的半控制体各为 h^2/8。
    vol = np.empty(N + 1)
    vol[0] = h * h / 8.0
    vol[1:N] = xi[1:N] * h
    vol[N] = h / 2.0 - h * h / 8.0

    C = np.full(N + 1, cm.C0)
    T = np.full(N + 1, cm.T0)

    nstep = int(round(t_end / dt))
    trep = np.asarray(t_report, dtype=float)
    # 报告时刻按"步号 = round(t/dt)"定位，落在半步上的时刻取不到值，其行保持 nan。
    idx = {int(round(t / dt)): k for k, t in enumerate(trep)}
    Ctab = np.full((len(trep), N + 1), np.nan)
    Ttab = np.full((len(trep), N + 1), np.nan)

    def R_of(t):
        return float(cm.R_of(t)) if shrink else cm.R0

    # 半径与环境一律取步末 t_new：后向 Euler 的解就落在 t_new，系数与状态取同一时刻。
    # fv_theta.simulate 取的是步首，两者相差 O(dt)，故对比时必须各自做时间步外推。
    for n in range(1, nstep + 1):
        t_new = n * dt
        t_old = (n - 1) * dt
        R = R_of(t_new)
        Cair = float(cm.C_air(t_new))
        Tair = float(cm.T_air(t_new))

        # Picard 迭代：物性 D、k、rho、cp 依赖未知的 C、T，于是先把它们冻结在上一轮
        # 迭代值上（首轮用上一时间步的解）解出新的 C、T，再刷新物性，如此重复 n_picard 次。
        Cnew, Tnew = C.copy(), T.copy()
        for _ in range(n_picard):
            rho, cp, k, D = cm.props(prob, Cnew, Tnew)
            # ---- 水分三对角：隐式项整理成 (I - dt*L) C_new = b ----
            # 内点的两个面各占一个非对角元；轴心侧因 g_{-1/2} = 0 只剩一个面；
            # 表面节点的 Robin 通量把 hm 折进对角，把 hm*C_air 移到右端 b。
            A = np.zeros((3, N + 1))
            b = C.copy()
            for i in range(N + 1):
                diag = 1.0
                if i == 0:
                    gf = xif[0] * 0.5 * (D[0] + D[1]) / h
                    coef = dt * gf / (R ** 2 * vol[0])
                    diag += coef
                    A[0, i] = -coef
                elif i == N:
                    gf = xif[N - 1] * 0.5 * (D[N - 1] + D[N]) / h
                    coef_i = dt * gf / (R ** 2 * vol[N])
                    coef_b = dt * (R * cm.hm) / (R ** 2 * vol[N])
                    diag += coef_i + coef_b
                    A[2, i] = -coef_i          # 系数 C_{N-1} 在次对角
                    b[i] += coef_b * Cair
                else:
                    gf_p = xif[i] * 0.5 * (D[i] + D[i + 1]) / h
                    gf_m = xif[i - 1] * 0.5 * (D[i - 1] + D[i]) / h
                    cp_ = dt * gf_p / (R ** 2 * vol[i])
                    cm_ = dt * gf_m / (R ** 2 * vol[i])
                    diag += cp_ + cm_
                    A[0, i] = -cp_
                    A[2, i] = -cm_
                A[1, i] = diag
            Cnew = _thomas(A, b)
            # ---- 温度三对角 ----
            # 结构与水分完全相同，只把 D 换成 k、并在体积测度上再除以节点热容 rho*cp。
            A = np.zeros((3, N + 1))
            b = T.copy()
            for i in range(N + 1):
                diag = 1.0
                if i == 0:
                    gf = xif[0] * 0.5 * (k[0] + k[1]) / h
                    coef = dt * gf / (R ** 2 * rho[0] * cp[0] * vol[0])
                    diag += coef
                    A[0, i] = -coef
                elif i == N:
                    gf = xif[N - 1] * 0.5 * (k[N - 1] + k[N]) / h
                    coef_i = dt * gf / (R ** 2 * rho[N] * cp[N] * vol[N])
                    coef_b = dt * (R * cm.h) / (R ** 2 * rho[N] * cp[N] * vol[N])
                    diag += coef_i + coef_b
                    A[2, i] = -coef_i          # 系数 T_{N-1} 在次对角
                    b[i] += coef_b * Tair
                else:
                    gf_p = xif[i] * 0.5 * (k[i] + k[i + 1]) / h
                    gf_m = xif[i - 1] * 0.5 * (k[i - 1] + k[i]) / h
                    cp_ = dt * gf_p / (R ** 2 * rho[i] * cp[i] * vol[i])
                    cm_ = dt * gf_m / (R ** 2 * rho[i] * cp[i] * vol[i])
                    diag += cp_ + cm_
                    A[0, i] = -cp_
                    A[2, i] = -cm_
                A[1, i] = diag
            Tnew = _thomas(A, b)

        C, T = Cnew, Tnew
        # 报告时刻按步号命中，快照直接存该步末的整条剖面。
        if n in idx:
            Ctab[idx[n]] = C
            Ttab[idx[n]] = T
    return xi, Ctab, Ttab


def _thomas(A, b):
    """追赶法（Thomas 算法）解三对角系统 A x = b，返回解向量 x。

    A 为 (3, n) 的紧凑带：A[0] = 上对角（偏移 +1）、A[1] = 主对角、A[2] = 下对角
    （偏移 -1）。不选主元：本模块的矩阵由隐式扩散加正的对角项构成，主对角严格大于
    两个非对角元之和（严格对角占优），消元过程中不会出现零主元。
    """
    n = len(b)
    cp = np.zeros(n)
    dp = np.zeros(n)
    cp[0] = A[0, 0] / A[1, 0]
    dp[0] = b[0] / A[1, 0]
    for i in range(1, n):
        m = A[1, i] - A[2, i] * cp[i - 1]
        cp[i] = A[0, i] / m
        dp[i] = (b[i] - A[2, i] * dp[i - 1]) / m
    x = np.zeros(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def sample_xi(xi, vals, x_target):
    """把节点值 vals 线性插值到目标位置 x_target（无量纲坐标 i/N）。

    x_target 可为标量或数组，返回同形状结果；超出 xi 范围时取端点值（np.interp 默认）。
    """
    return np.interp(np.asarray(x_target, float), xi, vals)
