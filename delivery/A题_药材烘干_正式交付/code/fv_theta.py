# -*- coding: utf-8 -*-
# ==============================================================================
# fv_theta.py —— 对照格式：节点中心有限体积 + 定步长 theta 法（CN 加 Rannacher 启动）
# ==============================================================================
# 做什么：与 spectral.py 属于完全不同的离散族（局部守恒、三对角、二阶空间精度、
#         定步长时间推进），用来交叉验证主求解器，并考察"系数冻结的定步长格式"的
#         时间收敛阶。水分与温度共用同一套算子，只换 coef 与 kappa 两个系数。
# 输入：不读文件。几何、环境与物性取自 common.py（R0、hm、h、T_air、C_air、按题目
#       附录 2/3/4 取系数的 props、以及附件 2 的 R(t)）；网格 N 与步长计划 dt_plan
#       由调用方 run_all.py、verify_data.py 给出。
# 输出：不写文件。simulate() 返回 xi 与各报告时刻的 C/T 快照 [kg/kg]、[degC]，
#       find_tdry() 返回烘干时刻 [s]。
# 关键变量：
#   N       径向控制体数（节点共 N+1 个），verify_data.py 取 1600
#   dt_plan 步长计划 [(t_lim, dt), ...]：t < t_lim 时用该 dt [s]
#   theta   时间格式权重，0.5 即 Crank-Nicolson（二阶），无量纲
#   n_ram   Rannacher 启动步数，前 n_ram 步强制 theta = 1（后向 Euler）
#   hi      烘干时刻搜索上限 3.0e5 [s]；crit 判据 C_max < 0.15 [kg/kg]（干基）
# 格式（节点 xi_i = i/N；体积测度 vol_0 = h^2/8、vol_i = xi_i*h、vol_N = h/2 - h^2/8）：
#   dphi_i/dt = (1/(R^2*kappa_i*vol_i)) * (g_{i+1/2} - g_{i-1/2})
#   g_{i+1/2} = xi_{i+1/2}*coef_{i+1/2}*(phi_{i+1}-phi_i)/h（面系数取算术平均）
#   g_{-1/2} = 0（轴心对称）；g_{N+1/2} = -R*beta*(phi_N - phi_amb)（表面 Robin 通量）
#   水分：phi=C、coef=D、kappa=1、beta=hm、phi_amb=C_air；温度换成 coef=k、
#         kappa=rho*cp、beta=h、phi_amb=T_air。物性冻结在上一时间步的值。
# ==============================================================================
import numpy as np
from scipy.linalg import solve_banded
import common as cm


def simulate(N, prob, t_end, dt_plan, shrink=False, theta=0.5, n_ram=4,
             report_idx=None):
    """在给定步长计划下用 theta 法积分，返回 (xi, out_t, out_C, out_T)。

    N 为控制体数，t_end [s] 为终了时刻，prob、shrink 决定物性与是否启用 R(t)；
    theta、n_ram、dt_plan 见文件头。report_idx 是"必须记录的时刻"列表 [s]，升序；
    每步结束后把已越过的时刻快照进 out_*，快照值取自步末。

    与 fv_uniform.solve_fv 的口径差异（两者都是 O(dt) 精度，此处只作声明，不试图
    统一）：本函数把 R(t)、环境 T_air/C_air 取在**步首**，而 solve_fv 取在**步末**
    （t_new = n*dt）。差异量级为 O(dt)，故两者比较时必须各自做时间步外推，不能直接
    对比单步结果。
    """
    h = 1.0 / N
    xi = np.arange(N + 1) * h
    xif = (np.arange(N) + 0.5) * h
    # 节点中心格式的体积测度（已约去 2*pi）：节点 i 代表环形控制体
    # int_{xi_{i-1/2}}^{xi_{i+1/2}} xi dxi，故紧贴轴心与表面的半控制体各为 h^2/8。
    vol = np.empty(N + 1)
    vol[0] = h * h / 8.0
    vol[1:N] = xi[1:N] * h
    vol[N] = h / 2.0 - h * h / 8.0

    C = np.full(N + 1, cm.C0)
    T = np.full(N + 1, cm.T0)
    t = 0.0
    step = 0
    out_t, out_C, out_T = [], [], []
    ri = 0

    # 步长查表：取第一个满足 t < t_lim 的 dt，全不满足时用最后一项。调用方据此
    # 实现"先用大步长跑过预热段、再切小步长"的时间推进。
    def cur_dt(tt):
        for (t_lim, d) in dt_plan:
            if tt < t_lim:
                return d
        return dt_plan[-1][1]

    while t < t_end - 1e-9:
        dt = cur_dt(t)
        if t + dt > t_end:
            dt = t_end - t
        R = float(cm.R_of(t)) if shrink else cm.R0
        Cair = float(cm.C_air(t))
        Tair = float(cm.T_air(t))
        # Rannacher 启动：CN 对初值跳变敏感、会激发非物理的高频振荡，故前 n_ram 步
        # 改用后向 Euler（theta=1）先把快模态阻尼掉，再切回 theta=0.5 取二阶精度。
        th = 1.0 if step < n_ram else theta

        # 物性冻结在步首状态上、整步不再更新：每步于是只解两个常系数三对角系统，
        # 代价是时间精度降为一阶（这也是本格式被称为"系数冻结"的原因）。
        rho, cp, k, D = cm.props(prob, C, T)

        # 同一个 _operator 组装两个场：水分取 coef=D、容量 kappa=1、表面系数 hm、
        # 环境 C_air；温度取 coef=k、容量 kappa=rho*cp、表面系数 h、环境 T_air。
        Lc, cc = _operator(N, h, xi, xif, vol, D, np.ones(N + 1), R,
                           cm.hm, Cair)
        Lt, ct = _operator(N, h, xi, xif, vol, k, rho * cp, R,
                           cm.h, Tair)

        C = _theta_step(Lc, cc, C, dt, th)
        T = _theta_step(Lt, ct, T, dt, th)
        t += dt
        step += 1

        # 报告时刻可能落在步内部，故只要本步末时刻越过它就记录，快照取步末状态。
        while ri < len(report_idx) and t >= report_idx[ri] - 1e-9:
            out_t.append(t); out_C.append(C.copy()); out_T.append(T.copy())
            ri += 1
    return xi, out_t, out_C, out_T


def _operator(N, h, xi, xif, vol, coef, kappa, R, beta, amb):
    """组装三对角算子 (a, up, lo) 与边界常数向量 c，满足
    dphi_i/dt = a_i*phi_i + up_i*phi_{i+1} + lo_i*phi_{i-1} + c_i。

    coef 为扩散/导热系数（面上取算术平均），kappa 为容量（水分取 1、温度取 rho*cp），
    beta 为表面传递系数（水分 hm、温度 h），amb 为环境值（C_air 或 T_air），单位同 phi。
    i=0 侧由对称性 g_{-1/2} = 0 自然封闭，故第 0 行只有一个面；i=N 侧把 Robin 通量
    -R*beta*(phi_N - amb) 展开，未知量留在对角、环境值进常数项 c。
    """
    a = np.zeros(N + 1)     # 主对角
    up = np.zeros(N + 1)    # 上对角 (i -> i+1)
    lo = np.zeros(N + 1)    # 下对角 (i -> i-1)
    c = np.zeros(N + 1)

    face = 0.5 * (coef[:-1] + coef[1:])          # 面系数（算术平均）
    gf = xif * face / h                          # 内部面系数 xi_{i+1/2}*coef/h

    inv = 1.0 / (R ** 2 * kappa * vol)           # 每行都要除的体积热容因子
    # i = 0：轴心节点，只有外侧一个面（g_{-1/2} = 0），故只有一个非对角元
    a[0] = -inv[0] * gf[0]
    up[0] = inv[0] * gf[0]
    # i = 1..N-1：内部节点，左右各一个面，对角是两者的负和
    i = np.arange(1, N)
    a[i] = -inv[i] * (gf[i] + gf[i - 1])
    up[i] = inv[i] * gf[i]
    lo[i] = inv[i] * gf[i - 1]
    # i = N：表面节点，左侧一个面再加上 Robin 通量项 R*beta，环境值进常数项 c
    a[N] = -inv[N] * (gf[N - 1] + R * beta)
    lo[N] = inv[N] * gf[N - 1]
    c[N] = inv[N] * R * beta * amb
    return (a, up, lo), c


def _theta_step(op, c, phi, dt, th):
    """推进一个时间步，解 (I - dt*th*L) phi_new = phi + dt*((1-th)*L phi + c)。

    op 为 _operator 返回的 (a, up, lo)，th 为 theta（1 为后向 Euler、0.5 为 CN）。
    显式部分只做乘加，故右端不需要额外三对角求解；矩阵按 solve_banded 要求的
    (1, 1) 带宽装箱，不必构造稠密矩阵。
    """
    a, up, lo = op
    n = len(phi)
    lhs_a = 1.0 - dt * th * a
    lhs_up = -dt * th * up[:-1]
    lhs_lo = -dt * th * lo[1:]
    rhs = phi + dt * (1.0 - th) * (a * phi + up * np.r_[phi[1:], 0.0]
                                   + lo * np.r_[0.0, phi[:-1]] + c) + dt * th * c
    ab = np.zeros((3, n))
    ab[0, 1:] = lhs_up
    ab[1, :] = lhs_a
    ab[2, :-1] = lhs_lo
    return solve_banded((1, 1), ab, rhs)


def find_tdry(N, prob, dt_plan, shrink=False, hi=300000.0, crit=0.15):
    """从 t=0 单次前向推进，返回全场最大含水率首次下穿判据 crit 的时刻 [s]。

    N、prob、dt_plan、shrink 同 simulate；hi 为搜索上限 3.0e5 [s]，到上限仍未
    下穿则返回 nan；crit 为烘干判据 0.15 [kg/kg]（干基含水率）。判据取 max_r C
    而不是平均值：轴心处水分最后才干，最大值达标即整个截面达标。
    """
    h = 1.0 / N
    xi = np.arange(N + 1) * h
    xif = (np.arange(N) + 0.5) * h
    # 体积测度同 simulate：节点中心格式下由 int_{xi-1/2}^{xi+1/2} xi dxi 得到。
    vol = np.empty(N + 1)
    vol[0] = h * h / 8.0
    vol[1:N] = xi[1:N] * h
    vol[N] = h / 2.0 - h * h / 8.0

    C = np.full(N + 1, cm.C0)
    T = np.full(N + 1, cm.T0)
    t = 0.0
    step = 0
    Cprev, tprev = None, None
    while t < hi - 1e-9:
        dt = dt_plan[-1][1]
        for (t_lim, d) in dt_plan:
            if t < t_lim:
                dt = d
                break
        # 时间离散口径与 simulate 不同：半径取步末（与新解同点），环境取步中点值。
        R = float(cm.R_of(t + dt)) if shrink else cm.R0   # 取步末半径
        Cair = float(cm.C_air(t + 0.5 * dt))              # 环境取步中点值
        Tair = float(cm.T_air(t + 0.5 * dt))
        # 前 4 步 Rannacher 启动（与 simulate 的默认 n_ram=4 同口径）
        th = 1.0 if step < 4 else 0.5
        rho, cp, k, D = cm.props(prob, C, T)
        Lc, cc = _operator(N, h, xi, xif, vol, D, np.ones(N + 1), R, cm.hm, Cair)
        Lt, ct = _operator(N, h, xi, xif, vol, k, rho * cp, R, cm.h, Tair)
        Cnew = _theta_step(Lc, cc, C, dt, th)
        T = _theta_step(Lt, ct, T, dt, th)
        # 记下本步起点 (tprev, Cprev)：判据下穿时靠这两点做插值求交。
        cprev = float(np.max(C))
        C = Cnew
        tprev, Cprev = t, cprev
        t += dt
        step += 1
        cnew = float(np.max(C))
        if cnew < crit:
            # Cmax 在 [tprev, t] 上线性插值求交
            if Cprev is None or Cprev == cnew:
                return t
            return tprev + (crit - Cprev) / (cnew - Cprev) * (t - tprev)
    return float("nan")
