# -*- coding: utf-8 -*-
"""对照格式：经典节点中心有限体积 + 定步长 theta 法（CN + Rannacher 启动）。

本模块与 spectral.py 属于完全不同的离散族（局部守恒、三对角、二阶精度），
用于交叉验证主求解器；同时用于考察"定步长显式系数冻结"格式的时间收敛阶。

节点 xi_i = i/N；体积测度 vol_0=h^2/8, vol_i=xi_i*h, vol_N=h/2-h^2/8；
dphi_i/dt = (1/(R^2*kappa_i*vol_i)) (g_{i+1/2}-g_{i-1/2})
g_{i+1/2} = xi_{i+1/2} * coef_{i+1/2} * (phi_{i+1}-phi_i)/h   (coef 取算术平均)
g_{-1/2}=0;  g_{N+1/2} = -R*beta*(phi_N - phi_amb)
水分: phi=C, coef=D, kappa=1,     beta=hm, phi_amb=C_air
温度: phi=T, coef=k, kappa=rho*cp, beta=h,  phi_amb=T_air
时间: theta 法, theta=0.5(CN), 前 4 步 theta=1 (Rannacher)；物性冻结于上一时间步。
"""
import numpy as np
from scipy.linalg import solve_banded
import common as cm


def simulate(N, prob, t_end, dt_plan, shrink=False, theta=0.5, n_ram=4,
             t_report=None, report_idx=None):
    h = 1.0 / N
    xi = np.arange(N + 1) * h
    xif = (np.arange(N) + 0.5) * h
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
        th = 1.0 if step < n_ram else theta

        rho, cp, k, D = cm.props(prob, C, T)

        Lc, cc = _operator(N, h, xi, xif, vol, D, np.ones(N + 1), R,
                           cm.hm, Cair)
        Lt, ct = _operator(N, h, xi, xif, vol, k, rho * cp, R,
                           cm.h, Tair)

        C = _theta_step(Lc, cc, C, dt, th)
        T = _theta_step(Lt, ct, T, dt, th)
        t += dt
        step += 1

        while ri < len(report_idx) and t >= report_idx[ri] - 1e-9:
            out_t.append(t); out_C.append(C.copy()); out_T.append(T.copy())
            ri += 1
    return xi, out_t, out_C, out_T


def _operator(N, h, xi, xif, vol, coef, kappa, R, beta, amb):
    """构造 L（三对角, 以 banded 形式返回）与边界常数向量 c。"""
    a = np.zeros(N + 1)     # 主对角
    up = np.zeros(N + 1)    # 上对角 (i -> i+1)
    lo = np.zeros(N + 1)    # 下对角 (i -> i-1)
    c = np.zeros(N + 1)

    face = 0.5 * (coef[:-1] + coef[1:])          # 面系数（算术平均）
    gf = xif * face / h                          # 内部面系数 xi_{i+1/2}*coef/h

    inv = 1.0 / (R ** 2 * kappa * vol)
    # i = 0
    a[0] = -inv[0] * gf[0]
    up[0] = inv[0] * gf[0]
    # i = 1..N-1
    i = np.arange(1, N)
    a[i] = -inv[i] * (gf[i] + gf[i - 1])
    up[i] = inv[i] * gf[i]
    lo[i] = inv[i] * gf[i - 1]
    # i = N
    a[N] = -inv[N] * (gf[N - 1] + R * beta)
    lo[N] = inv[N] * gf[N - 1]
    c[N] = inv[N] * R * beta * amb
    return (a, up, lo), c


def _theta_step(op, c, phi, dt, th):
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
    """单次前向推进，在 max_r C 首次下穿判据处线性插值求交。"""
    h = 1.0 / N
    xi = np.arange(N + 1) * h
    xif = (np.arange(N) + 0.5) * h
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
        R = float(cm.R_of(t + dt)) if shrink else cm.R0   # 取步末半径
        Cair = float(cm.C_air(t + 0.5 * dt))              # 环境取步中点值
        Tair = float(cm.T_air(t + 0.5 * dt))
        th = 1.0 if step < 4 else 0.5
        rho, cp, k, D = cm.props(prob, C, T)
        Lc, cc = _operator(N, h, xi, xif, vol, D, np.ones(N + 1), R, cm.hm, Cair)
        Lt, ct = _operator(N, h, xi, xif, vol, k, rho * cp, R, cm.h, Tair)
        Cnew = _theta_step(Lc, cc, C, dt, th)
        T = _theta_step(Lt, ct, T, dt, th)
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
