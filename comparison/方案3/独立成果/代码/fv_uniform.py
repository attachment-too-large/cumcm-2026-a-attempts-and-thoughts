# -*- coding: utf-8 -*-
"""独立参照解：均匀 r 网格 + 守恒型有限体积 + 后向 Euler + Picard 迭代。

刻意与谱配置法在三个方面都不同：
  * 坐标：均匀 r 网格（ξ=i/N），不用 u=ξ^2；
  * 离散：守恒型有限体积（面通量、体积测度），非谱微分；
  * 时间：等步长后向 Euler（一阶），非自适应 BDF。
因此两法一致 => 解的是同一个 PDE，且各自离散误差均已被压到目标精度以下。
"""
import numpy as np
import common as cm


def solve_fv(N, dt, t_end, prob, t_report, shrink=False, n_picard=3):
    """返回 (t_report 数组, C_tab[Nr?], T_tab)。C,T 形状 (len(t_report), N+1)。"""
    h = 1.0 / N
    xi = np.arange(N + 1) * h
    xif = (np.arange(N) + 0.5) * h           # 内部面 ξ_{i+1/2}, i=0..N-1

    vol = np.empty(N + 1)
    vol[0] = h * h / 8.0
    vol[1:N] = xi[1:N] * h
    vol[N] = h / 2.0 - h * h / 8.0

    C = np.full(N + 1, cm.C0)
    T = np.full(N + 1, cm.T0)

    nstep = int(round(t_end / dt))
    trep = np.asarray(t_report, dtype=float)
    idx = {int(round(t / dt)): k for k, t in enumerate(trep)}
    Ctab = np.full((len(trep), N + 1), np.nan)
    Ttab = np.full((len(trep), N + 1), np.nan)

    def R_of(t):
        return float(cm.R_of(t)) if shrink else cm.R0

    for n in range(1, nstep + 1):
        t_new = n * dt
        t_old = (n - 1) * dt
        R = R_of(t_new)
        Cair = float(cm.C_air(t_new))
        Tair = float(cm.T_air(t_new))

        Cnew, Tnew = C.copy(), T.copy()
        for _ in range(n_picard):
            rho, cp, k, D = cm.props(prob, Cnew, Tnew)
            # ---- 水分三对角 ----
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
        if n in idx:
            Ctab[idx[n]] = C
            Ttab[idx[n]] = T
    return xi, Ctab, Ttab


def _thomas(A, b):
    """三对角求解：A 为 (3, n)，A[0]=上对角(偏移+1), A[1]=主对角, A[2]=下对角(偏移-1)。"""
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
    """线性插值到 ξ 目标位置。"""
    return np.interp(np.asarray(x_target, float), xi, vals)
