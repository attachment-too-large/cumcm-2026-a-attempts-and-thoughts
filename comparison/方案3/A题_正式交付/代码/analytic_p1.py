# -*- coding: utf-8 -*-
"""问题 1 温度场的闭式解（Bessel-Duhamel 特征函数展开）。

模型（常物性、轴对称一维径向）：
    rho*cp*dT/dt = (1/r) d/dr ( r*k*dT/dr ),  0<r<R0
    dT/dr(0,t)=0,   -k*dT/dr(R0,t) = h*(T(R0,t) - T_air(t)),   T(r,0)=28
烘房温度取一阶惯性式  T_air(t) = Tset - A*exp(-t/tauT),  A = Tset-28。

推导（详见报告第 5 章）：令 w = T - T_air(t)，则
    dw/dt = alpha*Lap(w) - (A/tauT)*exp(-t/tauT)
    -k*dw/dr(R0,t) = h*w(R0,t)          （齐次 Robin）
    w(r,0) = 0
作准静态分解 w = W(r)*exp(-t/tauT) + w_hat：
    Lap(W) + kappa^2 W = kappa^2 A,  kappa^2 = 1/(alpha*tauT),  -k W'(R0)=h W(R0)
    => W(r) = A + B*J0(kappa*r),  B = h*A/(k*kappa*J1(kappa*R0) - h*J0(kappa*R0))
    w_hat 满足齐次热方程 + 齐次 Robin + w_hat(r,0) = -W(r)
    => w_hat = sum_n b_n J0(beta_n r/R0) exp(-gamma_n t)
       beta_n: beta*J1(beta) = Bi*J0(beta),  Bi = h*R0/k,  gamma_n = alpha*beta_n^2/R0^2
       b_n = A*kappa^2*N_n/(beta_n^2/R0^2 - kappa^2),
       N_n = 2*J1(beta_n)/(beta_n*(J0(beta_n)^2+J1(beta_n)^2))

最终  T(r,t) = T_air(t) + W(r)exp(-t/tauT) + sum_n b_n J0(beta_n r/R0) exp(-gamma_n t)
"""
import numpy as np
from scipy.special import j0, j1
import common as cm

# 问题 1 的常数物性（附录 2）
RHO, CP, K = 820.0, 2600.0, 0.36
R = cm.R0
ALPHA = K / (RHO * CP)
BI = cm.h * R / K
A = cm.Tset - cm.T0
KAPPA = np.sqrt(1.0 / (ALPHA * cm.tauT))
NROOT = 200
BET = cm.bessel_roots(BI, NROOT)
GAM = ALPHA * BET ** 2 / R ** 2
_B = cm.h * A / (K * KAPPA * j1(KAPPA * R) - cm.h * j0(KAPPA * R))
_NN = 2 * j1(BET) / (BET * (j0(BET) ** 2 + j1(BET) ** 2))
BN = A * KAPPA ** 2 * _NN / (BET ** 2 / R ** 2 - KAPPA ** 2)


def T_exact(r, t, Tset=None, tauT=None, nterm=NROOT):
    """闭式解 T(r,t) [degC]。可传入不同的环境标定值做敏感性分析。"""
    Ts = cm.Tset if Tset is None else Tset
    tau = cm.tauT if tauT is None else tauT
    a = Ts - cm.T0
    kap = np.sqrt(1.0 / (ALPHA * tau))
    B = cm.h * a / (K * kap * j1(kap * R) - cm.h * j0(kap * R))
    bn = a * kap ** 2 * _NN / (BET ** 2 / R ** 2 - kap ** 2)
    r = np.atleast_1d(np.asarray(r, dtype=float))
    val = np.full(r.shape, Ts - (Ts - cm.T0) * np.exp(-t / tau))
    val += (a + B * j0(kap * r)) * np.exp(-t / tau)
    for n in range(nterm):
        val += bn[n] * j0(BET[n] * r / R) * np.exp(-GAM[n] * t)
    return val


def _selftest():
    out = []
    out.append("alpha = %.6e m^2/s ; Bi = %.6f ; kappa = %.6f 1/m ; kappa*R = %.6f"
               % (ALPHA, BI, KAPPA, KAPPA * R))
    out.append("beta_1..4 = %s" % np.array2string(BET[:4], precision=10))
    out.append("特征方程残差 max|beta*J1-Bi*J0| = %.3e"
               % np.max(np.abs(BET * j1(BET) - BI * j0(BET))))
    out.append("b_1..4   = %s" % np.array2string(BN[:4], precision=6))
    out.append("gamma_n*t(1800) n=1..4 = %s" % np.array2string(GAM[:4] * 1800, precision=3))
    rchk = np.linspace(0, R, 41)
    out.append("t=0 全场极差（截断残差，应趋于 0）= %.3e" % np.ptp(T_exact(rchk, 0.0)))
    out.append("  注：t>0 时 n>=2 项已被 exp(-gamma_n t) 压到 e^-13 以下，")
    out.append("      故 t=1800 s 的结果对 nterm 不敏感（见下）。")
    for nt in [1, 2, 3, 4, 8]:
        out.append("nterm=%d  T(0,1800)=%.10f  T(R,1800)=%.10f"
                   % (nt, T_exact(0.0, 1800.0, nterm=nt)[0],
                      T_exact(R, 1800.0, nterm=nt)[0]))
    return "\n".join(out)


if __name__ == "__main__":
    print(_selftest())
