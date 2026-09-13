# -*- coding: utf-8 -*-
"""自检：谱配置法的空间收敛性（问题1 温度场，与解析解逐点比较）。

若误差随 N 呈指数下降 => 谱精度成立；若停在某个平台 => 格式只有代数阶，需改格式。
"""
import numpy as np
from scipy.special import j0, j1
import common as cm
from spectral import HerbSolver

# ---------- 解析解（问题1 温度，仅 T 场，与 C 解耦） ----------
rho, cp, k = 820.0, 2600.0, 0.36
alpha = k / (rho * cp)
R = cm.R0
Bi = cm.h * R / k
A = cm.Tset - cm.T0
kappa2 = 1.0 / (alpha * cm.tauT)
kappa = np.sqrt(kappa2)
BET = cm.bessel_roots(Bi, 8)
B = cm.h * A / (k * kappa * j1(kappa * R) - cm.h * j0(kappa * R))
GAM = alpha * BET ** 2 / R ** 2
NN = 2 * j1(BET) / (BET * (j0(BET) ** 2 + j1(BET) ** 2))
BN = A * kappa2 * NN / (BET ** 2 / R ** 2 - kappa2)


def T_exact(r, t):
    r = np.atleast_1d(np.asarray(r, float))
    val = np.full(r.shape, float(cm.T_air(t)))
    val += (A + B * j0(kappa * r)) * np.exp(-t / cm.tauT)
    for n in range(8):
        val += BN[n] * j0(BET[n] * r / R) * np.exp(-GAM[n] * t)
    return val


out = []
out.append("== 谱配置法 vs 解析解（问题1 温度场，t=1800 s）==")
out.append("N      T(0)误差      T(R)误差      全场max误差    条件数估计")
rchk = np.linspace(0, R, 201)
for N in [8, 12, 16, 20, 24, 32, 48, 64, 96, 128]:
    s = HerbSolver(N, prob=1)
    sol = s.solve(1800.0, t_eval=[1800.0], rtol=1e-12, atol=1e-14)
    y = sol.y[:, -1]
    # 只用 T 场：把 C 场噪声排除，重算纯 T 方程更干净 —— 这里 C 初值恒定，D 不变，
    # 故耦合项不影响 T（附录2 下 T 与 C 完全解耦）。
    Tn = y[s.Np:]
    e0 = Tn[0] - T_exact(0.0, 1800.0)[0]
    eR = Tn[-1] - T_exact(R, 1800.0)[0]
    Ti = s._bary(Tn, (rchk / R) ** 2)
    emax = np.max(np.abs(Ti - T_exact(rchk, 1800.0)))
    out.append("N=%-4d %+11.3e %+11.3e %11.3e" % (N, e0, eR, emax))

# 时间步长收敛性（固定 N=64）
out.append("")
out.append("== 时间积分精度（N=64）==")
s = HerbSolver(64, prob=1)
for rt, at in [(1e-8, 1e-11), (1e-10, 1e-13), (1e-12, 1e-15)]:
    sol = s.solve(1800.0, t_eval=[1800.0], rtol=rt, atol=at)
    y = sol.y[:, -1]
    out.append("rtol=%.0e atol=%.0e  T(0)-exact=%+.3e  T(R)-exact=%+.3e  nfev=%d"
               % (rt, at, y[s.Np] - T_exact(0.0, 1800.0)[0],
                  y[-1] - T_exact(R, 1800.0)[0], sol.nfev))

open(r"C:\Users\qing1\Desktop\数A\new_approach\step3_spectral_conv.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
