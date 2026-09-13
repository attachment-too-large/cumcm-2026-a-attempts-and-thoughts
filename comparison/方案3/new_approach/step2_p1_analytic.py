# -*- coding: utf-8 -*-
"""问题1 温度场：解析 Bessel-Duhamel 特征函数解（新解法 A）。

模型（与题目一致，常物性、轴对称一维径向）：
    rho*cp*dT/dt = (1/r) d/dr ( r*k*dT/dr ),  0<r<R0
    dT/dr(0,t)=0,   -k*dT/dr(R0,t) = h*(T(R0,t) - T_air(t))
    T(r,0) = 28

令 w = T - T_air(t)，T_air(t)=Tset - A*exp(-t/tauT)，A=Tset-28。则
    dw/dt = alpha*Lap(w) - (A/tauT)*exp(-t/tauT)
    -k*dw/dr(R0,t) = h*w(R0,t)         （齐次 Robin）
    w(r,0) = 0

准静态分解 w = W(r)*exp(-t/tauT) + w_hat：
    Lap(W) + kappa^2 W = kappa^2 A,  kappa^2 = 1/(alpha*tauT)
    -k W'(R0) = h W(R0)
  => W(r) = A + B*J0(kappa*r),  B = h*A/(k*kappa*J1(kappa*R0) - h*J0(kappa*R0))
    w_hat: 齐次热方程 + 齐次 Robin + w_hat(r,0) = -W(r)
  => w_hat = sum_n b_n J0(beta_n r/R0) exp(-gamma_n t)
     beta_n: beta*J1(beta)=Bi*J0(beta),  Bi=h*R0/k,  gamma_n = alpha*beta_n^2/R0^2
     b_n = A*kappa^2*N_n/(beta_n^2/R0^2 - kappa^2),
     N_n = 2*J1(beta_n)/(beta_n*(J0(beta_n)^2+J1(beta_n)^2))

最终 T(r,t) = T_air(t) + W(r)exp(-t/tauT) + sum_n b_n J0(beta_n r/R0) exp(-gamma_n t)
"""
import numpy as np
from scipy.special import j0, j1
import common as cm

out = []
P = 1
rho, cp, k = 820.0, 2600.0, 0.36
alpha = k / (rho * cp)
R = cm.R0
Bi = h = cm.h * R / k
A = cm.Tset - cm.T0
kappa2 = 1.0 / (alpha * cm.tauT)
kappa = np.sqrt(kappa2)
gam0 = 1.0 / cm.tauT

out.append("== 问题1 温度场解析解 参数自检 ==")
out.append("alpha = k/(rho*cp) = %.6e m^2/s   (论文值 1.6886e-7)" % alpha)
out.append("Bi = hR/k = %.6f" % Bi)
out.append("A = Tset-28 = %.6f  tauT = %.6f s" % (A, cm.tauT))
out.append("kappa = %.6f 1/m, kappa*R = %.6f, kappa^2 = %.4f" % (kappa, kappa * R, kappa2))
out.append("R^2/alpha = %.2f s (热特征时间)" % (R ** 2 / alpha))

# ---- 根 ----
NROOT = 40
bet = cm.bessel_roots(Bi, NROOT)
# 自检：特征方程残差
res = bet * j1(bet) - Bi * j0(bet)
out.append("特征方程残差 max|beta*J1-Bi*J0| = %.3e ; beta_1..4 = %s"
           % (np.max(np.abs(res)), np.array2string(bet[:4], precision=10)))

# ---- W(r) ----
den = k * kappa * j1(kappa * R) - cm.h * j0(kappa * R)
B = cm.h * A / den
W0 = A + B                       # r=0
WR = A + B * j0(kappa * R)

# 自检1：W 是否满足 Robin 条件与 ODE
lapW_R = -kappa2 * B * j0(kappa * R)
ode_res = alpha * lapW_R + WR / cm.tauT - A / cm.tauT
out.append("W 的 ODE 残差 alpha*LapW + W/tauT - A/tauT @r=R : %.3e" % ode_res)
out.append("W 的 Robin 残差 -k*W'(R) - h*W(R) : %.3e"
           % (-k * (-kappa * B * j1(kappa * R)) - cm.h * WR))
out.append("B = %.8f   W(0) = %.8f   W(R) = %.8f" % (B, W0, WR))

# ---- b_n ----
gam = alpha * bet ** 2 / R ** 2
Nn = 2 * j1(bet) / (bet * (j0(bet) ** 2 + j1(bet) ** 2))
bn = A * kappa2 * Nn / (bet ** 2 / R ** 2 - kappa2)

# 自检2：sum b_n J0(beta_n r/R) 是否等于 -W(r)（取足够多项，含指数权重=1）
rchk = np.linspace(0, R, 21)
rec = sum(bn[n] * j0(bet[n] * rchk / R) for n in range(NROOT))
Wchk = A + B * j0(kappa * rchk)
out.append("分解自检 max|sum b_n J0 - (-W)| = %.3e  (前%d项)" % (np.max(np.abs(rec + Wchk)), NROOT))
out.append("b_1..4 = %s" % np.array2string(bn[:4], precision=6))
out.append("gamma_n t=1800 for n=1..4 = %s" % np.array2string(gam[:4] * 1800, precision=3))


def T_p1(r, t, nterm=NROOT):
    r = np.atleast_1d(np.asarray(r, dtype=float))
    val = np.full(r.shape, float(cm.T_air(t)))
    val += (A + B * j0(kappa * r)) * np.exp(-t / cm.tauT)
    for n in range(nterm):
        val += bn[n] * j0(bet[n] * r / R) * np.exp(-gam[n] * t)
    return val


# ---- 自检3：初始条件 ----
out.append("t=0 时 T(r) 极差 max-min = %.3e （应为0，即恒为28）"
           % (np.max(T_p1(np.linspace(0, R, 41), 0.0)) - np.min(T_p1(np.linspace(0, R, 41), 0.0))))

# ---- 截断项数收敛性 ----
out.append("")
out.append("== 截断收敛性（t=1800 s）==")
ref = None
for nt in [1, 2, 3, 4, 6, 10, 20, 40]:
    c = T_p1(0.0, 1800.0, nt)[0]
    s = T_p1(R, 1800.0, nt)[0]
    out.append("nterm=%2d  T(0)=%.10f  T(R)=%.10f" % (nt, c, s))
out.append("解析解 nterm=40: T(0,1800)=%.6f  T(R,1800)=%.6f" % (T_p1(0.0, 1800.0, 40)[0], T_p1(R, 1800.0, 40)[0]))
out.append("目标: 中心 34.0425, 表面 37.1989;  偏差 = %+.4f / %+.4f"
           % (T_p1(0.0, 1800.0, 40)[0] - 34.0425, T_p1(R, 1800.0, 40)[0] - 37.1989))

# ---- 题目要求的表1 ----
out.append("")
out.append("== 表1 解析解（°C）==")
rs = [0.0, 0.005, 0.010, 0.015, 0.020]
out.append("时间/s  " + "  ".join("r=%.1fcm" % (x * 100) for x in rs))
for t in [100, 300, 600, 900, 1200, 1500, 1800]:
    out.append("%6d  " % t + "  ".join("%8.4f" % T_p1(x, float(t))[0] for x in rs))

open(r"C:\Users\qing1\Desktop\数A\new_approach\step2_p1_analytic.txt", "w", encoding="utf-8").write("\n".join(out))
print("done")
