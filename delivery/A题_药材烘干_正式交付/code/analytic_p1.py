# -*- coding: utf-8 -*-
# ==============================================================================
# analytic_p1.py —— 问题 1 温度场的闭式解（Bessel-Duhamel 特征函数展开）
# ==============================================================================
# 做什么：常物性、轴对称一维径向导热问题的解析解，用来独立校核谱配置法（spectral.py）
#         在问题 1 下的温度场并给出空间收敛阶的参照。解由准静态项与齐次瞬态级数两段
#         组成；烘房温度按一阶惯性式 T_air(t) = Tset - A*exp(-t/tauT) 变化。
# 输入：不读文件。常数物性取题目附录 2（直接写在下面）；R0、h、Tset、tauT、T0 取自
#       common.py，其中的环境式由附件 1 的 241 组实测数据标定。
# 输出：不写文件。T_exact(r, t) 供 run_all.py、make_figures.py 做误差对照；
#       直接运行本模块会打印自检值（特征值残差、t=0 截断残差、nterm 敏感性）。
# 关键变量：
#   RHO,CP,K  常数物性 820.0 [kg/m3]、2600.0 [J/(kg K)]、0.36 [W/(m K)]
#   ALPHA     热扩散率 k/(rho*cp) = 1.6886e-7 [m2/s]
#   BI        传热 Biot 数 h*R0/k = 1.3889（无量纲）
#   A, KAPPA  惯性式幅值 Tset - T0 = 22.212 [degC]；准静态项波数 56.163 [1/m]
#   BET,GAM   特征值 beta_n 与前 200 个衰减率 gamma_n = alpha*beta_n^2/R0^2 [1/s]
#   BN        瞬态级数的 200 个系数 b_n（由初值展开定出）[degC]
# 方程与解法（0 < r < R0 的常物性径向导热）：
#   rho*cp*dT/dt = (1/r) d/dr ( r*k*dT/dr )，dT/dr(0,t) = 0，T(r,0) = 28 degC，
#   -k*dT/dr(R0,t) = h*(T(R0,t) - T_air(t))。
#   令 w = T - T_air(t) 使边界齐次，再作准静态分解 w = W(r)exp(-t/tauT) + w_hat：
#     W(r) = A + B*J0(kappa*r)，B = h*A/(k*kappa*J1(kappa*R0) - h*J0(kappa*R0))
#     w_hat = sum_n b_n*J0(beta_n*r/R0)*exp(-gamma_n*t)，beta*J1(beta) = Bi*J0(beta)
# ==============================================================================
import numpy as np
from scipy.special import j0, j1
import common as cm

# ---------------- 常数物性与派生量（题目附录 2） ----------------
# 问题 1 的温度场与水分场解耦，故 rho、cp、k 取常数，不随 C、T 更新。
RHO, CP, K = 820.0, 2600.0, 0.36      # [kg/m3] [J/(kg K)] [W/(m K)]
R = cm.R0                             # 药材半径 R0 [m]
ALPHA = K / (RHO * CP)                # 热扩散率 alpha = k/(rho*cp) [m2/s]
BI = cm.h * R / K                     # 传热 Biot 数 h*R0/k（无量纲）
A = cm.Tset - cm.T0                   # 惯性式幅值 A = Tset - T0 [degC]
# 准静态项的波数 kappa [1/m]：把 w = W(r)exp(-t/tauT) 代回方程后，指数因子上
# 剩下的常数项给出 kappa^2 = 1/(alpha*tauT)，量纲上正是 1/热扩散长度。
KAPPA = np.sqrt(1.0 / (ALPHA * cm.tauT))

# ---------------- 特征函数展开的系数（取 NROOT 项） ----------------
# 特征值 beta_n 由 Robin 条件分离变量得到，即 beta*J1(beta) = Bi*J0(beta) 的正根。
NROOT = 200
BET = cm.bessel_roots(BI, NROOT)
GAM = ALPHA * BET ** 2 / R ** 2       # 模态衰减率 gamma_n = alpha*beta_n^2/R0^2 [1/s]
# _NN 是归一化内积 N_n = 2*J1/(beta*(J0^2+J1^2))（J0、J1 在 beta_n 处取值），
# BN 则是瞬态级数系数 b_n：把初值 w_hat(r,0) = -W(r) 按 J0(beta_n r/R0) 展开得到，
# 其分母 beta_n^2/R0^2 - kappa^2 是本征波数与准静态波数的平方差：kappa 恰等于某个
# beta_n/R0 时该项发散，本题二者相差较远，各 b_n 都是有限值。
_NN = 2 * j1(BET) / (BET * (j0(BET) ** 2 + j1(BET) ** 2))
BN = A * KAPPA ** 2 * _NN / (BET ** 2 / R ** 2 - KAPPA ** 2)


def T_exact(r, t, Tset=None, tauT=None, nterm=NROOT):
    """闭式解 T(r,t) [degC]，r [m]、t [s]；返回与 r 同形状的数组。

    Tset、tauT 不给就取 common.py 的标定值，传别的值可做环境参数敏感性分析；
    nterm 为瞬态级数的截断项数，默认 NROOT = 200（n >= 2 的模态衰减很快）。
    级数形式：T = T_air(t) + (A + B*J0(kappa*r))*exp(-t/tauT)
                        + sum_n b_n*J0(beta_n*r/R0)*exp(-gamma_n*t)
    """
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
    """返回自检报告文本：特征值方程残差、t=0 处的全场极差与 nterm 敏感性表。

    t=0 的极差就是瞬态级数在初值处的截断残差（理论上为 0），nterm 加密后该值
    与 T(0,1800)、T(R0,1800) 是否变化，说明截断项数取得够不够。供 run_all.py 打印。
    """
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
