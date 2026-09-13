# -*- coding: utf-8 -*-
"""问题1：完整求解（谱配置法），输出表1/表2 与关键值。"""
import numpy as np
import common as cm
from spectral import HerbSolver

out = []
rows_r = [0.0, 0.005, 0.010, 0.015, 0.020]
t_tab = [100, 300, 600, 900, 1200, 1500, 1800]
t_all = np.arange(1, 1801, 1, dtype=float)

out.append("== 问题1 谱配置法：网格收敛性（t=1800 s）==")
res = {}
for N in [24, 32, 48, 64, 96, 128]:
    s = HerbSolver(N, prob=1)
    sol = s.solve(1800.0, t_eval=[1800.0], rtol=1e-11, atol=1e-14)
    y = sol.y[:, -1]
    Cn, Tn = s.interp_u(y, (np.array(rows_r) / cm.R0) ** 2)
    out.append("N=%-4d C(0)=%.8f C(R)=%.8f | T(0)=%.8f T(R)=%.8f  nfev=%d"
               % (N, Cn[0], Cn[-1], Tn[0], Tn[-1], sol.nfev))
    res[N] = (Cn, Tn)
out.append("目标: C(R,1800)=1.5109 ; T(0,1800)=34.0425 ; T(R,1800)=37.1989")

# ---- 生产设置 ----
N = 64
s = HerbSolver(N, prob=1)
sol = s.solve(1800.0, t_eval=t_all, rtol=1e-11, atol=1e-14)
out.append("")
out.append("生产设置 N=%d, rtol=1e-11, nfev=%d, 时间点数=%d" % (N, sol.nfev, sol.y.shape[1]))

U = (np.array(rows_r) / cm.R0) ** 2
tabT = np.zeros((len(t_tab), len(rows_r)))
tabC = np.zeros((len(t_tab), len(rows_r)))
for i, t in enumerate(t_tab):
    j = int(np.argmin(np.abs(sol.t - t)))
    Cn, Tn = s.interp_u(sol.y[:, j], U)
    tabT[i] = Tn
    tabC[i] = Cn

out.append("")
out.append("== 表1 药材温度 / degC ==")
out.append("t/s    " + "".join("%10s" % ("r=%.1f" % (x * 100)) for x in rows_r))
for i, t in enumerate(t_tab):
    out.append("%-6d " % t + "".join("%10.4f" % v for v in tabT[i]))
out.append("== 论文表1 对照 ==")
ref_T = [[28.0001, 28.0005, 28.0054, 28.0425, 28.2219],
         [28.0500, 28.0776, 28.1857, 28.4550, 29.0261],
         [28.5456, 28.6438, 28.9618, 29.5667, 30.5576],
         [29.5568, 29.7074, 30.1740, 30.9976, 32.2359],
         [30.9069, 31.0882, 31.6397, 32.5824, 33.9433],
         [32.4397, 32.6353, 33.2248, 34.2150, 35.6117],
         [34.0425, 34.2411, 34.8362, 35.8249, 37.1989]]
for i, t in enumerate(t_tab):
    out.append("%-6d " % t + "".join("%10.4f" % v for v in ref_T[i]))

out.append("")
out.append("== 表2 水分浓度 / (kg/kg) ==")
out.append("t/s    " + "".join("%10s" % ("r=%.1f" % (x * 100)) for x in rows_r))
for i, t in enumerate(t_tab):
    out.append("%-6d " % t + "".join("%10.4f" % v for v in tabC[i]))
out.append("== 论文表2 对照 ==")
ref_C = [[2.5500, 2.5500, 2.5500, 2.5500, 2.2471],
         [2.5500, 2.5500, 2.5500, 2.5492, 2.0510],
         [2.5500, 2.5500, 2.5500, 2.5353, 1.8773],
         [2.5500, 2.5500, 2.5497, 2.5046, 1.7552],
         [2.5500, 2.5500, 2.5482, 2.4646, 1.6592],
         [2.5500, 2.5499, 2.5445, 2.4207, 1.5794],
         [2.5500, 2.5497, 2.5383, 2.3755, 1.5109]]
for i, t in enumerate(t_tab):
    out.append("%-6d " % t + "".join("%10.4f" % v for v in ref_C[i]))

out.append("")
out.append("== 关键值对照 ==")
out.append("表面水分 C(R,1800) = %.6f -> %.4f   (目标 1.5109)" % (tabC[-1][-1], tabC[-1][-1]))
out.append("中心温度 T(0,1800) = %.6f -> %.4f   (目标 34.0425)" % (tabT[-1][0], tabT[-1][0]))
out.append("表面温度 T(R,1800) = %.6f -> %.4f   (目标 37.1989)" % (tabT[-1][-1], tabT[-1][-1]))

np.save(r"C:\Users\qing1\Desktop\数A\new_approach\p1_state.npy",
        np.vstack([sol.t, sol.y]))
open(r"C:\Users\qing1\Desktop\数A\new_approach\step4_p1_full.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
