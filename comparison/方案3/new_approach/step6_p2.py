# -*- coding: utf-8 -*-
"""问题2：整段烘干过程（预热平衡 + 恒温干燥），谱配置法。"""
import numpy as np
import common as cm
from spectral import HerbSolver

out = []
rows_r = [0.0, 0.005, 0.010, 0.015, 0.020]
t_tab = [1800.0, 3600.0, 5400.0, 7200.0, 9000.0, 10800.0]
U = (np.array(rows_r) / cm.R0) ** 2

out.append("== 问题2 谱配置法网格收敛性（t=10800 s = 3 h）==")
for N in [24, 32, 48, 64, 96]:
    s = HerbSolver(N, prob=2)
    at = np.concatenate([np.full(s.Np, 1e-12), np.full(s.Np, 1e-10)])
    sol = s.solve(10800.0, t_eval=[10800.0], rtol=1e-11, atol=at)
    y = sol.y[:, -1]
    Cn, Tn = s.interp_u(y, U)
    out.append("N=%-4d T(0)=%.8f C(0)=%.8f C(R)=%.8f | T(R)=%.8f  nfev=%d"
               % (N, Tn[0], Cn[0], Cn[-1], Tn[-1], sol.nfev))
out.append("目标:   T(0)=49.9051  C(0)=1.7673  C(R)=1.0083")

# 生产设置
N = 64
s = HerbSolver(N, prob=2)
at = np.concatenate([np.full(s.Np, 1e-12), np.full(s.Np, 1e-10)])
sol = s.solve(10800.0, t_eval=t_tab, rtol=1e-11, atol=at)
out.append("")
out.append("生产 N=%d rtol=1e-11 nfev=%d" % (N, sol.nfev))

tabT = np.zeros((len(t_tab), len(rows_r)))
tabC = np.zeros((len(t_tab), len(rows_r)))
for i in range(len(t_tab)):
    Cn, Tn = s.interp_u(sol.y[:, i], U)
    tabT[i] = Tn
    tabC[i] = Cn

out.append("")
out.append("== 表3 温度 / degC（新解法） ==")
out.append("t/h    " + "".join("%10s" % ("r=%.1f" % (x * 100)) for x in rows_r))
for i, t in enumerate(t_tab):
    out.append("%-6.1f " % (t / 3600) + "".join("%10.4f" % v for v in tabT[i]))
out.append("== 论文表3 ==")
ref_T = [[32.5656, 32.7627, 33.3570, 34.3573, 35.7831],
         [40.4220, 40.5818, 41.0548, 41.8230, 42.8629],
         [45.5851, 45.6719, 45.9272, 46.3366, 46.8815],
         [48.2253, 48.2650, 48.3814, 48.5675, 48.8138],
         [49.4129, 49.4293, 49.4775, 49.5546, 49.6565],
         [49.9051, 49.9115, 49.9302, 49.9602, 49.9999]]
for i, t in enumerate(t_tab):
    out.append("%-6.1f " % (t / 3600) + "".join("%10.4f" % v for v in ref_T[i]))

out.append("")
out.append("== 表4 水分浓度 / (kg/kg)（新解法） ==")
out.append("t/h    " + "".join("%10s" % ("r=%.1f" % (x * 100)) for x in rows_r))
for i, t in enumerate(t_tab):
    out.append("%-6.1f " % (t / 3600) + "".join("%10.4f" % v for v in tabC[i]))
out.append("== 论文表4 ==")
ref_C = [[2.5499, 2.5488, 2.5247, 2.3237, 1.6538],
         [2.5248, 2.4935, 2.3562, 2.0227, 1.4709],
         [2.3860, 2.3258, 2.1350, 1.8020, 1.3444],
         [2.1731, 2.1105, 1.9250, 1.6256, 1.2282],
         [1.9589, 1.9025, 1.7369, 1.4716, 1.1157],
         [1.7673, 1.7174, 1.5705, 1.3332, 1.0083]]
for i, t in enumerate(t_tab):
    out.append("%-6.1f " % (t / 3600) + "".join("%10.4f" % v for v in ref_C[i]))

out.append("")
out.append("== 关键值 ==")
out.append("T(0,3h) = %.6f -> %.4f   (目标 49.9051)" % (tabT[-1][0], tabT[-1][0]))
out.append("C(0,3h) = %.6f -> %.4f   (目标 1.7673)" % (tabC[-1][0], tabC[-1][0]))
out.append("C(R,3h) = %.6f -> %.4f   (目标 1.0083)" % (tabC[-1][-1], tabC[-1][-1]))

np.savez(r"C:\Users\qing1\Desktop\数A\new_approach\p2_state.npz",
         t=sol.t, y=sol.y, U=U, rows_r=np.array(rows_r), tabT=tabT, tabC=tabC)
open(r"C:\Users\qing1\Desktop\数A\new_approach\step6_p2.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
