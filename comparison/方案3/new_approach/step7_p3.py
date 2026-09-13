# -*- coding: utf-8 -*-
"""问题3：烘干时间 t_dry（处处水分浓度 < 0.15 kg/kg）。

新做法（与 paper 的 Richardson 外推不同）：用自适应 BDF 的稠密输出 +
终端事件（Brent 求交）直接定位 C(0,t)=0.15 的首次穿越时刻，并做网格收敛性
研究确认谱解已收敛到连续极限。
"""
import numpy as np
import common as cm
from spectral import HerbSolver

CRIT = 0.15
out = []
rows_r = [0.0, 0.005, 0.010, 0.015, 0.020]
U = (np.array(rows_r) / cm.R0) ** 2
T_END = 300000.0


def run(N, rtol=1e-10):
    s = HerbSolver(N, prob=3)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])

    def ev(t, y):
        return np.max(y[:s.Np]) - CRIT
    ev.terminal = True
    ev.direction = -1
    sol = s.solve(T_END, t_eval=None, rtol=rtol, atol=at, events=ev)
    return s, sol


out.append("== 问题3 烘干时间：网格收敛性 ==")
for N in [24, 32, 48, 64, 96, 128]:
    s, sol = run(N)
    td = sol.t_events[0]
    td = float(td[0]) if len(td) else np.nan
    yend = sol.y[:, -1]
    out.append("N=%-4d tdry=%12.3f s = %8.4f h   Cmax_end=%.8f  nfev=%d"
               % (N, td, td / 3600.0, np.max(yend[:s.Np]), sol.nfev))
out.append("目标 tdry = 205559 s = 57.0997 h")

# 生产：N=64
N = 64
s, sol = run(N, rtol=1e-11)
tdry = float(sol.t_events[0][0])
out.append("")
out.append("生产 N=%d: tdry = %.3f s = %.4f h = %.4f d  (nfev=%d)"
           % (N, tdry, tdry / 3600.0, tdry / 86400.0, sol.nfev))

# 判据自检：终点处全场最大值是否在中心
yend = sol.y[:, -1]
Cend = yend[:s.Np]
out.append("终点 Cmax=%.8f 位于 u=%.4g (0=中心,1=表面); C(0)=%.8f C(R)=%.8f"
           % (np.max(Cend), s.u[int(np.argmax(Cend))], Cend[0], Cend[-1]))

# 表5：每 6 h
hours = [6, 12, 18, 24, 30, 36, 42, 48, 54]
tt = np.array([h * 3600.0 for h in hours])
at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])
sol2 = s.solve(tdry, t_eval=list(tt) + [tdry], rtol=1e-11, atol=at)
tab = np.zeros((len(hours) + 1, len(rows_r)))
for i in range(len(hours) + 1):
    Cn, _ = s.interp_u(sol2.y[:, i], U)
    tab[i] = Cn

out.append("")
out.append("== 表5 水分浓度（新解法） ==")
out.append("t/h    " + "".join("%10s" % ("r=%.1f" % (x * 100)) for x in rows_r))
for i, h in enumerate(hours):
    out.append("%-6d " % h + "".join("%10.4f" % v for v in tab[i]))
out.append("结束   " + "".join("%10.4f" % v for v in tab[-1]))
out.append("== 论文表5 ==")
ref = [[1.0149, 0.9863, 0.9003, 0.7548, 0.5343],
       [0.4546, 0.4417, 0.4019, 0.3291, 0.1651],
       [0.2978, 0.2905, 0.2678, 0.2247, 0.0853],
       [0.2370, 0.2320, 0.2160, 0.1850, 0.0673],
       [0.2052, 0.2012, 0.1886, 0.1636, 0.0607],
       [0.1853, 0.1820, 0.1713, 0.1500, 0.0575],
       [0.1716, 0.1686, 0.1593, 0.1404, 0.0557],
       [0.1614, 0.1588, 0.1503, 0.1331, 0.0546],
       [0.1535, 0.1510, 0.1433, 0.1274, 0.0539],
       [0.1500, 0.1477, 0.1402, 0.1249, 0.0536]]
for i, h in enumerate(hours):
    out.append("%-6d " % h + "".join("%10.4f" % v for v in ref[i]))
out.append("结束   " + "".join("%10.4f" % v for v in ref[-1]))

out.append("")
out.append("tdry = %.3f s -> %.0f s (目标 205559)  偏差 %+.1f s (%.4f%%)"
           % (tdry, round(tdry), tdry - 205559.0, (tdry - 205559.0) / 205559.0 * 100))
np.savez(r"C:\Users\qing1\Desktop\数A\new_approach\p3_result.npz",
         tdry=tdry, tab=tab, hours=np.array(hours), sol_t=sol2.t, sol_y=sol2.y)
open(r"C:\Users\qing1\Desktop\数A\new_approach\step7_p3.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
