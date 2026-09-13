# -*- coding: utf-8 -*-
"""问题4：考虑尺寸收缩（附件2 的 R(t) + 附录4 物性）。

新解法在物质坐标 xi=r/R(t) 上工作，自变量取 u=xi^2，Laplace 算子
    (1/r) d/dr ( r D d/dr ) = (4/R(t)^2) d/du ( u D d/du )
（无对流项，u=0 正则）。输出时按当前物理距离 r_j 反解 xi_j=r_j/R(t_i) 采样，
r_j > R(t_i) 处该位置不存在（留空）。
"""
import numpy as np
import common as cm
from spectral import HerbSolver

CRIT = 0.15
out = []

# ---- R(t) 自检 ----
out.append("== 附件2 半径数据自检 ==")
out.append("R(0)=%.4f cm  R(259200)=%.4f cm  R/R0=%.4f"
           % (cm.R_of(0) * 100, cm.R_of(259200) * 100, cm.R_of(259200) / cm.R0))
out.append("R(182764.7 s)=%.4f cm" % (cm.R_of(182764.7) * 100))
seg = cm.R_of(np.array([0.0, 900.0, 1800.0, 182700.0, 183600.0]) * 1.0)
out.append("R(t) 线性插值 抽样: t=0,900,1800,182700,183600 -> %s cm"
           % np.array2string(seg * 100, precision=4))


def run(N, rtol=1e-10):
    s = HerbSolver(N, prob=4, shrink=True)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])

    def ev(t, y):
        return np.max(y[:s.Np]) - CRIT
    ev.terminal = True
    ev.direction = -1
    sol = s.solve(300000.0, rtol=rtol, atol=at, events=ev)
    return s, sol


out.append("")
out.append("== 问题4 烘干时间：网格收敛性 ==")
for N in [24, 32, 48, 64, 96, 128]:
    s, sol = run(N)
    td = float(sol.t_events[0][0])
    out.append("N=%-4d tdry = %11.2f s = %8.4f h = %.4f d  nfev=%d"
               % (N, td, td / 3600, td / 86400, sol.nfev))
out.append("目标 tdry = 182764.7 s = 50.7680 h")

# ---- 生产解 ----
N = 64
s, sol = run(N, rtol=1e-11)
tdry = float(sol.t_events[0][0])
out.append("")
out.append("生产 N=%d: tdry = %.2f s = %.4f h = %.4f d (nfev=%d)"
           % (N, tdry, tdry / 3600, tdry / 86400, sol.nfev))
yend = sol.y[:, -1]
out.append("终点 Cmax=%.8f 位于 u=%.4g ; C(0)=%.8f C(表面)=%.8f"
           % (np.max(yend[:s.Np]), s.u[int(np.argmax(yend[:s.Np]))],
              yend[0], yend[s.N]))
out.append("终点 R(tdry)=%.4f cm" % (cm.R_of(tdry) * 100))

# ---- 表6 ----
hours = [6, 12, 18, 24, 30, 36, 42, 48]
tt = np.array([h * 3600.0 for h in hours] + [tdry])
at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])
sol2 = s.solve(tdry, t_eval=list(tt), rtol=1e-11, atol=at)
rcols = [0.0, 0.005, 0.010, 0.015]      # 0,0.5,1.0,1.5 cm ; 末列为当前表面
tab = np.full((len(tt), len(rcols) + 1), np.nan)
Rcol = np.zeros(len(tt))
for i, t in enumerate(tt):
    Rt = s.R_of(t)
    Rcol[i] = Rt
    yi = sol2.y[:, i]
    for j, r in enumerate(rcols):
        if r <= Rt + 1e-12:
            Cn, _ = s.interp_u(yi, (r / Rt) ** 2)
            tab[i, j] = Cn[0]
    tab[i, -1] = yi[s.N]                # 当前表面 C(xi=1)

out.append("")
out.append("== 表6 水分浓度（新解法）/ (kg/kg) ==")
out.append("t/h    R/cm   r=0.0    r=0.5    r=1.0    r=1.5    表面")
for i, h in enumerate(hours):
    row = "%-6d %.3f " % (h, Rcol[i] * 100)
    for j in range(len(rcols)):
        row += "%9s" % ("—" if np.isnan(tab[i, j]) else "%.4f" % tab[i, j])
    row += " %8.4f" % tab[i, -1]
    out.append(row)
row = "结束   %.3f " % (Rcol[-1] * 100)
for j in range(len(rcols)):
    row += "%9s" % ("—" if np.isnan(tab[-1, j]) else "%.4f" % tab[-1, j])
row += " %8.4f" % tab[-1, -1]
out.append(row)
out.append("== 论文表6 ==")
ref = [["6", 1.374, 1.7157, 1.5345, 1.0215, None, 0.4228],
       ["12", 1.248, 0.7335, 0.6512, 0.4069, None, 0.1677],
       ["18", 1.214, 0.4062, 0.3665, 0.2388, None, 0.0899],
       ["24", 1.204, 0.2835, 0.2596, 0.1783, None, 0.0681],
       ["30", 1.201, 0.2251, 0.2083, 0.1488, None, 0.0603],
       ["36", 1.200, 0.1919, 0.1788, 0.1313, None, 0.0568],
       ["42", 1.200, 0.1705, 0.1597, 0.1196, None, 0.0550],
       ["48", 1.200, 0.1555, 0.1462, 0.1113, None, 0.0539],
       ["结束50.7678h", 1.200, 0.1500, 0.1413, 0.1081, None, 0.0535]]
for r_ in ref:
    row = "%-6s %.3f " % (r_[0], r_[1])
    for j in range(4):
        row += "%9s" % ("—" if r_[2 + j] is None else "%.4f" % r_[2 + j])
    row += " %8.4f" % r_[6]
    out.append(row)

out.append("")
out.append("tdry = %.2f s -> %.1f s (目标 182764.7)  偏差 %+.1f s (%.4f%%)"
           % (tdry, tdry, tdry - 182764.7, (tdry - 182764.7) / 182764.7 * 100))
np.savez(r"C:\Users\qing1\Desktop\数A\new_approach\p4_result.npz",
         tdry=tdry, tab=tab, hours=np.array(hours), Rcol=Rcol,
         sol_t=sol2.t, sol_y=sol2.y)
open(r"C:\Users\qing1\Desktop\数A\new_approach\step9_p4.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
