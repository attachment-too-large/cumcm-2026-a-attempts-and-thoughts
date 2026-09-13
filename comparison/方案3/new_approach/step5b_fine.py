# -*- coding: utf-8 -*-
"""高精度对照：谱解 vs 有限体积(时间 Richardson) 在若干"第4位小数临界"点。"""
import numpy as np
import common as cm
from spectral import HerbSolver
from ref_fv import solve_fv, sample_xi

out = []
pts = [(1200.0, 0.020), (1800.0, 0.015), (1800.0, 0.020), (900.0, 0.010)]

out.append("== 谱配置法（提高 N）==")
for N in [48, 64, 96, 128, 160]:
    s = HerbSolver(N, prob=1)
    at = np.concatenate([np.full(s.Np, 1e-14), np.full(s.Np, 1e-12)])
    tt = sorted(set(p[0] for p in pts))
    sol = s.solve(1800.0, t_eval=tt, rtol=1e-12, atol=at)
    line = "N=%-4d " % N
    for (t, r) in pts:
        j = tt.index(t)
        Cn, Tn = s.interp_u(sol.y[:, j], (r / cm.R0) ** 2)
        line += " t=%.0f r=%.3f C=%.8f |" % (t, r, Cn[0])
    out.append(line)

out.append("")
out.append("== 有限体积 N=800，时间步收敛 + Richardson ==")
res = {}
for dt in [2.0, 1.0, 0.5, 0.25, 0.125]:
    xi, Ctab, Ttab = solve_fv(800, dt, 1800.0, prob=1, t_report=[900, 1200, 1800])
    v = {}
    for (t, r) in pts:
        i = [900, 1200, 1800].index(int(t))
        v[(t, r)] = float(np.atleast_1d(sample_xi(xi, Ctab[i], r / cm.R0))[0])
    res[dt] = v
    out.append("dt=%-6.3f " % dt + " ".join("t=%.0f r=%.3f C=%.8f |" % (t, r, v[(t, r)])
                                            for (t, r) in pts))

out.append("")
out.append("== Richardson 外推 (dt->0, 一阶) ==")
for (t, r) in pts:
    a = res[0.5][(t, r)]
    b = res[0.25][(t, r)]
    c = res[0.125][(t, r)]
    ext = c - (b - c)               # 假定一阶
    out.append("t=%.0f r=%.3f : dt=0.5 %.8f  dt=0.25 %.8f  dt=0.125 %.8f  -> 外推 %.8f (->%.4f)"
               % (t, r, a, b, c, ext, ext))
open(r"C:\Users\qing1\Desktop\数A\new_approach\step5b_fine.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
