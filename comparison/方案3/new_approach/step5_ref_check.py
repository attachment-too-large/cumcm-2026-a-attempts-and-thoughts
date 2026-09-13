# -*- coding: utf-8 -*-
"""问题1：有限体积参照解 vs 解析解(T) / 谱解(C)。"""
import numpy as np
import common as cm
from ref_fv import solve_fv, sample_xi

TARGET_C = 1.5109
TARGET_T0 = 34.0425
TARGET_TR = 37.1989
t_report = [100, 300, 600, 900, 1200, 1500, 1800]
rows_r = [0.0, 0.005, 0.010, 0.015, 0.020]

out = []
out.append("== 问题1 有限体积参照解收敛性 ==")
out.append("N     dt      T(0,1800)    T(R,1800)    C(0,1800)    C(R,1800)   C(1.5cm)   C(R,1200)")
out.append("目标值         %.4f      %.4f      (2.5500)     %.4f      2.3755      1.6592"
           % (TARGET_T0, TARGET_TR, TARGET_C))
for N, dt in [(200, 2.0), (200, 1.0), (200, 0.5), (400, 2.0), (400, 1.0),
              (400, 0.5), (800, 1.0), (800, 0.5), (800, 0.25)]:
    xi, Ctab, Ttab = solve_fv(N, dt, 1800.0, prob=1, t_report=t_report)
    xj = np.array(rows_r) / cm.R0
    vC0 = sample_xi(xi, Ctab[-1], xj)
    vT0 = sample_xi(xi, Ttab[-1], xj)
    C_R_1200 = sample_xi(xi, Ctab[4], xj)[-1]
    out.append("N=%-4d dt=%-5.2f %11.6f %12.6f %12.6f %13.6f %9.4f %10.4f"
               % (N, dt, vT0[0], vT0[-1], vC0[0], vC0[-1], vC0[3], C_R_1200))

# 精算一组并输出完整表2 供对照
xi, Ctab, Ttab = solve_fv(800, 0.25, 1800.0, prob=1, t_report=t_report)
xj = np.array(rows_r) / cm.R0
out.append("")
out.append("== 表2 有限体积 N=800 dt=0.25 ==")
out.append("t/s    " + "".join("%10s" % ("r=%.1f" % (x * 100)) for x in rows_r))
for i, t in enumerate(t_report):
    v = sample_xi(xi, Ctab[i], xj)
    out.append("%-6d " % t + "".join("%10.4f" % z for z in v))
out.append("== 表2 谱配置法 N=64 对照 ==")
ref = [[2.5500, 2.5500, 2.5500, 2.5500, 2.2471],
       [2.5500, 2.5500, 2.5500, 2.5492, 2.0510],
       [2.5500, 2.5500, 2.5500, 2.5353, 1.8773],
       [2.5500, 2.5500, 2.5497, 2.5046, 1.7552],
       [2.5500, 2.5500, 2.5482, 2.4646, 1.6591],
       [2.5500, 2.5499, 2.5445, 2.4207, 1.5794],
       [2.5500, 2.5497, 2.5383, 2.3756, 1.5109]]
for i, t in enumerate(t_report):
    out.append("%-6d " % t + "".join("%10.4f" % z for z in ref[i]))

open(r"C:\Users\qing1\Desktop\数A\new_approach\step5_ref_check.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
