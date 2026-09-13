# -*- coding: utf-8 -*-
"""重算七张教学插图中实际画出的数据，供图解析文档引用。

只读，不写任何数据文件，结果打印到屏幕。
"""
import sys
import numpy as np

CODEDIR = r"C:\Users\qing1\Desktop\bestway\药材烘干问题_正则化谱方法\代码"
sys.path.insert(0, CODEDIR)
import common as cm
from spectral import HerbSolver

RCOLS = np.array([0.0, 0.005, 0.010, 0.015, 0.020])   # 米
U = (RCOLS / cm.R0) ** 2

print("=" * 70)
print("图 5（问题 1 剖面）实际画出的数据")
print("=" * 70)
s1 = HerbSolver(64, prob=1)
t1 = np.arange(0.0, 1800.1, 60.0)
at1 = np.concatenate([np.full(s1.Np, 1e-13), np.full(s1.Np, 1e-11)])
sol1 = s1.solve(1800.0, t_eval=list(t1), rtol=1e-11, atol=at1)
print("时刻点数：%d，从 %g 到 %g 秒，间隔 60 秒" % (len(t1), t1[0], t1[-1]))
print("")
print("%-8s %14s %14s %14s %14s" % ("t/s", "T(0)", "T(R)", "C(0)", "C(R)"))
for tt in [100, 600, 1200, 1800]:
    i = int(np.argmin(np.abs(t1 - tt)))
    C, T = s1.interp_u(sol1.y[:, i], U)
    print("%-8d %14.6f %14.6f %14.6f %14.6f" % (tt, T[0], T[-1], C[0], C[-1]))

print("")
print("图中相邻两条剖面曲线的温差（用于判断画出的层次）：")
for tt in [100, 600, 1200, 1800]:
    i = int(np.argmin(np.abs(t1 - tt)))
    C, T = s1.interp_u(sol1.y[:, i], U)
    print("  t=%4d s  中心到表面的温差 %.4f 摄氏度，表面水分降幅 %.4f" %
          (tt, T[-1] - T[0], 2.55 - C[-1]))

print("")
print("=" * 70)
print("图 6（问题 3 干燥曲线）实际画出的数据")
print("=" * 70)
s3 = HerbSolver(48, prob=3)
at3 = np.concatenate([np.full(s3.Np, 1e-13), np.full(s3.Np, 1e-11)])


def ev3(t, y):
    return np.max(y[:s3.Np]) - 0.15


ev3.terminal = True
ev3.direction = -1
solE = s3.solve(300000.0, rtol=1e-10, atol=at3, events=ev3)
tdry = float(solE.t_events[0][0])
print("烘干时长 tdry = %.3f 秒 = %.4f 小时" % (tdry, tdry / 3600.0))
print("终点时全场最大值 %.8f，位于半径 %.5f 米" %
      (np.max(solE.y[:s3.Np, -1]), s3.r_nodes(tdry)[int(np.argmax(solE.y[:s3.Np, -1]))]))

tv = np.linspace(0.0, tdry * 1.06, 1500)
sol3 = s3.solve(tdry * 1.06, t_eval=list(tv), rtol=1e-10, atol=at3)
C3 = np.array([s3.interp_u(sol3.y[:, i], U)[0] for i in range(len(tv))])
hh = tv / 3600.0

print("")
print("横轴范围：0 到 %.3f 小时，共 %d 个点" % (hh[-1], len(hh)))
print("")
print("%-8s %12s %12s %12s %12s %12s" % ("t/h", "C(0)", "C(0.5)", "C(1.0)", "C(1.5)", "C(2.0)"))
for tt in [6, 12, 24, 36, 48, 57]:
    i = int(np.argmin(np.abs(hh - tt)))
    print("%-8d %12.5f %12.5f %12.5f %12.5f %12.5f" %
          (tt, C3[i, 0], C3[i, 1], C3[i, 2], C3[i, 3], C3[i, 4]))

i_end = int(np.argmin(np.abs(hh - tdry / 3600.0)))
print("%-8s %12.5f %12.5f %12.5f %12.5f %12.5f" %
      ("tdry", C3[i_end, 0], C3[i_end, 1], C3[i_end, 2], C3[i_end, 3], C3[i_end, 4]))

# 面积加权截面平均
r_nodes = s3.r_nodes(0.0)
w = r_nodes.copy()
mean = np.array([np.sum(C3[i] * w) / np.sum(w) for i in range(len(tv))])
idx = np.where(mean < 0.15)[0]
print("")
if idx.size:
    i_cross = idx[0]
    print("截面面积加权平均首次低于 0.15 的时刻：%.3f 小时" % hh[i_cross])
    print("该时刻中心值 %.5f，表面值 %.5f" % (C3[i_cross, 0], C3[i_cross, 4]))
    print("中心达标比平均达标晚 %.3f 小时，即晚 %.1f%%" %
          (tdry / 3600.0 - hh[i_cross],
           (tdry / 3600.0 - hh[i_cross]) / (tdry / 3600.0) * 100))

print("")
print("=" * 70)
print("图 7(a)（谱方法收敛）实际画出的数据")
print("=" * 70)
import analytic_p1 as an
print("%-6s %16s" % ("N", "全场最大误差/K"))
for N in [4, 6, 8, 12, 16, 24, 32, 64]:
    s = HerbSolver(N, prob=1)
    at = np.concatenate([np.full(s.Np, 1e-14), np.full(s.Np, 1e-12)])
    sol = s.solve(1800.0, t_eval=[1800.0], rtol=1e-12, atol=at)
    rc = np.linspace(0, cm.R0, 201)
    Ti = s._bary(sol.y[s.Np:, -1], (rc / cm.R0) ** 2)
    Te = an.T_exact(rc, 1800.0)
    print("%-6d %16.3e" % (N, np.max(np.abs(Ti - Te))))

print("")
print("=" * 70)
print("图 7(b)（有限体积收敛）画出的数据，来自专家报告第 9.3 节")
print("=" * 70)
for N, v in [(800, 205550.98), (1600, 205567.49), (3200, 205572.43), (6400, 205573.77)]:
    print("  N=%-5d tdry=%10.2f 秒，与 205574.2 秒之差 %6.2f 秒" % (N, v, 205574.2 - v))
