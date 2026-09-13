# -*- coding: utf-8 -*-
"""取图 6 需要的读数（节点参数 64，与画图脚本一致）。"""
import sys
import numpy as np

CODEDIR = r"C:\Users\qing1\Desktop\bestway\药材烘干问题_正则化谱方法\代码"
sys.path.insert(0, CODEDIR)
import common as cm
from spectral import HerbSolver

RCOLS = np.array([0.0, 0.005, 0.010, 0.015, 0.020])
U = (RCOLS / cm.R0) ** 2

s = HerbSolver(64, prob=3)
at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])


def ev(t, y):
    return np.max(y[:s.Np]) - 0.15


ev.terminal = True
ev.direction = -1
solE = s.solve(300000.0, rtol=1e-10, atol=at, events=ev)
tdry = float(solE.t_events[0][0])
tv = np.linspace(0.0, tdry * 1.06, 1500)
sol = s.solve(tdry * 1.06, t_eval=list(tv), rtol=1e-10, atol=at)
C3 = np.array([s.interp_u(sol.y[:, i], U)[0] for i in range(len(tv))])
hh = tv / 3600.0

print("tdry = %.3f s = %.4f h" % (tdry, tdry / 3600.0))
print("横轴范围 0 到 %.3f 小时，共 %d 个点" % (hh[-1], len(hh)))
print("")
print("%-8s %10s %10s %10s %10s %10s" % ("t/h", "0.0cm", "0.5cm", "1.0cm", "1.5cm", "2.0cm"))
for tt in [6, 12, 24, 36, 48]:
    i = int(np.argmin(np.abs(hh - tt)))
    print("%-8d %10.5f %10.5f %10.5f %10.5f %10.5f" %
          (tt, C3[i, 0], C3[i, 1], C3[i, 2], C3[i, 3], C3[i, 4]))
i = int(np.argmin(np.abs(hh - tdry / 3600.0)))
print("%-8s %10.5f %10.5f %10.5f %10.5f %10.5f" %
      ("tdry", C3[i, 0], C3[i, 1], C3[i, 2], C3[i, 3], C3[i, 4]))

print("")
print("放大区间 48 到 60 小时内的采样点数：%d" % int(np.sum((hh > 48) & (hh < 60))))
