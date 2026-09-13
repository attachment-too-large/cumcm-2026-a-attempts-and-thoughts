# -*- coding: utf-8 -*-
"""确认 35.1 与 35.8 的差别来自积分点数，不是物理差异。

做法：用同一套谱解，分别插值到 21 个均匀半径点与 2001 个均匀半径点，
各自按面积加权做梯形积分，比较两者的截面平均值与达标时刻。
"""
import sys
import numpy as np

CODEDIR = r"C:\Users\qing1\Desktop\bestway\药材烘干问题_正则化谱方法\代码"
sys.path.insert(0, CODEDIR)
import common as cm
from spectral import HerbSolver


def bary_matrix(u_nodes, u_targets):
    N = len(u_nodes) - 1
    w = (-1.0) ** np.arange(N + 1)
    w[0] *= 0.5
    w[N] *= 0.5
    diff = u_targets[:, None] - u_nodes[None, :]
    safe = np.where(np.abs(diff) < 1e-14, 1.0, diff)
    num = (w[None, :] / safe)
    return num / num.sum(axis=1, keepdims=True)


s = HerbSolver(64, prob=3)
at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])


def ev(t, y):
    return np.max(y[:s.Np]) - 0.15


ev.terminal = True
ev.direction = -1
solE = s.solve(300000.0, rtol=1e-10, atol=at, events=ev)
tdry = float(solE.t_events[0][0])

# 用与结果文件相同的时间网格：每 60 秒一点
tt = np.arange(60.0, np.floor(tdry / 60.0) * 60.0 + 1.0, 60.0)
sol = s.solve(tdry, t_eval=list(tt), rtol=1e-11, atol=at)
Call = sol.y[:s.Np, :]
hh = tt / 3600.0
R = cm.R0


def mean_on(npts):
    r = np.linspace(0.0, R, npts)
    B = bary_matrix(s.u, (r / R) ** 2)
    Cf = B @ Call
    return np.trapezoid(Cf * r[:, None], r, axis=0) * 2.0 / R ** 2


for npts in [21, 41, 101, 401, 2001]:
    m = mean_on(npts)
    idx = np.where(m < 0.15)[0]
    tt_cross = hh[idx[0]] if idx.size else float("nan")
    print("半径点数 %-5d 达标时刻 %8.4f 小时    t=36 h 时的平均值 %.6f" %
          (npts, tt_cross, m[int(np.argmin(np.abs(hh - 36)))]))

print("")
print("结论：点数越多越接近真实积分。21 点对应交付文件的列数，")
print("      它与收敛值的差别就是梯形积分的离散误差。")

# 顺带看看 36 小时处剖面在表面的陡峭程度
i36 = int(np.argmin(np.abs(hh - 36)))
r21 = np.linspace(0.0, R, 21)
B = bary_matrix(s.u, (r21 / R) ** 2)
c21 = B @ Call[:, i36]
print("")
print("t=36 h 时 21 个采样点上的含水率（厘米，数值）：")
for j in range(0, 21, 2):
    print("   r=%.1f cm  C=%.6f" % (r21[j] * 100, c21[j]))
print("   最后两点之差 %.6f，说明表面附近变化最快，粗网格会低估积分" %
      (c21[-2] - c21[-1]))
