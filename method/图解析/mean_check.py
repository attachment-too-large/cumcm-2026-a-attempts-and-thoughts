# -*- coding: utf-8 -*-
"""用三种算法计算截面平均含水率，确认哪一种才是真正的面积加权平均。

真正的定义： Cbar = (2/R^2) * 积分 C(r) r dr，从 0 到 R。
换元 u 等于半径比平方以后，dr 与 r 里的根号相消，得到
    Cbar = 积分 C(u) du，从 0 到 1。
所以截面平均等于 C 在 u 上的普通积分，与半径无关。

三种算法：
  甲  直接把节点值按节点半径加权求平均（节点分布不均匀，是错的）
  乙  用 Clenshaw-Curtis 权重在切比雪夫节点上精确积分（多项式意义下精确）
  丙  插值到均匀 r 网格后按面积加权做梯形积分（与乙互为独立验证）
"""
import sys
import numpy as np

CODEDIR = r"C:\Users\qing1\Desktop\bestway\药材烘干问题_正则化谱方法\代码"
sys.path.insert(0, CODEDIR)
import common as cm
from spectral import HerbSolver


def bary_matrix(u_nodes, u_targets):
    """返回重心插值矩阵，形状 (目标点数, 节点数)。"""
    N = len(u_nodes) - 1
    w = (-1.0) ** np.arange(N + 1)
    w[0] *= 0.5
    w[N] *= 0.5
    diff = u_targets[:, None] - u_nodes[None, :]
    num = w[None, :] * np.where(np.abs(diff) < 1e-14, 1.0, 1.0 / diff)
    den = w[None, :] * np.where(np.abs(diff) < 1e-14, 0.0, 1.0 / diff)
    return num / den.sum(axis=1, keepdims=True)


def cc_weights(N):
    """切比雪夫高斯洛巴托节点上的 Clenshaw-Curtis 积分权重，对应区间 [-1,1]。"""
    j = np.arange(N + 1)
    c = np.ones(N + 1)
    c[1:N] = 2.0
    w = np.zeros(N + 1)
    K = N // 2
    for jj in range(N + 1):
        s = 0.0
        for k in range(1, K + 1):
            b = 1.0 if 2 * k == N else 2.0
            s += b / (4.0 * k * k - 1.0) * np.cos(2.0 * k * jj * np.pi / N)
        w[jj] = c[jj] / N * (1.0 - s)
    return w


s = HerbSolver(64, prob=3)
at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])


def ev(t, y):
    return np.max(y[:s.Np]) - 0.15


ev.terminal = True
ev.direction = -1
solE = s.solve(300000.0, rtol=1e-10, atol=at, events=ev)
tdry = float(solE.t_events[0][0])
tv = np.linspace(0.0, tdry, 600)
sol = s.solve(tdry, t_eval=list(tv), rtol=1e-10, atol=at)
hh = tv / 3600.0
Call = sol.y[:s.Np, :]
u = s.u
r_nodes = np.sqrt(u) * cm.R0

# 甲：按节点半径直接加权
mA = (Call * r_nodes[:, None]).sum(axis=0) / r_nodes.sum()

# 乙：Clenshaw-Curtis 权重，Cbar = 0.5 * 积分[-1,1]
wcc = cc_weights(s.N)
mB = 0.5 * (wcc @ Call)

# 丙：均匀 r 网格上的面积加权梯形积分
rf = np.linspace(0.0, cm.R0, 601)
B = bary_matrix(u, (rf / cm.R0) ** 2)
Cf = B @ Call                                    # 形状 (r 点数, 时刻数)
mC = (np.trapezoid(Cf * rf[:, None], rf, axis=0)) * 2.0 / cm.R0 ** 2

print("烘干时长 %.3f 秒 = %.4f 小时" % (tdry, tdry / 3600.0))
print("Clenshaw-Curtis 权重之和应为 2：%.12f" % wcc.sum())
print("")
print("%-6s %14s %14s %14s" % ("t/h", "甲 节点加权", "乙 谱积分", "丙 r 积分"))
for tt in [6, 12, 24, 30, 36, 48, 57]:
    i = int(np.argmin(np.abs(hh - tt)))
    print("%-6d %14.6f %14.6f %14.6f" % (tt, mA[i], mB[i], mC[i]))

print("")
for name, m in [("甲 节点加权", mA), ("乙 谱积分", mB), ("丙 r 积分", mC)]:
    idx = np.where(m < 0.15)[0]
    if idx.size:
        i = idx[0]
        print("%-12s 第一次低于 0.15 在 %.3f 小时" % (name, hh[i]))
    else:
        print("%-12s 区间内未低于 0.15" % name)

print("")
print("乙与丙最大差别 %.3e（两者独立，说明积分正确）" % np.max(np.abs(mB - mC)))
print("甲与乙最大差别 %.3e（说明节点直接加权确实有偏差）" % np.max(np.abs(mA - mB)))
