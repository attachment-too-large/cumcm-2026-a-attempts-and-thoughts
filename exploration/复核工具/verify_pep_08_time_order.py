# -*- coding: utf-8 -*-
"""干净的收敛实验：一阶(BE) vs 二阶(CN+Rannacher) 的时间收敛阶。

改进点：正确构造节点网格 r_i = R·(i+0.5)/N，并与 Bessel 解析解直接比较。
"""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import indep_solver as M
import indep_solver_v2 as V2
from indep_solver import FV
from indep_solver_v2 import FV2

R0, k, rho, cp, h = 0.02, 0.36, 820.0, 2600.0, 25.0
alpha = k / (rho * cp)
Bi = h * R0 / k
Tair, T0 = 60.0, 28.0
M.Tinf = V2.Tinf = lambda t: Tair
M.Cinf = V2.Cinf = lambda t: 0.0

zz = jn_zeros(0, 60)
mu = []
edges = [0.0] + list(zz)
for a, b in zip(edges[:-1], edges[1:]):
    lo, hi = a + 1e-10, b - 1e-10
    f = lambda m: m * j1(m) - Bi * j0(m)
    if f(lo) * f(hi) < 0:
        mu.append(brentq(f, lo, hi))
mu = np.array(mu)


def bessel_T(r, t):
    r = np.atleast_1d(r)
    out = np.full_like(r, Tair, dtype=float)
    for m in mu[:40]:
        An = 2 * j1(m) / (m * (j0(m) ** 2 + j1(m) ** 2))
        out = out - An * j0(m * r / R0) * np.exp(-m ** 2 * alpha * t / R0 ** 2) * (Tair - T0)
    return out


N = 800
rn = R0 * (np.arange(N) + 0.5) / N
ex = bessel_T(rn, 3600.0)


def solve(cls, dt, **kw):
    b = cls('p1', False, N=N, **kw)
    T = np.full(N, T0); C = np.full(N, 2.55); t = 0.0
    while t < 3600.0 - 1e-9:
        hs = min(dt, 3600.0 - t)
        T, C = b.step(T, C, t, hs)
        t += hs
    return float(np.max(np.abs(T - ex))), b


print('圆柱对流边界解析解对照（N=800，节点与解析解同点比较）')
print('  %-8s %16s %8s %16s %8s' % ('dt/s', 'BE 误差/K', '阶', 'CN+Ran 误差/K', '阶'))
pb = pc = None
for dt in (8.0, 4.0, 2.0, 1.0, 0.5, 0.25):
    eB, _ = solve(FV, dt)
    eC, _ = solve(FV2, dt, face='harm')
    print('  %-8.2f %16.4e %8s %16.4e %8s'
          % (dt, eB, '' if pb is None else '%.2f' % np.log2(pb / eB),
             eC, '' if pc is None else '%.2f' % np.log2(pc / eC)))
    pb, pc = eB, eC
print()
print('说明：CN+Ran 的误差很快触到空间离散平台；BE 需更小 dt 才能同级。')
print('论文 10.1 节报告的同类实验为：Δt 4→0.25 s 误差降 195 倍（二阶理论 256 倍）。')
