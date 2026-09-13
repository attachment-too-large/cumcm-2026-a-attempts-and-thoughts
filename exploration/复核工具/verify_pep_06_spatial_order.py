# -*- coding: utf-8 -*-
"""空间收敛阶复核（避开时间误差平台）：常物性导热 vs Bessel 解析解。"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import brentq
sys.path.insert(0, r'C:\Users\qing1\Desktop\A题\复核工具')
import indep_solver as M
from indep_solver import FV

R0, k, rho, cp, h = 0.02, 0.36, 820.0, 2600.0, 25.0
alpha = k / (rho * cp)
Bi = h * R0 / k
Tair, T0 = 60.0, 28.0

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


class CE:
    def T(self, t):
        return Tair
    def C(self, t):
        return 0.0


_oT, _oC = M.Tinf, M.Cinf
M.Tinf, M.Cinf = CE().T, CE().C
r_ref = np.linspace(0, R0, 401)
print('Bessel 自检: t=0 中心 = %.5f（应 28）' % bessel_T(0.0, 0.0)[0])
ex = bessel_T(r_ref, 3600.0)
print()
print('空间收敛（dt=0.05 s，时间误差约 1e-5 K）：')
print('  %6s %16s %8s' % ('N', 'max|err| / K', '观测阶'))
prev = None
for N in (100, 200, 400, 800):
    s = FV('p1', False, N=N)
    T = np.full(N, T0); C = np.full(N, 2.55); t = 0.0; dt = 0.05
    while t < 3600.0 - 1e-9:
        hs = min(dt, 3600.0 - t)
        T, C = s.step(T, C, t, hs)
        t += hs
    num = np.interp(r_ref, R0 * s.xi, T)
    err = float(np.max(np.abs(num - ex)))
    print('  %6d %16.4e %8s' % (N, err, '' if prev is None else '%.2f' % np.log2(prev / err)))
    prev = err
print()
print('时间收敛（N=800，与 dt=0.0625 s 的参考解比较）：')
s = FV('p1', False, N=800)
T = np.full(800, T0); C = np.full(800, 2.55); t = 0.0
while t < 3600.0 - 1e-9:
    hs = min(0.0625, 3600.0 - t); T, C = s.step(T, C, t, hs); t += hs
ref = T.copy()
prev = None
for dt in (4.0, 2.0, 1.0, 0.5, 0.25):
    s = FV('p1', False, N=800)
    T = np.full(800, T0); C = np.full(800, 2.55); t = 0.0
    while t < 3600.0 - 1e-9:
        hs = min(dt, 3600.0 - t); T, C = s.step(T, C, t, hs); t += hs
    err = float(np.max(np.abs(T - ref)))
    print('  dt=%6.4f  max|err| = %.4e K  观测阶 %s'
          % (dt, err, '' if prev is None else '%.2f' % np.log2(prev / err)))
    prev = err
M.Tinf, M.Cinf = _oT, _oC
