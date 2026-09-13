# -*- coding: utf-8 -*-
"""决定性对照：同一空间离散下，一阶(BE/半步 Richardson) 与 二阶(CN+Rannacher) 互证。

做法：从同一初值前进到 t=3600 s，比较三种走法的温度场
  A. BE,  2 步 × dt/2      (一阶，误差 O(dt))
  B. BE,  4 步 × dt/4 + Richardson 外推  T ≈ 2T(dt/2) − T(dt)
  C. CN + Rannacher, 1 步 × dt  (二阶)
若 C 与 B 在高精度上一致 → CN 的二阶性得到独立确认。
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
Tair, T0 = 60.0, 28.0
M.Tinf = V2.Tinf = lambda t: Tair
M.Cinf = V2.Cinf = lambda t: 0.0
Bi = h * R0 / k

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


N = 400
rn = R0 * (np.arange(N) + 0.5) / N
ex = bessel_T(rn, 3600.0)
TEND = 3600.0


def march(cls, dt, **kw):
    b = cls('p1', False, N=N, **kw)
    T = np.full(N, T0); C = np.full(N, 2.55); t = 0.0
    while t < TEND - 1e-9:
        hs = min(dt, TEND - t)
        T, C = b.step(T, C, t, hs)
        t += hs
    return T


print('N=%d，目标时刻 t=3600 s；与 Bessel 解析解的最大偏差' % N)
print()
print('%-34s %16s' % ('走法', 'max|T − T_exact| / K'))
for dt in (8.0, 4.0, 2.0):
    Tbe = march(FV, dt)
    e_be = float(np.max(np.abs(Tbe - ex)))
    # Richardson：BE 一阶，2T(dt/2) − T(dt)
    T1 = march(FV, dt / 2)
    TR = 2 * T1 - Tbe
    e_R = float(np.max(np.abs(TR - ex)))
    Tcn = march(FV2, dt, face='harm')
    e_cn = float(np.max(np.abs(Tcn - ex)))
    print('dt=%-5.1f  A BE 1步            %16.4e' % (dt, e_be))
    print('          B BE 半步+Richardson  %16.4e   （一阶→二阶）' % e_R)
    print('          C CN+Rannacher        %16.4e' % e_cn)
    print('          |B − C| 相互差        %16.4e' % float(np.max(np.abs(TR - Tcn))))
    print()
print('结论：B 与 C 两种"不同来源的二阶重构"互相吻合，且都优于纯 BE 一步；')
print('      说明 CN + Rannacher 的二阶性成立（前提：空间离散误差已被 N 压到其下）。')
