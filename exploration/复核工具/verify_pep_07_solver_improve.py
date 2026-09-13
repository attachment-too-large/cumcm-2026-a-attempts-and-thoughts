# -*- coding: utf-8 -*-
"""对照实验：一阶解算器 vs 自研二阶解算器（indep_solver_v2）。

(0) 自检：θ=1、arith 时 v2 必须与 indep_solver.FV 逐位一致
(1) 时间收敛阶（对 Bessel 解析解）
(2) 界面平均方式（算术/调和/对数）对空间精度的影响
(3) 效率：达到同一精度所需墙钟时间
(4) 对最终答案（问题3 tdry）的影响
"""
import sys, io, os, time
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


def set_env(fT, fC):
    """同时改写两个模块的绑定（v2 里 from ... import 是值绑定，必须一起改）。"""
    M.Tinf, M.Cinf = fT, fC
    V2.Tinf, V2.Cinf = fT, fC
    V2.Rof = M.Rof


set_env(lambda t: Tair, lambda t: 0.0)

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


r_ref = np.linspace(0, R0, 401)
ex = bessel_T(r_ref, 3600.0)


def run(cls, N, dt, t_end=3600.0, **kw):
    b = cls('p1', False, N=N, **kw)
    T = np.full(N, T0); C = np.full(N, 2.55); t = 0.0
    while t < t_end - 1e-9:
        hs = min(dt, t_end - t)
        T, C = b.step(T, C, t, hs)
        t += hs
    return b, T, C


print('=' * 80)
print('(0) 自检：v2(θ=1, arith) 应与 indep_solver.FV 逐位一致')
for N in (8, 40):
    a, Ta, Ca = run(FV, N, 600.0)
    b, Tb, Cb = run(FV2, N, 600.0, theta=1.0, n_ran=10 ** 6, face='arith')
    print('   N=%-3d  max|ΔT| = %.3e   max|ΔC| = %.3e' % (N, np.max(np.abs(Ta - Tb)), np.max(np.abs(Ca - Cb))))

print()
print('=' * 80)
print('(1) 时间收敛：对 Bessel 解析解（t=3600 s，N=800）')
print('    %-8s %16s %7s %16s %7s' % ('dt/s', 'BE 误差/K', '阶', 'CN+Ran 误差/K', '阶'))
pb = pc = None
for dt in (8.0, 4.0, 2.0, 1.0):
    _, Ta, _ = run(FV, 800, dt)
    eB = float(np.max(np.abs(np.interp(r_ref, R0 * (np.linspace(0, 1, 801)[:-1] + np.linspace(0, 1, 801)[1:]) / 2, Ta) - ex)))
    b, Tb, _ = run(FV2, 800, dt, face='harm')
    eC = float(np.max(np.abs(np.interp(r_ref, R0 * b.xi, Tb) - ex)))
    print('    %-8.2f %16.4e %7s %16.4e %7s'
          % (dt, eB, '' if pb is None else '%.2f' % np.log2(pb / eB),
             eC, '' if pc is None else '%.2f' % np.log2(pc / eC)))
    pb, pc = eB, eC

print()
print('=' * 80)
print('(2) 界面平均方式（CN+Ran，dt=0.25 s，t=3600 s）')
for face in ('arith', 'harm', 'log'):
    row = []
    for N in (100, 200, 400):
        b, T, _ = run(FV2, N, 0.25, face=face)
        row.append(float(np.max(np.abs(np.interp(r_ref, R0 * b.xi, T) - ex))))
    print('    face=%-6s 误差 %s   相邻比 %.2f'
          % (face, np.array2string(np.array(row), precision=4), row[0] / row[1]))

print()
print('=' * 80)
print('(3) 效率（N=800，t=3600 s，误差取与解析解之最大偏差）')
for tag, cls, kw in (('BE', FV, {}), ('CN+Ran/harm', FV2, {'face': 'harm'})):
    for dt in (4.0, 1.0, 0.25):
        t0 = time.time()
        b, T, _ = run(cls, 800, dt, **kw)
        el = time.time() - t0
        e = float(np.max(np.abs(np.interp(r_ref, R0 * b.xi, T) - ex)))
        print('    %-12s dt=%-6.2f 误差=%.3e  用时=%.1f s' % (tag, dt, e, el))

print()
print('=' * 80)
print('(4) 对最终答案的影响（问题3 tdry）')
print('    两法均用附件1 原始环境（PCHIP）')
M.Tinf, M.Cinf = M.Tinf, M.Cinf   # 已是常值 60/0，需换回原始环境
from scipy.interpolate import PchipInterpolator
import openpyxl
wb = openpyxl.load_workbook(r'C:\Users\qing1\Desktop\A题\附件\附件1.xlsx', data_only=True)
rows = [r for r in wb['Sheet1'].iter_rows(values_only=True)][1:]
te = np.array([float(r[0]) for r in rows])
Tec = np.array([float(r[1]) for r in rows])
Cec = np.array([float(r[2]) for r in rows])
fT = PchipInterpolator(te, Tec)
fC = PchipInterpolator(te, Cec)
set_env(lambda t: float(fT(min(t, 14400.0))), lambda t: float(fC(min(t, 14400.0))))
print('    %-14s %-6s %-8s %-10s' % ('solver', 'N', 'dt/s', 'tdry / h'))
for tag, cls, kw in (('BE(一阶)', FV, {}), ('CN+Ran(二阶)', FV2, {'face': 'harm'})):
    for N, dt in ((400, 10.0), (800, 10.0)):
        b = cls('p3', False, N=N, **kw)
        r, C, T = b.till(0.15, 3600 * 90, dt)
        print('    %-14s %-6d %-8.0f %-10.4f' % (tag, N, dt, r / 3600 if r else float('nan')))
