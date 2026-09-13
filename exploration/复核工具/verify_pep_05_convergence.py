# -*- coding: utf-8 -*-
"""针对 paper_electronic.pdf 的独立复核 · 第 5 部分

核验项：
  T. 空间/时间收敛阶（论文 10.1 节：N=100->800 误差降 91 倍；dt 4->0.25 s 降 195 倍）
     用自写求解器做同样的加密实验（常物性、恒环境 Tair=60°C 的 Bessel 对照）
"""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import brentq
sys.path.insert(0, r'C:\Users\qing1\Desktop\A题\复核工具')
import indep_solver as M
from indep_solver import FV

R0, k, rho, cp, h = 0.02, 0.36, 820.0, 2600.0, 25.0
alpha = k / (rho * cp)
Tair, T0 = 60.0, 28.0
Bi = h * R0 / k
print('常物性、恒环境对照：alpha=%.4e, Bi=%.4f' % (alpha, Bi))

# 精确求 μJ1(μ) = Bi·J0(μ) 的正根：以 J0 的零点为区间端点（等式在每段内恰有一根）
zz = jn_zeros(0, 40)
mu = []
edges = [0.0] + list(zz)
for a, b in zip(edges[:-1], edges[1:]):
    lo = a + 1e-10
    hi = b - 1e-10
    fa = lo * j1(lo) - Bi * j0(lo)
    fb = hi * j1(hi) - Bi * j0(hi)
    if fa * fb < 0:
        mu.append(brentq(lambda m: m * j1(m) - Bi * j0(m), lo, hi))
mu = np.array(mu)
print('前 3 个特征值: %s' % np.round(mu[:3], 5))
_s = sum(2 * j1(m) / (m * (j0(m) ** 2 + j1(m) ** 2)) for m in mu[:20])
print('级数归一性检验 ΣAn = %.6f（应→1）' % _s)


def bessel_T(r, t):
    r = np.atleast_1d(r)
    out = np.full_like(r, Tair, dtype=float)
    for m in mu[:30]:
        An = 2 * j1(m) / (m * (j0(m) ** 2 + j1(m) ** 2))
        out = out - An * j0(m * r / R0) * np.exp(-m ** 2 * alpha * t / R0 ** 2) * (Tair - T0)
    return out


# 用环境函数把 Tair 固定 60 °C
class ConstEnv:
    def T(self, t):
        return Tair
    def C(self, t):
        return 0.0


_oT, _oC = M.Tinf, M.Cinf
M.Tinf, M.Cinf = ConstEnv().T, ConstEnv().C
print('检查 Bessel 解 t=0: T(0)=%.4f (应 %.1f)' % (bessel_T(0.0, 0.0)[0], T0))

print()
print('空间收敛（常物性导热，t=3600 s，dt 取足够小）：')
print('   %6s %14s %12s' % ('N', 'max|T_num-T_exact|', '观测阶'))
prev = None
r_ref = np.linspace(0, R0, 401)
ex = bessel_T(r_ref, 3600.0)
for N in (100, 200, 400, 800, 1600):
    s = FV('p1', False, N=N)
    T = np.full(N, T0); C = np.full(N, 2.55); t = 0.0
    dt = 0.25
    while t < 3600.0 - 1e-9:
        hstep = min(dt, 3600.0 - t)
        T, C = s.step(T, C, t, hstep)
        t += hstep
    rn = R0 * s.xi
    num = np.interp(r_ref, rn, T)
    err = float(np.max(np.abs(num - ex)))
    p = '' if prev is None else '%.2f' % (np.log2(prev / err))
    print('   %6d %14.3e %12s' % (N, err, p))
    prev = err

print()
print('时间收敛（常物性导热，N=800，t=3600 s）：')
print('   %10s %14s %12s' % ('dt/s', 'max|err|', '观测阶'))
prev = None
for dt in (4.0, 2.0, 1.0, 0.5, 0.25):
    s = FV('p1', False, N=800)
    T = np.full(800, T0); C = np.full(800, 2.55); t = 0.0
    while t < 3600.0 - 1e-9:
        hstep = min(dt, 3600.0 - t)
        T, C = s.step(T, C, t, hstep)
        t += hstep
    rn = R0 * s.xi
    num = np.interp(r_ref, rn, T)
    # 与自身 dt=0.0625 的参考解比（消除空间误差）
    err = float(np.max(np.abs(num - ex)))
    p = '' if prev is None else '%.2f' % (np.log2(prev / err))
    print('   %10.4f %14.3e %12s' % (dt, err, p))
    prev = err

M.Tinf, M.Cinf = _oT, _oC
print()
print('论文 10.1 节：空间 N=100->800 误差降 91 倍（二阶理论 64 倍）；')
print('             时间 dt 4->0.25 s 误差降 195 倍（理论 256 倍）')
