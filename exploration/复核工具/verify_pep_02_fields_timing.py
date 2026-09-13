# -*- coding: utf-8 -*-
"""针对 paper_electronic.pdf 的独立复核 · 第 2 部分：用自写求解器复算。

核验项：
  E. 问题1 表1/表2（1800 s）
  F. 问题2 表3/表4（3 h）
  G. 问题3 tdry（拟合环境 & 原始环境）
  H. 问题4 tdry（含收缩）
  I. 不收缩对照（附录4 物性、R 固定 2 cm）
  J. 截面平均的正确值（面积加权） vs 论文口径
"""
import sys, io, os, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import openpyxl
from scipy.optimize import curve_fit
from scipy.interpolate import PchipInterpolator

ROOT = r'C:\Users\qing1\Desktop\A题'
ATT = os.path.join(ROOT, '附件')
sys.path.insert(0, os.path.join(ROOT, '复核工具'))
import indep_solver as M
from indep_solver import FV

# ---- 环境数据与拟合 ----
wb = openpyxl.load_workbook(os.path.join(ATT, '附件1.xlsx'), data_only=True)
rows = [r for r in wb['Sheet1'].iter_rows(values_only=True)][1:]
t_env = np.array([float(r[0]) for r in rows])
Ta = np.array([float(r[1]) for r in rows])
Ca = np.array([float(r[2]) for r in rows])
f1 = lambda tt, A, tau: A - (A - 28.0) * np.exp(-tt / tau)
g1 = lambda tt, A, tau: A - (A - 0.01963) * np.exp(-tt / tau)
pT, _ = curve_fit(f1, t_env, Ta, p0=[50.2, 1877.0], maxfev=20000)
pC, _ = curve_fit(g1, t_env, Ca, p0=[0.0509, 2776.0], maxfev=20000)
print('环境拟合: Tset=%.4f tau=%.1f | Cset=%.6f tau=%.1f' % (pT[0], pT[1], pC[0], pC[1]))


class RawEnv:
    """附件1 原始数据 + PCHIP（t>14400 取末点常值）"""
    def __init__(self):
        self.fT = PchipInterpolator(t_env, Ta)
        self.fC = PchipInterpolator(t_env, Ca)
    def T(self, tt):
        return float(self.fT(tt)) if tt <= 14400 else float(Ta[-1])
    def C(self, tt):
        return float(self.fC(tt)) if tt <= 14400 else float(Ca[-1])


class FitEnv:
    """论文式(14)：0-4h 一阶惯性曲线，t>4h 恒定"""
    def T(self, tt):
        return f1(min(tt, 14400.0), *pT)
    def C(self, tt):
        return g1(min(tt, 14400.0), *pC)


crit = ['0', '0.5', '1', '1.5', '2']
A = 0.02

# ---------------- E 问题1 ----------------
print()
print('=' * 78)
print('E. 问题1：独立解 vs 论文表1/表2（1800 s，环境=拟合曲线，附录2 物性）')
_oT, _oC = M.Tinf, M.Cinf
M.Tinf, M.Cinf = FitEnv().T, FitEnv().C


def run_eval(mode, shrink, t_end, dt, times, N=800, epsi=1e-9):
    s = FV(mode, shrink, N=N)
    T = np.full(N, M.T0); C = np.full(N, M.C0); tt = 0.0
    out = {}
    want = sorted(times)
    wi = 0
    while tt < t_end - 1e-12 and wi < len(want):
        nxt = want[wi]
        h = min(dt, max(1e-6, nxt - tt), t_end - tt)
        T, C = s.step(T, C, tt, h)
        tt += h
        if abs(tt - nxt) < 1e-6:
            out[nxt] = (T.copy(), C.copy())
            wi += 1
    return s, out


s1, o1 = run_eval('p1', False, 1800.0, 1.0, [100, 300, 600, 900, 1200, 1500, 1800])
rn = A * s1.xi
paper_T1 = {100: [28.0001, 28.0005, 28.0054, 28.0425, 28.2219],
            300: [28.0500, 28.0776, 28.1857, 28.4550, 29.0261],
            600: [28.5456, 28.6438, 28.9618, 29.5667, 30.5576],
            900: [29.5568, 29.7074, 30.1740, 30.9976, 32.2359],
            1200: [30.9069, 31.0882, 31.6397, 32.5824, 33.9433],
            1500: [32.4397, 32.6353, 33.2248, 34.2150, 35.6117],
            1800: [34.0425, 34.2411, 34.8362, 35.8249, 37.1989]}
paper_C1 = {100: [2.5500, 2.5500, 2.5500, 2.5500, 2.2471],
            300: [2.5500, 2.5500, 2.5500, 2.5492, 2.0510],
            600: [2.5500, 2.5500, 2.5500, 2.5353, 1.8773],
            900: [2.5500, 2.5500, 2.5497, 2.5046, 1.7552],
            1200: [2.5500, 2.5500, 2.5482, 2.4646, 1.6592],
            1500: [2.5500, 2.5499, 2.5445, 2.4207, 1.5794],
            1800: [2.5500, 2.5497, 2.5383, 2.3755, 1.5109]}
rc = [0.0, 0.005, 0.01, 0.015, 0.02]
print('%-6s %-30s %-30s' % ('t/s', '温度 我/论文（最大差）', '浓度 我/论文（最大差）'))
wT = wC = 0.0
for k in sorted(paper_T1):
    Tv, Cv = o1[float(k)]
    Ti = [float(np.interp(x, rn, Tv)) for x in rc]
    Ci = [float(np.interp(x, rn, Cv)) for x in rc]
    dT = max(abs(a - b) for a, b in zip(Ti, paper_T1[k]))
    dC = max(abs(a - b) for a, b in zip(Ci, paper_C1[k]))
    wT = max(wT, dT); wC = max(wC, dC)
    print('%-6d %.4f / %.4f  差 %.4f      %.4f / %.4f  差 %.4f'
          % (k, Ti[0], paper_T1[k][0], abs(Ti[0] - paper_T1[k][0]),
             Ci[4], paper_C1[k][4], abs(Ci[4] - paper_C1[k][4])))
print('   全表最大差: 温度 %.4f °C, 浓度 %.4f kg/kg' % (wT, wC))

# 截面平均正确值
Tv, Cv = o1[1800.0]
wv = 2 * rn / A ** 2
print('   1800 s 截面平均（我的解，21 点面积加权）= %.4f  [论文 2.2907]' % np.average(Cv, weights=wv))
print('   1800 s 截面平均（等权 21 点）           = %.4f' % Cv.mean())
m_d = 820 / (1 + 2.55) * np.pi * A ** 2 * 0.25
print('   按面积加权算失水 = %.4f g   [论文 18.71]' % (m_d * (2.55 - np.average(Cv, weights=wv)) * 1e3))

# ---------------- F 问题2 ----------------
print()
print('=' * 78)
print('F. 问题2：独立解 vs 论文表3/表4（3 h，拟合环境，附录3 物性）')
s2, o2 = run_eval('p3', False, 10800.0, 2.0,
                  [1800, 3600, 5400, 7200, 9000, 10800], N=800)
paper_T3 = {1800: [32.5656, 32.7627, 33.3570, 34.3573, 35.7831],
            3600: [40.4220, 40.5818, 41.0548, 41.8230, 42.8629],
            5400: [45.5851, 45.6719, 45.9272, 46.3366, 46.8815],
            7200: [48.2253, 48.2650, 48.3814, 48.5675, 48.8138],
            9000: [49.4129, 49.4293, 49.4775, 49.5546, 49.6565],
            10800: [49.9051, 49.9115, 49.9302, 49.9602, 49.9999]}
paper_C3 = {1800: [2.5499, 2.5488, 2.5247, 2.3237, 1.6538],
            3600: [2.5248, 2.4935, 2.3562, 2.0227, 1.4709],
            5400: [2.3860, 2.3258, 2.1350, 1.8020, 1.3444],
            7200: [2.1731, 2.1105, 1.9250, 1.6256, 1.2282],
            9000: [1.9589, 1.9025, 1.7369, 1.4716, 1.1157],
            10800: [1.7673, 1.7174, 1.5705, 1.3332, 1.0083]}
rn2 = A * s2.xi
wT2 = wC2 = 0.0
for k in sorted(paper_T3):
    Tv, Cv = o2[float(k)]
    Ti = [float(np.interp(x, rn2, Tv)) for x in rc]
    Ci = [float(np.interp(x, rn2, Cv)) for x in rc]
    dT = max(abs(a - b) for a, b in zip(Ti, paper_T3[k]))
    dC = max(abs(a - b) for a, b in zip(Ci, paper_C3[k]))
    wT2 = max(wT2, dT); wC2 = max(wC2, dC)
    print('   t=%5d 中心T 我%.4f/论文%.4f 差%.4f | 表面C 我%.4f/论文%.4f 差%.4f'
          % (k, Ti[0], paper_T3[k][0], abs(Ti[0] - paper_T3[k][0]),
             Ci[4], paper_C3[k][4], abs(Ci[4] - paper_C3[k][4])))
print('   全表最大差: 温度 %.4f °C, 浓度 %.4f kg/kg' % (wT2, wC2))

# ---------------- G/H/I 时长 ----------------
print()
print('=' * 78)
print('G/H/I. 烘干时长')


def td(mode, shrink, N, dt, tmax_h=200.0, env=None):
    if env is not None:
        M.Tinf, M.Cinf = env.T, env.C
    s = FV(mode, shrink, N=N)
    r, C, T = s.till(0.15, 3600 * tmax_h, dt)
    return (r / 3600) if r else None


print('   问题3（附录3, R 固定）：')
for N, dt in ((400, 10.0), (800, 10.0)):
    print('     拟合环境 N=%d dt=%.0fs -> %.4f h   [论文 57.108 / 外推 57.099]'
          % (N, dt, td('p3', False, N, dt, env=FitEnv())))
for N, dt in ((400, 10.0), (800, 10.0)):
    print('     原始环境 N=%d dt=%.0fs -> %.4f h   [论文 10.7 节 +0.118 h 对照]'
          % (N, dt, td('p3', False, N, dt, env=RawEnv())))
print('   问题4（附录4, 含收缩）：')
for N, dt in ((400, 10.0), (800, 10.0)):
    print('     拟合环境 N=%d dt=%.0fs -> %.4f h   [论文 50.769 / 外推 50.768]'
          % (N, dt, td('p4', True, N, dt, env=FitEnv())))
print('   问题4 不收缩对照（附录4 物性, R≡2 cm）：')
for N, dt in ((400, 10.0), (800, 10.0)):
    print('     拟合环境 N=%d dt=%.0fs -> %.4f h   [论文 129.01]'
          % (N, dt, td('p4', False, N, dt, env=FitEnv(), tmax_h=300.0)))
M.Tinf, M.Cinf = _oT, _oC
