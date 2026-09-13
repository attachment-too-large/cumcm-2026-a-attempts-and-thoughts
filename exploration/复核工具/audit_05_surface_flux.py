# -*- coding: utf-8 -*-
"""核实：独立求解器是否严格满足表面质量收支 J_surf = km(C_s - C_inf)。"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, r'C:\Users\qing1\Desktop\A题')
import numpy as np
import indep_solver as M
from indep_solver import FV, Cinf

M.km = 8e-7
N = 800
s = FV('p3', False, N=N)
T = np.full(N, M.T0); C = np.full(N, M.C0); t = 0.0
dt = 60.0
A = 0.02 * 1.0                      # 2πH 归一后的表面积
cum_mass, cum_flux = 0.0, 0.0
mass0 = float(np.sum(C * (0.02 ** 2 * s.wV)))
hist = []
while t < 3600 * 60:
    Cs = C[-1]
    Js = M.km * (Cs - Cinf(t))       # 表面通量 kg/(m2 s)
    cum_flux += Js * dt * A
    T, C = s.step(T, C, t, dt)
    t += dt
    mass = float(np.sum(C * (0.02 ** 2 * s.wV)))
    cum_mass = mass0 - mass
    if abs(t % 3600) < 1e-9:
        hist.append((t / 3600, Cs, Js, cum_mass, cum_flux))
print('           t/h    C_surf      J_surf      累计失水(体)   累计表面通量   比值')
for h in hist:
    print('%14.1f  %9.5f  %11.3e  %12.5e  %12.5e  %7.4f'
          % (h[0], h[1], h[2], h[3], h[4], h[3] / h[4] if h[4] else float('nan')))
print()
print('起始总水量 = %.6e ; 60 h 后 = %.6e' % (mass0, mass))
print('需要排出的总水量（降到 0.15）≈ %.6e' % (mass0 - 0.15 * (0.02 ** 2 * s.wV).sum()))
