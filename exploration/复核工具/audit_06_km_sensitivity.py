# -*- coding: utf-8 -*-
"""复核 km 敏感性是否真实（排除网格/时间步假象）。"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, r'C:\Users\qing1\Desktop\A题')
import numpy as np
import indep_solver as M
from indep_solver import FV

print('%-10s %-6s %-6s %-12s %-10s' % ('km', 'N', 'dt', 't_dry / h', '相对'))
for km in (8e-7, 8e-6, 8e-5):
    M.km = km
    for N, dt in ((400, 60.0), (800, 60.0), (800, 15.0)):
        s = FV('p3', False, N=N)
        tt, C, T = s.till(0.15, 3600 * 300, dt)
        print('%-10.1e %-6d %-6.0f %-12.4f' % (km, N, dt, tt / 3600))
