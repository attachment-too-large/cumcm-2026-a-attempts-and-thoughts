# -*- coding: utf-8 -*-
"""M1 容量因子偏差的完整时间序列证据（无需推断）。

对每个输出时刻，用同一个状态 C(t) 比较两种"单元水量"：
  代码口径 : Σ_i C_i · M · Δs_i                  (M = πR²/I)   —— 代码时间推进实际所用的容量
  物理口径 : Σ_i C_i · πR²·ΔJ_i / v0             —— 真实 ∫ ρ_w dV（v0 = 纯干物质比容）
若两者始终相差常数倍 v0/I → 证明容量项少乘了 v0/I（在 geometric 模式下 v≡v0 故恰好为 1，
所以问题 1—4 基线不受影响）。
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, r'C:\Users\qing1\Desktop\A题\cumcm2026-a-herb-drying-main\outputs\code')
import numpy as np
import herb_v2 as v2

t_oven, T_inf, C_inf = v2.hm.load_oven()
env = v2.hm.Environment(t_oven, T_inf, C_inf)
t_rad, R_rad = v2.hm.load_radius()
rad = v2.hm.RadiusLaw(t_rad, R_rad)
p4 = v2.PropsV2('p4')
C0 = 2.55
v_d0 = 1.0 / (p4.rho(C0) / (1.0 + C0))

for mapping in ('geometric', 'local'):
    s = v2.LagSolver(p4, rad, env, N=400, mapping=mapping, cap_mode='mass')
    T = np.full(401, 28.0); C = np.full(401, C0); t = 0.0
    print('=' * 72)
    print('mapping =', mapping, '  (v0 = %.6e)' % v_d0)
    print('%8s %14s %14s %10s %10s' % ('t/h', 'Σ C·MΔs (代码)', 'Σ C·πR²ΔJ/v0', '比值', 'v0/I'))
    for k in range(3600 // 60 * 30):
        T, C = s.step(T, C, t, 60.0)
        t += 60.0
        if (k + 1) % (3600 // 60 * 5) == 0:
            R = float(rad.R_of(t))
            r, J, I, v = s.mapping_arrays(C, R)
            M = np.pi * R ** 2 / I
            m_code = float(np.sum(C * M * s.dsm))
            dJ = np.diff(J)
            m_phys = float(np.sum(0.5 * (C[:-1] + C[1:]) * np.pi * R ** 2 * dJ / v_d0))
            print('%8.1f %14.6e %14.6e %10.4f %10.4f'
                  % (t / 3600, m_code, m_phys, m_code / m_phys, v_d0 / I))
