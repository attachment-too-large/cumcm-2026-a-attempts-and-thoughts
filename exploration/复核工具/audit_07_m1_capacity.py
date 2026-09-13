# -*- coding: utf-8 -*-
"""复核 M1（local mapping）的守恒因子：直接复跑题库 herb_v2.LagSolver。

对比三组：
  A. geometric（仿射收缩）           —— 已有基线
  B. local（题库实现）
  C. local + 修正的容量因子 v0/I
"""
import sys, io, os, json, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, r'C:\Users\qing1\Desktop\A题\cumcm2026-a-herb-drying-main\outputs\code')
import numpy as np
import herb_v2 as v2

t_oven, T_inf, C_inf = v2.hm.load_oven()
env = v2.hm.Environment(t_oven, T_inf, C_inf)
t_rad, R_rad = v2.hm.load_radius()
rad = v2.hm.RadiusLaw(t_rad, R_rad)

te = np.arange(0, 120 * 3600 + 1, 60.0)


def dry_time(t, C, thr=0.15):
    Cmax = C.max(axis=1)
    idx = np.where(Cmax < thr)[0]
    if idx.size == 0:
        return None
    i = idx[0]
    return float(t[i - 1] + (thr - Cmax[i - 1]) * (t[i] - t[i - 1]) / (Cmax[i] - Cmax[i - 1]))


class LagFixed(v2.LagSolver):
    """在容量项上补回 M1 推导中要求的比容因子：Δm_i = M*(v0/I)*Δs_i。"""

    def _coeffs(self, T, C, R):
        gC, gT, gC_out, gT_out, gL, capT = v2.LagSolver._coeffs(self, T, C, R)
        r, J, I, v = self.mapping_arrays(C, R)
        self._dsm_eff = self.dsm * (self.v0 / I)
        return gC, gT, gC_out, gT_out, gL, capT

    def step(self, T_n, C_n, t, dt, maxiter=50, tol=1e-11, theta=None):
        th = self.theta if theta is None else theta
        R = float(self.radius.R_of(t + 0.5 * dt))
        Tinf_o = float(self.env.T_inf(t)); Tinf_n = float(self.env.T_inf(t + dt))
        Cinf_o = float(self.env.C_inf(t)); Cinf_n = float(self.env.C_inf(t + dt))
        T_g, C_g = T_n.copy(), C_n.copy()
        for it in range(maxiter):
            Tm = 0.5 * (T_n + T_g); Cm = 0.5 * (C_n + C_g)
            gC, gT, gC_out, gT_out, gL, capT = self._coeffs(Tm, Cm, R)
            cap = getattr(self, '_dsm_eff', self.dsm)
            C_new = self._solve(C_n, cap.copy(), gC, gC_out, Cinf_o, Cinf_n, dt, th)
            extra = 0.0
            T_new = self._solve(T_n, cap * capT, gT, gT_out, Tinf_o, Tinf_n, dt, th, extra=extra)
            d = max(float(np.max(np.abs(T_new - T_g))), float(np.max(np.abs(C_new - C_g))))
            T_g, C_g = T_new, C_new
            if d < tol:
                break
        return T_g, C_g


results = {}
for tag, cls, mapping in (('A geometric', v2.LagSolver, 'geometric'),
                          ('B local(题库)', v2.LagSolver, 'local'),
                          ('C local(修正容量)', LagFixed, 'local')):
    t0 = time.time()
    s = cls(v2.PropsV2('p4'), rad, env, N=400, mapping=mapping, cap_mode='mass')
    r = s.run(120 * 3600.0, v2.schedule_prod, t_eval=te)
    td = dry_time(r['t'], r['C'])
    results[tag] = td / 3600
    print('%-18s t_dry = %.4f h   (C_max末=%.8f, 用时 %.0fs)'
          % (tag, td / 3600, r['C'].max(axis=1)[np.where(r['C'].max(axis=1) < 0.15)[0][0]],
             time.time() - t0), flush=True)

print()
print('题库公布值：几何基线 50.8119 h；local(M1) 47.6208 h')
print('复核：', json.dumps({k: round(v, 4) for k, v in results.items()}, ensure_ascii=False))
print('A 与 B 之比 = %.4f；A 与 C 之比 = %.4f'
      % (results['A geometric'] / results['B local(题库)'],
         results['A geometric'] / results['C local(修正容量)']))
