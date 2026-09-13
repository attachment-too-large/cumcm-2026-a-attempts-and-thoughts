# -*- coding: utf-8 -*-
"""能量审计（量纲正确版）：统一到"单位长度"基
   干物质线密度 m_d' = pi*R0^2 / v0 ,  v0 = (1+C0)/rho(C0)  [m^3/kg干]
   Q_conv = 2*pi*R*h*int(Tinf-Ts)dt                [J/m]
   dE     = m_d' * sum(dsm*capT*T)|0^T             [J/m]
   E_evap = m_d' * mean(C0-C) * L                  [J/m]
"""
import numpy as np
from hb_core import (Ambient, load_att1, load_att2, Radius, props_q1, props_q23,
                     props_q4, R0, LZ, H_CONV, HM_MASS, C0)
from sem_core import SEM1D, Model

L_VAP = 2.4e6


def run_audit(tag, props, radius, t_end, mesh=None, rtol=1e-8, atol=1e-10, cap_mode='fixed_vol'):
    S = SEM1D(**(mesh or dict(breaks=(0., .6, .9, 1.), n=16)))
    t, T, C = load_att1()
    amb = Ambient(t, T, C)
    m = Model(S, props, amb, radius)
    tel = np.unique(np.concatenate([np.linspace(0, t_end, 500), [t_end]]))
    sol = m.solve(t_end, t_eval=tel, rtol=rtol, atol=atol, dense=True)
    sig, w = S.global_quad()          # int_0^1 (.) sigma dsigma 权重
    v0 = (1 + C0) / float(props['rho'](np.array([C0]))[0])
    md_line = np.pi * R0 ** 2 / v0    # 干物质线密度 kg/m
    Qconv = 0.0
    for k in range(len(sol.t) - 1):
        tm = 0.5 * (sol.t[k] + sol.t[k + 1]); dt = sol.t[k + 1] - sol.t[k]
        ym = sol.sol(tm)
        Tinf, _ = amb(tm)
        R = float(radius(tm)) if radius is not None else R0
        Qconv += 2 * np.pi * R * H_CONV * (Tinf - ym[S.surf_dof]) * dt

    def energy(y):
        Tq = S.interp(y[:S.ndof], sig)
        Cq = S.interp(y[S.ndof:], sig)
        if cap_mode == 'fixed_vol':
            capT = props['rho'](Cq) * props['cp'](Cq) * v0      # J/(kg干 K)
        else:
            capT = (1 + Cq) * props['cp'](Cq)
        E = float(np.sum(w * capT * Tq)) / w.sum() * 1.0        # 平均 J/kg干
        return E, Cq
    E0, Cq0 = energy(sol.y[:, 0])
    E1, Cq1 = energy(sol.y[:, -1])
    dE = md_line * (E1 - E0)
    dC = float(np.sum(w * (Cq0 - Cq1)) / w.sum())
    E_evap = md_line * dC * L_VAP
    print('  %-14s Q_conv=%9.3e J/m  dE=%9.3e J/m  Cbar降=%.4f kg/kg  m_evap=%.4f kg/m  E_lat=%9.3e J/m'
          % (tag, Qconv, dE, dC, md_line * dC, E_evap))
    print('                 干物质线密度=%.4f kg/m   (Q_conv-dE)/E_lat = %.4f  -> 供热仅为潜热需求的 1/%.1f'
          % (md_line, (Qconv - dE) / E_evap, E_evap / max(Qconv - dE, 1e-30)))
    return dict(Qconv=Qconv, dE=dE, dC=dC, E_evap=E_evap, md=md_line)


if __name__ == '__main__':
    t2, R2 = load_att2()
    print('== 能量审计（统一到单位长度基, 量纲一致）==')
    run_audit('问题1 0-1800s', props_q1(), None, 1800.0)
    run_audit('问题3 0-60h', props_q23(), None, 216000.0, rtol=1e-7, atol=1e-9)
    run_audit('问题4 0-60h', props_q4(), Radius(fixed=False, t=t2, R=R2), 216000.0,
              rtol=1e-7, atol=1e-9, cap_mode='mass')
