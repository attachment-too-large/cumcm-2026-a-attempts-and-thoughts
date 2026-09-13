# -*- coding: utf-8 -*-
"""SEM 求解器：自收敛 + 与独立有限体积对拍 + 守恒对账"""
import time
import numpy as np
from hb_core import Ambient, props_q1, R0, T0_C, C0, H_CONV, HM_MASS
from sem_core import SEM1D, Model
from selftest import fv_solve, fv_eval

np.set_printoptions(precision=6, suppress=True)
SIG = np.array([0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0])


def run_sem(breaks, n, amb, props, t_end=1800.0, ts=(60., 300., 900., 1800.)):
    S = SEM1D(breaks, n)
    m = Model(S, props, amb)
    t0 = time.time()
    sol = m.solve(t_end, t_eval=[0.0] + list(ts), rtol=1e-11, atol=1e-13, dense=True)
    return m, sol, time.time() - t0


if __name__ == '__main__':
    amb = Ambient(np.array([0.0, 1e9]), np.array([60.0, 60.0]), np.array([0.05, 0.05]))
    print('--- (1) 自收敛（基准=4单元 n=24，ndof=97，t=1800 s）---')
    mR, sR, _ = run_sem((0., .5, .8, .95, 1.), 24, amb, props_q1())
    Tr = mR.S.interp(sR.y[:mR.ndof, -1], SIG)
    Cr = mR.S.interp(sR.y[mR.ndof:, -1], SIG)
    print('  基准 T =', Tr)
    print('  基准 C =', Cr)
    caselist = [((0., 1.), 16), ((0., 1.), 24), ((0., 1.), 32), ((0., 1.), 40),
                ((0., .6, .9, 1.), 12), ((0., .6, .9, 1.), 16), ((0., .6, .9, 1.), 20),
                ((0., .6, .9, 1.), 24), ((0., .5, .8, .95, 1.), 16)]
    for breaks, n in caselist:
        m, sol, dt = run_sem(breaks, n, amb, props_q1())
        T = m.S.interp(sol.y[:m.ndof, -1], SIG)
        C = m.S.interp(sol.y[m.ndof:, -1], SIG)
        print('  %-24s n=%2d ndof=%3d %5.2fs  dT=%9.2e  dC=%9.2e   Cs=%.8f'
              % (str(breaks), n, m.ndof, dt, np.max(np.abs(T - Tr)),
                 np.max(np.abs(C - Cr)), C[-1]))

    print('--- (2) 与独立有限体积解对拍（隐式Euler + Richardson 外推）---')
    def fv_extrap(ncell, dt):
        T1, C1, sc, _ = fv_solve(ncell, dt, 1800.0, props_q1(), amb, store_t=[1800.0])
        T2, C2, sc, _ = fv_solve(ncell, dt / 2, 1800.0, props_q1(), amb, store_t=[1800.0])
        return 2 * T2 - T1, 2 * C2 - C1, sc
    for nc, dt in [(600, 0.5), (1200, 0.25)]:
        t0 = time.time()
        Tf, Cf, sc = fv_extrap(nc, dt)
        fvT = fv_eval(Tf, sc, SIG)
        fvC = fv_eval(Cf, sc, SIG)
        print('  FV nc=%4d dt=%.2f->%.3f  %.1fs  Ts=%.6f Cs=%.6f'
              % (nc, dt, dt / 2, time.time() - t0, fvT[-1], fvC[-1]))
        print('   FV  T =', fvT)
        print('   FV  C =', fvC)
        print('   SEM-FV 对比: dT=%.2e dC=%.2e' % (np.max(np.abs(Tr - fvT)), np.max(np.abs(Cr - fvC))))

    print('--- (3) 守恒恒等式（离散层面应精确到舍入）---')
    import numpy as _np
    for tt in (0.0, 1.0, 100.0, 900.0, 1800.0):
        y = sR.sol(tt)
        yd = mR.rhs(tt, y)
        T = y[:mR.ndof]
        C = y[mR.ndof:]
        Tq = mR.S.gather(T)
        Cq = mR.S.gather(C)
        dTq = mR.S.gather(yd[:mR.ndof])
        dCq = mR.S.gather(yd[mR.ndof:])
        w = [el['w'] for el in mR.S.elems]
        rcp = [mR.pr['rho'](c) * mR.pr['cp'](c) for c in Cq]
        Erate = sum(_np.sum(wi * ri * ti) for wi, ri, ti in zip(w, rcp, dTq))
        Wrate = sum(_np.sum(wi * ti) for wi, ti in zip(w, dCq))
        Tinf, Cinf = amb(tt)
        ref_E = -(H_CONV / R0) * (T[mR.S.surf_dof] - Tinf)
        ref_W = -(HM_MASS / R0) * (C[mR.S.surf_dof] - Cinf)
        print('  t=%7.1f  dE/dt=%.6f (理论 %.6f, 绝对差 %.2e)   dW/dt=%.3e (理论 %.3e, 绝对差 %.2e)'
              % (tt, Erate, ref_E, abs(Erate - ref_E), Wrate, ref_W, abs(Wrate - ref_W)))
