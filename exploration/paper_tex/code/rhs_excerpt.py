# -*- coding: utf-8 -*-
"""Excerpt 2: right-hand side, time integration and the drying-time event
(from sem_core.py / run_all.py)."""
import numpy as np
from scipy.integrate import solve_ivp

R0, T0C, C0 = 0.02, 28.0, 2.55
H_CONV, HM_MASS = 25.0, 8.0e-7
C_TARGET = 0.15


class Model:
    """State y = [T nodal values (ndof); C nodal values (ndof)].
    mapping='affine' : material coordinate sigma = r/R(t) for uniform shrinkage."""

    def __init__(self, sem, props, amb, radius=None, mapping='affine'):
        self.S, self.pr, self.amb, self.rad = sem, props, amb, radius
        self.mapping, self.ndof = mapping, sem.ndof
        M0 = sem.assemble([np.ones_like(el['sig_q']) for el in sem.elems], 'L')
        from scipy.linalg import lu_factor
        self.M0_lu = lu_factor(M0)

    def fields(self, y):
        S = self.S
        T, C = y[:S.ndof], y[S.ndof:]
        return T, C, S.gather(T), S.gather(C)

    def rhs(self, t, y):
        S = self.S
        T, C, Tq, Cq = self.fields(y)
        R = float(self.rad(t))
        Tinf, Cinf = self.amb(t)
        # coefficients at quadrature points
        rcp = [self.pr['rho'](c) * self.pr['cp'](c) for c in Cq]
        kq = [self.pr['k'](c) for c in Cq]
        Dq = [self.pr['D'](c, tt) for c, tt in zip(Cq, Tq)]
        MT = S.assemble(rcp, 'L')          # int rho*cp*phi_i*phi_j*sigma dsigma
        AT = S.assemble(kq, 'DL')          # int k*phi_i'*phi_j'*sigma dsigma
        AC = S.assemble(Dq, 'DL')
        e_s = np.zeros(S.ndof)
        e_s[S.surf_dof] = 1.0              # surface indicator (phi_j(1) = 1)
        Ts, Cs = T[S.surf_dof], C[S.surf_dof]
        invR2 = 1.0 / (R * R)
        ra = -invR2 * (AT @ T) - (H_CONV / R) * (Ts - Tinf) * e_s
        rc = -invR2 * (AC @ C) - (HM_MASS / R) * (Cs - Cinf) * e_s
        from scipy.linalg import lu_solve
        dT = np.linalg.solve(MT, ra)       # variable thermal inertia
        dC = lu_solve(self.M0_lu, rc)      # mass matrix is constant and diagonal
        return np.concatenate([dT, dC])

    def solve(self, t_end, t_eval=None, rtol=1e-10, atol=1e-12,
              events=None, dense=False, method='BDF'):
        y0 = np.concatenate([np.full(self.ndof, T0C), np.full(self.ndof, C0)])
        sol = solve_ivp(self.rhs, (0.0, t_end), y0, method=method,
                        t_eval=t_eval, rtol=rtol, atol=atol,
                        events=events, dense_output=dense)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol


def drying_time(model, t_max=4.0e5, rtol=1e-9, atol=1e-11):
    """Root finding on the continuous event g(t) = C(0,t) - 0.15."""
    S = model.S
    cdof = int(np.argmin(S.sig_nodes))     # node at sigma = 0

    def ev(t, y):
        return y[S.ndof + cdof] - C_TARGET
    ev.terminal, ev.direction = True, -1
    sol = model.solve(t_max, t_eval=[0.0], rtol=rtol, atol=atol,
                      events=ev, dense=True)
    if sol.t_events[0].size == 0:
        raise RuntimeError('target not reached')
    return float(sol.t_events[0][0])
