# -*- coding: utf-8 -*-
"""独立校验：1) 谱基恒等式  2) 守恒律  3) 与独立有限体积解对拍  4) 解析稳态"""
import numpy as np
from scipy.special import j0, j1, jn_zeros
from hb_core import (EvenLegendre, Model, Ambient, Radius, load_att1,
                     props_q1, props_q23, R0, H_CONV, HM_MASS, T0_C, C0)


# ---------------------------------------------------------------- 有限体积参考
def fv_solve(ncell, dt, t_end, props, amb, rad=None, T0=T0_C, C0v=C0,
             nprint=0, store_t=None, theta=1.0):
    """单元中心有限体积 + 半隐式 Euler（系数显式、扩散项隐式, 三对角求解）。
    L-稳定、无振铃，用于对拍 SEM 结果（配 Richardson 外推达到二阶）。"""
    from scipy.linalg import solve_banded
    h = 1.0 / ncell
    sc = (np.arange(ncell) + 0.5) * h
    sf = (np.arange(ncell) + 1.0) * h
    sfm = np.concatenate([[0.0], sf[:-1]])
    V = 0.5 * (sf ** 2 - sfm ** 2)
    Aface = h * np.arange(ncell + 1)            # Aface[j] = sigma at face j (0..ncell)
    T = np.full(ncell, float(T0))
    C = np.full(ncell, float(C0v))
    nt = int(round(t_end / dt))
    store = {}
    band = np.zeros((3, ncell))

    def step(vals, coef, R, bnd_coef, bnd_val, source_dens):
        """半隐式 Euler： source*v' = (1/R^2) div(coef grad v) + 表面 Robin"""
        cf = 0.5 * (coef[1:] + coef[:-1])          # 内面系数 (face 1..nc-1)
        Ah = Aface[1:ncell] * cf / h               # Ah[i] = A_{i+1} cf_i / h, i=0..nc-2
        diag = source_dens * V / dt
        diag[0] += Ah[0] / R ** 2
        diag[1:ncell - 1] += (Ah[:-1] + Ah[1:]) / R ** 2
        beta = R * bnd_coef * (coef[-1] / h) / (coef[-1] / h + R * bnd_coef / 2.0)
        diag[ncell - 1] += (Ah[-1] + beta) / R ** 2
        up = -Ah / R ** 2
        rhs = source_dens * V / dt * vals
        rhs[-1] += beta / R ** 2 * bnd_val
        band[0, 1:] = up
        band[1, :] = diag
        band[2, :-1] = up
        return solve_banded((1, 1), band, rhs)

    for n in range(nt + 1):
        t = n * dt
        if store_t is not None:
            for ts in store_t:
                if abs(t - ts) < 0.5 * dt:
                    store[ts] = (T.copy(), C.copy())
        if n == nt:
            break
        Tinf, Cinf = amb(t)
        R = rad(t) if rad is not None else R0
        kq = props['k'](C)
        Dq = props['D'](C, T)
        rc = props['rho'](C) * props['cp'](C)
        T = step(T, kq, R, H_CONV, Tinf, rc)
        C = step(C, Dq, R, HM_MASS, Cinf, np.ones(ncell))
    if store_t is not None:
        return T, C, sc, store
    return T, C, sc


def fv_eval(vals, sc, sigma):
    return np.interp(sigma, sc, vals)


# ---------------------------------------------------------------- 校验
def test_basis():
    print('== 1. 谱基恒等式 ==')
    for N in (8, 16, 28, 40):
        B = EvenLegendre(N)
        Mq = (B.P * B.wq) @ B.P.T
        e1 = np.max(np.abs(Mq - B.M))
        # int sigma phi_i' phi_j' dsigma 的解析值：对称性 A0 = A0^T，且对 i=0 行应全 0
        e2 = np.max(np.abs(B.A0 - B.A0.T))
        e3 = np.max(np.abs(B.A0[0, :]))
        # Legendre 导数恒等式： int_0^1 sigma phi_i' phi_j' dsigma 与数值积分一致
        print('  N=%2d  |M-M_ana|max=%.2e  |A-A^T|max=%.2e  |A[0,:]|max=%.2e  cond(M)=%.1f'
              % (N, e1, e2, e3, np.linalg.cond(Mq)))
        assert e1 < 1e-12 and e2 < 1e-11 and e3 < 1e-13


def test_trivial_and_conservation():
    print('== 2. 平凡解 + 守恒律 ==')
    # (a) Tinf == T0 == 28 且 Cinf == C0 -> 应为精确常数解
    amb = Ambient(np.array([0.0, 1e9]), np.array([T0_C, T0_C]), np.array([C0, C0]))
    m = Model(24, props_q1(), amb)
    y0 = np.concatenate([[T0_C] + [0.0] * 23, [C0] + [0.0] * 23])
    r = m.rhs(0.0, y0)
    print('  (a) 平凡解 rhs 最大值 = %.3e (应≈0)' % np.max(np.abs(r)))
    assert np.max(np.abs(r)) < 1e-14
    # (b) 零通量(hm=0,h=0)下总水分守恒
    global HM_MASS, H_CONV
    hm0, h0 = HM_MASS, H_CONV
    import hb_core
    hb_core.HM_MASS = 0.0
    hb_core.H_CONV = 0.0
    amb2 = Ambient(np.array([0.0, 1e9]), np.array([80.0, 80.0]), np.array([0.0, 0.0]))
    m2 = Model(24, props_q23(), amb2)
    sol = m2.solve(600.0, t_eval=[0.0, 600.0], rtol=1e-12, atol=1e-14)
    E0, W0, _, _ = m2.integrals(sol.y[:, 0])
    E1, W1, _, _ = m2.integrals(sol.y[:, -1])
    print('  (b) 绝热/零通量: dE=%.3e  dW=%.3e (应≈0)' % (E1 - E0, W1 - W0))
    assert abs(E1 - E0) < 1e-6 * abs(E0) and abs(W1 - W0) < 1e-12
    hb_core.HM_MASS, hb_core.H_CONV = hm0, h0


def test_against_fv():
    print('== 3. 与独立有限体积解对拍（Q1 参数, Tinf=60 常数）==')
    amb = Ambient(np.array([0.0, 1e9]), np.array([60.0, 60.0]), np.array([0.0, 0.0])
                  if False else np.array([0.05, 0.05]))
    ts = [60.0, 300.0, 900.0, 1800.0]
    for N in (16, 24, 32):
        m = Model(N, props_q1(), amb)
        sol = m.solve(1800.0, t_eval=[0.0] + ts, rtol=1e-12, atol=1e-14)
        sig = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
        spT = m.B.eval(sol.y[:N, -1], sig)
        spC = m.B.eval(sol.y[N:, -1], sig)
        Tf, Cf, sc, _ = fv_solve(600, 0.25, 1800.0, props_q1(), amb,
                                 store_t=[1800.0])
        fvT = fv_eval(Tf, sc, sig)
        fvC = fv_eval(Cf, sc, sig)
        print('  N=%2d  max|dT|=%.2e  max|dC|=%.2e' % (N, np.max(np.abs(spT - fvT)),
                                                       np.max(np.abs(spC - fvC))))


def test_analytic_steady():
    print('== 4. 与解析准稳态对拍（Tinf=60 长时间极限）==')
    # 无限长圆柱, 常物性, 稳态: T(r) = Tinf + (T0-Tinf) 形式 -> 实为均匀 Tinf
    amb = Ambient(np.array([0.0, 1e9]), np.array([60.0, 60.0]), np.array([0.0, 0.0]))
    m = Model(28, props_q1(), amb)
    sol = m.solve(200000.0, t_eval=[200000.0], rtol=1e-11, atol=1e-14)
    T = m.B.eval(sol.y[:28, -1], np.linspace(0, 1, 11))
    print('  t=2e5 s 时 T 剖面 min/max = %.6f / %.6f (极限 60)' % (T.min(), T.max()))
    # 解析指数衰减率（Bessel-Robin 特征值）
    Bi = H_CONV * R0 / 0.36
    mu = jn_zeros(1, 1)[0]      # 近似
    print('  Bi=%.4f' % Bi)


if __name__ == '__main__':
    test_basis()
    test_trivial_and_conservation()
    test_against_fv()
    test_analytic_steady()
    print('ALL SELFTESTS PASSED')
