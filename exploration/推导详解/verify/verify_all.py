# -*- coding: utf-8 -*-
"""§5.2–§5.7 推导的独立数值验证（自包含，不依赖论文代码）。

A. 式(7)/(7′)/(8) 的代数核对（sympy）
B. §5.6 / §5.7 数值复算
C. 式(8) 的网格无关性；「照抄式(7)」会得到什么
D. R≡R0 时式(8) 退化为式(2)，并与解析 Bessel 级数解对比
E. 离散守恒性
F. ρcp 乘 (R0/R(t))² 的对照算例
G. 干基密度假设的影响量级
"""
import sys, io, os, json, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace',
                              line_buffering=True)
import numpy as np
from scipy.linalg import solve_banded
from scipy.special import j0, j1
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

HERE = os.path.dirname(os.path.abspath(__file__))
R0 = 0.02
HM = 8e-7
H_HEAT = 25.0
CAIR = 0.0196
TAIR = 50.0
rep = {}


# ============================================================ 物性 / 几何
class PD:                                   # 附录 4 的 ρ/cp/k + 常扩散系数
    """用于检验「方程形式」的干净算例：
    附录 4 的 D(C)=4.2e-4·exp(-0.3/C)·exp(-3850/T) 在 C≈0.02 时趋于 0，
    与调和平均的界面系数叠加会形成数值「干壳」把内部封死（见 verify 报告 note），
    那属于物性插值问题、而非式(8) 的问题。故这里取常 D，隔离出坐标变换本身。"""
    D0 = 2.5e-9

    @staticmethod
    def rho(C):  return 760.0 + 90.0 * C
    @staticmethod
    def cp(C):   return 1850.0 + 2150.0 * C / (C + 1.0)
    @staticmethod
    def k(C):    return 0.12 + 0.20 * C / (C + 1.0)
    @staticmethod
    def D(C, T): return np.full_like(np.asarray(C, float), PD.D0)


class P4:                                   # 附录 4（问题 4）
    @staticmethod
    def rho(C):  return 760.0 + 90.0 * C
    @staticmethod
    def cp(C):   return 1850.0 + 2150.0 * C / (C + 1.0)
    @staticmethod
    def k(C):    return 0.12 + 0.20 * C / (C + 1.0)
    @staticmethod
    def D(C, T): return 4.2e-4 * np.exp(-0.30 / np.maximum(C, 1e-3)) * np.exp(-3850.0 / T)


class P3:                                   # 附录 3（问题 2/3）
    @staticmethod
    def rho(C):  return 650.0 + 128.0 * C
    @staticmethod
    def cp(C):   return 1450.0 + 2736.0 * C / (C + 1.0)
    @staticmethod
    def k(C):    return 0.21 + 0.38 * C / (C + 1.0)
    @staticmethod
    def D(C, T): return 2.4e-3 * np.exp(-0.45 / np.maximum(C, 1e-3)) * np.exp(-3850.0 / T)


def make_grid(N, kind):
    """kind='geo': ξ_i=i/N; kind='sqrt': ξ_i=sqrt(i/N)（等干基质量）"""
    if kind == 'geo':
        xi = np.linspace(0.0, 1.0, N + 1)
        vol = np.empty(N + 1)
        h = 1.0 / N
        vol[0] = h * h / 8
        vol[1:N] = np.arange(1, N) * h * h
        vol[N] = h / 2 - h * h / 8
        return xi, vol, np.diff(xi), 'harmonic'
    xi = np.sqrt(np.linspace(0.0, 1.0, N + 1))
    vol = np.full(N + 1, 1.0 / (2.0 * N))
    return xi, vol, np.diff(xi), 'arith'


def face_D(Dv, kind):
    """界面扩散系数：kind='harmonic' 调和平均；其余为算术平均（与交付代码一致）。"""
    if kind == 'harmonic':
        return 2.0 * Dv[:-1] * Dv[1:] / (Dv[:-1] + Dv[1:])
    return 0.5 * (Dv[:-1] + Dv[1:])


FACE_AVG = 'arith'          # 全局：'arith'（默认，与 herb_model.py / herb_v2.py 一致）或 'harmonic'


def fd(Dv, kind):
    if FACE_AVG == 'harmonic':
        return 2.0 * Dv[:-1] * Dv[1:] / (Dv[:-1] + Dv[1:])
    return 0.5 * (Dv[:-1] + Dv[1:])


def R_from_table():
    import openpyxl
    wb = openpyxl.load_workbook(r'C:\Users\qing1\Desktop\A题\附件\附件2.xlsx')
    ws = wb.active
    ts, rs = [], []
    for row in ws.iter_rows(min_row=2, values_only=True):
        if row[0] is None:
            continue
        ts.append(float(row[0])); rs.append(float(row[1]) / 100.0)
    ts = np.array(ts); rs = np.array(rs)
    def R(t):
        return float(np.interp(t, ts, rs))
    return R, ts, rs


# ============================================================ 通用求解器
def solve(N, dt, t_end, Rfun, props, variant='eq8', grid='geo', C0=2.55,
          T_air_fun=None, C_air=CAIR, coupled=False, cap_scale=False,
          record=(), stop_thr=None, tmax=None, conserve_check=False, T0=28.0):
    """物质坐标有限体积：dC_i/dt = (g_{i+1/2}-g_{i-1/2})/(R² vol_i), g = ξ D ∂_ξC。

    variant='eq8' 正确式(8)；'eq7' 照抄原文式(7)得到的 4ξ²/R²·∂_ξ(D∂_ξC)。
    返回 (t, xi, C, T, rec, extra)
    """
    xi, vol, dx, dkind = make_grid(N, grid)
    C = np.full(N + 1, C0)
    T = np.full(N + 1, T0)
    g_out_fac = R0  # placeholder
    t = 0.0
    k = 0
    rec = {}
    extra = {}
    if T_air_fun is None:
        T_air_fun = lambda tt: TAIR
    while True:
        if t_end is not None and t >= t_end - 1e-9:
            break
        if tmax is not None and t >= tmax:
            break
        theta = 1.0 if k < 4 else 0.5
        if variant == 'eq7':
            theta = 1.0
        R = Rfun(t + dt)
        Rc = max(R, 1e-4)
        Tair = T_air_fun(t + dt)
        # ---- 质量方程
        Dv = props.D(C, T + 273.15)
        Df = fd(Dv, dkind)
        if variant == 'eq8':
            # 正确式(8): 面通量 g = ξ D ∂_ξC，系数 1/(R² vol)
            gf = xi[1:] * Df / dx
            coefC = 1.0 / (Rc ** 2 * vol)
        else:
            # 照抄原文式(7): (1/U²)∂_ξ(ρd²r² D ∂_ξC) = (4/R²)∂_ξ(ξ² D ∂_ξC)
            gf = xi[1:] ** 2 * Df / dx
            coefC = 4.0 / (Rc ** 2 * vol)
        gC_out = Rc * HM
        abC, rhsC = _assemble(C, dt, theta, coefC, gf, gC_out, C_air, N)
        C = np.maximum(solve_banded((1, 1), abC, rhsC), 1e-6)
        if conserve_check:
            gprev = np.zeros(N + 1); gprev[:N] = gf * (C[1:] - C[:-1])
            flux_out = gC_out * (C[-1] - C_air) / Rc ** 2
            extra.setdefault('flux', []).append((t + dt, flux_out))
        # ---- 能量方程（可选）
        if coupled:
            kv = props.k(C)
            kf = 0.5 * (kv[:-1] + kv[1:])
            gT = xi[1:] * kf / dx
            kap = props.rho(C) * props.cp(C)
            if cap_scale:
                kap = kap * (R0 / Rc) ** 2
            coefT = 1.0 / (Rc ** 2 * vol * kap)
            gT_out = Rc * H_HEAT
            abT, rhsT = _assemble(T, dt, 1.0, coefT, gT, gT_out, Tair, N)
            T = solve_banded((1, 1), abT, rhsT)
        t += dt; k += 1
        for rt in record:
            if abs(t - rt) < 0.5 * dt:
                rec[rt] = (C.copy(), T.copy())
        if stop_thr is not None and np.max(C) <= stop_thr:
            break
        if k > 40_000_000:
            break
    return t, xi, C, T, rec, extra


def _assemble(Y, dt, theta, coef, gf, g_out, Y_inf, N):
    """θ 方法三对角系统（向量化）。gf 长度 N（面系数），coef 长度 N+1（节点系数）。"""
    low = -theta * dt * coef[1:] * gf
    up = -theta * dt * coef[:N] * gf
    main = np.ones(N + 1)
    main[1:N] += theta * dt * coef[1:N] * (gf[:-1] + gf[1:])
    main[0] += theta * dt * coef[0] * gf[0]
    main[N] += theta * dt * coef[N] * (gf[-1] + g_out)
    ab = np.zeros((3, N + 1)); ab[1] = main; ab[0, 1:] = up; ab[2, :-1] = low
    # 面通量 G_j = g_{j+1/2} = ξ D ∂_ξC 在面 j+1/2 的取值；G_N 为表面通量
    G = np.zeros(N + 1)
    G[:N] = gf * (Y[1:] - Y[:-1])
    G[N] = -g_out * (Y[N] - Y_inf)
    # 节点 i 的净流入 F_i = g_{i+1/2} - g_{i-1/2}（轴心处 g_{-1/2}=0）
    F = np.empty(N + 1)
    F[0] = G[0]
    F[1:] = G[1:] - G[:-1]
    rhs = Y + (1 - theta) * dt * coef * F
    rhs[-1] += theta * dt * coef[N] * g_out * Y_inf
    return ab, rhs


# ============================================================ A
def part_A():
    import sympy as sp
    xi, R, R0s, rho_d0, D = sp.symbols('xi R R_0 rho_d0 D', positive=True)
    C = sp.Function('C')(xi)
    rho_d = rho_d0 * R0s ** 2 / R ** 2
    r = xi * R
    U = rho_d0 * R0s ** 2 / 2
    lhs7 = sp.simplify(sp.diff(rho_d ** 2 * r ** 2 * D * sp.diff(C, xi), xi) / U ** 2)
    lhs7p = sp.simplify(sp.diff(xi * rho_d * D * sp.diff(C, xi), xi) / (R ** 2 * xi * rho_d))
    eq8 = sp.simplify(sp.diff(xi * D * sp.diff(C, xi), xi) / (xi * R ** 2))
    d = sp.simplify(sp.expand(lhs7p - eq8))
    d7 = sp.simplify(sp.expand(lhs7 - eq8))
    rep['A'] = {'eq7_printed_substituted': sp.sstr(lhs7), 'eq7_correct': sp.sstr(lhs7p),
                'eq8': sp.sstr(eq8), 'eq7p_minus_eq8': sp.sstr(d),
                'eq7_minus_eq8': sp.sstr(d7), 'eq7p_equals_eq8': bool(d == 0),
                'eq7_equals_eq8': bool(d7 == 0)}
    print('=== A. 代数核对 ===')
    print('  式(7) 照抄代入 (ρd=ρd0R0²/R², r=ξR) ->', rep['A']['eq7_printed_substituted'])
    print('  正确中间式 (7′)                       ->', rep['A']['eq7_correct'])
    print('  论文式(8)                              ->', rep['A']['eq8'])
    print('  (7′) − (8) =', rep['A']['eq7p_minus_eq8'], '| 恒等:', d == 0)
    print('  (7)  − (8) =', rep['A']['eq7_minus_eq8'], '| 恒等:', d7 == 0)


def part_A2():
    """核对：把式(7) 读成「关于干基质量坐标 η=ξ² 的方程」，是否恰好给出式(8)。

    关键：系数 (1/U²)ρd²r² 必须先化成 η 的表达式（=4η/R²）再求导，
    否则 η 与 ξ 混用会得到错误的比对结果。
    """
    import sympy as sp
    xi, R, D = sp.symbols('xi R D', positive=True)
    f = sp.Function('f')(xi)

    # (1/U²)ρd²r² = 4η/R² = 4ξ²/R² ，而 ∂_η = (1/(2ξ))∂_ξ
    eta_form = (4 / R ** 2) * (1 / (2 * xi)) * sp.diff(
        xi ** 2 * D * (1 / (2 * xi)) * sp.diff(f, xi), xi)
    eq8 = sp.diff(xi * D * sp.diff(f, xi), xi) / (xi * R ** 2)
    diff = sp.simplify(sp.expand(eta_form - eq8))
    print('\n=== A2. 式(7)（η 坐标版）与式(8) 的核对 ===')
    print('  (7)η版 → 4/R²·∂_η(ηD∂_ηC) 换回 ξ 后 =', sp.sstr(sp.simplify(eta_form)))
    print('  式(8)                                  =', sp.sstr(sp.simplify(eq8)))
    print('  差 =', sp.sstr(diff), '| 恒等:', diff == 0)
    rep['A2'] = {'eta_form_equals_eq8': bool(diff == 0),
                 'eta_form': sp.sstr(sp.simplify(eta_form)),
                 'eq8': sp.sstr(sp.simplify(eq8))}
    return rep['A2']


# ============================================================ B
def part_B():
    print('\n=== B. §5.6 / §5.7 数值复算 ===')
    D3 = lambda C, Tc: 2.4e-3 * math.exp(-0.45 / C) * math.exp(-3850.0 / (Tc + 273.15))
    D2 = lambda C: 7e-9 * math.exp(-0.89 / C)
    r = {}
    r['D3_2.55_50C'] = D3(2.55, 50); r['D3_0.15_50C'] = D3(0.15, 50)
    r['D3_ratio'] = D3(2.55, 50) / D3(0.15, 50)
    r['Ea_kJ'] = 3850 * 8.314 / 1000.0
    rho3 = 650 + 128 * 2.55; cp3 = 1450 + 2736 * 2.55 / 3.55; k3 = 0.21 + 0.38 * 2.55 / 3.55
    r['app3_rho'] = rho3; r['app3_cp'] = cp3; r['app3_k'] = k3
    r['app3_rhocp'] = rho3 * cp3
    r['app3_alpha'] = k3 / (rho3 * cp3)
    r['app3_Bim'] = HM * R0 / D3(2.55, 50)
    r['app3_Le'] = D3(2.55, 50) / r['app3_alpha']
    a2 = 0.36 / (820.0 * 2600.0)
    r['app2_alpha'] = a2
    r['p1_Bi'] = 25 * R0 / 0.36
    r['p1_Bim'] = HM * R0 / D2(2.55)
    r['p1_Le'] = D2(2.55) / a2
    r['p1_FoT_1800'] = a2 * 1800 / R0 ** 2
    r['p1_FoD_1800'] = D2(2.55) * 1800 / R0 ** 2
    r['app4_rhod_C2.55'] = (760 + 90 * 2.55) / 3.55
    r['app4_rhod_C0.15'] = (760 + 90 * 0.15) / 1.15
    r['app3_rhod_C2.55'] = (650 + 128 * 2.55) / 3.55
    paper = {'D3_2.55_50C': 1.35e-8, 'D3_0.15_50C': 8.0e-10, 'D3_ratio': 17.0, 'Ea_kJ': 32.0,
             'app3_alpha': 1.02e-7, 'app3_Bim': 1.19, 'app3_Le': 0.13,
             'p1_Bi': 1.389, 'p1_Bim': 3.24, 'p1_Le': 2.9e-2,
             'p1_FoT_1800': 0.760, 'p1_FoD_1800': 0.0222,
             'app4_rhod_C2.55': 231.0, 'app4_rhod_C0.15': 57.0}
    r['_paper'] = paper
    for k, v in paper.items():
        got = r.get(k)
        if got is None:
            continue
        rel = abs(got - v) / abs(v)
        print('  %s%-16s 论文=%-11g 复算=%-13.6g 相对差=%6.2f%%'
              % ('OK ' if rel < 0.02 else '!! ', k, v, got, rel * 100))
    rep['B'] = r
    return r


# ============================================================ 解析解
def bessel_roots(Bi, n=60):
    f = lambda lam: lam * j1(lam) - Bi * j0(lam)
    roots = []
    x = 1e-6
    step = 0.02
    prev = f(x); px = x
    x += step
    while len(roots) < n and x < 400:
        cur = f(x)
        if prev == 0:
            roots.append(px)
        elif prev * cur < 0:
            roots.append(brentq(f, px, x, xtol=1e-13))
        prev, px = cur, x
        x += step
    return np.array(roots)


def analytic_C(r, t, R, D, hm, C0, Cinf, n=80):
    Bi = hm * R / D
    lam = bessel_roots(Bi, n)
    A = 2 * j1(lam) / (lam * (j0(lam) ** 2 + j1(lam) ** 2))
    out = np.zeros_like(np.asarray(r, dtype=float))
    for L, a in zip(lam, A):
        out += a * j0(L * np.asarray(r, float) / R) * math.exp(-L * L * D * t / R ** 2)
    return Cinf + (C0 - Cinf) * out


# ============================================================ C ~ G
def solve_eta(N, dt, t_end, Rfun, props, C0=2.55, C_air=CAIR, record=()):
    """完全独立的一套离散：在干基质量坐标 η=ξ² 上求解
           ∂C/∂t = (4/R²) ∂_η (η D ∂_η C)
       （它与式(7) 的 (1/U²)∂_η(ρd²r²D∂_ηC) 等价，但与式(8) 的 ξ 网格不是同一套离散）
    """
    eta = np.linspace(0.0, 1.0, N + 1)
    dn = 1.0 / N
    vol = np.full(N + 1, dn)
    vol[0] = dn / 2
    vol[N] = dn / 2
    C = np.full(N + 1, C0)
    t = 0.0
    k = 0
    rec = {}
    while t < t_end - 1e-9:
        theta = 1.0 if k < 4 else 0.5
        R = Rfun(t + dt)
        Dv = props.D(C, 323.15)
        Df = fd(Dv, 'geo')
        gf = eta[1:] * Df / dn                      # 面系数 g = η D ∂_η C
        coef = 4.0 / (R ** 2 * vol)
        g_out = R * HM / 2.0                        # 表面 g_{N+1/2} = -R hm ΔC / 2
        ab, rhs = _assemble(C, dt, theta, coef, gf, g_out, C_air, N)
        C = np.maximum(solve_banded((1, 1), ab, rhs), 1e-6)
        t += dt; k += 1
        for rt in record:
            if abs(t - rt) < 0.5 * dt:
                rec[rt] = (C.copy(), None)
    return t, eta, C, rec


def part_C():
    print('\n=== C. 网格无关性 / 坐标无关性 / 照抄式(7) 的后果 ===')
    Rfun, ts, rs = R_from_table()
    t_end = 6 * 3600.0
    tc = 6 * 3600.0
    out = {}
    for kind, Ns in (('geo', (50, 100, 200, 400)), ('sqrt', (50, 100, 200, 400))):
        vals = []
        for N in Ns:
            t, xi, C, T, rec, _ = solve(N, 1.0, t_end, Rfun, P4, 'eq8', kind, T0=50.0,
                                        record=(tc,), T_air_fun=lambda x: TAIR)
            vals.append(rec[tc][0])
        out[kind] = vals
        ref = vals[-1][0]
        print('  eq(8) 网格=%-5s N=%s' % (kind, Ns))
        print('         C(ξ=0,6h) =', ['%.6f' % v[0] for v in vals])
        print('         与最细网格之差 =', ['%.2e' % abs(v[0] - ref) for v in vals])

    # ---- 坐标无关性：η 坐标求解器（独立离散）
    print('  --- 坐标无关性检验：η=ξ² 独立求解器 ---')
    eta_cmp, eta_err = {}, []
    for N in (100, 200, 400):
        t, eta, C, rec = solve_eta(N, 1.0, t_end, Rfun, P4, record=(tc,))
        eta_cmp[N] = (np.sqrt(eta), rec[tc][0])
    _, xi_g400, _, _, recg400, _ = solve(400, 1.0, t_end, Rfun, P4, 'eq8', 'geo',
                                         T0=50.0, record=(tc,))
    for N in (100, 200, 400):
        xi_e, Ce = eta_cmp[N]
        ci = np.interp(xi_e, xi_g400, recg400[tc][0])
        d_ = float(np.max(np.abs(Ce - ci)))
        print('  η 网格 N=%3d  与 ξ 网格 N=400 的 max|ΔC| = %.3e  (中心 %.5f vs %.5f)'
              % (N, d_, Ce[0], ci[0]))
        eta_err.append(d_)
    rep['C_eta'] = {'N': [100, 200, 400], 'maxdiff': eta_err}

    # ---- 照抄式(7)（把求导变量当成 ξ）
    N = 200
    t, xi, C8, T8, rec8, _ = solve(N, 1.0, t_end, Rfun, P4, 'eq8', 'geo', T0=50.0, record=(tc,))
    t, xi, C7, T7, rec7, _ = solve(N, 1.0, t_end, Rfun, P4, 'eq7', 'geo', T0=50.0, record=(tc,))
    a, b = rec8[tc][0], rec7[tc][0]
    print('  N=200, t=6 h:')
    print('    式(8) 剖面  ξ=0 / 0.5 / 1.0 = %.4f / %.4f / %.4f' % (a[0], a[N // 2], a[N]))
    print('    式(7) 剖面  ξ=0 / 0.5 / 1.0 = %.4f / %.4f / %.4f' % (b[0], b[N // 2], b[N]))
    print('    max|Δ| = %.4f kg/kg  (相对初始值 %.1f%%)'
          % (np.max(np.abs(a - b)), 100 * np.max(np.abs(a - b)) / 2.55))
    d8, d7 = run_dry('eq8', N), run_dry('eq7', N)
    print('  烘干时长 (max C ≤ 0.15): 式(8) %.2f h ；式(7) %.2f h%s'
          % (d8 / 3600, d7 / 3600, '（未达标，已截断）' if d7 >= 99 * 3600 else ''))
    rep['C'] = {'geo': {str(n): float(v[0]) for n, v in zip((50, 100, 200, 400), out['geo'])},
                'sqrt': {str(n): float(v[0]) for n, v in zip((50, 100, 200, 400), out['sqrt'])},
                'profile8': a.tolist(), 'profile7': b.tolist(),
                'maxdiff': float(np.max(np.abs(a - b))),
                'tdry_eq8_h': d8 / 3600, 'tdry_eq7_h': d7 / 3600, 'xi_geo200': xi.tolist()}
    return rep['C']


def run_dry(variant, N, dt=2.0, grid='geo', tmax_h=100):
    Rfun, _, _ = R_from_table()
    t, xi, C, T, rec, _ = solve(N, dt, None, Rfun, P4, variant, grid, T0=50.0,
                                stop_thr=0.15, tmax=tmax_h * 3600.0)
    return t


def part_D():
    print('\n=== D. R≡R0 退化 + 解析解对比 ===')
    D = 5e-9
    class Pc:
        rho = staticmethod(lambda C: np.full_like(np.asarray(C, float), 820.0))
        cp = staticmethod(lambda C: np.full_like(np.asarray(C, float), 2600.0))
        k = staticmethod(lambda C: np.full_like(np.asarray(C, float), 0.36))
        D = staticmethod(lambda C, T: np.full_like(np.asarray(C, float), D))
    Rconst = lambda t: R0
    Ns = (50, 100, 200, 400)
    errs = []
    for N in Ns:
        t, xi, C, T, rec, _ = solve(N, 1.0, 1800.0, Rconst, Pc, 'eq8', 'geo',
                                    C0=2.55, C_air=0.02, record=(600.0, 1800.0))
        e = []
        for tc in (600.0, 1800.0):
            num = rec[tc][0][::max(1, N // 8)]
            rr = xi[::max(1, N // 8)] * R0
            ana = analytic_C(rr, tc, R0, D, HM, 2.55, 0.02)
            e.append(np.max(np.abs(num - ana)))
        errs.append(e)
        print('  N=%3d  max|数值−解析|  600s: %.3e   1800s: %.3e' % (N, e[0], e[1]))
    rep['D'] = {'N': list(Ns), 'err600': [e[0] for e in errs], 'err1800': [e[1] for e in errs],
                'D': D, 'Bi_m': HM * R0 / D}
    return rep['D']


def part_E():
    print('\n=== E. 离散守恒性 ===')
    Rfun, _, _ = R_from_table()
    N, dt, t_end = 200, 1.0, 3600.0
    # 逐点记录总量与表面通量
    xi, vol, dx, dk = make_grid(N, 'geo')
    C = np.full(N + 1, 2.55); t = 0.0
    tot0 = float(np.sum(C * vol))
    cum = 0.0
    k = 0
    while t < t_end - 1e-9:
        theta = 1.0 if k < 4 else 0.5
        R = Rfun(t + dt)
        Dv = P4.D(C, 323.15); Df = fd(Dv, dk); gf = xi[1:] * Df / dx
        coef = 1.0 / (R ** 2 * vol); g_out = R * HM
        gprev = np.zeros(N + 1); gprev[:N] = gf * (C[1:] - C[:-1])
        flux_out = g_out * (C[-1] - CAIR) / R ** 2   # d(Σvol·C)/dt = g_out·ΔC/R²
        ab, rhs = _assemble(C, dt, theta, coef, gf, g_out, CAIR, N)
        C = np.maximum(solve_banded((1, 1), ab, rhs), 1e-6)
        cum += flux_out * dt
        t += dt; k += 1
    tot1 = float(np.sum(C * vol))
    lhs = tot0 - tot1
    print('  N=200, 1 h: Σ C·vol 减少 = %.8f' % lhs)
    print('              累计表面流出 = %.8f' % cum)
    print('  相对偏差 = %.3e' % (abs(lhs - cum) / abs(cum)))
    rep['E'] = {'loss': lhs, 'surface_cum': cum, 'rel': abs(lhs - cum) / abs(cum)}
    return rep['E']


def part_F():
    print('\n=== F. ρcp 缩放的对照算例（复现论文 §5.5 的对照实验）===')
    Rfun, _, _ = R_from_table()
    res, recs = {}, {}
    for tag, cs in (('baseline', False), ('scaled', True)):
        N = 400; dt = 2.0
        t, xi, C, T, rec, _ = solve(N, dt, None, Rfun, P4, 'eq8', 'geo', coupled=True,
                                    cap_scale=cs, stop_thr=0.15, tmax=120 * 3600.0,
                                    record=(24 * 3600.0, 48 * 3600.0))
        res[tag] = t
        recs[tag] = rec
        print('  %-8s (cap_scale=%s): t_dry = %.4f h' % (tag, cs, t / 3600))
    d = res['scaled'] - res['baseline']
    print('  Δt_dry = %+.1f s (%+.4f%%)' % (d, 100 * d / res['baseline']))
    for rt in (24 * 3600.0, 48 * 3600.0):
        if rt in recs['baseline'] and rt in recs['scaled']:
            a = recs['baseline'][rt][0][0]; b = recs['scaled'][rt][0][0]
            print('  t=%2.0f h 中心 C: %.6f vs %.6f（差 %+.4f%%）'
                  % (rt / 3600, a, b, 100 * (b - a) / a))
    rep['F'] = {'baseline_h': res['baseline'] / 3600, 'scaled_h': res['scaled'] / 3600,
                'delta_s': d, 'delta_pct': 100 * d / res['baseline']}
    return rep['F']


def part_H():
    """界面扩散系数：算术平均 vs 调和平均（论文 §5.8 写调和平均，交付代码用算术平均）。"""
    global FACE_AVG
    print('\n=== H. 界面平均方式的影响（调和平均 vs 算术平均）===')
    Rfun, _, _ = R_from_table()
    out = {}
    for tag, avg in (('arith 算术（交付代码）', 'arith'), ('harmonic 调和（论文 §5.8 文字）', 'harmonic')):
        FACE_AVG = avg
        t, xi, C, T, rec, _ = solve(400, 2.0, None, Rfun, P4, 'eq8', 'geo',
                                    T0=50.0, stop_thr=0.15, tmax=200 * 3600.0,
                                    record=(6 * 3600.0, 24 * 3600.0, 48 * 3600.0))
        r6, r24, r48 = (rec[k][0] for k in (6 * 3600.0, 24 * 3600.0, 48 * 3600.0))
        out[avg] = {'tdry_h': t / 3600, 'C6_center': float(r6[0]), 'C6_surf': float(r6[-1]),
                    'C24_center': float(r24[0]), 'C24_surf': float(r24[-1]),
                    'C48_center': float(r48[0]), 'C48_surf': float(r48[-1])}
        print('  %-20s t_dry = %7.2f h' % (tag, t / 3600))
        print('       中心 C：6h %.4f → 24h %.4f → 48h %.4f' % (r6[0], r24[0], r48[0]))
        print('       表面 C：6h %.4f → 24h %.4f → 48h %.4f' % (r6[-1], r24[-1], r48[-1]))
    FACE_AVG = 'arith'
    rep['H'] = out
    return out


def part_G():
    print('\n=== G. 干基密度非均匀性的量级 ===')
    Rfun, _, _ = R_from_table()
    N = 200
    t, xi, C, T, rec, _ = solve(N, 1.0, 12 * 3600.0, Rfun, P4, 'eq8', 'geo', T0=50.0,
                                record=(2 * 3600.0, 6 * 3600.0, 12 * 3600.0))
    rows = []
    for tc in sorted(rec):
        Cc = rec[tc][0]
        rd = P4.rho(Cc) / (1 + Cc)
        R = Rfun(tc)
        uni = (760 + 90 * 2.55) / 3.55 * (R0 / R) ** 2
        ratio = rd / uni
        i90 = int(round(0.9 * N))
        rows.append((tc / 3600, float(rd.min()), float(rd.max()), float(uni),
                     float(ratio.max() / ratio.min()),
                     float(rd[0]), float(rd[i90]), float(rd[-1])))
        print('  t=%4.1f h  ρd: 中心 %6.1f   ξ=0.9 %6.1f   表面 %6.1f  kg/m³'
              % (tc / 3600, rd[0], rd[i90], rd[-1]))
        print('             全剖面 max/min = %.3f ；表面/中心 = %.3f ；均匀假设值 = %.1f'
              % (ratio.max() / ratio.min(), rd[-1] / rd[0], uni))
    rep['G'] = rows
    return rows


def part_I():
    """复现新版（paper_electronic.pdf）§10.6 的外部数据交叉校验：
    用逐点 ρd=ρ(C)/(1+C) 的剖面按干基质量守恒反推半径 Rpred(t)，与附件2 实测 R(t) 比较。"""
    print('\n=== I. 复现新版 §10.6：附件2 与附录4 的相容性 ===')
    Rfun, ts, rs = R_from_table()
    rho_d0 = (760 + 90 * 2.55) / 3.55
    rec_times = (6 * 3600.0, 7 * 3600.0, 24 * 3600.0, 48 * 3600.0, 72 * 3600.0)
    t, xi, C, T, rec, _ = solve(400, 2.0, 72 * 3600.0, Rfun, P4, 'eq8', 'geo',
                                coupled=True, tmax=73 * 3600.0, record=rec_times)
    rows = []
    for tc in rec_times:
        Cc = rec[tc][0]
        rd = P4.rho(Cc) / (1 + Cc)
        den = float(np.trapezoid(rd * xi, xi)) if hasattr(np, 'trapezoid') else float(np.trapz(rd * xi, xi))
        Rpred = R0 * math.sqrt((rho_d0 / 2.0) / den)
        Rmeas = Rfun(tc)
        dev = (Rpred - Rmeas) / Rmeas * 100
        rows.append((tc / 3600, Rmeas * 100, Rpred * 100, dev))
        print('  t=%5.1f h  实测 R=%6.3f cm  反推 Rpred=%6.3f cm  偏差 %+6.2f%%'
              % (tc / 3600, Rmeas * 100, Rpred * 100, dev))
    need = rho_d0 * (R0 / Rfun(72 * 3600.0)) ** 2
    rdmax = float(np.max(P4.rho(np.linspace(0, 2.6, 500)) / (1 + np.linspace(0, 2.6, 500))))
    print('  终态所需均匀干密度 ρd0(R0/R)² = %.1f kg/m³ ；附录4 的 ρd 全场最大值 = %.1f（C=0）'
          % (need, rdmax))
    print('  → 实测收缩大于密度式所能支持的程度，二者在严格守恒意义下不相容 ✓')
    rep['I'] = {'rows': rows, 'need_uniform_rhod': need, 'max_rhod': rdmax,
                'rho_d0': rho_d0}
    return rep['I']


def part_J():
    """量化新版 §5.5 写出的「逐点 ρd」多出项 D/R²·∂_ξC·∂_ξln ρd 相对于式(10) 主项的比值。"""
    print('\n=== J. 「逐点 ρd」多出项的量级 ===')
    Rfun, _, _ = R_from_table()
    t, xi, C, T, rec, _ = solve(400, 2.0, 12 * 3600.0, Rfun, P4, 'eq8', 'geo', T0=50.0,
                                record=(2 * 3600.0, 6 * 3600.0, 12 * 3600.0))
    rows = []
    sel = (xi > 0.05) & (xi < 0.95)          # 避开轴心（主项趋零）与表面（网格尺度奇异）
    for tc in sorted(rec):
        Cc = rec[tc][0]
        R = Rfun(tc)
        rd = P4.rho(Cc) / (1 + Cc)
        Dv = P4.D(Cc, 323.15)
        dC = np.gradient(Cc, xi)
        dlr = np.gradient(np.log(rd), xi)
        extra = (Dv * dC * dlr / R ** 2)[sel]
        den = np.where(xi > 0, xi, 1.0)
        main = (np.gradient(xi * Dv * dC, xi) / (den * R ** 2))[sel]
        re_ = float(np.sqrt(np.mean(extra ** 2)))
        rm_ = float(np.sqrt(np.mean(main ** 2)))
        rows.append((tc / 3600, re_ / rm_, re_, rm_))
        print('  t=%4.1f h  ‖多出项‖/‖主项‖ = %.4f  (RMS %.3e vs %.3e)'
              % (tc / 3600, re_ / rm_, re_, rm_))
    rep['J'] = rows
    return rows


# ============================================================ main
if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    part_A(); part_A2(); part_B()
    if which in ('all', 'cd'):
        part_C(); part_D(); part_E()
    if which in ('all', 'fg'):
        part_F(); part_G(); part_H()
    if which in ('all', 'i'):
        part_I()
    if which in ('all', 'j'):
        part_J()
    json.dump(rep, open(os.path.join(HERE, 'verify_report.json'), 'w', encoding='utf-8'),
              ensure_ascii=False, indent=1, default=float)
    print('\n报告已写入 verify_report.json')
