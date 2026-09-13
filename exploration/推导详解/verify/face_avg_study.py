# -*- coding: utf-8 -*-
"""判决性研究：界面扩散系数用调和平均还是算术平均？

内容
  Part 1  理论验证：单界面精确电导 vs 调和平均 vs 算术平均（对真实剖面逐面计算）
  Part 2  网格收敛：两种平均方式下 t_dry 随 N 的变化（N = 200/400/800/1600）
  Part 3  与论文报告值（50.79 h）的比对

严格按论文 §5.8 的格式复刻：
  节点 ξ_i = i·h, h = 1/N
  体积测度 vol_0 = h²/8, vol_i = ξ_i·h, vol_N = h/2 − h²/8
  面位置 ξ_{i+1/2} = (i+1/2)h，面通量 g = ξ_{i+1/2}·D_{i+1/2}·(C_{i+1}−C_i)/h
  表面 g_{N+1/2} = −R·h_m·(C_N − C_air)，轴心 g_{−1/2} = 0
  θ 法（θ=1/2 + 前 4 步 Rannacher），物性取上一步显式值
"""
import sys, io, os, json, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace',
                              line_buffering=True)
import numpy as np
from scipy.linalg import solve_banded
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
R0 = 0.02
HM = 8e-7
H_HEAT = 25.0
CAIR = 0.0196
TAIR = 50.0
rep = {}


# ---------------------------------------------------------------- 物性 / 数据
def rho4(C):  return 760.0 + 90.0 * C
def cp4(C):   return 1850.0 + 2150.0 * C / (C + 1.0)
def k4(C):    return 0.12 + 0.20 * C / (C + 1.0)
def D4(C, T): return 4.2e-4 * np.exp(-0.30 / np.maximum(C, 1e-3)) * np.exp(-3850.0 / T)


def R_from_table():
    wb = openpyxl.load_workbook(r'C:\Users\qing1\Desktop\A题\附件\附件2.xlsx')
    ws = wb.active
    ts, rs = [], []
    for row in ws.iter_rows(min_row=2, values_only=True):
        if row[0] is None:
            continue
        ts.append(float(row[0])); rs.append(float(row[1]) / 100.0)
    ts, rs = np.array(ts), np.array(rs)
    return (lambda t: float(np.interp(t, ts, rs))), ts, rs


def face_D(Dv, kind):
    if kind == 'harmonic':
        return 2.0 * Dv[:-1] * Dv[1:] / np.maximum(Dv[:-1] + Dv[1:], 1e-30)
    return 0.5 * (Dv[:-1] + Dv[1:])


# ---------------------------------------------------------------- 论文 §5.8 格式
def make_grid(N):
    h = 1.0 / N
    xi = np.linspace(0.0, 1.0, N + 1)
    vol = np.empty(N + 1)
    vol[0] = h * h / 8.0
    vol[1:N] = xi[1:N] * h
    vol[N] = h / 2.0 - h * h / 8.0
    xi_f = (np.arange(N) + 0.5) * h              # 面位置 ξ_{i+1/2}
    return xi, vol, xi_f, h


def theta_solve(Y, dt, theta, coef, gf, g_out, Y_inf, N):
    """θ 法三对角系统；coef_i = 1/(R²vol_i)，gf_i = ξ_{i+1/2}D_{i+1/2}/h。"""
    main = np.ones(N + 1)
    low = np.zeros(N); up = np.zeros(N)
    low[:] = -theta * dt * coef[1:] * gf
    up[:] = -theta * dt * coef[:N] * gf
    main[1:N] += theta * dt * coef[1:N] * (gf[:-1] + gf[1:])
    main[0] += theta * dt * coef[0] * gf[0]
    main[N] += theta * dt * coef[N] * (gf[-1] + g_out)
    ab = np.zeros((3, N + 1)); ab[1] = main; ab[0, 1:] = up; ab[2, :-1] = low
    G = np.zeros(N + 1)
    G[:N] = gf * (Y[1:] - Y[:-1])
    G[N] = -g_out * (Y[N] - Y_inf)
    F = np.empty(N + 1)
    F[0] = G[0]
    F[1:] = G[1:] - G[:-1]
    rhs = Y + (1 - theta) * dt * coef * F
    rhs[-1] += theta * dt * coef[N] * g_out * Y_inf
    return ab, rhs


def run(N, avg, dt=2.0, tmax_h=220.0, T0=50.0, record=(), stop=0.15, coupled=True):
    """按论文 §5.8 格式跑问题 4（含收缩）。返回 (t_dry 或截断时刻, xi, C, T, rec)。"""
    Rfun, _, _ = R_from_table()
    xi, vol, xi_f, h = make_grid(N)
    C = np.full(N + 1, 2.55)
    T = np.full(N + 1, T0)
    t = 0.0; k = 0; rec = {}
    tmax = tmax_h * 3600.0
    while t < tmax - 1e-9:
        theta = 1.0 if k < 4 else 0.5
        R = max(Rfun(t + dt), 1e-4)
        Dv = D4(C, T + 273.15)
        Df = face_D(Dv, avg)
        gf = xi_f * Df / h                      # ξ_{i+1/2} D_{i+1/2} / h
        coefC = 1.0 / (R ** 2 * vol)
        gC_out = R * HM
        ab, rhs = theta_solve(C, dt, theta, coefC, gf, gC_out, CAIR, N)
        C = np.maximum(solve_banded((1, 1), ab, rhs), 1e-6)
        if coupled:
            kv = k4(C)
            kf = 0.5 * (kv[:-1] + kv[1:])       # 导热界面系数（两种方案都用算术）
            gT = xi_f * kf / h
            kap = rho4(C) * cp4(C)
            coefT = 1.0 / (R ** 2 * vol * kap)
            abT, rhsT = theta_solve(T, dt, 1.0, coefT, gT, R * H_HEAT, TAIR, N)
            T = solve_banded((1, 1), abT, rhsT)
        t += dt; k += 1
        for rt in record:
            if abs(t - rt) < 0.5 * dt:
                rec[rt] = (C.copy(), T.copy())
        if stop is not None and np.max(C) <= stop:
            break
    return t, xi, C, T, rec


# ---------------------------------------------------------------- Part 1
def exact_face_conductance(C, xi, xi_f, h, T=323.15, nq=64):
    """逐面计算：精确电导 G_exact = 1/⟨1/D⟩，调和平均，算术平均。"""
    # 在细网格上重构 D(ξ)
    n = len(xi) - 1
    xs = np.linspace(0.0, 1.0, n * nq + 1)
    Cs = np.interp(xs, xi, C)
    Ds = D4(Cs, T)
    Ge = np.empty(n); Gh = np.empty(n); Ga = np.empty(n)
    for i in range(n):
        seg = Ds[i * nq:(i + 1) * nq + 1]
        Ge[i] = 1.0 / np.mean(1.0 / seg)
    Dv = D4(C, T)
    Gh = 2.0 * Dv[:-1] * Dv[1:] / np.maximum(Dv[:-1] + Dv[1:], 1e-30)
    Ga = 0.5 * (Dv[:-1] + Dv[1:])
    return Ge, Gh, Ga


def part1():
    print('=' * 78)
    print('Part 1  单界面精确电导 vs 调和平均 vs 算术平均')
    print('=' * 78)
    Rfun, _, _ = R_from_table()
    N = 400
    t, xi, C, T, rec = run(N, 'arith', dt=2.0, tmax_h=72.0,
                           record=(6 * 3600.0, 24 * 3600.0, 48 * 3600.0))
    h = 1.0 / N
    xi_f = (np.arange(N) + 0.5) * h
    rows = []
    print('%-8s %-14s %-12s %-12s %-12s %-11s' %
          ('t/h', '面（自表面数）', 'G_exact', 'G_调和', 'G_算术', '调和/算术'))
    for tc in sorted(rec):
        Cc = rec[tc][0]
        Ge, Gh, Ga = exact_face_conductance(Cc, xi, xi_f, h)
        for back in (1, 5, 20):                       # 用等价 ξ 索引定位近表面面
            j = N - back
            rows.append((tc / 3600, back, float(Ge[j]), float(Gh[j]), float(Ga[j])))
            print('%-8.1f 倒数第 %-7d %.4e   %.4e   %.4e   %.4f' %
                  (tc / 3600, back, Ge[j], Gh[j], Ga[j], Gh[j] / Ga[j]))
        # 全局相对偏差
        m = Ge > 1e-14
        eh = np.max(np.abs(Gh[m] - Ge[m]) / Ge[m])
        ea = np.max(np.abs(Ga[m] - Ge[m]) / Ge[m])
        print('        → 与精确值的最大相对偏差：调和 %.1f%% ，算术 %.1f%%' %
              (100 * eh, 100 * ea))
    rep['part1'] = rows
    return rows


# ---------------------------------------------------------------- Part 2
def part2():
    print()
    print('=' * 78)
    print('Part 2  网格收敛：t_dry 随 N 的变化（两种平均方式）')
    print('=' * 78)
    rec_t = (24 * 3600.0, 48 * 3600.0, 72 * 3600.0)
    out = {}
    for avg in ('arith', 'harmonic'):
        out[avg] = []
        print('--- 界面系数：%s ---' % ('算术平均' if avg == 'arith' else '调和平均'))
        print('%-6s %-12s %-12s %-12s %-12s' %
              ('N', 't_dry/h', 'C(0,24h)', 'C(0,48h)', 'C(0,72h)'))
        for N in (200, 400, 800, 1600):
            dt = 2.0 if N <= 400 else 4.0
            t, xi, C, T, rec = run(N, avg, dt=dt, tmax_h=220.0, record=rec_t)

            def g(rt):
                v = rec.get(rt)
                return float(v[0][0]) if v is not None else float('nan')
            row = {'N': N, 't_dry_h': t / 3600, 'C24': g(rec_t[0]), 'C48': g(rec_t[1]),
                   'C72': g(rec_t[2]),
                   'surface24': (float(rec[rec_t[0]][0][-1]) if rec_t[0] in rec else float('nan'))}
            out[avg].append(row)
            print('%-6d %-12.2f %-12.6f %-12.6f %-12.6f' %
                  (N, row['t_dry_h'], row['C24'], row['C48'], row['C72']))
        print('    （t_dry 达到 220 h 上界表示该网格下中心未降到 0.15）')
    rep['part2'] = out
    return out


# ---------------------------------------------------------------- Part 3
def part3():
    """阻力量级分析：真实干壳厚度、单元阻力、整条半径的内部阻力。"""
    print()
    print('=' * 78)
    print('Part 3  阻力量级分析（谁在控制通量？）')
    print('=' * 78)
    Rfun, _, _ = R_from_table()
    N = 400
    t, xi, C, T, rec = run(N, 'arith', dt=2.0, tmax_h=72.0,
                           record=(6 * 3600.0, 24 * 3600.0, 48 * 3600.0))
    rows = []
    print('%-6s %-11s %-12s %-12s %-12s %-12s' %
          ('t/h', 'Cs', '干壳厚度/m', '干壳阻力', '内部阻力', '比值'))
    for tc in sorted(rec):
        Cc = rec[tc][0]
        R = Rfun(tc)
        Cs = Cc[-1]
        rhod = rho4(Cs) / (1 + Cs)
        Js = rhod * HM * (Cs - CAIR)                 # kg/(m²·s)
        dCdr = Js / (rhod * D4(np.array([Cs]), 323.15)[0])   # kg/kg per m
        delta = (Cs - CAIR) / dCdr                   # m
        D_skin = D4(np.array([(Cs + CAIR) / 2]), 323.15)[0]
        R_skin = delta / D_skin
        R_int = R / D4(np.array([Cc[0]]), 323.15)[0]
        rows.append((tc / 3600, float(Cs), float(delta), float(R_skin), float(R_int),
                     float(R_skin / R_int)))
        print('%-6.1f %-11.4f %-12.3e %-12.3e %-12.3e %-12.4f' %
              (tc / 3600, Cs, delta, R_skin, R_int, R_skin / R_int))

    # 单元阻力（两种平均）在 N=400 与 N=1600 下
    print()
    print('%-8s %-6s %-12s %-12s %-12s %-12s' %
          ('t/h', 'N', 'D_调和', 'D_算术', '单元阻力(调和)', '单元阻力(算术)'))
    for tc in sorted(rec):
        Cc = rec[tc][0]
        R = Rfun(tc)
        for NN in (400, 1600):
            xi2 = np.linspace(0, 1, NN + 1)
            C2 = np.interp(xi2, xi, Cc)
            Dv = D4(C2, 323.15)
            Dh = face_D(Dv, 'harmonic')[-1]
            Da = face_D(Dv, 'arith')[-1]
            h_phys = R / NN
            print('%-8.1f %-6d %-12.3e %-12.3e %-12.3e %-12.3e' %
                  (tc / 3600, NN, Dh, Da, h_phys / Dh, h_phys / Da))
    rep['part3'] = rows
    return rows


# ---------------------------------------------------------------- main
if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('all', '1'):
        part1()
    if which in ('all', '2'):
        part2()
    if which in ('all', '3'):
        part3()
    json.dump(rep, open(os.path.join(HERE, 'face_avg_study.json'), 'w', encoding='utf-8'),
              ensure_ascii=False, indent=1, default=float)
    print('\n结果已写入 face_avg_study.json')
