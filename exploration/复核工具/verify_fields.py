# -*- coding: utf-8 -*-
"""独立复核 C：用独立求解器算出的浓度/温度场 vs 题库 result3 / result2 的存储值。"""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import openpyxl
from indep_solver import FV, Rof, T0, C0

OUT = r'C:\Users\qing1\Desktop\A题\cumcm2026-a-herb-drying-main\outputs\result_files'
N = 800


def read_sheet(path, sheet):
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    rows = list(wb[sheet].iter_rows(values_only=True))
    wb.close()
    hdr = list(rows[0])
    t = np.array([float(r[0]) for r in rows[1:]])
    ncol = len(hdr) - 1
    mat = np.full((len(rows) - 1, ncol), np.nan)
    for i, r in enumerate(rows[1:]):
        for j in range(min(ncol, len(r) - 1)):
            v = r[j + 1]
            if v is not None:
                mat[i, j] = float(v)
    return hdr, t, mat


def colidx(hdr, c):
    return hdr.index(c) if c in hdr else hdr.index(float(c))


hdr, t_file, C_file = read_sheet(os.path.join(OUT, 'result3.xlsx'), '水分浓度')
cols = ['0', '0.5', '1', '1.5', '2']
ci = [hdr.index(int(float(c))) - 1 for c in cols]        # 表头第 0 列是标题行
print('表头长度 =', len(hdr), ' 数据列数 =', C_file.shape[1], ' 列索引 =', dict(zip(cols, ci)))
targets = [3600 * k for k in (6, 12, 18, 24, 36, 48, 54)]

s = FV('p3', False, N=N)
T = np.full(N, T0); C = np.full(N, C0); t = 0.0
dt = 60.0
watch = set(int(x) for x in targets)
snap = {}
while t < max(targets) - 1e-9:
    T, C = s.step(T, C, t, dt)
    t += dt
    if int(t) in watch:
        snap[int(t)] = (T.copy(), C.copy())
print('独立求解完成 N=%d，输出时刻 %s' % (N, sorted(snap)))
r_nodes = 0.02 * s.xi

print()
print('%-7s %-7s %-12s %-12s %-10s' % ('时刻/h', '距离/cm', 'result3存储', '独立复核', '相对差'))
worst = 0.0
for tt in targets:
    T_, C_ = snap[int(tt)]
    for c, j in zip(cols, ci):
        v_file = C_file[int(np.where(np.isclose(t_file, tt))[0][0]), j]
        v_ind = float(np.interp(float(c) / 100.0, r_nodes, C_))
        rel = abs(v_ind - v_file) / max(abs(v_file), 1e-12)
        worst = max(worst, rel)
        print('%-7d %-7s %-12.4f %-12.4f %-10.2e' % (tt / 3600, c, v_file, v_ind, rel))
print('浓度场最大相对差 = %.2e' % worst)

hdrT, tT, Tf = read_sheet(os.path.join(OUT, 'result2_3h.xlsx'), '温度')
ciT = [hdrT.index(int(float(c))) - 1 for c in cols]
print()
print('温度场对照（result2_3h，附录3 物性，0—3 h）')
for hh in (0.5, 1.0, 2.0, 3.0):
    tt = hh * 3600
    s2 = FV('p3', False, N=N)
    T2 = np.full(N, T0); C2 = np.full(N, C0); t2 = 0.0
    while t2 < tt - 1e-9:
        d = min(5.0, tt - t2)
        T2, C2 = s2.step(T2, C2, t2, d)
        t2 += d
    i = int(np.where(np.isclose(tT, tt))[0][0])
    vals_f = [Tf[i, j] for j in ciT]
    vals_i = [float(np.interp(float(c) / 100.0, r_nodes, T2)) for c in cols]
    dif = max(abs(a - b) for a, b in zip(vals_f, vals_i))
    print('  t=%.1fh 文件  =%s' % (hh, np.round(vals_f, 4)))
    print('         独立  =%s  最大绝对差=%.4f °C' % (np.round(vals_i, 4), dif))
