# -*- coding: utf-8 -*-
"""审计脚本 3：result4 空单元格几何模式 + 4 位小数舍入后的达标裕量。"""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import numpy as np
import openpyxl

OUT = r'C:\Users\qing1\Desktop\A题\cumcm2026-a-herb-drying-main\outputs\result_files'


def load(f, sh):
    wb = openpyxl.load_workbook(os.path.join(OUT, f), data_only=True, read_only=True)
    ws = wb[sh]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    return rows


print('=' * 74)
print('A. result4 空单元格是否恰为 "到中心距离 > R(t)" 的单元')
print('=' * 74)
rows = load('result4.xlsx', '水分浓度')
hdr = rows[0]
cols = np.array([float(c) for c in hdr[1:-1]])           # 0..2.0 cm
t = np.array([float(r[0]) for r in rows[1:]])
D = np.array([[np.nan if v is None else float(v) for v in r[1:-1]] for r in rows[1:]])
ii = np.isnan(D)
wrong = 0
for i in range(len(t)):
    for j, c in enumerate(cols):
        if ii[i, j]:
            wrong += 1
print(f'  空单元总数 = {ii.sum()}（应为 0）' if wrong else '  无空单元？')
# 与表面位置表对照
rowsR = load('result4.xlsx', '表面位置')
tR = np.array([float(r[0]) for r in rowsR[1:]])
RR = np.array([float(r[1]) for r in rowsR[1:]])
print(f'  表面位置表: {len(tR)} 行, t {tR[0]} .. {tR[-1]}, R {RR[0]} .. {RR[-1]} cm')
print(f'  时间轴一致: {np.allclose(t, tR)}')
bad = 0
for i in range(len(t)):
    mask = cols > RR[i] + 1e-9
    filled = ~ii[i]
    if np.any(filled & mask):
        bad += 1
        if bad <= 5:
            print(f'    !! t={t[i]:.1f} R={RR[i]:.4f} 但 {cols[filled & mask]} 列已填值')
print(f'  "R 之外仍填值" 的行数 = {bad}')
# 表面列与最后有效列
i0 = 0
print(f'  t=0 行 2.0cm 列值 = {D[0, -1]}, 表面列值 = {rows[1][-1]}')
i_end = len(t) - 1
print(f'  末行 t={t[i_end]:.4f} R={RR[i_end]:.4f}cm: 2.0cm列={D[i_end, -1]}, '
      f'1.2cm列={D[i_end, cols.tolist().index(1.2)]}, 表面列={rows[i_end+1][-1]}')

print()
print('=' * 74)
print('B. result3 末几行中心点数值（阈值 0.15 的舍入达标性）')
print('=' * 74)
rows3 = load('result3.xlsx', '水分浓度')
t3 = np.array([float(r[0]) for r in rows3[1:]])
C3 = np.array([[np.nan if v is None else float(v) for v in r[1:]] for r in rows3[1:]])
n = len(t3)
for k in range(n - 6, n):
    row = C3[k]
    print(f'  t={t3[k]:12.4f}s ({t3[k]/3600:8.4f}h)  C_center={row[0]:.8f} '
          f' C_max={np.nanmax(row):.8f}  C_surf={row[-1]:.6f}')
print(f'  末行是否 60s 整倍数: {t3[-1] % 60 == 0}')
print(f'  上一行(t={t3[-2]})与末行间隔: {t3[-1]-t3[-2]:.4f} s')
print(f'  倒数第二个 60s 整倍数时刻 t={t3[-2]:.1f}s 时 C_max={np.nanmax(C3[-2]):.8f}')

print()
print('=' * 74)
print('C. "严格达标" 判定：4 位小数舍入后是否仍 < 0.15')
print('=' * 74)
for tag, arr in [('result3 末行', C3[-1]), ('result4 末行', None)]:
    if arr is None:
        continue
    r = np.round(arr, 4)
    print(f'  {tag}: max(round(C,4)) = {r.max():.4f}; 是否 <0.15: {r.max() < 0.15}; '
          f'是否 <=0.15: {r.max() <= 0.15}')
rows4 = load('result4.xlsx', '水分浓度')
C4 = np.array([[np.nan if v is None else float(v) for v in r[1:]] for r in rows4[1:]])
r4 = np.round(np.nan_to_num(C4, nan=-1), 4)
print(f'  result4 末行: max(round(C,4)) = {r4[-1].max():.4f}; 是否 <0.15: {r4[-1].max() < 0.15}')
rows4l = load('result4_local.xlsx', '水分浓度')
C4l = np.array([[np.nan if v is None else float(v) for v in r[1:]] for r in rows4l[1:]])
r4l = np.round(np.nan_to_num(C4l, nan=-1), 4)
print(f'  result4_local 末行: max(round(C,4)) = {r4l[-1].max():.4f}; 是否 <0.15: {r4l[-1].max() < 0.15}')

print()
print('=' * 74)
print('D. 阈值穿越：整秒/整分栅格上真正 "< 0.15" 的首个时刻')
print('=' * 74)
ok3 = np.nanmax(C3, axis=1) < 0.15
print(f'  result3 栅格内首个 C_max<0.15 的行: idx={int(np.argmax(ok3))}, '
      f't={t3[int(np.argmax(ok3))]:.1f}s, 是否末行={int(np.argmax(ok3))==n-1}')
ok4 = np.nanmax(np.nan_to_num(C4, nan=-1), axis=1) < 0.15
print(f'  result4 栅格内首个 C_max<0.15 的行: idx={int(np.argmax(ok4))}, '
      f't={t4v:.1f}s' if False else f'  result4 首个达标行 idx={int(np.argmax(ok4))}, '
      f't={np.array([float(r[0]) for r in rows4[1:]])[int(np.argmax(ok4))]:.1f}s')
