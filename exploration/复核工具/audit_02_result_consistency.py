# -*- coding: utf-8 -*-
"""审计脚本 2：结果文件内部一致性 + 合规性（阈值/单调性/网格/表格命中）。"""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import numpy as np
import openpyxl

ROOT = r'C:\Users\qing1\Desktop\A题'
OUT = os.path.join(ROOT, 'cumcm2026-a-herb-drying-main', 'outputs', 'result_files')


def load(path, sheet):
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    ws = wb[sheet] if sheet else wb.worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    hdr = rows[0]
    t = np.array([float(r[0]) for r in rows[1:]])
    data = np.array([[np.nan if v is None else float(v) for v in r[1:]] for r in rows[1:]])
    return hdr, t, data


print('=' * 78)
print('A. 表头 / 行数 / 时间栅格 合规性')
print('=' * 78)
specs = [('result1.xlsx', '温度', 0, 1800, 1),
         ('result1.xlsx', '水分浓度', 0, 1800, 1),
         ('result2.xlsx', '温度', 0, 205729, 1),
         ('result2.xlsx', '水分浓度', 0, 205729, 1),
         ('result3.xlsx', '水分浓度', 0, None, 60),
         ('result4.xlsx', '水分浓度', 0, None, 60),
         ('result4_local.xlsx', '水分浓度', 0, None, 60)]
cache = {}
for f, sh, t0, t1, dt in specs:
    key = (f, sh)
    if key not in cache:
        cache[key] = load(os.path.join(OUT, f), sh)
    hdr, t, data = cache[key]
    ncol = len(hdr) - 1
    d = np.diff(t)
    regular = np.allclose(d[:-1], dt, atol=1e-6)
    n_reg = int(np.sum(np.isclose(t % dt, 0, atol=1e-6)))
    print(f'-- {f}[{sh}]: {len(t)} 时间点, {ncol} 数据列, 表头={hdr[0]}')
    print(f'   列坐标: {hdr[1:]}')
    print(f'   t: {t[0]:.4f} .. {t[-1]:.4f}; 步长是否恒为 {dt}s: {regular};'
          f' 落在 {dt}s 整倍数上的点数: {n_reg}')

print()
print('=' * 78)
print('B. result3 / result4 达标性核查（是否真的 "各处 < 0.15"）')
print('=' * 78)
for f, sh in [('result3.xlsx', '水分浓度'), ('result4.xlsx', '水分浓度'),
              ('result4_local.xlsx', '水分浓度')]:
    hdr, t, data = cache[(f, sh)]
    # 排除最后的 "药材表面" 列做中心/全场判断，也单独看它
    ncol = len(hdr) - 1
    has_surf = (hdr[-1] == '药材表面')
    body = data[:, :ncol - 1] if has_surf else data
    surf = data[:, -1] if has_surf else None
    cmax = np.nanmax(body, axis=1)
    # 早停：首个满足 < 0.15 的时刻（含表面列一并判断）
    if surf is not None:
        allc = np.nanmax(np.concatenate([body, surf[:, None]], axis=1), axis=1)
    else:
        allc = cmax
    ok = allc < 0.15
    idx = int(np.argmax(ok)) if ok.any() else -1
    print(f'-- {f}: 最后一行 t={t[-1]:.4f}')
    print(f'   末行 max(C)（不含表面列）= {cmax[-1]:.8f}; 含表面列 = {allc[-1]:.8f}')
    if idx >= 0:
        print(f'   首个 "全部 <0.15" 的行: t={t[idx]:.4f} (第 {idx+2} 行), '
              f'max={allc[idx]:.8f}, 前一行为 {allc[idx-1]:.8f}')
        print(f'   达标行是否就是末行? {idx == len(t)-1}')
    # 单调性（中心点 & 全场）
    dec = np.all(np.diff(body[:, 0]) <= 1e-9)
    print(f'   中心点 C 单调不增: {dec}')
    # 是否存在 NaN（材料外留空）
    n_nan = int(np.isnan(body).sum())
    print(f'   材料外留空单元数: {n_nan} ({n_nan/body.size*100:.2f}%)')
    if surf is not None:
        print(f'   表面列与 2.0cm 列是否相同? '
              f'{np.allclose(surf, data[:, -2], equal_nan=True)}')

print()
print('=' * 78)
print('C. result2 的前 3h 与 result2_3h / result1 交叉一致')
print('=' * 78)
# result2_3h 与 result1 在 0..1800 应一致（同一模型不同阶段物性? 实际上问题1用附录2,
# 问题2用附录3 -> 二者应当不同；这里检查 result2 与 result2_3h 的重叠段）
h2, t2, d2 = cache[('result2.xlsx', '水分浓度')]
h23, t23, d23 = load(os.path.join(OUT, 'result2_3h.xlsx'), '水分浓度')
m = min(len(t2), len(t23))
print(f'result2 与 result2_3h 前 {m} 行时间一致: {np.allclose(t2[:m], t23[:m])}')
print(f'二者最大差: {np.nanmax(np.abs(d2[:m]-d23[:m])):.3e}')

# 温度
h2t, t2t, d2t = cache[('result2.xlsx', '温度')]
h23t, t23t, d23t = load(os.path.join(OUT, 'result2_3h.xlsx'), '温度')
print(f'温度表最大差: {np.nanmax(np.abs(d2t[:m]-d23t[:m])):.3e}')

print()
print('=' * 78)
print('D. 表 1/2/3/4/5/6 数值命中检查（与文档表格比对）')
print('=' * 78)
h1T, t1T, d1T = cache[('result1.xlsx', '温度')]
h1C, t1C, d1C = cache[('result1.xlsx', '水分浓度')]


def at(t, val, times, cols):
    i = int(np.argmin(np.abs(times - val)))
    assert abs(times[i] - val) < 1e-6, (times[i], val)
    j = int(np.argmin(np.abs(np.array(cols) - cols[0])))  # placeholder
    return i


cols_req = [0, 0.5, 1, 1.5, 2]
col_index = {c: list(h1T[1:]).index(c) for c in cols_req}
print('表1（温度，文档值 vs result1）')
doc_T1 = {100: [28.0000, 28.0003, 28.0038, 28.0319, 28.1795],
          300: [28.0406, 28.0634, 28.1513, 28.3673, 28.8486],
          600: [28.4533, 28.5361, 28.8039, 29.3155, 30.1641],
          900: [29.3243, 29.4583, 29.8757, 30.6164, 31.7295],
          1200: [30.5429, 30.7099, 31.2225, 32.1130, 33.4294],
          1500: [31.9960, 32.1870, 32.7667, 33.7476, 35.1211]}
for tt, vals in doc_T1.items():
    i = int(np.where(np.isclose(t1T, tt, atol=1e-6))[0][0])
    got = [d1T[i, col_index[c]] for c in cols_req]
    dif = max(abs(a - b) for a, b in zip(got, vals))
    print(f'   t={tt:5d} 文件={np.round(got,4)} 文档={vals} 最大差={dif:.2e}')

print('表2（水分浓度）')
doc_C1 = {100: [2.5500, 2.5500, 2.5500, 2.5500, 2.2471],
          300: [2.5500, 2.5500, 2.5500, 2.5492, 2.0509],
          600: [2.5500, 2.5500, 2.5500, 2.5353, 1.8770],
          900: [2.5500, 2.5500, 2.5497, 2.5045, 1.7547],
          1200: [2.5500, 2.5500, 2.5482, 2.4646, 1.6586],
          1500: [2.5500, 2.5499, 2.5445, 2.4206, 1.5787]}
for tt, vals in doc_C1.items():
    i = int(np.where(np.isclose(t1C, tt, atol=1e-6))[0][0])
    got = [d1C[i, col_index[c]] for c in cols_req]
    dif = max(abs(a - b) for a, b in zip(got, vals))
    print(f'   t={tt:5d} 文件={np.round(got,4)} 文档={vals} 最大差={dif:.2e}')

h3, t3, d3 = cache[('result3.xlsx', '水分浓度')]
ci3 = {c: list(h3[1:]).index(c) for c in cols_req}
print('表5（问题3）文档 vs result3')
doc_C5 = {6: [1.0156, 0.9868, 0.9006, 0.7546, 0.5335],
          12: [0.4548, 0.4418, 0.4019, 0.3289, 0.1640],
          18: [0.2979, 0.2906, 0.2679, 0.2247, 0.0842],
          24: [0.2371, 0.2321, 0.2161, 0.1851, 0.0663],
          30: [0.2053, 0.2013, 0.1887, 0.1637, 0.0597],
          36: [0.1854, 0.1821, 0.1714, 0.1500, 0.0565],
          42: [0.1716, 0.1687, 0.1593, 0.1404, 0.0547],
          48: [0.1614, 0.1588, 0.1503, 0.1331, 0.0536],
          54: [0.1535, 0.1511, 0.1433, 0.1274, 0.0528]}
for hh, vals in doc_C5.items():
    i = int(np.where(np.isclose(t3, hh * 3600, atol=1e-6))[0][0])
    got = [d3[i, ci3[c]] for c in cols_req]
    dif = max(abs(a - b) for a, b in zip(got, vals))
    print(f'   t={hh:2d}h 文件={np.round(got,4)} 文档={vals} 最大差={dif:.2e}')
i_end = len(t3) - 1
got = [d3[i_end, ci3[c]] for c in cols_req]
print(f'   末行 t={t3[i_end]:.4f}s ({t3[i_end]/3600:.4f}h) 文件={np.round(got,6)}')
