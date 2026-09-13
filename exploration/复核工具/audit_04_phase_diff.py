# -*- coding: utf-8 -*-
"""审计脚本 5：问题1（附录2 常物性）与问题2（附录3 变物性）在时间重叠段的差异量级。"""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import openpyxl

OUT = r'C:\Users\qing1\Desktop\A题\cumcm2026-a-herb-drying-main\outputs\result_files'


def read_sheet(path, sheet):
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    rows = list(wb[sheet].iter_rows(values_only=True))
    wb.close()
    hdr = list(rows[0])
    t = np.array([float(r[0]) for r in rows[1:]])
    ncol = len(hdr) - 1
    mat = np.full((len(rows) - 1, ncol), np.nan)
    for i, r in enumerate(rows[1:]):
        for j in range(ncol):
            v = r[j + 1] if j + 1 < len(r) else None
            if v is not None:
                mat[i, j] = float(v)
    return hdr, t, mat


for name in ['温度', '水分浓度']:
    h1, t1, A = read_sheet(os.path.join(OUT, 'result1.xlsx'), name)
    h2, t2, B = read_sheet(os.path.join(OUT, 'result2_3h.xlsx'), name)
    # 对齐 0..1800
    k = int(np.where(np.isclose(t2, 1800))[0][0]) + 1
    a = A[:k]
    b = B[:k]
    d = np.abs(a - b)
    print('== %s（0—1800 s 重叠段，问题1 附录2 vs 问题2 附录3）' % name)
    print('   最大绝对差 = %.4f，出现在 t=%d s、距中心 %.1f cm'
          % (np.nanmax(d), t1[int(np.nanargmax(d) // d.shape[1])],
             float(h1[1:][int(np.nanargmax(d) % d.shape[1])])))
    for tt in (600, 1200, 1800):
        i = int(np.where(np.isclose(t1, tt))[0][0])
        print('   t=%4d s: 问题1 中心=%.4f 表面=%.4f | 问题2 中心=%.4f 表面=%.4f'
              % (tt, a[i, 0], a[i, -1], b[i, 0], b[i, -1]))
    print('   中心点最大差 = %.4f，表面最大差 = %.4f'
          % (np.nanmax(np.abs(a[:, 0] - b[:, 0])), np.nanmax(np.abs(a[:, -1] - b[:, -1]))))
    print('   1800 s 时两阶段的中心差 = %.4f（问题1 %.4f vs 问题2 %.4f）'
          % (a[k - 1, 0] - b[k - 1, 0], a[k - 1, 0], b[k - 1, 0]))
