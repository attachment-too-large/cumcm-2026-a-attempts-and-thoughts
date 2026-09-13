# -*- coding: utf-8 -*-
"""加密网格：把两种界面平均方式的外推极限钉死。"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import face_avg_study as S      # 该模块已自行配置好 stdout 编码

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
for avg in ('harmonic', 'arith'):
    out[avg] = []
    print('--- %s ---' % ('算术平均' if avg == 'arith' else '调和平均'))
    for N in (3200, 6400):
        t, xi, C, T, rec = S.run(N, avg, dt=8.0, tmax_h=220.0,
                                 record=(24 * 3600., 48 * 3600., 72 * 3600.))
        row = {'N': N, 't_dry_h': t / 3600.}
        for k, rt in (('C24', 24 * 3600.), ('C48', 48 * 3600.), ('C72', 72 * 3600.)):
            row[k] = float(rec[rt][0][0]) if rt in rec else float('nan')
        out[avg].append(row)
        print('N=%-6d t_dry=%-8.2f C24=%.6f C48=%.6f C72=%.6f'
              % (N, row['t_dry_h'], row['C24'], row['C48'], row['C72']))
json.dump(out, open(os.path.join(HERE, 'face_avg_ext.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)
print('done')
