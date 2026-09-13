# -*- coding: utf-8 -*-
"""E2b  表面层分辨的等价检验：末单元细化 10 倍后，中期整个剖面是否变化

末段(C≈0.05)表面层仅约 3.8 um，任何工程网格都无法分辨；
但表面通量由 Robin 边界条件直接给出（不依赖近表面梯度），
因此只需检验“解本身”对末单元分辨是否敏感。
"""
import io
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hb_core import Ambient, load_att1, props_q23
from sem_core import SEM1D, Model

MESH_A = dict(breaks=(0.0, 0.6, 0.9, 1.0), n=18)           # 生产网格
MESH_B = dict(breaks=(0.0, 0.6, 0.9, 0.99, 1.0), n=12)     # 末单元宽度缩小 10 倍
MESH_C = dict(breaks=(0.0, 0.5, 0.8, 0.95, 0.99, 1.0), n=12)
CHK = np.array([0.0, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0])
TE = [600.0, 3600.0, 50000.0]
OUT = []


def P(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)
    io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '_e2b.txt'),
            'w', encoding='utf-8').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    t, T, C = load_att1()
    amb = Ambient(t, T, C)
    res = {}
    for tag, mesh in [('生产 (0,.6,.9,1) n=18', MESH_A),
                      ('末单元细化10倍 (0,.6,.9,.99,1) n=12', MESH_B),
                      ('单元更细 (0,.5,.8,.95,.99,1) n=12', MESH_C)]:
        s = SEM1D(**mesh)
        m = Model(s, props_q23(), amb, None, 'affine')
        t0 = time.time()
        sol = m.solve(50000.0, t_eval=TE, rtol=1e-9, atol=1e-11)
        vals = {te: s.interp(sol.y[s.ndof:, k], CHK) for k, te in enumerate(TE)}
        res[tag] = vals
        P('%-38s ndof=%3d  (%.0fs)  Csurf(3600)=%.6f  Ccenter(50000)=%.6f'
          % (tag, s.ndof, time.time() - t0, vals[3600.0][-1], vals[50000.0][0]))
    P('')
    P('以生产网格为基准的偏差:')
    base = res['生产 (0,.6,.9,1) n=18']
    for tag in list(res)[1:]:
        for te in TE:
            d = np.max(np.abs(res[tag][te] - base[te]))
            P('  %-38s t=%6.0f  max|dC| = %.3e' % (tag, te, d))
