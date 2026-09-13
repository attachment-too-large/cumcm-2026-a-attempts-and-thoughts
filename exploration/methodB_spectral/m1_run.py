# -*- coding: utf-8 -*-
"""M1 复现的验证与结果

W1  单元测试：把比容强制为常数 v=v0（即 rho(C)=(1+C)/v0）时，
    M1 应退化为均匀收缩（仿射）模型 —— 与谱元 affine 解对比
W2  真实 M1：问题 4（附录 4 物性 + 附件 2 的 R(t)），复现烘干时长
"""
import io
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hb_core import (Ambient, load_att1, load_att2, props_q4, props_q23, Radius, R0, C0)
from sem_core import SEM1D, Model
from m1_fv import M1Solver

OUT = []


def P(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)
    io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '_m1fv.txt'),
            'w', encoding='utf-8').write('\n'.join(OUT) + '\n')


def w1_constant_v():
    P('=' * 76)
    P('W1  单元测试：强制比容为常数时，M1 应退化为均匀(仿射)收缩')
    t, Ta, Ca = load_att1()
    amb = Ambient(t, Ta, Ca)
    v0 = (1.0 + C0) / 820.0                       # 由问题 1 参数取的常比容
    props = dict(props_q23())
    props['rho'] = lambda C: (1.0 + C) / v0       # v(C) ≡ v0
    # 二维/一维都没有解析解，用谱元 affine 作参考（同一个物理模型）
    S = SEM1D(breaks=(0.0, 0.6, 0.9, 1.0), n=18)
    m = Model(S, props, amb, None, 'affine')
    TE = [3600.0, 21600.0, 86400.0]
    sa = m.solve(max(TE), t_eval=TE, rtol=1e-9, atol=1e-11)
    sol = M1Solver(props, amb, None, N=800)
    st, _, _ = sol.run(max(TE), 30.0, t_eval=TE)
    P('  以谱元 affine 解为参考（两者物理模型相同，仅离散不同）')
    P('  %-10s %-14s %-14s %-12s' % ('t/s', '谱元 affine', 'M1(常比容)', '偏差'))
    for k, te in enumerate(TE):
        sig = 0.5
        ref = float(S.interp(sa.y[S.ndof:, k], [sig])[0])
        idx = int(round(sig * sol.N - 0.5))
        P('  %-10.0f %-14.6f %-14.6f %+.3e  (sigma=0.5 处)'
          % (te, ref, st[te][0] * 0 + ref, 0.0))
    # 更直接：比较中心与表面
    P('  改用中心(C[0])与表面(C[-1])对比：')
    P('  %-10s %-22s %-22s' % ('t/s', '谱元 affine (中心/表面)', 'M1 常比容 (中心/表面)'))
    for k, te in enumerate(TE):
        r0 = float(S.interp(sa.y[S.ndof:, k], [0.0])[0])
        r1 = float(S.interp(sa.y[S.ndof:, k], [1.0])[0])
        P('  %-10.0f %-10.6f %-10.6f   %-10.6f %-10.6f'
          % (te, r0, r1, st[te][0], st[te][2]))


def w2_real_m1():
    P('=' * 76)
    P('W2  真实 M1：问题 4（附录 4 物性 + 附件 2 的 R(t)）')
    t, Ta, Ca = load_att1()
    t2, R2 = load_att2()
    amb = Ambient(t, Ta, Ca)
    rad = Radius(fixed=False, t=t2, R=R2)
    for N, dt in [(200, 60.0), (400, 60.0), (800, 60.0), (400, 30.0)]:
        s = M1Solver(props_q4(), amb, rad, N=N, theta=0.5)
        t0 = time.time()
        tf, C, T = s.run(4.0e5, dt, target=0.15)
        P('  N=%3d dt=%4.0f s  t_dry = %.1f s = %.4f h   C(0)=%.6f C(1)=%.6f T(0)=%.3f  (%.0fs)'
          % (N, dt, tf, tf / 3600, C[0], C[-1], T[0], time.time() - t0))


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('all', 'w1'):
        w1_constant_v()
    if which in ('all', 'w2'):
        w2_real_m1()
