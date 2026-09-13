# -*- coding: utf-8 -*-
"""二维轴对称有限体积：验证与端面效应量化

V1  一维极限对拍：长柱 + z 边界绝热 -> 应与一维谱元解一致（径向剖面）
V2  网格与时间步收敛
V3  真实尺寸（L=25 cm，含端面）-> 中心 C(t) 与 t_dry，量化端面影响
"""
import io
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hb_core import (Ambient, load_att1, props_q23, R0, LZ)
from sem_core import SEM1D, Model
import fv2d as F

OUT = []
MESH1D = dict(breaks=(0.0, 0.6, 0.9, 1.0), n=18)


def P(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)
    io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '_fv2d.txt'),
            'w', encoding='utf-8').write('\n'.join(OUT) + '\n')


def v1_one_dim():
    P('=' * 76)
    P('V1  一维极限对拍（长柱 L=2 m + z 边界绝热，应等于一维谱元解）')
    t, Ta, Ca = load_att1()
    amb = Ambient(t, Ta, Ca)
    # 一维谱元参考
    S = SEM1D(**MESH1D)
    m = Model(S, props_q23(), amb, None, 'affine')
    TE = [600.0, 3600.0, 21600.0]
    sol = m.solve(max(TE), t_eval=TE, rtol=1e-9, atol=1e-11)
    # 二维（一维极限）
    s2 = F.Solver2D(props_q23(), amb, None, Nr=80, Nz=4, L=2.0, bnd_z=False)
    t0 = time.time()
    st = s2.solve(max(TE), 15.0, t_eval=TE, field=True)
    P('  二维耗时 %.0f s' % (time.time() - t0))
    P('  %-8s %-14s %-14s %-12s' % ('t/s', '一维(谱元)', '二维(z绝热)', '最大偏差'))
    for te in TE:
        prof = st[te][4]
        rc = s2.g['rc']
        sig = rc / R0
        ref = S.interp(sol.y[S.ndof:, TE.index(te)], sig)
        d = np.max(np.abs(prof - ref))
        P('  %-8.0f %-14.6f %-14.6f %.2e' % (te, ref[0], prof[0], d))
        P('           端面处(中心=%.6f) 剖面最大偏差 %.3e' % (ref[0], d))
    return


def v2_convergence():
    P('=' * 76)
    P('V2  网格与时间步收敛（真实尺寸 L=25 cm，含端面，积分到 21600 s）')
    t, Ta, Ca = load_att1()
    amb = Ambient(t, Ta, Ca)
    for Nr, Nz, dt in [(40, 25, 60.0), (60, 40, 60.0), (80, 50, 60.0),
                       (80, 50, 30.0), (120, 75, 60.0)]:
        s2 = F.Solver2D(props_q23(), amb, None, Nr=Nr, Nz=Nz, L=LZ, bnd_z=True)
        t0 = time.time()
        st = s2.solve(21600.0, dt, t_eval=[21600.0])
        c0 = st[21600.0][0]
        P('  Nr=%3d Nz=%3d dt=%4.0f s  C_center(6h)=%.6f  (%.0fs)'
          % (Nr, Nz, dt, c0, time.time() - t0))


def v3_endface():
    P('=' * 76)
    P('V3  端面效应：二维(含端面) 与 一维(无限长) 的中心 C(t) 与烘干时长对比')
    from sem_core import drying_time
    t, Ta, Ca = load_att1()
    amb = Ambient(t, Ta, Ca)
    S = SEM1D(**MESH1D)
    m = Model(S, props_q23(), amb, None, 'affine')
    t0 = time.time()
    tf1, _ = drying_time(m, t_max=400000., rtol=1e-8, atol=1e-10)
    P('  一维（无限长圆柱）t_dry = %.1f s = %.4f h  (%.0fs)' % (tf1, tf1 / 3600, time.time() - t0))
    TCHK = [6 * 3600, 24 * 3600, 48 * 3600, tf1]
    s1 = m.solve(tf1, t_eval=TCHK, rtol=1e-8, atol=1e-10)
    c1s = [float(S.interp(s1.y[S.ndof:, k], [0.0])[0]) for k in range(len(TCHK))]
    # 二维
    s2 = F.Solver2D(props_q23(), amb, None, Nr=60, Nz=40, L=LZ, bnd_z=True)
    P('  二维网格 Nr=60 Nz=40，dt=60 s，积分至一维 t_dry ...')
    t0 = time.time()
    st = s2.solve(tf1, 60.0, t_eval=TCHK)
    P('  二维耗时 %.0f s' % (time.time() - t0))
    P('  %-10s %-14s %-14s %-10s' % ('t/h', '一维 C_center', '二维 C_center', '相对差'))
    for k, te in enumerate(TCHK):
        c2 = st[te][0]
        P('  %-10.2f %-14.6f %-14.6f %+.3e'
          % (te / 3600, c1s[k], c2, (c2 - c1s[k]) / max(c1s[k], 1e-12)))
    tf2, _, _ = s2.drying_time(60.0, t_max=400000.)
    P('  二维（含端面）t_dry = %.1f s = %.4f h' % (tf2, tf2 / 3600))
    P('  端面使 t_dry 变化 %+.3f%%（相对一维）' % ((tf2 - tf1) / tf1 * 100))
    P('  中心单元与端面的距离 z=%.2f cm，药材半长 %.2f cm' % (LZ / 2, LZ / 2))


def v4_tdry_convergence():
    P('=' * 76)
    P('V4  二维烘干时长的网格与时间步收敛（含端面）')
    t, Ta, Ca = load_att1()
    amb = Ambient(t, Ta, Ca)
    for Nr, Nz, dt in [(40, 25, 60.0), (60, 40, 60.0), (80, 50, 60.0), (60, 40, 30.0)]:
        s2 = F.Solver2D(props_q23(), amb, None, Nr=Nr, Nz=Nz, L=LZ, bnd_z=True)
        t0 = time.time()
        tf, _, _ = s2.drying_time(dt, t_max=400000.)
        P('  Nr=%3d Nz=%3d dt=%4.0f s  t_dry = %.1f s = %.4f h  (%.0fs)'
          % (Nr, Nz, dt, tf, tf / 3600, time.time() - t0))


def v5_isolate():
    P('=' * 76)
    P('V5  严格隔离端面效应：同一网格、同一时间步，只切换 z=L/2 的边界条件')
    P('    (a) bnd_z=True  -> 端面暴露于热风（真实情形）')
    P('    (b) bnd_z=False -> 端面绝热（等价于无限长柱）')
    P('    两者之差即端面效应的纯贡献，不含离散误差的干扰。')
    t, Ta, Ca = load_att1()
    amb = Ambient(t, Ta, Ca)
    TCHK = [6 * 3600, 24 * 3600, 48 * 3600]
    for Nr, Nz in [(40, 25), (60, 40), (80, 50)]:
        res = {}
        for bz in (True, False):
            s2 = F.Solver2D(props_q23(), amb, None, Nr=Nr, Nz=Nz, L=LZ, bnd_z=bz)
            res[bz] = s2.solve(48 * 3600.0, 60.0, t_eval=TCHK)
        P('  网格 Nr=%3d Nz=%2d  (dt=60 s)' % (Nr, Nz))
        for te in TCHK:
            a, b = res[True][te][0], res[False][te][0]
            P('    t=%5.0f h :  含端面 %.6f   绝热端面 %.6f   端面效应 %+.3e (%.3f%%)'
              % (te / 3600, a, b, a - b, (a - b) / b * 100))
        s2t = F.Solver2D(props_q23(), amb, None, Nr=Nr, Nz=Nz, L=LZ, bnd_z=True)
        s2f = F.Solver2D(props_q23(), amb, None, Nr=Nr, Nz=Nz, L=LZ, bnd_z=False)
        tt, _, _ = s2t.drying_time(60.0)
        tf, _, _ = s2f.drying_time(60.0)
        P('    t_dry :      含端面 %.1f s   绝热端面 %.1f s   差 %+.1f s (%.3f%%)'
          % (tt, tf, tt - tf, (tt - tf) / tf * 100))


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('all', 'v1'):
        v1_one_dim()
    if which in ('all', 'v2'):
        v2_convergence()
    if which in ('all', 'v3'):
        v3_endface()
    if which in ('all', 'v4'):
        v4_tdry_convergence()
    if which in ('all', 'v5'):
        v5_isolate()
