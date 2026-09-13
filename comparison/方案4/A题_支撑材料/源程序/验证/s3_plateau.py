# -*- coding: utf-8 -*-
"""S3  恒温段设定值取“平台均值”而非“末端值”的影响"""
import io
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hb_core import Ambient, load_att1, props_q23
from sem_core import SEM1D, Model, drying_time

MESH = dict(breaks=(0.0, 0.6, 0.9, 1.0), n=18)
OUT = []


def P(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)
    io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '_s3.txt'),
            'w', encoding='utf-8').write('\n'.join(OUT) + '\n')


def go(amb, label):
    s = SEM1D(**MESH)
    m = Model(s, props_q23(), amb, None, 'affine')
    t0 = time.time()
    tf, _ = drying_time(m, t_max=400000., rtol=1e-8, atol=1e-10)
    P('  %-30s t_dry = %.1f s = %.4f h  (%.0fs)' % (label, tf, tf / 3600, time.time() - t0))
    return tf


if __name__ == '__main__':
    t, T, C = load_att1()
    base = go(Ambient(t, T, C), '末端值 50.165 ℃（基准）')
    # 平台均值：t>=7200 s 的平均
    m = t >= 7200
    Tm, Cm = float(T[m].mean()), float(C[m].mean())
    P('  平台均值: T=%.4f ℃  C=%.5f kg/kg' % (Tm, Cm))
    T2, C2 = T.copy(), C.copy()
    T2[t > 14400] = Tm
    C2[t > 14400] = Cm
    # 直接把 14400 s 之后的数据点整体替换为平台均值
    idx = t >= 14400
    T2[idx] = Tm
    C2[idx] = Cm
    alt = go(Ambient(t, T2, C2), '平台均值 %.3f ℃' % Tm)
    P('  相对基准 %+.3f%%' % ((alt - base) / base * 100))
