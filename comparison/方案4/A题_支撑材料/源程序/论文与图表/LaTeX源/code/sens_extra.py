# -*- coding: utf-8 -*-
"""补充可信度实验：
S1  烘房环境插值方式对 t_dry 的影响（线性 / 保形三次 PCHIP / Savitzky-Golay 平滑）
S2  对流传质系数 h_m 的灵敏度（×0.8 / ×1.2）
"""
import io
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hb_core import Ambient, load_att1, props_q23, HM_MASS
import hb_core
import sem_core                      # 注意: sem_core 用 from ... import 绑定了常量副本,
from sem_core import SEM1D, Model, drying_time   # 必须同时改 sem_core 的模块全局量

MESH = dict(breaks=(0.0, 0.6, 0.9, 1.0), n=18)
OUT = []


def P(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)
    io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sens.txt'),
            'w', encoding='utf-8').write('\n'.join(OUT) + '\n')


def set_hm(v):
    """同时更新 hb_core 与 sem_core 两个模块里的 HM_MASS"""
    hb_core.HM_MASS = v
    sem_core.HM_MASS = v


def run(amb, label, hm=None):
    old = sem_core.HM_MASS
    if hm is not None:
        set_hm(hm)
    # 断言生效（防止静默失效）
    assert abs(sem_core.HM_MASS - (hm if hm is not None else old)) < 1e-18
    s = SEM1D(**MESH)
    m = Model(s, props_q23(), amb, None, 'affine')
    t0 = time.time()
    tf, _ = drying_time(m, t_max=400000., rtol=1e-8, atol=1e-10)
    set_hm(old)
    P('  %-34s t_dry = %9.1f s = %.4f h   (%.0fs)  [hm=%.3e]'
      % (label, tf, tf / 3600, time.time() - t0, hm if hm is not None else old))
    return tf


if __name__ == '__main__':
    t, T, C = load_att1()
    P('=' * 70)
    P('S1 环境插值方式对问题3烘干时长的影响')
    base = run(Ambient(t, T, C, mode='linear'), '基准：分段线性')
    pch = run(Ambient(t, T, C, mode='pchip'), '保形三次 PCHIP')
    sg = run(Ambient(t, T, C, mode='linear', smooth=41), 'Savitzky-Golay(41,2) 平滑')
    P('    PCHIP 相对基准 %+.3f%%   平滑相对基准 %+.3f%%   散布 %.3f%%'
      % ((pch - base) / base * 100, (sg - base) / base * 100,
         (max(base, pch, sg) - min(base, pch, sg)) / base * 100))
    P('')
    P('=' * 70)
    P('S2 对流传质系数 h_m 的灵敏度（附录2 值 = %.2e m/s）' % HM_MASS)
    a = run(Ambient(t, T, C), 'h_m x 0.8', hm=0.8 * HM_MASS)
    b = run(Ambient(t, T, C), 'h_m x 1.0（基准）', hm=HM_MASS)
    c = run(Ambient(t, T, C), 'h_m x 1.2', hm=1.2 * HM_MASS)
    P('    h_m 变化 -20%%/+20%% 对应 t_dry 变化 %+.3f%% / %+.3f%%（弹性 %.2f）'
      % ((a - b) / b * 100, (c - b) / b * 100, ((c - a) / b * 100) / 40.0 * 100))
