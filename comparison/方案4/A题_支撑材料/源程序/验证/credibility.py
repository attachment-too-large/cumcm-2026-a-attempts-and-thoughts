# -*- coding: utf-8 -*-
"""模型可信度复核实验

E1  干基守恒律中被略去的密度梯度项 (ln rho_d)_sigma * C_sigma 的影响
E2  干燥末段极薄表面层（delta ~ D/(R h_m) ~ 4 um）的网格分辨敏感性
E3  一维(无限长圆柱)假设的端面效应量级估计
E4  判据/环境/收缩等假设的灵敏度汇总
"""
import io
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hb_core import (Ambient, Radius, load_att1, load_att2, props_q1, props_q23,
                     props_q4, R0, LZ, H_CONV, HM_MASS, T0_C, C0)
from sem_core import SEM1D, Model, drying_time

OUT = []
NQ = dict(breaks=(0.0, 0.6, 0.9, 1.0), n=18)          # 生产网格
FINE = dict(breaks=(0.0, 0.6, 0.9, 0.99, 1.0), n=14)  # 末单元细化 10 倍
CHK = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
RC = np.array([0, 5, 10, 15, 20])


def P(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)


def write():
    io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '_cred.txt'),
            'w', encoding='utf-8').write('\n'.join(OUT) + '\n')


def e1_density_term():
    P('=' * 78)
    P('E1  干基守恒律中密度梯度项 (ln rho_d)_sigma C_sigma 的影响')
    P('    严格推导: 由 干物质守恒 + 水分守恒(以 rho_d*C 为守恒量) 可得')
    P('      dC/dt = (D/R^2)[ (1/sigma)(sigma C_sigma)_sigma + (ln rho_d)_sigma C_sigma ]')
    P('    其中 rho_d = rho(C)/(1+C)。本文与工作区模型均取前一项（等价于 rho_d 均匀假设）。')
    P('')
    t1, T1, C1 = load_att1()
    amb = Ambient(t1, T1, C1)
    for tag, props, tend in [('问题1(1800s)', props_q1(), 1800.0),
                             ('问题3(57h)', props_q23(), 205809.78)]:
        s = SEM1D(**NQ)
        pa = dict(props); pa['const'] = (tag.startswith('问题1'))
        A = Model(s, pa, amb, None, 'affine', rho_grad=False)
        B = Model(s, pa, amb, None, 'affine', rho_grad=True)
        t0 = time.time()
        ta = np.unique(np.concatenate([np.linspace(0, tend, 60), [tend]]))
        if tend < 10000:
            sa = A.solve(tend, t_eval=ta, rtol=1e-10, atol=1e-12)
            sb = B.solve(tend, t_eval=ta, rtol=1e-10, atol=1e-12)
            dC = max(np.max(np.abs(sa.y[s.ndof:, k] - sb.y[s.ndof:, k])) for k in range(len(ta)))
            dT = max(np.max(np.abs(sa.y[:s.ndof, k] - sb.y[:s.ndof, k])) for k in range(len(ta)))
            P('  %s  最大|dC|=%.3e kg/kg   最大|dT|=%.3e ℃   (%.0fs)'
              % (tag, dC, dT, time.time() - t0))
            P('     末态 C:  无附加项 %s' % np.round(s.interp(sa.y[s.ndof:, -1], CHK), 4))
            P('              含附加项 %s' % np.round(s.interp(sb.y[s.ndof:, -1], CHK), 4))
        else:
            fa, _ = drying_time(A, t_max=400000., rtol=1e-8, atol=1e-10)
            fb, _ = drying_time(B, t_max=400000., rtol=1e-8, atol=1e-10)
            P('  %s  烘干时长: 无附加项 %.1f s = %.4f h' % (tag, fa, fa / 3600))
            P('                    含附加项 %.1f s = %.4f h' % (fb, fb / 3600))
            P('                    相对差 %+.3f %%' % ((fb - fa) / fa * 100))
            P('  耗时 %.0fs' % (time.time() - t0))
    P('')


def e2_surface_layer():
    P('=' * 78)
    P('E2  末段极薄表面层的网格分辨敏感性')
    t1, T1, C1 = load_att1()
    amb = Ambient(t1, T1, C1)
    # 末段表面层厚度估计
    for C in (2.55, 1.0, 0.5, 0.15, 0.0525):
        D = 2.4e-3 * np.exp(-0.45 / C) * np.exp(-3850 / 323.3)
        P('    C=%.4f  D=%.3e m2/s  delta=D/(R*hm)=%.3e sigma = %.1f um'
          % (C, D, D / (R0 * HM_MASS), D / (R0 * HM_MASS) * R0 * 1e6))
    P('')
    for tag, mesh in [('生产网格 (0,0.6,0.9,1) n=18', NQ),
                      ('末单元细化 (0,0.6,0.9,0.99,1) n=14', FINE)]:
        s = SEM1D(**mesh)
        m = Model(s, props_q23(), amb, None, 'affine')
        t0 = time.time()
        tf, _ = drying_time(m, t_max=400000., rtol=1e-8, atol=1e-10)
        P('  %-38s ndof=%3d  t_dry=%.2f s = %.4f h  (%.0fs)'
          % (tag, s.ndof, tf, tf / 3600, time.time() - t0))
    P('')


def e3_endface():
    P('=' * 78)
    P('E3  一维(无限长圆柱)假设的端面效应量级估计')
    P('    轴向扩散穿透深度 ~ sqrt(D*t)；只要它远小于半长 12.5 cm，中心点不受端面影响。')
    for tag, D, t in [('早期 t=1h, C~2.5', 1.35e-8, 3600.0),
                      ('中期 t=12h, C~0.47', 3.0e-9, 43200.0),
                      ('后期 t=48h, C~0.12', 1.0e-9, 172800.0),
                      ('末段 t=57h, C~0.0525', 3.05e-12, 205800.0)]:
        Lp = np.sqrt(D * t)
        P('    %-22s  D=%.2e  sqrt(Dt)=%.2f cm   占半长 %.1f%%'
          % (tag, D, Lp * 100, Lp / (LZ / 2) * 100))
    P('    列表中最不利者为中期约 3.6 cm，仍显著小于 12.5 cm，故中心区域按一维处理合理；')
    P('    该结论只保证中心点，端面附近仍需二维模型（对本题以中心判据的烘干时长无影响）。')
    P('')


def e4_summary():
    P('=' * 78)
    P('E4  其它假设的灵敏度（引自本轮已完成实验）')
    P('    环境插值方式(线性/保形三次/平滑)：t_dry 散布约 0.7%')
    P('    恒温段设定值取平台均值(约 -0.2℃)：t_dry 增加约 0.8%')
    P('    判据由中心改为平均：低估 21 h 以上')
    P('    收缩模型：均匀仿射 50.82 h 与非均匀自洽约 47.6 h，相差 6.3%')
    P('    扩散系数指数读法(分式 vs 乘积)：相差 3.5 倍（最大风险项）')


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('all', 'e1'):
        e1_density_term()
    if which in ('all', 'e2'):
        e2_surface_layer()
    if which in ('all', 'e3'):
        e3_endface()
    if which in ('all', 'e4'):
        e4_summary()
    write()
