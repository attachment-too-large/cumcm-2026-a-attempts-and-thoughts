# -*- coding: utf-8 -*-
"""生成 §10.7 的剖面图 -> figs/fig_profile.png"""
import io, os, sys, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['savefig.dpi'] = 200
plt.rcParams['figure.facecolor'] = 'white'

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'verify'))
import face_avg_study as S

FIG = os.path.join(HERE, 'figs')
tsel = (6 * 3600.0, 24 * 3600.0, 48 * 3600.0)
t, xi, C, T, rec = S.run(400, 'arith', dt=2.0, tmax_h=60.0, record=tsel)

fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.4))
col = ['#2e75b6', '#e0a83c', '#c0392b']
a = ax[0]
for c, tc in zip(col, tsel):
    a.plot(xi, rec[tc][0], color=c, lw=2.2, label='$t=%.0f$ h' % (tc / 3600))
a.axhline(0.15, color='#2e8b57', ls=':', lw=1.6)
a.text(0.60, 0.165, '判据 $C=0.15$', color='#2e8b57', fontsize=8.4)
a.set_xlabel(r'物质坐标 $\xi=r/R(t)$', fontsize=9.5)
a.set_ylabel(r'含水率 $C$ (kg/kg)', fontsize=9.5)
a.set_title('(a) 含水率剖面：表面附近形成陡峭边界层', fontsize=10)
a.legend(fontsize=8.6, loc='lower left')
a.grid(alpha=0.25)

b = ax[1]
c = np.linspace(0.01, 2.6, 400)
rd = (760 + 90 * c) / (1 + c)
b.plot(c, rd, color='#2e75b6', lw=2.2)
b.set_xlabel(r'含水率 $C$ (kg/kg)', fontsize=9.5)
b.set_ylabel(r'干基密度 $\rho_d$ (kg/m³)', fontsize=9.5)
b.set_title(r'(b) $\rho_d=\rho(C)/(1+C)$：278.7 → 672.6 kg/m³', fontsize=10)
b.set_ylim(200, 800)
b.grid(alpha=0.25)
b.annotate(r'$C=2.55$：278.7', xy=(2.55, 278.7), xytext=(1.25, 330),
           arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.0),
           fontsize=8.6, color='#c0392b')
b.annotate(r'$C=0.15$：672.6', xy=(0.15, 672.6), xytext=(0.60, 745),
           arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.0),
           fontsize=8.6, color='#c0392b')

fig.tight_layout()
fig.savefig(os.path.join(FIG, 'fig_profile.png'), bbox_inches='tight')
plt.close(fig)
print('figs/fig_profile.png 已生成')
