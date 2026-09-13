# -*- coding: utf-8 -*-
"""界面平均方式判决图 -> figs/fig_faceavg.png"""
import io, os, sys
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
FIG = os.path.join(HERE, 'figs')

NS = [200, 400, 800, 1600, 3200, 6400]
ARITH_TDRY = [50.43, 50.49, 50.52, 50.54, 50.55, 50.55]
HARM_TDRY = [220.0] * 6
HARM_C72 = [0.300752, 0.274814, 0.254442, 0.238181, 0.221119, 0.207959]

fig, ax = plt.subplots(1, 2, figsize=(9.8, 3.7))

# ---------------- (a) t_dry vs N
a = ax[0]
a.semilogx(NS, [50.7691] * 6, ':', color='#2e8b57', lw=2.2, zorder=1,
           label='论文报告值 50.7691 h')
a.semilogx(NS, ARITH_TDRY, 'o-', color='#2e75b6', lw=2, ms=6, zorder=3,
           label='算术平均（交付代码）')
a.semilogx(NS, HARM_TDRY, 's--', color='#c0392b', lw=2, ms=6, zorder=3,
           label='调和平均（正文文字）')
a.annotate('四档网格全部跑不到 0.15\n（220 h 上界截断）', xy=(1600, 216), xytext=(215, 150),
           fontsize=8.6, color='#c0392b', va='center',
           arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.0))
a.annotate('', xy=(6400, 52.6), xytext=(2600, 92),
           arrowprops=dict(arrowstyle='->', color='#2e75b6', lw=1.0))
a.text(2700, 96, '加密 32 倍仅变 0.12 h\n（收敛到 50.55 h）', fontsize=8.6, color='#2e75b6')
a.set_xlabel('网格数 $N$（对数轴）', fontsize=9.5)
a.set_ylabel(r'烘干时长 $t_{dry}$ / h', fontsize=9.5)
a.set_title('(a) 烘干时长随网格的变化', fontsize=10)
a.set_ylim(25, 260); a.set_xlim(170, 8000)
a.legend(fontsize=8.2, loc='upper right', framealpha=0.95)
a.grid(alpha=0.25, which='both')

# ---------------- (b) C(0,72h) vs N
b = ax[1]
b.semilogx(NS, HARM_C72, 's--', color='#c0392b', lw=2, ms=6,
           label='调和平均：$C(0,72\\,\\mathrm{h})$')
b.axhline(0.15, color='#2e8b57', ls=':', lw=2.0)
b.text(180, 0.132, '达标判据 $C=0.15$', fontsize=8.6, color='#2e8b57', va='bottom')
b.semilogx(NS, [0.1548, 0.1549, 0.1550, 0.1550, 0.1550, 0.1550], 'o-', color='#2e75b6',
           lw=2, ms=5, label='算术平均：$C(0,72\\,\\mathrm{h})$')
b.annotate('$N$ 加密 32 倍\n仍停在 0.208', xy=(5600, 0.2065), xytext=(950, 0.268),
           fontsize=8.6, color='#c0392b', va='center',
           arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.0))
b.set_xlabel('网格数 $N$（对数轴）', fontsize=9.5)
b.set_ylabel(r'中心含水率 $C$ (kg/kg)', fontsize=9.5)
b.set_title('(b) 72 h 时的中心含水率', fontsize=10)
b.set_ylim(0.125, 0.335); b.set_xlim(170, 8000)
b.legend(fontsize=8.2, loc='upper right', framealpha=0.95)
b.grid(alpha=0.25, which='both')

fig.tight_layout()
fig.savefig(os.path.join(FIG, 'fig_faceavg.png'), bbox_inches='tight')
plt.close(fig)
print('figs/fig_faceavg.png 已生成')
