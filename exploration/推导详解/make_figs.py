# -*- coding: utf-8 -*-
"""生成详解正文所需的插图 -> figs/*.png"""
import sys, io, os, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, FancyArrowPatch, Rectangle, Ellipse
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['savefig.dpi'] = 200
plt.rcParams['figure.facecolor'] = 'white'

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, 'figs')
os.makedirs(FIG, exist_ok=True)
REP = os.path.join(HERE, 'verify', 'verify_report.json')


# ---------------------------------------------------------------- 图 1 控制体
def fig_shell():
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.5))
    a = ax[0]
    a.set_aspect('equal'); a.axis('off')
    a.add_patch(Rectangle((-0.9, -1.3), 1.8, 2.6, fc='#f2f6fa', ec='#9db6cc', lw=0.8))
    a.add_patch(Ellipse((0, 1.3), 1.8, 0.55, fc='#dce7f2', ec='#7c9cba', lw=1.0))
    a.add_patch(Ellipse((0, -1.3), 1.8, 0.55, fc='#dce7f2', ec='#7c9cba', lw=1.0))
    a.add_patch(Ellipse((0, 0.55), 1.8, 0.55, fc='#ffffff', ec='#7c9cba', lw=0.8, ls='--'))
    a.annotate('', xy=(0.9, -1.3), xytext=(-0.9, -1.3),
               arrowprops=dict(arrowstyle='<->', color='#c0392b', lw=1.1))
    a.text(0, -1.62, r'$2R_0=4$ cm', ha='center', color='#c0392b', fontsize=9)
    a.annotate('', xy=(0.98, -1.3), xytext=(0.98, 1.3),
               arrowprops=dict(arrowstyle='<->', color='#c0392b', lw=1.1))
    a.text(1.06, 0, r'$L=25$ cm', va='center', rotation=90, color='#c0392b', fontsize=9)
    a.text(0, -1.95, '(a) 药材：圆柱 $L/R_0=12.5$', ha='center', fontsize=9.5)
    a.set_xlim(-1.5, 1.6); a.set_ylim(-2.2, 1.75)

    b = ax[1]
    b.set_aspect('equal'); b.axis('off')
    for rr, fc in ((1.0, '#dce7f2'), (0.62, '#eef4fa'), (0.30, '#ffffff')):
        b.add_patch(Circle((0, 0), rr, fc=fc, ec='#7c9cba', lw=0.9))
    b.add_patch(Wedge((0, 0), 1.0, 52, 128, width=0.16, fc='#f6c85f', ec='#b8860b', lw=1.1))
    b.annotate('', xy=(0, 0), xytext=(0.92, 0),
               arrowprops=dict(arrowstyle='->', color='#1f4e79', lw=1.3))
    b.text(0.45, 0.075, r'$r$', color='#1f4e79', fontsize=11)
    b.annotate('', xy=(0, 0), xytext=(0.62, 0.79),
               arrowprops=dict(arrowstyle='->', color='#2e8b57', lw=1.3))
    b.text(0.44, 0.90, r'$r+\mathrm{d}r$', color='#2e8b57', fontsize=9.5)
    b.text(-2.35, 0.30, '控制体：\n环形薄壳\n$V=2\\pi rL\\,\\mathrm{d}r$', fontsize=8.6,
           color='#b8860b', ha='left', va='center', linespacing=1.6)
    b.text(-2.35, -0.62, '内面通量 $q_r(r)$\n外面通量 $q_r(r+\\mathrm{d}r)$',
           fontsize=7.8, color='#5a6b7a', ha='left', va='center', linespacing=1.6)
    b.text(0, -1.42, '(b) 横截面与控制体', ha='center', fontsize=9.5)
    b.set_xlim(-2.5, 1.35); b.set_ylim(-1.65, 1.35)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig_shell.png'), bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- 图 2 物质坐标
def fig_matcoord():
    R_tab = np.array([[0, 2.0], [3600, 1.794], [21600, 1.352], [43200, 1.268],
                      [86400, 1.226], [259200, 1.198]])
    xi_list = [0.0, 0.25, 0.5, 0.75, 1.0]
    fig, ax = plt.subplots(figsize=(9.4, 3.8))
    cmap = plt.get_cmap('viridis')
    tsel = [0, 21600, 43200, 259200]
    labels = ['$t=0$', '$t=6$ h', '$t=12$ h', '$t=72$ h']
    for k, (t, lb) in enumerate(zip(tsel, labels)):
        R = np.interp(t, R_tab[:, 0], R_tab[:, 1])
        y = 3 - k
        ax.plot([0, 2.02], [y, y], color='#c8d3de', lw=7, solid_capstyle='round', zorder=1)
        ax.plot([0, R], [y, y], color='#2e75b6', lw=7, solid_capstyle='round', zorder=2)
        for xi in xi_list:
            r = xi * R
            ax.plot([r], [y], 'o', ms=6, color=cmap(xi * 0.85), zorder=3,
                    mec='white', mew=0.8)
            if k == 0:
                ax.text(r, y + 0.28, r'$\xi=%.2f$' % xi, ha='center', fontsize=7.4,
                        color=cmap(xi * 0.85))
        ax.text(2.12, y, lb, va='center', fontsize=9.5)
        ax.text(-0.06, y, r'$R=%.3f$ cm' % R, va='center', ha='right', fontsize=8.4,
                color='#2e75b6')
        for xi in xi_list[1:]:
            r0 = xi * np.interp(0, R_tab[:, 0], R_tab[:, 1])
            r1 = xi * R
            if k > 0:
                ax.annotate('', xy=(r1, y + 0.10), xytext=(r0, 3 - 0 + 0.10),
                            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=0.8,
                                            ls=':', alpha=0.75))
    ax.set_xlim(-0.55, 2.75); ax.set_ylim(-0.55, 3.75)
    ax.set_yticks([]); ax.set_xlabel('物理半径 $r$ (cm)', fontsize=9.5)
    ax.spines[['left', 'right', 'top']].set_visible(False)
    ax.set_title('物质坐标：同一材料点（同色圆点）的 $\\xi$ 不变，物理位置 $r=\\xi R(t)$ 随收缩内移',
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig_matcoord.png'), bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- 图 3 ρd 非均匀
def fig_rhod():
    if not os.path.exists(REP):
        print('缺 verify_report.json，跳过 fig_rhod'); return
    rep = json.load(open(REP, encoding='utf-8'))
    reps = rep.get('C', {})
    xi = np.array(reps.get('xi_geo200', []))
    p8 = np.array(reps.get('profile8', []))
    p7 = np.array(reps.get('profile7', []))
    if not len(xi):
        return
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.4))
    a = ax[0]
    a.plot(xi, p8, '-', color='#2e75b6', lw=2.0, label='式 (8)（正确）')
    a.plot(xi, p7, '--', color='#c0392b', lw=1.8, label='式 (7) 按字面（对 $\\xi$ 求导）')
    a.axhline(0.15, color='#2e8b57', lw=1.0, ls=':')
    a.text(0.02, 0.18, '干燥判据 0.15', color='#2e8b57', fontsize=8.4)
    a.set_xlabel(r'物质坐标 $\xi=r/R(t)$', fontsize=9.5)
    a.set_ylabel(r'含水率 $C$ (kg/kg)', fontsize=9.5)
    a.set_title('(a) 6 h 时刻的剖面：$\\xi$ 网格 $N=200$', fontsize=10)
    a.legend(fontsize=8.4); a.grid(alpha=0.25)

    b = ax[1]
    c = np.linspace(0.05, 2.6, 300)
    rd4 = (760 + 90 * c) / (1 + c)
    b.plot(c, rd4, color='#2e75b6', lw=2)
    b.set_xlabel(r'含水率 $C$ (kg/kg)', fontsize=9.5)
    b.set_ylabel(r'干基密度 $\rho_d$ (kg/m³)', fontsize=9.5)
    b.set_title('(b) $\\rho_d=\\rho(C)/(1+C)$：从 278.7 变到 672.6  kg/m³', fontsize=10, pad=8)
    b.set_ylim(190, 830)
    b.grid(alpha=0.25)
    b.annotate(r'$C=2.55$: 278.7', xy=(2.55, 278.7), xytext=(1.30, 340),
               arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.0), fontsize=8.6,
               color='#c0392b')
    b.annotate(r'$C=0.15$: 672.6', xy=(0.15, 672.6), xytext=(0.62, 745),
               arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.0), fontsize=8.6,
               color='#c0392b')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig_eq7_rhod.png'), bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- 图 4 解析解
def fig_analytic():
    if not os.path.exists(REP):
        print('缺 verify_report.json，跳过 fig_analytic'); return
    rep = json.load(open(REP, encoding='utf-8'))
    d = rep.get('D', {})
    if not d:
        return
    N = d['N']
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.4))
    a = ax[0]
    a.loglog(N, d['err600'], 'o-', color='#2e75b6', label='$t=600$ s')
    a.loglog(N, d['err1800'], 's-', color='#2e8b57', label='$t=1800$ s')
    ref = np.array(d['err1800'])[0] * (np.array(N) / N[0]) ** -2.0
    a.loglog(N, ref, 'k:', lw=1.2, label=r'$\propto N^{-2}$（二阶收敛）')
    a.set_xlabel('网格数 $N$', fontsize=9.5)
    a.set_ylabel(r'$\max|C_{\mathrm{num}}-C_{\mathrm{ana}}|$', fontsize=9.5)
    a.set_title('(a) 与解析解（Bessel 级数）的误差', fontsize=10)
    a.legend(fontsize=8.2); a.grid(alpha=0.25, which='both')

    b = ax[1]
    xi = np.array(rep.get('C', {}).get('xi_geo200', []))
    p8 = np.array(rep.get('C', {}).get('profile8', []))
    if len(xi):
        b.plot(xi, p8, color='#2e75b6', lw=2)
        b.fill_between(xi, 0.15, p8, color='#2e75b6', alpha=0.12)
    b.axhline(0.15, color='#2e8b57', ls=':', lw=1.0)
    b.set_xlabel(r'$\xi$', fontsize=9.5); b.set_ylabel(r'$C$ (kg/kg)', fontsize=9.5)
    b.set_title('(b) 收缩算例 6 h 剖面（$R/R_0=0.68$）', fontsize=10)
    b.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig_verify.png'), bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- 图 5 时间尺度
def fig_scales():
    items = [
        ('环境温度 $\\tau_T$', 1877.5, '#e0a83c'),
        ('径向导热 $R_0^2/\\alpha$（问题1）', 2368.8, '#c0392b'),
        ('径向导热 $R_0^2/\\alpha$（问题3）', 2761.6, '#e07b39'),
        ('水分扩散 $R_0^2/D$（问题3）', 29694., '#2e8b57'),
        ('水分扩散 $R_0^2/D$（问题1）', 81041., '#2e75b6'),
        ('轴向导热 $L^2/\\alpha$', 370140., '#7f8c8d'),
    ]
    fig, ax = plt.subplots(figsize=(9.0, 3.4))
    y = np.arange(len(items))[::-1]
    for yy, (lb, v, c) in zip(y, items):
        ax.barh(yy, v, color=c, alpha=0.85, height=0.6)
        if v < 3600:
            txt = '%.0f s (%.0f min)' % (v, v / 60)
        elif v < 86400:
            txt = '%.2g s (%.1f h)' % (v, v / 3600)
        else:
            txt = '%.2g s (%.1f d)' % (v, v / 86400)
        ax.text(v * 1.15, yy, txt, va='center', fontsize=8.6)
    ax.set_yticks(y); ax.set_yticklabels([i[0] for i in items], fontsize=9)
    ax.set_xscale('log'); ax.set_xlim(1e3, 3e6)
    ax.set_xlabel('特征时间 (s，对数轴)', fontsize=9.5)
    ax.set_title('各过程的特征时间：传热比传质快约 1.5 个数量级', fontsize=10)
    ax.grid(alpha=0.25, axis='x')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig_scales.png'), bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- 图 6 物性
def fig_props():
    C = np.linspace(0, 3, 300)
    w = C / (C + 1)
    fig, ax = plt.subplots(1, 3, figsize=(9.6, 2.9))
    ax[0].plot(C, 650 + 128 * C, color='#2e75b6', lw=2, label='附录 3')
    ax[0].plot(C, 760 + 90 * C, color='#c0392b', lw=2, label='附录 4')
    ax[0].set_xlabel(r'$C$ (kg/kg)'); ax[0].set_ylabel(r'$\rho$ (kg/m³)')
    ax[0].set_title(r'$\rho$ 随含水率线性增长', fontsize=9.6)
    ax[1].plot(w, 1450 + 2736 * w, color='#2e75b6', lw=2, label='附录 3')
    ax[1].plot(w, 1850 + 2150 * w, color='#c0392b', lw=2, label='附录 4')
    ax[1].axhline(4186, color='#2e8b57', ls=':', lw=1.2)
    ax[1].text(0.03, 4290, '水的 $c_p=4186$', fontsize=7.6, color='#2e8b57')
    ax[1].set_xlabel(r'湿基含水率 $w=C/(C+1)$'); ax[1].set_ylabel(r'$c_p$ (J/(kg·K))')
    ax[1].set_title(r'$c_p$ 是质量分数混合律', fontsize=9.6)
    ax[2].semilogy(C, 2.4e-3 * np.exp(-0.45 / np.maximum(C, 1e-3)) * np.exp(-3850 / 323.15),
                   color='#2e75b6', lw=2, label='附录 3')
    ax[2].semilogy(C, 4.2e-4 * np.exp(-0.30 / np.maximum(C, 1e-3)) * np.exp(-3850 / 323.15),
                   color='#c0392b', lw=2, label='附录 4')
    ax[2].set_xlabel(r'$C$ (kg/kg)'); ax[2].set_ylabel(r'$D$ (m²/s)')
    ax[2].set_title(r'$D$ 随含水率急剧下降（$T=50$ °C）', fontsize=9.6)
    for a in ax:
        a.grid(alpha=0.25); a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig_props.png'), bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    fig_shell(); fig_matcoord(); fig_props()
    fig_rhod(); fig_analytic(); fig_scales()
    print('figs ->', FIG)
    for f in sorted(os.listdir(FIG)):
        if f.endswith('.png') and not f.startswith('_'):
            print('  %-22s %6.1f KB' % (f, os.path.getsize(os.path.join(FIG, f)) / 1024))
