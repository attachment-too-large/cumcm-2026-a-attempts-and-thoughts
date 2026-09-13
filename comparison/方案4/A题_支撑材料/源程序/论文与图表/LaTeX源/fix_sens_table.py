# -*- coding: utf-8 -*-
"""更新 tab_sens 表：全部改用本文自验数据"""
import os
import io

HERE = os.path.dirname(os.path.abspath(__file__))
rows = [
    ['基准', '本文模型 + 分段线性插值', '57.169 h', '---'],
    ['数值离散', '自由度 43 / 55 / 57', '57.1694 / 57.1696 / 57.1695 h', r'$<10^{-6}$'],
    ['插值方式', '保形三次 PCHIP', '57.169 h', r'$0.00\%$'],
    ['环境平滑', 'Savitzky--Golay (41,2)', '57.425 h', r'$+0.45\%$'],
    ['守恒形式', r'水分方程取 $\rho_d=\rho(C)/(1+C)$', '55.266 h', r'$-3.33\%$'],
    ['恒温段取值', r'取平台均值 49.969$\,^\circ$C', '57.527 h', r'$+0.63\%$'],
    ['传质系数', r'$h_m\times0.8$ / $\times1.2$', '58.867 / 56.130 h', r'$+2.97\%$ / $-1.82\%$'],
    ['收缩模型', '非均匀自洽收缩（仅问题 4，本文已复现）', r'$\approx47.61$ h', r'$-6.3\%$'],
    ['几何假设', '二维轴对称含端面（同网格隔离实验）', r'$-7.7$ s', r'$-0.004\%$'],
    ['判据口径', '改用平均含水率判定', r'$\approx35.5$ h', r'$-38\%$'],
    ['系数读法', r'把 $\ee^{-0.89/C}$ 读成 $\ee^{-0.89\,C}$（乘积）', r'$\approx16.4$ h', r'$-71\%$'],
]
hdr = ['项目', '做法', r'$t_{dry}$（问题 3）', '相对影响']
s = [r'\begin{tabular}{@{}llll@{}}', r'\toprule',
     ' & '.join(hdr) + r' \\', r'\midrule']
for r in rows:
    s.append(' & '.join(r) + r' \\')
s += [r'\bottomrule', r'\end{tabular}', '']
io.open(os.path.join(HERE, 'tabs', 'tab_sens.tex'), 'w', encoding='utf-8').write('\n'.join(s))
print('tab_sens.tex updated')
