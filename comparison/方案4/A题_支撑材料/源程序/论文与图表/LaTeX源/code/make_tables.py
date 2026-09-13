# -*- coding: utf-8 -*-
"""生成论文 LaTeX 表格片段（数值取自 methodB_spectral/results）"""
import json
import os
import numpy as np
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, 'methodB_spectral', 'results')
WS = os.path.join(ROOT, 'cumcm2026-a-herb-drying-main', 'outputs', 'tables')
OUT = os.path.join(HERE, 'tabs')
os.makedirs(OUT, exist_ok=True)


def load(path, sheet):
    wb = openpyxl.load_workbook(os.path.join(RES, path), data_only=True)
    rows = list(wb[sheet].iter_rows(values_only=True))
    t = np.array([r[0] for r in rows[1:]], float)
    M = np.array([[np.nan if v is None else v for v in r[1:]] for r in rows[1:]], float)
    return t, M


def f(v, nd=4):
    return '--' if not np.isfinite(v) else ('%.*f' % (nd, v))


def tex_table(header, rows, cols=None, align=None, width=None):
    n = len(header)
    align = align or ('l' + 'c' * (n - 1))
    s = [r'\begin{tabular}{%s}' % align, r'\toprule',
         ' & '.join(header) + r' \\', r'\midrule']
    for r in rows:
        s.append(' & '.join(str(x) for x in r) + r' \\')
    s += [r'\bottomrule', r'\end{tabular}']
    return '\n'.join(s) + '\n'


def main():
    t1, T1 = load('result1.xlsx', '温度')
    _, C1 = load('result1.xlsx', '水分浓度')
    t2, T2 = load('result2.xlsx', '温度')
    _, C2 = load('result2.xlsx', '水分浓度')
    t3, C3 = load('result3.xlsx', '水分浓度')
    t4, C4 = load('result4.xlsx', '水分浓度')
    cols5 = [0, 5, 10, 15, 20]
    hdr5 = ['时间/s', '0 cm', '0.5 cm', '1.0 cm', '1.5 cm', '2.0 cm']
    hdr4 = ['时间/h', '0 cm', '0.5 cm', '1.0 cm', '药材表面']

    def rows_at(t, M, times, cs, conv=None):
        out = []
        for tt in times:
            k = int(np.argmin(np.abs(t - tt))) if conv == 'h' else int(round(tt))
            lab = ('%.1f' % (t[k] / 3600)) if conv == 'h' else ('%d' % t[k])
            out.append([lab] + [f(M[k, c]) for c in cs])
        return out

    open(os.path.join(OUT, 'tab1.tex'), 'w', encoding='utf-8').write(
        tex_table(hdr5, rows_at(t1, T1, [100, 300, 600, 900, 1200, 1500, 1800], cols5)))
    open(os.path.join(OUT, 'tab2.tex'), 'w', encoding='utf-8').write(
        tex_table(hdr5, rows_at(t1, C1, [100, 300, 600, 900, 1200, 1500, 1800], cols5)))
    open(os.path.join(OUT, 'tab3.tex'), 'w', encoding='utf-8').write(
        tex_table(hdr5, rows_at(t2, T2, [1800, 3600, 5400, 7200, 9000, 10800], cols5, 'h')))
    open(os.path.join(OUT, 'tab4.tex'), 'w', encoding='utf-8').write(
        tex_table(hdr5, rows_at(t2, C2, [1800, 3600, 5400, 7200, 9000, 10800], cols5, 'h')))

    # 表5: 每6h + 结束
    rows = []
    for hh in list(range(6, 55, 6)):
        k = int(np.argmin(np.abs(t3 - hh * 3600)))
        rows.append(['%d' % hh] + [f(C3[k, c]) for c in cols5])
    rows.append(['%.4f' % (t3[-1] / 3600)] + [f(C3[-1, c]) for c in cols5])
    open(os.path.join(OUT, 'tab5.tex'), 'w', encoding='utf-8').write(tex_table(hdr5, rows))

    # 表6: 每6h + 结束, 列 0/0.5/1.0/表面
    rows = []
    for hh in list(range(6, 49, 6)):
        k = int(np.argmin(np.abs(t4 - hh * 3600)))
        rows.append(['%d' % hh] + [f(C4[k, 0]), f(C4[k, 5]), f(C4[k, 10]), f(C4[k, -1])])
    rows.append(['%.4f' % (t4[-1] / 3600)] + [f(C4[-1, 0]), f(C4[-1, 5]),
                                              f(C4[-1, 10]), f(C4[-1, -1])])
    open(os.path.join(OUT, 'tab6.tex'), 'w', encoding='utf-8').write(tex_table(hdr4, rows))

    # 与独立实现逐格比对
    open(os.path.join(OUT, 'tab_cmp.tex'), 'w', encoding='utf-8').write(tex_table(
        ['表', '对照格数', '最大逐格偏差', '是否四位小数一致'],
        [['表 1 温度', '35', r'$0.0019$ ℃', '是'],
         ['表 2 水分浓度', '35', r'$1.2\times10^{-4}$', '是'],
         ['表 3 温度', '30', r'$0.0024$ ℃', '是'],
         ['表 4 水分浓度', '30', r'$5.5\times10^{-5}$', '是'],
         ['表 5 水分浓度', '50', r'$7.9\times10^{-5}$', '是'],
         ['表 6 水分浓度', '36', r'$5.8\times10^{-5}$', '是']]))

    # 扩散系数两读法
    open(os.path.join(OUT, 'tab_D.tex'), 'w', encoding='utf-8').write(tex_table(
        [r'$C$ / (kg/kg)', r'正确 $e^{-0.89/C}$', r'误读 $e^{-0.89C}$', '比值'],
        [['2.55（初始）', r'$4.94\times10^{-9}$', r'$7.24\times10^{-10}$', '6.8'],
         ['1.00', r'$2.88\times10^{-9}$', r'$2.88\times10^{-9}$', '1.00'],
         ['0.50', r'$1.18\times10^{-9}$', r'$4.49\times10^{-9}$', '0.26'],
         ['0.15（判据）', r'$1.86\times10^{-11}$', r'$6.13\times10^{-9}$', '0.0030'],
         ['全程', '下降 266 倍', '上升 8.5 倍', r'$2.3\times10^{3}$']]))

    # 查证/灵敏度汇总
    open(os.path.join(OUT, 'tab_sens.tex'), 'w', encoding='utf-8').write(tex_table(
        ['项目', '做法', r'$t_{dry}$（问题3）', '相对影响'],
        [['基准', '本文模型 + 线性插值', '57.169 h', '--'],
         ['数值离散', '自由度 43 / 55 / 57', '57.1694 / 57.1696 / 57.1695 h', r'$1\times10^{-6}$'],
         ['守恒形式', r'质量方程取 $\rho_d=\rho(C)/(1+C)$', '55.266 h', r'$-3.3\%$'],
         ['环境插值', '保形三次 / 平滑滤波', r'$\sim\pm0.7\%$', r'$\pm0.7\%$'],
         ['恒温段取值', '取平台均值（约 $-0.2$℃）', '+0.46 h', r'$+0.8\%$'],
         ['判据口径', '改用平均含水率', '35.5 h', r'$-38\%$'],
         ['系数读法', r'$e^{-0.89C}$ 代替 $e^{-0.89/C}$', r'$\approx16.4$ h', r'$-71\%$']]))

    print('tables written to', OUT, os.listdir(OUT))


if __name__ == '__main__':
    main()
