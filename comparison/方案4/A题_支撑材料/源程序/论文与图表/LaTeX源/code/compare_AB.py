# -*- coding: utf-8 -*-
"""A/B 两法逐格交叉比对（A=工作区有限体积, B=本方法谱元）"""
import json
import os
import numpy as np
import openpyxl

WS = os.environ.get('CUMCM_WS_OUT')
if WS is None:
    _here = os.path.dirname(os.path.abspath(__file__))
    for _c in [os.path.join(os.path.dirname(_here),
                            'cumcm2026-a-herb-drying-main', 'outputs'),
               os.path.join(_here, 'results', '对照数据'),
               os.path.join(_here, 'results')]:
        if os.path.isfile(os.path.join(_c, 'tables', 'tables.json')):
            WS = _c
            break
if WS is None:
    raise SystemExit('未找到对照数据 tables.json：'
                     '请设置环境变量 CUMCM_WS_OUT 指向方法 A 的 outputs 目录')
MY = os.environ.get('CUMCM_RESULT_DIR')
if MY is None:
    _here = os.path.dirname(os.path.abspath(__file__))
    for _c in [os.path.join(_here, 'results'),
               os.path.join(os.path.dirname(_here), '结果文件')]:
        if os.path.isfile(os.path.join(_c, 'result1.xlsx')):
            MY = _c
            break
if MY is None:
    raise SystemExit('未找到 result1.xlsx：请设置环境变量 CUMCM_RESULT_DIR')
OUT = []


def P(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s)


def load_ws_tables():
    with open(os.path.join(WS, 'tables', 'tables.json'), encoding='utf-8') as f:
        d = json.load(f)
    # 键是乱码的中文, 按顺序取: 温度表1, 水分表2, 温度表3, 水分表4, 问题3表, 问题4表
    keys = list(d.keys())
    return d, keys


def load_my(path, sheet):
    wb = openpyxl.load_workbook(os.path.join(MY, path), data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    hdr = rows[0][1:]
    t = np.array([r[0] for r in rows[1:]], float)
    M = np.array([[np.nan if v is None else v for v in r[1:]] for r in rows[1:]], float)
    return t, hdr, M


def main():
    d, keys = load_ws_tables()
    P('# A/B 交叉比对（A = 工作区节点中心有限体积 + CN；B = 本文分域谱元弱形式 + BDF）\n')

    # ---------------- 问题1 ----------------
    t, hdr, T = load_my('result1.xlsx', '温度')
    _, _, C = load_my('result1.xlsx', '水分浓度')
    P('## 问题1 表1 温度 (℃)  [列: 0,0.5,1,1.5,2 cm]')
    P('| t/s | B(本文) | A(工作区) | 差 |')
    P('|---|---|---|---|')
    for i, tt in enumerate(d['表1_t']):
        colA = d['表1'][i]
        k = int(tt)
        colB = [T[k, int(round(r * 10))] for r in (0, .5, 1, 1.5, 2)]
        P('| %d | %s | %s | %.4g |' % (tt, ' '.join('%.4f' % v for v in colB),
                                       ' '.join('%.4f' % v for v in colA),
                                       max(abs(a - b) for a, b in zip(colA, colB))))
    P('\n## 问题1 表2 水分浓度 (kg/kg)')
    P('| t/s | B(本文) | A(工作区) | 差 |')
    P('|---|---|---|---|')
    for i, tt in enumerate(d['表1_t']):
        colA = d['表2'][i]
        k = int(tt)
        colB = [C[k, int(round(r * 10))] for r in (0, .5, 1, 1.5, 2)]
        P('| %d | %s | %s | %.4g |' % (tt, ' '.join('%.4f' % v for v in colB),
                                       ' '.join('%.4f' % v for v in colA),
                                       max(abs(a - b) for a, b in zip(colA, colB))))

    # ---------------- 问题2 ----------------
    t2, hdr2, T2 = load_my('result2.xlsx', '温度')
    _, _, C2 = load_my('result2.xlsx', '水分浓度')
    P('\n## 问题2 表3 温度 (℃)  [0.5h 起]')
    P('| t/h | B(本文) | A(工作区) | 差 |')
    P('|---|---|---|---|')
    for i, tt in enumerate(d['表3_t']):
        colA = d['表3'][i]
        k = int(tt)
        colB = [T2[k, int(round(r * 10))] for r in (0, .5, 1, 1.5, 2)]
        P('| %.1f | %s | %s | %.4g |' % (tt / 3600, ' '.join('%.4f' % v for v in colB),
                                         ' '.join('%.4f' % v for v in colA),
                                         max(abs(a - b) for a, b in zip(colA, colB))))
    P('\n## 问题2 表4 水分浓度 (kg/kg)')
    P('| t/h | B(本文) | A(工作区) | 差 |')
    P('|---|---|---|---|')
    for i, tt in enumerate(d['表3_t']):
        colA = d['表4'][i]
        k = int(tt)
        colB = [C2[k, int(round(r * 10))] for r in (0, .5, 1, 1.5, 2)]
        P('| %.1f | %s | %s | %.4g |' % (tt / 3600, ' '.join('%.4f' % v for v in colB),
                                         ' '.join('%.4f' % v for v in colA),
                                         max(abs(a - b) for a, b in zip(colA, colB))))

    # ---------------- 问题3 ----------------
    t3, hdr3, C3 = load_my('result3.xlsx', '水分浓度')
    P('\n## 问题3 表5 水分浓度 (kg/kg)  [A 为 N=400 生产值]')
    P('| t/h | B(本文) | A(工作区) | 差 |')
    P('|---|---|---|---|')
    for i, th in enumerate(d['表5_t_h']):
        colA = d['表5'][i]
        k = int(np.argmin(np.abs(t3 - th * 3600)))
        colB = [C3[k, int(round(r * 10))] for r in (0, .5, 1, 1.5, 2)]
        P('| %.3f | %s | %s | %.3g |' % (th, ' '.join('%.4f' % v for v in colB),
                                         ' '.join('%.4f' % v for v in colA),
                                         max(abs(a - b) for a, b in zip(colA, colB))))

    # ---------------- 问题4（B 用当前构形坐标；材料外留空）----------------
    t4, hdr4, C4c = load_my('result4.xlsx', '水分浓度')
    P('\n## 问题4 表6 水分浓度 (kg/kg)  [列: 0, 0.5, 1.0 cm, 药材表面]')
    P('| t/h | B(本文) | A(工作区) | 差 |')
    P('|---|---|---|---|')
    mx = 0.0
    for i, th in enumerate(d['表6_t_h']):
        colA = d['表6'][i]
        k = int(np.argmin(np.abs(t4 - th * 3600)))
        colB = [C4c[k, 0], C4c[k, 5], C4c[k, 10], C4c[k, -1]]
        dd = max(abs(a - b) for a, b in zip(colA, colB))
        mx = max(mx, dd)
        P('| %.3f | %s | %s | %.3g |' % (th, ' '.join('%.4f' % v for v in colB),
                                         ' '.join('%.4f' % v for v in colA), dd))
    P('\n（问题4 最大逐格差 %.3g kg/kg）' % mx)
    P('\n## 烘干时长汇总')
    P('| 问题 | B(本文) | A(工作区) | 绝对差 | 相对差 |')
    P('|---|---|---|---|---|')
    P('| 问题3 | 205809.78 s = 57.1694 h | 205783.21 s = 57.1620 h | 26.6 s | 0.013% |')
    P('| 问题4(仿射) | 182954.56 s = 50.8207 h | 182951.26 s = 50.8198 h | 3.3 s | 0.0018% |')
    P('| 问题4(非均匀收缩) | 171236.5 s = 47.5657 h（N=800，外推 47.606 h） | '
      '171434.71 s = 47.6208 h | 198 s | 0.116% |')
    with open('_compare_AB.md', 'w', encoding='utf-8') as f:
        f.write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
