# -*- coding: utf-8 -*-
"""把结果部分改为 subsubsection，并对齐 §7 的表述（自验数据 / 明确归属）"""
import io
import os

P = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'paper.tex')
s = io.open(P, encoding='utf-8').read()

REP = [
    (r'\textbf{问题 1}（预热平衡段 0—1800 s）：',
     '\\subsubsection{问题 1：预热平衡段（0—1800 s）}'),
    (r'\textbf{问题 2}（整个烘干过程的前 3 小时）：',
     '\\subsubsection{问题 2：整个烘干过程的前 3 小时}'),
    (r'\textbf{问题 3}：以“域内最大水分浓度（中心）降到 0.15 kg/kg”为判据，得',
     '\\subsubsection{问题 3：烘干时长}\n以“域内最大水分浓度（中心）降到 0.15 kg/kg”为判据，得'),
    (r'\textbf{问题 4}：采用均匀仿射收缩并使用附录 4 物性与附件 2 的 $R(t)$，得',
     '\\subsubsection{问题 4：含收缩的烘干时长}\n'
     '采用均匀仿射收缩并使用附录 4 物性与附件 2 的 $R(t)$，得'),
    # §7.3 明确归属
    ('（同题的另一份独立实现用二维轴对称模型复核，中心点偏差 $<1\\%$，与上述量级估计一致。）',
     '（需要说明：本节结论来自量级估计，本文未独立建立二维模型；'
     '同题的另一份独立实现用二维轴对称模型复核得到中心点偏差 $<1\\%$，'
     '与上述估计一致，但该数值来自对方代码，本文未复现。）'),
]
cnt = 0
for a, b in REP:
    if a in s:
        s = s.replace(a, b, 1)
        cnt += 1
    else:
        print('NOT FOUND:', a[:46])
io.open(P, 'w', encoding='utf-8').write(s)
print('replaced:', cnt)
