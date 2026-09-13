# -*- coding: utf-8 -*-
"""生成两份 LaTeX 文档并编译（xelatex）。"""
import io
import os
import subprocess
import sys

_ORIG_STDOUT = sys.stdout
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import md2tex

HERE = os.path.dirname(os.path.abspath(__file__))
TEXDIR = os.path.join(HERE, 'latex')
REP = os.path.join(HERE, '报告')
os.makedirs(TEXDIR, exist_ok=True)

PREAMBLE_COMMON = r"""
\usepackage{amsmath,amssymb,bm}
\usepackage{booktabs,array,longtable,makecell}
\usepackage{graphicx}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{caption}
\usepackage{enumitem}
\definecolor{lstbg}{RGB}{248,248,248}
\definecolor{lstcmt}{RGB}{0,128,0}
\definecolor{lstkw}{RGB}{0,0,200}
\lstset{language=Python, backgroundcolor=\color{lstbg},
  basicstyle=\ttfamily\fontsize{7}{8.4}\selectfont,
  keywordstyle=\color{lstkw}\bfseries, commentstyle=\color{lstcmt},
  numbers=left, numberstyle=\tiny\color{gray}, numbersep=5pt, frame=single,
  rulecolor=\color{gray!45}, breaklines=true, columns=flexible,
  keepspaces=true, showstringspaces=false, extendedchars=true, tabsize=4}
\captionsetup{font={small,bf},labelsep=quad}
\setlength{\parskip}{2pt}
"""

# ---------------- 文档 1：逐步推导详解（cumcmthesis） ----------------
TITLE = '中药材热风烘干问题的耦合扩散建模与逐步推导'
TIHAO = 'A'
ABSTRACT = r"""
中药材热风烘干是「热湿耦合传递 + 尺寸收缩」的典型问题。本文以圆柱形药材（长 25 cm、
半径 2 cm）为对象，建立温度—干基含水率耦合的轴对称扩散方程组，并给出\textbf{从物理守恒
到可编程实现的完整推导}。

\textbf{建模方面}，以环形壳层为控制体，由能量守恒与 Fourier 定律推得
$\rho c_p\partial_tT=r^{-1}\partial_r(rk\partial_rT)$；把控制体绑定在干物质上，由水分守恒与
Fick 定律推得 $\partial_tC=r^{-1}\partial_r(rD\partial_rC)$。表面取 Robin 条件，中心取对称
条件，初始条件由题面给出。引入无量纲数（$Bi=1.389$、$Bi_m=3.24$、$\mathrm{Le}=0.029$）
先做量级判断：热扩散特征时间 2369 s、水分扩散特征时间 $8.1\times10^4$ s，两者相差 34 倍，
预告了「温度显著上升、水分仅表层变化」的结果形态。

\textbf{收缩处理方面}，引入物质坐标 $\xi=r/R(t)$，证明几何相似收缩下对流项与坐标变换项
精确抵消，移动边界问题化为固定域纯扩散问题。本文同时指出该推导中「干物质守恒」与
「$\rho_d$ 由局部含水率决定」之间的内在张力，并给出一条只需单条假设、无逻辑缺口的替代
推法，两者最终得到同一个方程。

\textbf{数值方法方面}，采用节点中心有限体积离散（表面节点恰好落在物理边界上）与
$\theta$ 法时间推进（取 $\theta=1/2$ 为 Crank--Nicolson，前 4 步改用隐式 Euler 作
Rannacher 启动以抑制初值跳跃引起的振荡）。本文给出该格式守恒性的严格证明（望远镜求和），
以及放大因子 $g=(1+z/2)/(1-z/2)\to-1$ 的推导，说明为何必须做 Rannacher 启动。

\textbf{结果}：问题 1 在 1800 s 内温度整体升至 $33.6\sim36.8$ ℃、水分仅表层由 2.55 降至
1.5109；问题 2 在 3 h 末温度趋于均匀（中心 49.85 ℃）；问题 3 的烘干时长为
$205\,559\sim205\,810$ s（取决于环境数据处理口径）；问题 4 的收缩对照实验表明
\textbf{收缩使烘干时间缩短 60.6\%}，远大于物性表差异的影响。

\textbf{关键词}：热湿耦合扩散；干基含水率；Robin 边界；物质坐标；有限体积法；
Rannacher 启动；尺寸收缩
"""
doc1 = r'''% !TeX program = xelatex
\PassOptionsToPackage{hidelinks}{hyperref}
\documentclass[withoutpreface,bwprint]{cumcmthesis}
''' + PREAMBLE_COMMON + r'''
\title{''' + TITLE + r'''}
\tihao{''' + TIHAO + r'''}
\baominghao{}\schoolname{}\membera{}\memberb{}\memberc{}\supervisor{}
\yearinput{2026}\monthinput{9}\dayinput{}
\begin{document}
\maketitle
\begin{abstract}
''' + ABSTRACT + r'''
\end{abstract}
\input{body_detail.tex}
\end{document}
'''

# ---------------- 文档 2：对比检验报告（ctexart） ----------------
doc2 = r'''% !TeX program = xelatex
\documentclass[12pt,a4paper]{ctexart}
\usepackage[a4paper,top=2.4cm,bottom=2.4cm,left=2.3cm,right=2.3cm]{geometry}
\usepackage[hidelinks]{hyperref}
''' + PREAMBLE_COMMON + r'''
\title{\bfseries 中药材烘干问题：两套独立方案的对比检验报告\\[4pt]
\large 方案 X（节点中心有限体积 + 一阶惯性环境）\\ vs \ 方案 B（分级谱元弱形式 + 分段线性环境）}
\author{}\date{2026 年 9 月}
\begin{document}
\maketitle
\input{body_compare.tex}
\end{document}
'''

JOBS = [
    ('detail', doc1,
     ['00_导读与方程推导.md', '01_收缩与数值方法.md', '02_环境求解与验证.md'],
     '中药材烘干_行方案_逐步推导详解_LaTeX版'),
    ('compare', doc2,
     ['10_对比报告.md'],
     '中药材烘干_两方案对比检验报告_LaTeX版'),
]

SUMMARY = []
for tag, tex, mds, outname in JOBS:
    # 拼 Markdown 并转换
    body_md = '\n\n'.join(io.open(os.path.join(REP, f), encoding='utf-8').read() for f in mds)
    io.open(os.path.join(TEXDIR, 'body_%s.md' % tag), 'w', encoding='utf-8', newline='\n').write(body_md)
    body_tex = md2tex.convert(body_md)
    io.open(os.path.join(TEXDIR, 'body_%s.tex' % tag), 'w', encoding='utf-8', newline='\n').write(body_tex)
    io.open(os.path.join(TEXDIR, '%s.tex' % outname), 'w', encoding='utf-8', newline='\n').write(tex)
    # 编译两遍
    for k in (1, 2):
        with io.open(os.path.join(TEXDIR, '_%s_pass%d.txt' % (tag, k)), 'w', encoding='utf-8') as fh:
            subprocess.run(['xelatex', '-interaction=nonstopmode', '%s.tex' % outname],
                           cwd=TEXDIR, stdout=fh, stderr=subprocess.STDOUT, timeout=1200)
    pdf = os.path.join(TEXDIR, '%s.pdf' % outname)
    log = os.path.join(TEXDIR, '%s.log' % outname)
    nerr = 0
    if os.path.exists(log):
        for line in io.open(log, encoding='utf-8', errors='replace'):
            if line.startswith('!'):
                nerr += 1
    size = ('%.0f KB' % (os.path.getsize(pdf) / 1024)) if os.path.exists(pdf) else '0'
    SUMMARY.append('%-34s %-5s %-9s 错误=%d' % (outname, 'OK' if os.path.exists(pdf) else 'FAIL', size, nerr))

io.open(os.path.join(TEXDIR, '_build_summary.txt'), 'w', encoding='utf-8').write('\n'.join(SUMMARY) + '\n')
print('\n'.join(SUMMARY))
