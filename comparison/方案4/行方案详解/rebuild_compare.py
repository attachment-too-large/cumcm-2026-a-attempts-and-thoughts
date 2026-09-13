# -*- coding: utf-8 -*-
"""只重建对比检验报告的 LaTeX 版。"""
import io
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import md2tex

HERE = os.path.dirname(os.path.abspath(__file__))
TEXDIR = os.path.join(HERE, 'latex')
NAME = '中药材烘干_两方案对比检验报告_LaTeX版'

md = io.open(os.path.join(HERE, '报告', '10_对比报告.md'), encoding='utf-8').read()
io.open(os.path.join(TEXDIR, 'body_compare.tex'), 'w', encoding='utf-8',
        newline='\n').write(md2tex.convert(md))

TEX = r'''% !TeX program = xelatex
\documentclass[12pt,a4paper]{ctexart}
\usepackage[a4paper,top=2.4cm,bottom=2.4cm,left=2.3cm,right=2.3cm]{geometry}
\usepackage[hidelinks]{hyperref}
\usepackage{amsmath,amssymb,bm}
\usepackage{booktabs,array,longtable,makecell}
\usepackage{graphicx}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{caption}
\usepackage{enumitem}
\definecolor{lstbg}{RGB}{248,248,248}
\lstset{language=Python, backgroundcolor=\color{lstbg},
  basicstyle=\ttfamily\fontsize{7}{8.4}\selectfont,
  numbers=left, numberstyle=\tiny\color{gray}, numbersep=5pt, frame=single,
  breaklines=true, columns=flexible, keepspaces=true,
  showstringspaces=false, extendedchars=true, tabsize=4}
\captionsetup{font={small,bf},labelsep=quad}
\title{\bfseries 中药材烘干问题：两套独立方案的对比检验报告\\[4pt]
\large 方案 X（节点中心有限体积 $+$ 一阶惯性环境）\ \ vs\ \ 方案 B（分级谱元弱形式 $+$ 分段线性环境）}
\author{}\date{2026 年 9 月}
\begin{document}
\maketitle
\input{body_compare.tex}
\end{document}
'''
io.open(os.path.join(TEXDIR, NAME + '.tex'), 'w', encoding='utf-8',
        newline='\n').write(TEX)

for k in (1, 2):
    with io.open(os.path.join(TEXDIR, '_cmp_pass%d.txt' % k), 'w', encoding='utf-8') as fh:
        subprocess.run(['xelatex', '-interaction=nonstopmode', NAME + '.tex'],
                       cwd=TEXDIR, stdout=fh, stderr=subprocess.STDOUT, timeout=1200)

log = os.path.join(TEXDIR, NAME + '.log')
nerr = 0
if os.path.exists(log):
    nerr = sum(1 for l in io.open(log, encoding='utf-8', errors='replace') if l.startswith('!'))
pdf = os.path.join(TEXDIR, NAME + '.pdf')
import fitz
out = '错误数=%d  页数=%s  大小=%s' % (
    nerr, len(fitz.open(pdf)) if os.path.exists(pdf) else 'NA',
    ('%d KB' % (os.path.getsize(pdf) // 1024)) if os.path.exists(pdf) else 'NA')
io.open(os.path.join(TEXDIR, '_cmp_summary.txt'), 'w', encoding='utf-8').write(out + '\n')
