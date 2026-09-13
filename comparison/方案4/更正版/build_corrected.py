# -*- coding: utf-8 -*-
"""构建更正版报告 PDF（正文 + 附录：代码更正补丁）。"""
import io
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(r'C:\Users\qing1\Desktop\2026CUMCM_', '行方案详解'))
import md2tex

HERE = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.join(HERE, 'latex')
os.makedirs(TEX, exist_ok=True)

md = io.open(os.path.join(HERE, '报告', '00_更正版解答与自查报告.md'), encoding='utf-8').read()
io.open(os.path.join(TEX, 'body_correct.md'), 'w', encoding='utf-8', newline='\n').write(md)
io.open(os.path.join(TEX, 'body_correct.tex'), 'w', encoding='utf-8',
        newline='\n').write(md2tex.convert(md))

patch = io.open(os.path.join(HERE, '补丁', '代码更正补丁.md'), encoding='utf-8').read()
patch = patch.replace('# 代码更正补丁（含验证）', '\\section{代码更正补丁（含验证）}')
io.open(os.path.join(TEX, 'body_patch.tex'), 'w', encoding='utf-8',
        newline='\n').write(md2tex.convert(patch))

MAIN = r'''% !TeX program = xelatex
\documentclass[12pt,a4paper]{ctexart}
\usepackage[a4paper,top=2.4cm,bottom=2.4cm,left=2.3cm,right=2.3cm]{geometry}
\usepackage[hidelinks]{hyperref}
\usepackage{amsmath,amssymb,bm}
\usepackage{booktabs,array,longtable,makecell}
\usepackage{graphicx,listings,xcolor,caption,enumitem}
\definecolor{lstbg}{RGB}{248,248,248}
\lstset{language=Python,backgroundcolor=\color{lstbg},
  basicstyle=\ttfamily\fontsize{7}{8.4}\selectfont,
  numbers=left,numberstyle=\tiny\color{gray},numbersep=5pt,frame=single,
  breaklines=true,columns=flexible,keepspaces=true,
  showstringspaces=false,extendedchars=true,tabsize=4}
\captionsetup{font={small,bf},labelsep=quad}
\title{\bfseries 更正版解答与自查报告\\[6pt]
\large 2026 CUMCM A 题「药材的烘干问题」\\[2pt]
\normalsize 按逐式审阅结论做出的更正 · 正确性验证 · 与答案 X 的偏离归因及自我反思}
\author{}\date{2026 年 9 月}
\begin{document}
\maketitle
\input{body_correct.tex}
\appendix
\input{body_patch.tex}
\end{document}
'''
io.open(os.path.join(TEX, 'corrected_report.tex'), 'w', encoding='utf-8',
        newline='\n').write(MAIN)

for k in (1, 2):
    with io.open(os.path.join(TEX, '_pass%d.txt' % k), 'w', encoding='utf-8') as fh:
        subprocess.run(['xelatex', '-interaction=nonstopmode', 'corrected_report.tex'],
                       cwd=TEX, stdout=fh, stderr=subprocess.STDOUT, timeout=1200)

log = os.path.join(TEX, 'corrected_report.log')
nerr = sum(1 for l in io.open(log, encoding='utf-8', errors='replace') if l.startswith('!'))
import fitz
pdf = os.path.join(TEX, 'corrected_report.pdf')
res = '错误数=%d  页数=%s  大小=%s' % (
    nerr,
    len(fitz.open(pdf)) if os.path.exists(pdf) else 'NA',
    ('%d KB' % (os.path.getsize(pdf) // 1024)) if os.path.exists(pdf) else 'NA')
io.open(os.path.join(TEX, '_summary.txt'), 'w', encoding='utf-8').write(res + '\n')
