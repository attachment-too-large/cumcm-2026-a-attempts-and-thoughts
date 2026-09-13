# -*- coding: utf-8 -*-
"""按竞赛格式规范改造 paper.tex：
   ① 删除目录（规范第四条：正文“不要目录”）
   ② 附录放入全部完整、可运行的源程序（规范第五条）
"""
import glob
import io
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, 'methodB_spectral')
CODE = os.path.join(HERE, 'code')
os.makedirs(CODE, exist_ok=True)

CORE = ['hb_core.py', 'sem_core.py', 'run_all.py']
VERIFY = ['selftest.py', 'test_sem.py', 'compare_AB.py', 'credibility.py',
          'sens_extra.py', 's3_plateau.py', 'e2b_surface.py', 'audit_energy.py',
          'analyze_shrinkage.py', 'analyze_B.py',
          'fv2d.py', 'fv2d_run.py', 'm1_fv.py', 'm1_run.py']
TOOLS = ['make_figs.py', 'make_tables.py']

for f in CORE + VERIFY + TOOLS:
    if f in CORE or f in VERIFY:
        s = os.path.join(SRC, f)
    elif f == 'make_tables.py':
        s = os.path.join(HERE, f)
    else:
        s = os.path.join(ROOT, 'paper_tools', f)
    if os.path.exists(s):
        shutil.copy(s, os.path.join(CODE, f))
    else:
        print('MISSING', f)

# ---- 组装附录 LaTeX ----
parts = [r'\newpage',
         r'\section*{附录 A\quad 支撑材料文件列表}',
         r'\addcontentsline{toc}{section}{附录 A}', '',
         r'按竞赛规范第五条，本附录先列出支撑材料的文件列表，随后给出全部完整、'
         r'可运行的源程序代码。支撑材料压缩包 \texttt{A题\_支撑材料.zip} 的目录结构如下：', '',
         r'\begin{center}', r'\small',
         r'\begin{tabular}{@{}p{5.0cm}p{8.6cm}@{}}', r'\toprule',
         r'路径 & 内容 \\', r'\midrule',
         r'\texttt{源程序/核心/} & \texttt{hb\_core.py}（模型与物性）、'
         r'\texttt{sem\_core.py}（分级谱元弱形式求解器）、'
         r'\texttt{run\_all.py}（问题 1—4 生产计算与结果输出） \\',
         r'\texttt{源程序/验证/} & \texttt{selftest.py}、\texttt{test\_sem.py}（数值验证）、'
         r'\texttt{compare\_AB.py}（跨实现逐格比对）、\texttt{credibility.py}、'
         r'\texttt{sens\_extra.py}、\texttt{s3\_plateau.py}（可信度实验）、'
         r'\texttt{fv2d.py}、\texttt{fv2d\_run.py}（二维轴对称复核）、'
         r'\texttt{m1\_fv.py}、\texttt{m1\_run.py}（非均匀收缩模型复现）、'
         r'\texttt{audit\_energy.py}、\texttt{analyze\_B.py}、'
         r'\texttt{analyze\_shrinkage.py}（能量审计与收缩一致性） \\',
         r'\texttt{源程序/论文与图表/} & \texttt{make\_figs.py}、\texttt{make\_paper.py}、'
         r'\texttt{make\_tables.py}、\texttt{md2pdf.py} 及论文 \LaTeX{} 源 \\',
         r'\texttt{结果文件/} & \texttt{result1}—\texttt{result4}.xlsx、'
         r'\texttt{result4\_material\_frame.xlsx}、\texttt{summary\_B.json}、'
         r'逐格比对表 \\',
         r'\texttt{中间结果/} & 论文全部插图、公式原版面核对截图 \\',
         r'\texttt{实验日志/} & 可信度复核各实验的原始输出文本 \\',
         r'\texttt{赛题原始数据/} & 附件 1、附件 2、附件 3（主办方发布，便于复现，版权归主办方） \\',
         r'\bottomrule', r'\end{tabular}', r'\end{center}', '',
         r'\noindent 运行环境：Python 3.9 及以上，依赖 NumPy / SciPy / openpyxl；'
         r'论文排版为 XeLaTeX{} + ctex。复现步骤见包内 \texttt{README.md}。', '',
         r'\newpage',
         r'\section*{附录 B\quad 核心程序（模型、求解器与结果生成）}',
         r'\addcontentsline{toc}{section}{附录 B}', '',
         r'以下为本文建模与求解的全部核心源程序，可直接运行。', '']
for f in CORE:
    cap = {'hb_core.py': '模型与物性参数（hb\_core.py）',
           'sem_core.py': '分级谱元弱形式求解器（sem\_core.py）',
           'run_all.py': '问题 1—4 的生产计算与结果文件输出（run\_all.py）'}[f]
    parts += [r'\subsection*{%s}' % cap,
              r'\lstinputlisting[language=Python]{%s}' % ('code/' + f), '']

parts += [r'\newpage',
          r'\section*{附录 C\quad 验证与可信度复核程序}',
          r'\addcontentsline{toc}{section}{附录 C}', '',
          r'以下程序对应第 6 章（数值验证）与第 7 章（可信度复核）的全部实验，'
          r'其中 \texttt{fv2d.py} / \texttt{fv2d\_run.py} 为二维轴对称复核，'
          r'\texttt{m1\_fv.py} / \texttt{m1\_run.py} 为非均匀收缩模型的独立复现。', '']
for f in VERIFY:
    parts += [r'\subsection*{%s}' % f.replace('_', r'\_'),
              r'\lstinputlisting[language=Python]{code/%s}' % f, '']

parts += [r'\newpage',
          r'\section*{附录 D\quad 图表生成程序}',
          r'\addcontentsline{toc}{section}{附录 D}', '']
for f in TOOLS:
    parts += [r'\subsection*{%s}' % f.replace('_', r'\_'),
              r'\lstinputlisting[language=Python]{code/%s}' % f, '']

appendix = '\n'.join(parts)

P = os.path.join(HERE, 'paper.tex')
s = io.open(P, encoding='utf-8').read()

# ① 删目录
toc_old = '\\newpage\n\n% ==================== 目录（可选） ====================\n\\tableofcontents\n\\newpage\n'
if toc_old in s:
    s = s.replace(toc_old, '\\newpage\n')
    print('TOC removed')
else:
    s2 = s.replace('\\tableofcontents', '')
    if s2 != s:
        s = s2
        print('TOC removed (fallback)')
    else:
        print('TOC NOT FOUND')

# ② 替换附录
i = s.find('% ==================== 附录 ====================')
if i < 0:
    raise SystemExit('appendix marker not found')
s = s[:i] + '% ==================== 附录 ====================\n' + appendix + '\n\\end{document}\n'

io.open(P, 'w', encoding='utf-8').write(s)
print('appendix assembled,', len(appendix), 'chars')
