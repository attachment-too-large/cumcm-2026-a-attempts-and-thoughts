# -*- coding: utf-8 -*-
"""修正 paper.tex：表题手动编号（题目要求的表1—表6 必须对应正确）+ 补插图"""
import io
import os

P = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'paper.tex')
s = io.open(P, encoding='utf-8').read()

REP = [
    (r'\caption{题给的密度、比热容与导热系数}',
     r'\caption*{表 A\quad 题给的密度、比热容与导热系数}'),
    (r'\caption{扩散系数两种读法的对比（问题 1，单位 m$^2$/s）}',
     r'\caption*{表 B\quad 扩散系数两种读法的对比（问题 1，单位 m$^2$/s）}'),
    (r'\caption{30 分钟内药材的温度（单位：$^\circ$C）}',
     r'\caption*{表 1\quad 30 分钟内药材的温度（单位：$^\circ$C）}'),
    (r'\caption{30 分钟内药材的水分浓度（单位：kg/kg）}',
     r'\caption*{表 2\quad 30 分钟内药材的水分浓度（单位：kg/kg）}'),
    (r'\caption{3 小时内药材的温度（单位：$^\circ$C）}',
     r'\caption*{表 3\quad 3 小时内药材的温度（单位：$^\circ$C）}'),
    (r'\caption{3 小时内药材的水分浓度（单位：kg/kg）}',
     r'\caption*{表 4\quad 3 小时内药材的水分浓度（单位：kg/kg）}'),
    (r'\caption{药材烘干过程的水分浓度（单位：kg/kg）}' + '\n' + r'\label{tab:5}',
     r'\caption*{表 5\quad 药材烘干过程的水分浓度（单位：kg/kg）}'),
    (r'\caption{药材烘干过程的水分浓度（单位：kg/kg）}' + '\n' + r'\label{tab:6}',
     r'\caption*{表 6\quad 药材烘干过程的水分浓度（单位：kg/kg）}'),
    (r'\caption{与独立实现的逐格比对（共 216 个输出格子）}',
     r'\caption*{表 C\quad 与独立实现的逐格比对（共 216 个输出格子）}'),
    (r'\caption{问题 3 烘干时长的不确定度分解}',
     r'\caption*{表 D\quad 问题 3 烘干时长的不确定度分解}'),
    (r'见表 \ref{tab:props}', '见表 A'),
    (r'表 \ref{tab:Dread} 给出', '表 B 给出'),
    (r'结果见表 \ref{tab:cmp}', '结果见表 C'),
]
cnt = 0
for a, b in REP:
    if a in s:
        s = s.replace(a, b, 1)
        cnt += 1
    else:
        print('NOT FOUND:', a[:50])

# 去掉失效的 label（caption* 不再产生编号）
s = s.replace(r'\label{tab:5}', '').replace(r'\label{tab:6}', '')
s = s.replace(r'\label{tab:props}', '').replace(r'\label{tab:Dread}', '')
s = s.replace(r'\label{tab:cmp}', '').replace(r'\label{tab:sens}', '')

# 插入扩散系数对比图（放在 §5.4 判读段落前）
anchor = r'\noindent 这一差异无法通过任何收敛性检验'
FIG = ('\\begin{figure}[htbp]\n\\centering\n'
       '\\includegraphics[width=0.95\\textwidth]{figs/fig2_D.png}\n'
       '\\caption{扩散系数对含水率的依赖：正确读法（分式指数）与误读（乘积指数）}\n'
       '\\label{fig:D}\n\\end{figure}\n\n')
if anchor in s and 'figs/fig2_D.png' not in s:
    s = s.replace(anchor, FIG + anchor, 1)
    print('fig2_D inserted')

io.open(P, 'w', encoding='utf-8').write(s)
print('captions replaced:', cnt)
print('remaining \\caption{ (numbered):', s.count(r'\caption{'))
print('caption* count:', s.count(r'\caption*{'))
