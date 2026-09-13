# -*- coding: utf-8 -*-
"""把 2D 端面复核与 M1 复现结果写入 paper.tex，并把结论口径改为主答案+区间"""
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, 'paper.tex')
s = io.open(P, encoding='utf-8').read()

# ---------- ① §7.2(2) 收缩模型：改为本文已复现 ----------
old2 = s[s.find('\\textbf{（2）收缩模型。}'):s.find('\\textbf{（3）能量口径。}')]
new2 = r"""\textbf{（2）收缩模型。}
问题 4 采用均匀仿射收缩 $r=\sigma R(t)$。若允许用附录 4 的 $\rho(C)$ 补充\textbf{局部}
比容信息，可建立非均匀收缩自洽模型：以干物质坐标 $s$ 为自变量（每单位 $s$ 含相同干物质质量），
局部比容 $v(C)=(1+C)/\rho(C)$，$I=\int_0^1v\,\ud s$，$J(s)=\int_0^sv\,\ud\xi/I$，
物理半径 $r=R(t)\sqrt{J(s)}$，控制方程为
\begin{equation}
\frac{\partial C}{\partial t}=\frac{4I^2}{R^2}\frac{\partial}{\partial s}
\!\left(\frac{JD}{v^2}\frac{\partial C}{\partial s}\right),\qquad
(1+C)c_p\frac{\partial T}{\partial t}=\frac{4I^2}{R^2}\frac{\partial}{\partial s}
\!\left(\frac{Jk}{v}\frac{\partial T}{\partial s}\right).
\label{eq:m1}
\end{equation}
本文用干物质坐标下的三对角有限体积（系数含 $4I^2/R^2$、面权重 $J_f$、边界 $2Ih_m/(Rv_N)$）
独立实现了该模型，结果如下：

\begin{center}
\small
\begin{tabular}{@{}lcccc@{}}
\toprule
网格 / 时间步 & $N=200$ & $N=400$ & $N=800$ & $N=400$, $\ud t=30$ s \\
\midrule
$t_{dry}$ / h & 47.4423 & 47.5250 & 47.5657 & 47.5338 \\
\bottomrule
\end{tabular}
\end{center}

\noindent 按一阶 Richardson 外推得 $t_{dry}\to47.606$ h，与另一份独立实现给出的
$47.623$ h（$N=800$）相差 $0.04\%$，\textbf{该模型已被独立复现}。
另做了单元测试：把比容强制为常数后式 \eqref{eq:m1} 应退化为均匀仿射模型，
实测与谱元仿射解吻合到 $10^{-4}\sim10^{-3}$（仅剩有限体积自身的离散误差），
说明该实现正确。

因此问题 4 的答案对收缩建模假设敏感：均匀仿射 $50.82$ h 对非均匀自洽 $47.61$ h，
相差 $-6.3\%$，该不确定度已被两份独立实现确认，是问题 4 的主要不确定度来源。

"""
s = s.replace(old2, new2)

# ---------- ② §7.3 几何假设：改为二维数值验证 ----------
i3 = s.find('\\subsection{几何假设可靠性}')
i4 = s.find('\\subsection{数据与参数可靠性}')
new3 = r"""\subsection{几何假设可靠性：二维轴对称复核}

为检验“无限长圆柱（一维）”假设，本文建立了\textbf{二维轴对称有限体积模型}：
域为 $r\in[0,R]$、$z\in[0,L/2]$（$z=0$ 为对称面，等价于整根长 $L=25$ cm），
$r=R$ 与 $z=L/2$ 为对流边界，五点格式 + $\theta$ 格式 + 稀疏直接解。

需要强调\textbf{实验设计}：若直接把二维解与一维解比较，会把二维自身的离散误差
（在 $N_r=80$ 时 $t_{dry}$ 仍有约 $1.2\times10^3$ s 的偏差）误当作端面效应。
因此本文采用\textbf{同一网格、仅切换 $z=L/2$ 边界条件（对流 / 绝热）}的方式隔离端面效应：

\begin{center}
\small
\begin{tabular}{@{}lcccc@{}}
\toprule
网格（$\ud t=60$ s） & $\Delta C_{center}(6\ \text{h})$ & $\Delta C_{center}(24\ \text{h})$
& $\Delta C_{center}(48\ \text{h})$ & $\Delta t_{dry}$ \\
\midrule
$N_r=40,\ N_z=25$ & $-1.04\times10^{-5}$ & $-4.4\times10^{-6}$ & $-2.8\times10^{-6}$ & $-8.2$ s \\
$N_r=60,\ N_z=40$ & $-1.01\times10^{-5}$ & $-3.9\times10^{-6}$ & $-2.6\times10^{-6}$ & $-7.7$ s \\
$N_r=80,\ N_z=50$ & $-1.00\times10^{-5}$ & $-3.8\times10^{-6}$ & $-2.5\times10^{-6}$ & $-7.6$ s \\
\bottomrule
\end{tabular}
\end{center}

\noindent 结论：端面效应使 $t_{dry}$ 缩短约 $7.7$ s（$\mathbf{-0.004\%}$），
中心含水率差异不超过 $1\times10^{-5}$，且该差值在三套网格上稳定（$-8.2/-7.7/-7.6$ s）。
因此\textbf{一维无限长圆柱假设在本题判据下成立，误差小于 $0.01\%$}。
这与量级估计（轴向扩散穿透深度 $\sqrt{Dt}$ 最大约 4 cm，小于半长 12.5 cm）一致。

作为该二维求解器的可靠性检查：在“长柱 + 端面绝热”的一维极限下，
二维解与一维谱元解的径向剖面最大偏差为 $1.5\times10^{-4}$（600 s）、
$6.7\times10^{-5}$（3600 s）、$2.0\times10^{-5}$（21600 s），
且随网格加密单调减小（$N_r=40/80/120$ 下 6 h 中心值变化仅 $5\times10^{-6}$）。

"""
s = s[:i3] + new3 + s[i4:]

# ---------- ③ §7.5 结论：主答案 + 区间 ----------
s = s.replace(r"""  \item 在\textbf{标准读法}（式 \eqref{eq:gov}，$\rho_d$ 均匀）下，问题 3 的时长为
        $57.17$ h；在严格守恒形式 \eqref{eq:cons_form} 下为 $55.27$ h。
        两者差异 $3.3\%$，本文以 $57.2$ h 为基准值，并给出 $55.3$—$57.2$ h 的区间。""",
r"""  \item 在\textbf{标准读法}（式 \eqref{eq:gov}，$\rho_d$ 均匀）下，问题 3 的时长为
        $57.17$ h；在严格守恒形式 \eqref{eq:cons_form} 下为 $55.27$ h。
        两者差异 $3.3\%$。本文以 $\mathbf{57.2\ h}$ 作为基准答案
        （理由：与题面“固定几何”的前提自洽，且与独立实现逐格对拍一致），
        并给出 $\mathbf{55.3\text{—}57.2\ h}$ 作为模型形式不确定度区间。""")
s = s.replace(r"""  \item 问题 4 的结论对收缩建模假设敏感（仿射 $50.82$ h 对非均匀 $47.6$ h，$6.3\%$），
        该不确定度大于问题 3。""",
r"""  \item 问题 4 的结论对收缩建模假设敏感（仿射 $50.82$ h 对非均匀 $47.61$ h，$6.3\%$），
        该不确定度大于问题 3；两种收缩模型本文均已独立实现并复现。
  \item \textbf{几何假设已获数值验证}：二维轴对称复核表明端面效应仅约 $-0.004\%$，
        一维模型的误差小于 $0.01\%$，不再是不确定度来源。""")

io.open(P, 'w', encoding='utf-8').write(s)
print('patched; len =', len(s))
