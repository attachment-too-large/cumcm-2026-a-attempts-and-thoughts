# -*- coding: utf-8 -*-
# ==============================================================================
# make_controlled_comparison.py —— 受控对照：把"物性效应"与"收缩效应"拆开
# ==============================================================================
# 做什么：从问题 3 到问题 4，物性关系（附录 3 -> 附录 4）与几何（半径固定 -> 实测收缩）
#         同时改变，两者直接相减只能得到"净差异"，无法归因。本脚本补算一个受控算例
#         ——物性仍取附录 4，但半径固定为 R_0 不收缩——于是可以把总差异拆成两段纯效应：
#           A 到 B：只换物性（附录 3 -> 附录 4），半径都固定 -> 物性效应
#           B 到 C：只加收缩，物性都用附录 4             -> 收缩效应
#         A 到 C 即两问的净差异，画在同一张图上便于逐段对照。
# 输入：不读文件。算例由 spectral.HerbSolver 求解（其几何、环境标定与物性系数取自
#       common.py：R0、h、hm、T_air(t)、C_air(t) 与 props() 的附录 4 关系）；
#       基准值 TD_PROBLEM3、TD_PROBLEM4_SHRINK 是问题 3、4 的烘干时刻记录值。
# 输出：figure/fig20_p4_controlled.png   受控对照图，240 dpi，左右两面板
# 关键变量：
#   CRIT  烘干判据阈值 0.15 kg/kg（干基含水率），与其余脚本口径一致
#   n     谱节点数；正式出图取 64，另跑一次 48 做网格敏感性对照
#   t_dry 由事件求交得到的烘干时刻 [s]
# ==============================================================================

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from spectral import HerbSolver


# ---------------- 算例常数与基准值 ----------------
CRIT = 0.15                 # kg/kg（干基）烘干判据阈值
TD_PROBLEM3 = 205575.5      # s   问题 3：附录 3 物性、半径固定不收缩
TD_PROBLEM4_SHRINK = 182778.3   # s  问题 4：附录 4 物性、含实测收缩 R(t)


def fixed_radius_time(n=64):
    """补算"固定半径 + 附录 4 物性"这个受控算例，返回烘干时刻 t_dry [s]。

    只把 HerbSolver 的 shrink 置 False（半径恒为 R_0 = 2 cm），物性仍取 prob=4 的附录 4
    关系，于是它与问题 3 的差别就只剩物性一项。事件定义、容差与积分器（BDF）都与其余
    脚本一致，故算出的差值不是求解设置造成的。
    """
    solver = HerbSolver(n, prob=4, shrink=False)
    # atol 按状态分块：前 Np 个是水分浓度 C [kg/kg]（量级 1），后 Np 个是温度 T [degC]
    # （量级 3e2）。两块量级不同，同一个绝对容差对两者并不等价。
    atol = np.concatenate(
        [np.full(solver.Np, 1.0e-13), np.full(solver.Np, 1.0e-11)]
    )

    def event(t, y):
        """事件函数：y 前半段是水分浓度场 [kg/kg]，取其最大值（干燥最慢的中心）减判据。"""
        return np.max(y[: solver.Np]) - CRIT

    event.terminal = True       # 命中阈值立即停止积分，不必白算到 600000 s
    event.direction = -1        # 只认由大变小穿越，避免初值恰好等于阈值时误判
    # 600000 s（约 167 h）是足够宽的上界，事件必定在此之前触发
    sol = solver.solve(
        600000.0, rtol=1.0e-10, atol=atol, events=event, method="BDF"
    )
    if not sol.t_events[0].size:
        raise RuntimeError("Drying threshold was not reached.")
    return float(sol.t_events[0][0])


def draw(output_path, td_fixed):
    """画受控对照图并存到 output_path，返回 None。

    (a) 横轴为 A/B/C 三个算例、纵轴 烘干时长 t_dry/h，柱顶标出小时数：A = 问题 3
        （附录 3 物性、不收缩），B = 固定半径 + 附录 4 物性（本脚本补算，即入参
        td_fixed [s]），C = 问题 4（附录 4 物性、含实测收缩）。
    (b) 横轴 相对变化/%、纵轴三根水平柱，柱端标百分数：A 到 B 是纯物性效应、B 到 C 是
        纯收缩效应、A 到 C 是净差异。
    配色语义固定（两面板一致，便于横向对照）：#287271 深青 = A 问题 3 基准、
    #6C757D 灰 = B 固定半径附录 4、#D95F43 砖红 = C 问题 4 含收缩；(b) 中 A 到 B 用灰、
    B 到 C 用砖红、A 到 C 用深青。
    """
    plt.rcParams.update(
        {
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "DengXian"],
            "axes.unicode_minus": False,
            "font.size": 14,
            "axes.titlesize": 15,
            "axes.labelsize": 14,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.linewidth": 1.1,
            "savefig.dpi": 240,
        }
    )

    labels = ["A", "B", "C"]        # 代号：A=问题3，B=固定半径+附录4，C=问题4含收缩
    seconds = np.array([TD_PROBLEM3, td_fixed, TD_PROBLEM4_SHRINK])
    hours = seconds / 3600.0
    colors = ["#287271", "#6C757D", "#D95F43"]

    # 画布按 0.86\textwidth = 5.42 in 排版：取 8.2 in 宽时 14 pt 的正文字号印后只剩
    # 6.98 pt，低于 7 pt 的可读下限；缩到 8.0 in 可到 7.16 pt，故宽度定为 8.0 in。
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.0, 3.27))

    bars = ax1.bar(np.arange(3), hours, width=0.62, color=colors, edgecolor="white")
    ax1.bar_label(bars, labels=[f"{v:.2f} h" for v in hours], padding=4, fontsize=12)
    ax1.set_xticks(np.arange(3), labels)
    ax1.set_ylabel("烘干时长 / h")
    ax1.set_title("(a) 烘干时长")
    ax1.set_ylim(0, 145)
    ax1.grid(axis="y", alpha=0.28, linewidth=0.8)
    ax1.spines[["top", "right"]].set_visible(False)

    # 三段效应的相对变化 [%]：正号表示烘干时长变长、负号变短，A 到 C 是净差异
    effects = np.array([125.8, -60.6, -11.1])
    effect_labels = [
        "A 到 B  物性",
        "B 到 C  收缩",
        "A 到 C  净差",
    ]
    y = np.arange(3)
    effect_colors = [colors[1], colors[2], colors[0]]
    ebar = ax2.barh(y, effects, height=0.56, color=effect_colors)
    ax2.axvline(0, color="#333333", linewidth=1.2)
    ax2.set_yticks(y, effect_labels)
    ax2.invert_yaxis()              # 反转纵轴，让 "A 到 B" 一行显示在最上方
    ax2.set_xlim(-75, 140)
    ax2.set_xlabel("相对变化 / %")
    ax2.set_title("(b) 效应分解")
    ax2.grid(axis="x", alpha=0.28, linewidth=0.8)
    ax2.spines[["top", "right"]].set_visible(False)
    for rect, value in zip(ebar, effects):
        # 柱端百分数一律写在柱内：正向长柱靠右端内对齐，负向柱靠左端内对齐，
        # 绝对值小于 30 的短柱若也写白色就会压在色块外，故改用深灰并移出一点
        if value > 30:
            x, ha, color = value - 5, "right", "white"
        elif value < -30:
            x, ha, color = value + 5, "left", "white"
        else:
            x, ha, color = value - 4, "right", "#222222"
        ax2.text(
            x,
            rect.get_y() + rect.get_height() / 2,
            f"{value:+.1f}%",
            va="center",
            ha=ha,
            color=color,
            fontsize=12,
        )

    # 手工设定边距：wspace=0.34 是留给 (b) 三条较长的纵轴标签（"A 到 B  物性" 等），
    # 否则它们会顶到 (a) 的柱子上
    fig.subplots_adjust(left=0.085, right=0.985, top=0.88, bottom=0.19, wspace=0.34)
    fig.savefig(output_path, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)


def main():
    """跑两档网格的固定半径算例、画出受控对照图，并把三个数值打到屏幕上。

    两档网格只用来检查空间离散够不够细：rel_gap 应远小于要归因的那几个百分点，
    出图用 N=64 的结果。
    """
    td48 = fixed_radius_time(48)
    td64 = fixed_radius_time(64)
    rel_gap = abs(td64 - td48) / td64
    out = Path(__file__).resolve().parents[1] / "figure" / "fig20_p4_controlled.png"
    draw(out, td64)
    print(f"N=48 fixed-R Appendix-4: {td48:.3f} s")
    print(f"N=64 fixed-R Appendix-4: {td64:.3f} s")
    print(f"relative grid gap: {rel_gap:.3e}")
    print(f"figure: {out}")


if __name__ == "__main__":
    main()
