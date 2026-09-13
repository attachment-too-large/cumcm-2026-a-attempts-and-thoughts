# -*- coding: utf-8 -*-
# ==============================================================================
# make_figures.py —— 生成论文插图（本脚本产出 25 张 PNG，其中 24 张被正文引用）
# ==============================================================================
# 做什么：按论文出现顺序逐张画出 25 幅插图：前 4 幅是附件 1/2 的数据处理与环境标定
#         诊断，中段是问题 1~4 的径向剖面、时空云图与烘干终点，末段是三层验证、

#         收敛性与误差分析，最后 1 幅是技术路线流程图。每幅图画的是哪个物理量、横纵轴
#         是什么、数据取自哪个文件的哪个键、各子图 (a)(b)(c) 分别是什么，都写在该段
#         代码上方的注释里，按注释对照图即可。
# 输入：中间数据/solutions.npz（四问解；缺失或版本/指纹不符时由 fig_core.build_cache()
#       自动重算）；中间数据/verification.json（收敛性与容差数据，fig22、fig23 使用，
#       由 verify_data.py --recompute 或 run_all.py 生成）；附件 1/2 原始序列经

#       common.py 载入。
# 输出：figure/fig01_env_fit.png … figure/fig25_flow.png 共 25 张，240 dpi。另有受控对照图
#       fig20_p4_controlled.png 由 make_controlled_comparison.py 单独生成，全包合计 26 张；
#       正文共引用 25 张，本脚本这 25 张里只有备用图 fig20_p4_vs_p3.png 未被引用——它把物性
#       与几何同时改变后的净差异画在一起，差值不能当作纯收缩效应。
# 用法：python make_figures.py（先跑过 run_all.py 或 make_results.py，再直接出图）

# 关键变量：
#   RCOLS  剖面采样的 5 个物理半径 [m]：0、0.005、0.010、0.015、0.020，即 0~2.0 cm
#   CRIT   烘干判据阈值 0.15 kg/kg（干基含水率）
#   FLUX_T_END  表面传质通量曲线的终了时刻 72 h = 259200 s
#   S1~S4  缓存解与插值接口打包成的轻量对象，交给 _samples() 复用
# ==============================================================================
import os
from console_utf8 import fix_console
fix_console()
import numpy as np
from scipy.special import j0, j1
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe                 # 等值线数值的白描边
from matplotlib.patches import FancyBboxPatch      # 技术路线图与示意图的圆角方框
from fig_core import (save, load_cache, build_cache, FIGDIR, DATADIR, DPI)
import fig_core as FC
import verify_data as VD

RCOLS = np.array([0.0, 0.005, 0.010, 0.015, 0.020])   # m，即 0、0.5、1.0、1.5、2.0 cm
CRIT = 0.15                                           # kg/kg（干基）烘干判据阈值

# ---------------- 全篇语义色表 ----------------
# 25 张图共用一套颜色语义：同一种颜色在任何一页上永远指同一种东西。这是"翻到哪页都
# 显得专业"的底层原因，也是读者不会看错图的保证。本表只作用于 make_figures.py；
# make_controlled_comparison.py 独立运行、不 import 本模块，那里另有一份同值的副本，
# 改动时必须同步——否则同一含义在两个脚本出的图里会是两种颜色。
COL_TEMP = "#c0392b"    # 温度：暖红
COL_MOIST = "#1f6fb2"   # 水分浓度：冷蓝
COL_CRIT = "#2c3e50"    # 判据 / 阈值线：深靛（一律点线）
COL_PHASE = "#b9770e"   # 阶段分界、环境设定值：赭黄（一律虚线）
COL_MEAS = "#8a8f98"    # 实测散点：中灰（数据是陪衬，主结论才上色）
COL_SIDE = "#9aa0a6"    # 陪衬曲线：浅灰


def _glow(ax, x, y, color, lw=2.0, glow=3.6, alpha=0.20, zorder=3, **kw):
    """画一条"发光主线"：先铺一层同色、宽而淡的底线，再压一条实线，返回实线的 Line2D。

    目的是让主结论曲线在一堆陪衬线里第一眼就跳出来，而不必把线加粗到失真。
    底层用 alpha 而不是预混浅色，是为了在深浅两种底色上都保持同一个色相。
    注意两层线都进图例时会出现重复条目，故主线的 label 只挂在实线上。
    """
    ax.plot(x, y, color=color, lw=glow, alpha=alpha, solid_capstyle="round",
            zorder=zorder - 1)
    return ax.plot(x, y, color=color, lw=lw, zorder=zorder, **kw)[0]


def _cross_time(t, v, crit):
    """求 v 由大变小、首次跌破 crit 的时刻 [s]，返回 float；始终大于 crit 则返回 nan。

    用相邻两点线性求交，而不是直接取第一个低于阈值的采样点：采样间隔在问题 3 里是
    百秒量级，直接取点会把达标时刻系统性地推迟若干步，图上的标注与正文就对不上了。
    """
    idx = np.where(np.asarray(v) <= crit)[0]
    if idx.size == 0:
        return np.nan
    k = int(idx[0])
    if k == 0:
        return float(t[0])
    return float(t[k - 1] + (v[k - 1] - crit) / (v[k - 1] - v[k]) * (t[k] - t[k - 1]))


def _samples(s, sol, t_idx=None, rcols=RCOLS):
    """把解按 u=(r/R)^2 插值到 rcols 给出的物理半径上，返回 (C, T) 两个 (nt, nr) 数组。

    C 为干基含水率 [kg/kg]、T 为温度 [degC]；nt 为取样时刻数、nr = len(rcols)。
    t_idx 给出要取的时刻下标（默认全部）；s 需提供 R_of、interp_u 与 Np 属性，
    形状取自传入的 sol（sol.y 的第 1 维是时刻、第 2 维是状态）。

    问题 4 的半径随收缩变小，固定的物理距离 r_j 可能在某个时刻已落到药材之外
    （r_j > R(t_i)）。该位置此时不存在——按论文 S7.4 与 make_results.py 的约定填
    NaN，而不能把 u_t=(r_j/R(t))^2 > 1 送去插值：那是对求解域之外的做外推，浮点下
    还会因重心插值分母抵消到 0 而返回 inf。
    """
    idx = range(sol.y.shape[1]) if t_idx is None else t_idx
    idx = list(idx)
    C = np.full((len(idx), len(rcols)), np.nan)
    T = np.full((len(idx), len(rcols)), np.nan)
    for k, i in enumerate(idx):
        Rt = s.R_of(sol.t[i])
        inside = rcols <= Rt + 1e-12
        if not np.any(inside):
            continue
        Cn, Tn = s.interp_u(sol.y[:, i], (rcols[inside] / Rt) ** 2)
        C[k, inside] = Cn
        T[k, inside] = Tn
    return C, T


def main():
    """按论文引用顺序生成 25 张插图，返回 (t_dry3, t_dry4) 两个烘干时刻 [s]。

    绘图语句全部写在本函数体内，按下面的分节注释一段一图。开头先把四问解从缓存取
    出来，再用 mk() 打包成只需 interp_u/R_of 的轻量对象交给 _samples()，这样采样
    逻辑只写一份，四个问题共用。
    """
    build_cache()
    z, HerbSolver, cm = load_cache()
    s1 = HerbSolver(64, prob=1)
    s2 = HerbSolver(64, prob=2)
    s3 = HerbSolver(64, prob=3)
    s4 = HerbSolver(64, prob=4, shrink=True)
    t1, y1 = z["p1_t"], z["p1_y"]
    t2, y2 = z["p2_t"], z["p2_y"]
    t3, y3 = z["p3_t"], z["p3_y"]
    t4, y4 = z["p4_t"], z["p4_y"]
    td3 = float(z["p3_tdry"][0])
    td4 = float(z["p4_tdry"][0])
    print("tdry3=%.3f s  tdry4=%.3f s" % (td3, td4))

    class _S:  # 轻量容器，复用 _samples 接口
        pass

    def mk(s, t, y):
        """把（求解器 s、缓存时刻 t [s]、解 y）打包成 _S，属性名与 _samples 的调用一致。"""
        o = _S(); o.y = y; o.t = t; o.R_of = s.R_of
        o.interp_u = s.interp_u; o.Np = s.Np; o.N = s.N; o.u = s.u
        return o

    S1, S2, S3, S4 = mk(s1, t1, y1), mk(s2, t2, y2), mk(s3, t3, y3), mk(s4, t4, y4)

    # ---------------- 第 1 组：附件数据的处理与环境标定诊断（fig01~fig04） ----------------
    # 这一组只用附件本身与 common.py 的标定结果，不涉及四问的解，先画完再进入问题 1。
    print("[fig01] 环境序列与拟合")
    # ---------------- fig01：烘房环境这条边界条件本身与其一阶惯性拟合 ----------------
    # 数据取自 common.py 载入的附件 1 实测序列：cm._t_env [s]、cm._Ta_env [degC]、
    # cm._Ca_env [kg/kg]（0~4 h，每 60 s 一点）。横轴均为 时间/h。
    # (a) 纵轴 烘房温度/degC（灰点实测、红实线一阶惯性拟合、黑点线设定值 50.212 degC、
    # 绿虚线标出 4 h 的阶段分界）；(b) 纵轴 含湿量/(kg/kg)，同一套画法，设定值
    # 0.050908 kg/kg。可见 4 h 之后实测与拟合基本重合，故其后直接取设定值。

    # 加花：实测散点压成中灰、一阶惯性拟合线发光——这张图承担两个论点（数据处理合格、
    # 阶段分界取 4 h），主角是拟合线而不是散点。两幅各加一处箭头标注，把正文那句
    # "t 约 4 h 达到设定值并稳定"直接画在图上，读者不必回正文找依据。
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    tt = np.linspace(0, 14400, 400)
    Ta_fit = cm.Tset - (cm.Tset - cm.T_air0) * np.exp(-tt / cm.tauT)
    Ca_fit = cm.Cset - (cm.Cset - cm.C_air0) * np.exp(-tt / cm.tauC)
    for k, (meas, fit, setv, lab_set) in enumerate([
            (cm._Ta_env, Ta_fit, cm.Tset, r"$T_{\rm set}=50.212$"),
            (cm._Ca_env, Ca_fit, cm.Cset, r"$C_{\rm set}=0.050908$")]):
        ax[k].plot(cm._t_env / 3600, meas, "o", ms=2.6, alpha=.45,
                   color=COL_MEAS, label="附件1 实测")
        _glow(ax[k], tt / 3600, fit, COL_TEMP, lw=1.8, label="一阶惯性拟合")
        ax[k].axhline(setv, color="k", ls=":", lw=1, label=lab_set)
        ax[k].axvline(4.0, color=COL_PHASE, ls="--", lw=1.3, label="阶段分界 4 h")
        ax[k].set_xlabel("时间 / h")
        ax[k].legend(fontsize=8, loc="center right")
    ax[0].set_ylabel("烘房温度 / $^\\circ$C")
    ax[0].set_title("(a) 烘房温度：实测与拟合")
    ax[1].set_ylabel("含湿量 / (kg/kg)")
    ax[1].set_title("(b) 烘房含湿量：实测与拟合")
    # 标注文字贴着 4 h 分界线放在轴底：不加引线，避免引线横穿拟合曲线（装饰不许压数据）
    for k, setv in enumerate([cm.Tset, cm.Cset]):
        lo, hi = ax[k].get_ylim()
        ax[k].text(3.88, lo + (hi - lo) * 0.03, "4 h 起\n达设定值",
                   fontsize=9.5, color=COL_PHASE, ha="right", va="bottom")
    save(fig, "fig01_env_fit.png")

    print("[fig02] 拟合残差诊断")
    # ---------------- fig02：拟合残差的时间序列诊断 ----------------
    # 残差定义为 实测 − 拟合式，取附件 1 全序列。横轴均为 时间/h。
    # (a) 纵轴 温度残差/degC，序列 RMSE = 0.336 degC；(b) 纵轴 含湿量残差/(g/kg)
    # （原值乘 1e3；含湿量残差只有 1e-4 kg/kg 量级，不放大在图上就是一条直线），
    # RMSE = 0.808 g/kg。两幅都画了 y=0 基准线，用来看残差有没有系统性偏移。
    rT = cm._Ta_env - (cm.Tset - (cm.Tset - cm.T_air0) * np.exp(-cm._t_env / cm.tauT))
    rC = cm._Ca_env - (cm.Cset - (cm.Cset - cm.C_air0) * np.exp(-cm._t_env / cm.tauC))
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(cm._t_env / 3600, rT, lw=1)
    ax[0].axhline(0, color="k", lw=.8)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("残差 / $^\\circ$C")
    ax[0].set_title("(a) 温度拟合残差时序（RMSE=0.336）")
    ax[1].plot(cm._t_env / 3600, rC * 1e3, lw=1, color="C1")
    ax[1].axhline(0, color="k", lw=.8)
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("残差 / (g/kg)")
    ax[1].set_title("(b) 含湿量拟合残差时序（RMSE=0.808 g/kg）")
    save(fig, "fig02_env_resid_time.png")

    print("[fig03] 残差统计（直方图+Q-Q）")
    # ---------------- fig03：残差的正态性统计诊断（2x2） ----------------
    # 上行 (a)(b) 是温度残差，下行 (c)(d) 是含湿量残差（放大 1e3 成 g/kg，同 fig02）。
    # (a)(c) 横轴 残差/对应单位、纵轴 概率密度：柱状为直方图，红实线为同均值同标准差的
    # 正态密度；(b)(d) 是 Q-Q 图，横轴 理论分位、纵轴 样本分位。两幅用来看标定的残差
    # 能否按白噪声处理，若明显偏离正态，后面用 RMSE 描述拟合优度就要加一句说明。
    from scipy import stats
    # hspace 必须手动放宽：默认 0.2 时下排两幅子图的标题会与上排子图的
    # x 轴刻度数字叠在一起（实测在 9.2 in 画布 + 11 pt 标题下已经压字）。
    fig, ax = plt.subplots(2, 2, figsize=(9.2, 6.6))
    fig.subplots_adjust(hspace=0.45)
    for k, (r, nm, un) in enumerate([(rT, "温度", "$^\\circ$C"), (rC * 1e3, "含湿量", "g/kg")]):
        ax[k, 0].hist(r, bins=24, density=True, alpha=.7, color="C%d" % k)
        xs = np.linspace(r.min(), r.max(), 200)
        ax[k, 0].plot(xs, stats.norm.pdf(xs, r.mean(), r.std()), "r-", lw=1.6)
        ax[k, 0].set_xlabel("残差 / %s" % un); ax[k, 0].set_ylabel("概率密度")
        ax[k, 0].set_title("(%s) %s残差直方图与正态密度" % ("ac"[k], nm))
        stats.probplot(r, dist="norm", plot=ax[k, 1])
        ax[k, 1].set_title("(%s) %s残差 Q-Q 图" % ("bd"[k], nm))
        ax[k, 1].set_xlabel("理论分位"); ax[k, 1].set_ylabel("样本分位")
    save(fig, "fig03_env_resid_stat.png")

    print("[fig04] 附件2 半径历史")
    # ---------------- fig04：附件 2 的药材半径历史（问题 4 的收缩输入） ----------------
    # 数据取自 common.py 载入的附件 2：cm._t_R [s] 与 cm._R_hist [m]，横轴均为 时间/h。
    # (a) 纵轴 半径 R/cm；(b) 纵轴 逐区间收缩速率 |ΔR/Δt|/(cm/h)，跨 3 个数量级用半对数轴；
    # (c) 纵轴 无量纲收缩比 R/R_0，红虚线标出终点 0.599。
    #
    # 为什么 (b) 用逐区间差分而不是 np.gradient：附件 2 是每 30 min 一点的实测序列，

    # 半径按 0.001~0.127 cm 的**步进**下降（144 个区间里 92 个 ΔR 恰为 0），中心差分会把
    # 一次步进摊到相邻两个区间上、凭空造出"半格"的假台阶；逐区间 ΔR/Δt 才是这一步实际
    # 发生的速率。半径不变的区间速率为**严格 0**，半对数轴无法表示，故置 NaN 断开不画
    # （早先用 +1e-12 兜底，结果在图上画出一条贴着 1e-12 的长直线，读起来像"速率降到
    # 1e-12 cm/h"，其实那一段是"收缩停了"）。
    tR = cm._t_R; Rh = cm._R_hist
    d_rate = np.abs(np.diff(Rh) * 100) / (np.diff(tR) / 3600.0)     # cm/h，逐区间速率
    d_xm = (tR[:-1] + tR[1:]) / 2 / 3600.0                          # 区间中点 [h]
    d_rate = np.where(d_rate > 0, d_rate, np.nan)
    # 画布与布局：三个面板都带中文/公式纵轴标签，默认间距会把 (b)、(c) 的纵轴标签压进
    # 左邻面板的绘图区（实测压入 0.157 in，观感就是"文字叠在一起"），故改用 constrained
    # 布局自动给标签留位。画布宽 8.2 in 是按"存出 PNG 宽 ÷ 版心宽"反推的：印到 6.30 in
    # 版心时轴标签 7.96 pt、刻度 7.20 pt，都过 7 pt 可读线（10.4 in 时只有 7.56 / 6.84）。
    fig = plt.figure(figsize=(8.2, 3.05), layout="constrained")
    ax = fig.subplots(1, 3)
    ax[0].plot(tR / 3600, Rh * 100, "o-", ms=2.2, lw=1.2)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("$R$ / cm")
    ax[0].set_title("(a) 半径历史")
    ax[1].semilogy(d_xm, d_rate, "-o", ms=2.4, lw=1.2, color="C1")
    ax[1].set_ylim(5e-4, 5e-1)                 # 数据范围 1e-3~2.5e-1，上下各留约半格余量
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel(r"$|\mathrm{d}R/\mathrm{d}t|$ / (cm/h)")
    ax[1].set_title("(b) 收缩速率（半对数）")
    ax[2].plot(tR / 3600, Rh / cm.R0, lw=1.5, color="C2")
    ax[2].axhline(0.599, color="r", ls="--", lw=1.2, label=r"$R/R_0=0.599$")
    ax[2].set_xlabel("时间 / h"); ax[2].set_ylabel("$R/R_0$")
    ax[2].set_title("(c) 无量纲收缩比"); ax[2].legend(fontsize=8)
    save(fig, "fig04_radius.png")

    print("[fig05] 无量纲分析")
    # ---------------- fig05：无量纲数与特征时间尺度（问题 1 的物性口径） ----------------
    # 本图不画解，只用问题 1 的常数物性算无量纲数：热扩散率 alpha=k/(rho*cp)、
    # 质扩散率 D_0=D(C_0)=7e-9*exp(-0.89/2.55)。横轴/纵轴见下。
    # (a) 横轴 时间/s、纵轴 Fourier 数，画 Fo_T=alpha*t/R_0^2 与 Fo_D=D_0*t/R_0^2 两条线，
    # 说明热快质慢；(b) 横轴 特征时间/h（对数）、纵轴三个尺度 R_0^2/alpha、R_0^2/D_0、
    # L^2/alpha；(c) 关掉坐标轴，直接列出 Bi、Bi_m、Le 与两个 Fo 在 t=1800 s 的值。
    alpha = 0.36 / (820 * 2600)
    D0 = 7e-9 * np.exp(-0.89 / 2.55)
    # 画布由 10.4 in 缩到 9.6 in（字号不动）：印后字号按"存出 PNG 宽 ÷ 6.30 in 版心"算，
    # 10.4 in 时刻度只有 6.84 pt、低于 7 pt 的可读下限；缩到 9.6 in 后刻度约 7.5 pt、
    # 轴标签约 8.3 pt。宽度仍是 \textwidth，图的印刷高度基本不变。
    fig, ax = plt.subplots(1, 3, figsize=(9.6, 2.85))
    tts = np.linspace(60, 1800, 200)
    ax[0].plot(tts, alpha * tts / cm.R0 ** 2, lw=1.6, label=r"${\rm Fo}_T$")
    ax[0].plot(tts, D0 * tts / cm.R0 ** 2, lw=1.6, label=r"${\rm Fo}_D$")
    ax[0].set_xlabel("时间 / s"); ax[0].set_ylabel("Fourier 数")
    ax[0].set_title("(a) ${\\mathrm{Fo}}_D\\ll{\\mathrm{Fo}}_T$：热快质慢"); ax[0].legend()
    names = ["$R_0^2/\\alpha$", "$R_0^2/D$", "$L^2/\\alpha$"]
    vals = [cm.R0 ** 2 / alpha / 3600, cm.R0 ** 2 / D0 / 3600, 0.25 ** 2 / alpha / 3600]
    ax[1].barh(names, vals, color=["C0", "C1", "C2"])
    ax[1].set_xscale("log"); ax[1].set_xlabel("特征时间 / h")
    ax[1].set_title("(b) 特征时间尺度")
    # 横轴上限放宽到 600 h：最长那条 102.8 h 的数值标签原来按 v*1.15 摆放，右端会跑出
    # 坐标框 0.34 in（印后 24 pt），正是"文字压到图上"的典型。放宽上限后条长仍占 76%，
    # 而标签改按固定 4 pt 摆，不再随数量级变化。
    ax[1].set_xlim(0.35, 600)
    ax[1].tick_params(axis="y", pad=5.0)        # 类别标签与色条左端之间留出约 5 pt
    for i, v in enumerate(vals):
        ax[1].annotate("%.1f h" % v, xy=(v, i), xytext=(4, 0), textcoords="offset points",
                       va="center", ha="left", fontsize=8)
    ax[2].axis("off")
    # (c) 原先把 5 个量用等宽字体拼成一个整块，字体与全篇不一致、(c) 与 (a)(b) 像是两篇
    # 文章的图。改成两列排版：左列符号右对齐、右列数值左对齐，"=" 自然成一条竖线。
    rows = [(r"$\mathrm{Bi}$", 25 * cm.R0 / 0.36),
            (r"$\mathrm{Bi}_m$", 8e-7 * cm.R0 / D0),
            (r"$\mathrm{Le}$", D0 / alpha),
            (r"$\mathrm{Fo}_T(1800)$", alpha * 1800 / cm.R0 ** 2),
            (r"$\mathrm{Fo}_D(1800)$", D0 * 1800 / cm.R0 ** 2)]
    for k, (sym, val) in enumerate(rows):
        yy = 0.78 - 0.16 * k
        ax[2].text(0.47, yy, sym, ha="right", va="center", fontsize=10.5)
        ax[2].text(0.53, yy, "= %.4f" % val, ha="left", va="center", fontsize=10.5)
    ax[2].set_title("(c) 无量纲数与判据")
    save(fig, "fig05_dimensionless.png")

    # ---------------- 第 2 组：坐标变换示意（fig06） ----------------
    # 下面一幅图不画数据，只画离散节点本身，用来说明为什么换到 u=xi^2 求解。
    print("[fig06] 坐标变换示意")
    # ---------------- fig06：常规 r 坐标与正则化坐标 u=xi^2 的节点分布 ----------------
    # 画的是离散节点位置，不涉及物理解。节点由 N=16 的 Chebyshev-Gauss-Lobatto 点
    # x_j 映射而来：r_j=(1-x_j)*R_0/2、u_j=(1-x_j)/2。
    # (a) 横轴 r/cm、纵轴隐去：常规 r 坐标下节点在轴心处高度密集，r=0 处算子奇异、需要
    # 特殊处理；(b) 横轴 u=xi^2：同一批节点在 u 坐标下等距，u=0 处算子已正则；
    # (c) 横轴 u、纵轴 r/cm：画出映射曲线 r=R*sqrt(u) 及节点落在曲线上的位置。

    # 画布与布局：默认间距下 (c) 的纵轴标签会压进 (b) 的绘图区（实测 0.048 in），故用
    # constrained 布局自动给标签留位；画布宽 8.2 in 是按"存出 PNG 宽 ÷ 版心宽"反推的，
    # 印到 6.30 in 版心时轴标签约 8.1 pt、刻度约 7.3 pt。
    fig = plt.figure(figsize=(8.2, 3.05), layout="constrained")
    ax = fig.subplots(1, 3)
    N = 16
    xj = np.cos(np.pi * np.arange(N + 1) / N)
    rj = (1 - xj) / 2 * cm.R0 * 100
    aj = (1 - xj) / 2
    ax[0].plot(rj, np.zeros_like(rj), "o", ms=6, color="C3")
    ax[0].set_ylim(-1, 1); ax[0].set_yticks([])
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_title("(a) 常规 $r$ 坐标节点\n轴心处密集、需特殊处理")
    ax[1].plot(aj, np.zeros_like(aj), "o", ms=6, color="C0")
    ax[1].set_ylim(-1, 1); ax[1].set_yticks([])
    ax[1].set_xlabel("$u=\\xi^2$"); ax[1].set_title("(b) 正则化坐标 $u$ 节点\n$u=0$ 处算子正则")
    uu = np.linspace(0, 1, 200)
    ax[2].plot(uu, np.sqrt(uu) * cm.R0 * 100, lw=1.8)
    ax[2].plot(aj, np.sqrt(aj) * cm.R0 * 100, "o", ms=6, color="C0")
    ax[2].set_xlabel("$u$"); ax[2].set_ylabel("$r$ / cm")
    ax[2].set_title("(c) 映射 $r=R\\sqrt{u}$")
    save(fig, "fig06_u_transform.png")

    # ---------------- 第 3 组：问题 1 的解与离散精度验证（fig07~fig11） ----------------
    # 问题 1 的物性取附录 2 的常数，温度场与水分场解耦；fig07~fig09 画数值解本身，
    # fig10~fig11 拿 Bessel 闭式解作基准校核离散精度。
    C1, T1 = _samples(S1, S1)
    print("[fig07-08] 问题1 剖面 / 云图")
    # ---------------- fig07：问题 1 的温度场与水分浓度场的径向剖面 ----------------
    # 数据取自 中间数据/solutions.npz 的 p1_t [s] 与 p1_y（前 65 行是水分浓度 C、后 65 行
    # 是温度 T），经 _samples() 按 u=(r/R_0)^2 插值到 RCOLS 五个物理半径，另存 C1、T1。
    # 横轴均为 到中心距离 r/cm。(a) 纵轴 温度/degC；(b) 纵轴 水分浓度 C/(kg/kg)。
    # 两幅都画 t = 100、600、1200、1800 s 四个时刻，对比即可看出热与质传播的快慢差异。
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    for k, tt in enumerate([100, 600, 1200, 1800]):
        i = int(np.argmin(np.abs(t1 - tt)))
        ax[0].plot(RCOLS * 100, T1[i], "o-", ms=3.5, lw=1.4, label="t=%d s" % tt)
        ax[1].plot(RCOLS * 100, C1[i], "s-", ms=3.5, lw=1.4, label="t=%d s" % tt)
    ax[0].set_xlabel("到中心距离 $r$ / cm"); ax[0].set_ylabel("温度 / $^\\circ$C")
    ax[0].set_title("(a) 问题1 温度剖面演化"); ax[0].legend(fontsize=8)
    ax[1].set_xlabel("到中心距离 $r$ / cm"); ax[1].set_ylabel("水分浓度 / (kg/kg)")
    ax[1].set_title("(b) 问题1 水分浓度剖面演化"); ax[1].legend(fontsize=8)
    save(fig, "fig07_p1_profiles.png")

    # ---------------- fig08：问题 1 的时空云图（温度场与水分浓度场） ----------------
    # 数据同样来自 p1_t、p1_y，这里改用 220 个 u 节点密采样，再按 r=sqrt(u)*R_0 换算
    # 成物理半径 rr [cm]。横轴均为 r/cm、纵轴均为 时间/min。
    # (a) 温度时空云图，色标 温度/degC；(b) 水分浓度时空云图，色标 C/(kg/kg)。
    # 与 fig05(c) 的 Fo_T、Fo_D 对照可见温度场已经先趋于均匀，水分浓度场却基本未动。
    # 加花：4 张云图原先只有一块裸色块，看不出"结构在哪"。这里叠等值线 + 线上标数值

    # （深色描边白字，压在任何底色上都读得清），并在 (b) 右上角加一处内嵌放大，把
    # r>1.8 cm 那条约 2 mm 厚的表面薄层摊开——问题 1 的物理内容全挤在这条薄带里。
    # 内嵌面板的字号与主面板同级（9.5 pt，印后 7.2 pt）：不能按"内嵌就该小一号"缩，
    # 缩了印到纸上就低于 7 pt 的可读下限。等值线是插值结果，本图 u 方向 220 点、
    # 时间方向 181 点，足以支撑这几个水平，不会造出假细节。
    uu = np.linspace(0, 1, 220)
    rr = np.sqrt(uu) * cm.R0 * 100
    Tg = np.array([s1.interp_u(y1[:, i], uu)[1] for i in range(len(t1))])
    Cg = np.array([s1.interp_u(y1[:, i], uu)[0] for i in range(len(t1))])
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    # colorbar 是从 ax[0] 里"偷"走宽度的，不占用两幅子图之间的间隙；
    # 默认 wspace 下 colorbar 的标签会顶到右图的 y 轴标题上，故手动放宽。
    fig.subplots_adjust(wspace=0.42)
    tmin = t1 / 60
    m0 = ax[0].pcolormesh(rr, tmin, Tg, shading="auto", cmap="inferno")
    plt.colorbar(m0, ax=ax[0], label="温度 / $^\\circ$C")
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_ylabel("时间 / min")
    ax[0].set_title("(a) 问题1 温度时空云图")
    m1 = ax[1].pcolormesh(rr, tmin, Cg, shading="auto", cmap="viridis")
    plt.colorbar(m1, ax=ax[1], label="水分浓度 / (kg/kg)")
    ax[1].set_xlabel("$r$ / cm"); ax[1].set_ylabel("时间 / min")
    ax[1].set_title("(b) 问题1 水分浓度时空云图")
    # 云图上关网格：全局 axes.grid=True 会把虚线网格叠在色块上，盖住数据纹理
    ax[0].grid(False); ax[1].grid(False)
    # 等值线取物理上有含义的整数值，数值直接标在线上（深色字 + 白描边，任何底色上都读得清）
    lv_T = [30, 32, 34, 36]
    lv_C = [2.0, 2.2, 2.4, 2.5]
    for a, Z, lv in ((ax[0], Tg, lv_T), (ax[1], Cg, lv_C)):
        cs = a.contour(rr, tmin, Z, levels=lv, colors="white", linewidths=0.9)
        for t in a.clabel(cs, fmt="%g", fontsize=9, inline=True, inline_spacing=2):
            t.set_color("#1a1a1a")
            t.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white")])
    # 内嵌放大：放在 (b) 左上角那一整块"纯色区"上——那里本来没有信息，压掉不损失任何东西，
    # 却能把 r>1.75 cm 那条约 2 mm 厚的表面薄层摊开。不放 indicate_inset_zoom 的连接线：
    # 放大区在右边缘，引线会横穿整个色块。
    axin = ax[1].inset_axes([0.075, 0.47, 0.47, 0.41])
    axin.pcolormesh(rr, tmin, Cg, shading="auto", cmap="viridis")
    cs2 = axin.contour(rr, tmin, Cg, levels=lv_C, colors="white", linewidths=0.8)
    for t in axin.clabel(cs2, fmt="%g", fontsize=9, inline=True, inline_spacing=2):
        t.set_color("#1a1a1a")
        t.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white")])
    axin.set_xlim(1.75, 2.0)
    axin.tick_params(labelsize=9.5, labelleft=False)
    axin.set_title("表面薄层 $r>1.75$ cm", fontsize=9.5, pad=3)
    axin.grid(False)
    for s in axin.spines.values():
        s.set_color("white"); s.set_linewidth(1.4)
    save(fig, "fig08_p1_spacetime.png")

    print("[fig09] 问题1 表面薄层与 Robin 条件核对")
    # ---------------- fig09：问题 1 表面薄层的 Robin 边界条件核对 ----------------
    # 取最后一条剖面（t=1800 s）。(a) 横轴 r/cm——只放大到 1.6~2.0 cm 的表面薄层，纵轴
    # 水分浓度 C/(kg/kg)，黑点线标出表面 r=2.0 cm；(b) 纵轴 表面浓度梯度 |dC/dr|/(kg/kg/m)：
    # 左柱是剖面数值差分值，右柱是 Robin 条件 -h_m*(C_s-C_air)/D 的预测值，两者的相对差
    # 写在子图标题里。这是边界条件实现得对不对的一道直接校验。
    i18 = len(t1) - 1
    Cn, _ = s1.interp_u(y1[:, i18], uu)
    Ds = 7e-9 * np.exp(-0.89 / max(Cn[-1], 1e-3))
    grad_num = (Cn[-1] - Cn[-2]) / ((rr[-1] - rr[-2]) / 100)
    grad_robin = -8e-7 * (Cn[-1] - float(cm.C_air(1800.0))) / Ds
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(rr, Cn, lw=1.8)
    ax[0].axvline(2.0, color="k", ls=":", lw=1)
    ax[0].set_xlim(1.6, 2.0)
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_title("(a) 1800 s 近表面水分剖面（薄层）")
    # 两根柱值几乎相等（相对差 <1%），柱宽取 0.32：默认 0.8 会把绘图区填成两块实心
    # 色块，过细则两根之间空出一大片。同时把组中心间距从默认的 1.0 收到 0.5 并收紧
    # xlim，让数据占满画幅；柱顶标出数值，便于直接核对。
    _v9 = [abs(grad_num), abs(grad_robin)]
    _x9 = np.array([0.0, 0.5])
    ax[1].bar(_x9, _v9, width=0.32, color=["C0", "C1"])
    ax[1].set_xticks(_x9)
    ax[1].set_xticklabels(["数值差分\n(Robin 左端)", "Robin 条件\n右端"])
    ax[1].set_xlim(-0.30, 0.80)
    for _x, _v in zip(_x9, _v9):
        ax[1].text(_x, _v, "%.1f" % _v, ha="center", va="bottom", fontsize=10.5)
    ax[1].set_ylim(0, max(_v9) * 1.14)
    ax[1].set_ylabel("$|\\partial C/\\partial r|$ / (kg/kg/m)")
    ax[1].set_title("(b) 表面梯度：差分值 vs Robin 预测\n(相对差 %.1f%%)"
                    % (abs(abs(grad_num) - abs(grad_robin)) / abs(grad_robin) * 100))
    save(fig, "fig09_p1_surface.png")

    print("[fig10-11] Bessel 展开与验证")
    # ---------------- fig10：问题 1 闭式解（Bessel 级数）的收敛速度 ----------------
    # 数据取自 analytic_p1.py 在导入时算出的展开系数 an.BN 与各阶衰减率 an.GAM。
    # (a) 横轴 阶数 n、纵轴 |b_n|（半对数）：系数随阶数快速衰减；
    # (b) 横轴 阶数 n、纵轴 gamma_n*t（t=1800 s），黑虚线为 y=1：n>=2 的阶次 gamma_n*t
    # 已远大于 1，对应 e^{-gamma_n t} 可忽略，故闭式解只取前几阶即可当作基准解。
    import analytic_p1 as an
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    n = np.arange(1, 9)
    ax[0].semilogy(n, np.abs(an.BN[:8]), "o-", lw=1.5)
    ax[0].set_xlabel("阶数 $n$"); ax[0].set_ylabel("$|b_n|$")
    ax[0].set_title("(a) 展开系数 $b_n$ 的衰减")
    ax[1].plot(n, an.GAM[:8] * 1800, "s-", lw=1.5, color="C3")
    ax[1].axhline(1, color="k", ls="--", lw=1)
    ax[1].set_xlabel("阶数 $n$"); ax[1].set_ylabel("$\\gamma_n t$ (t=1800 s)")
    ax[1].set_title("(b) 各阶衰减因子：$n\\geq 2$ 项已被压至 $e^{-13}$ 以下")
    save(fig, "fig10_bessel_coef.png")

    # ---------------- fig11：谱配置法与闭式解的一致性验证 ----------------
    # 基准是 analytic_p1.T_exact 给出的闭式温度解，被校核的是 HerbSolver 在 t=1800 s 的
    # 解：用 _bary 把 u 节点上的解重心插值到 201 个 (r/R_0)^2 处。
    # (a) 横轴 谱节点数 N、纵轴 与闭式解的最大偏差/K（半对数）：N=8 即降到 1e-10，红虚线
    # 是 1e-10 的时间积分误差平台；(b) 横轴 r/cm、纵轴 偏差/(1e-10 K)：偏差沿半径的分布，
    # 量级已是纯舍入水平，说明空间离散本身没有引入系统性误差。
    rr2 = np.linspace(0, cm.R0, 201)
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    Ns = [4, 6, 8, 12, 16, 24, 32, 48, 64]
    err = []
    for NN in Ns:
        ss = HerbSolver(NN, prob=1)
        a_ = np.concatenate([np.full(ss.Np, 1e-14), np.full(ss.Np, 1e-12)])
        so = ss.solve(1800.0, t_eval=[1800.0], rtol=1e-12, atol=a_)
        Ti = ss._bary(so.y[ss.Np:, -1], (rr2 / cm.R0) ** 2)
        err.append(np.max(np.abs(Ti - an.T_exact(rr2, 1800.0))))
    ax[0].semilogy(Ns, err, "o-", lw=1.6)
    ax[0].axhline(1e-10, color="r", ls="--", lw=1, label="时间积分误差平台 $10^{-10}$")
    ax[0].set_xlabel("谱节点数 $N$"); ax[0].set_ylabel("与闭式解的最大偏差 / K")
    ax[0].set_title("(a) 谱收敛性：$N=8$ 即达 $10^{-10}$"); ax[0].legend(fontsize=8)
    ss = HerbSolver(64, prob=1)
    a_ = np.concatenate([np.full(ss.Np, 1e-14), np.full(ss.Np, 1e-12)])
    so = ss.solve(1800.0, t_eval=[1800.0], rtol=1e-12, atol=a_)
    Ti = ss._bary(so.y[ss.Np:, -1], (rr2 / cm.R0) ** 2)
    res = Ti - an.T_exact(rr2, 1800.0)
    ax[1].plot(rr2 * 100, res * 1e10, lw=1.4)
    ax[1].set_xlabel("$r$ / cm"); ax[1].set_ylabel("偏差 / ($10^{-10}$ K)")
    ax[1].set_title("(b) 偏差沿半径的分布（纯舍入量级）")
    save(fig, "fig11_bessel_validate.png")

    # ---------------- 第 4 组：问题 2 的变物性耦合解（fig12~fig14） ----------------
    # 问题 2 的物性改取附录 3（含 C、T 的强非线性关系），温度场与水分场双向耦合。
    C2, T2 = _samples(S2, S2)
    print("[fig12-13] 问题2")
    # ---------------- fig12：问题 2 的温度与水分浓度剖面演化 ----------------
    # 数据取自 solutions.npz 的 p2_t [s]、p2_y，经 _samples() 插值到 RCOLS 得 C2、T2。
    # 横轴均为 r/cm。(a) 纵轴 温度/degC；(b) 纵轴 C/(kg/kg)。两幅都画 0.5、1.0、1.5、
    # 2.0、3.0 h 五个时刻；(b) 中表面浓度早早掉到很低的水平而中心仍然很高，即干壳形成。
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    hs = []
    for tt in [0.5, 1.0, 1.5, 2.0, 3.0]:
        i = int(np.argmin(np.abs(t2 - tt * 3600)))
        hh, = ax[0].plot(RCOLS * 100, T2[i], "o-", ms=3.5, lw=1.4, label="%.1f h" % tt)
        ax[1].plot(RCOLS * 100, C2[i], "s-", ms=3.5, lw=1.4, label="%.1f h" % tt)
        hs.append(hh)
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_ylabel("温度 / $^\\circ$C")
    ax[0].set_title("(a) 问题2 温度剖面演化")
    ax[1].set_xlabel("$r$ / cm"); ax[1].set_ylabel("C / (kg/kg)")
    ax[1].set_title("(b) 问题2 水分剖面演化（干壳形成）")
    # 两幅子图的时刻序列完全相同，故共用一条图例放到底部：
    # 5 条曲线在面板内铺满整个纵向量程，面板内已无空位可放图例
    # （放左下角会盖住 0.5 h 曲线，而该曲线横跨整个版心）。
    fig.legend(hs, ["%.1f h" % tt for tt in [0.5, 1.0, 1.5, 2.0, 3.0]],
               loc="lower center", ncol=5, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.015))
    fig.subplots_adjust(bottom=0.24)
    save(fig, "fig12_p2_profiles.png")

    # ---------------- fig13：问题 2 的时空云图（0~3 h） ----------------
    # 数据同 fig12（p2_t、p2_y），密采样 220 个 u 节点。横轴均为 r/cm、纵轴均为 时间/h。
    # (a) 温度云图，色标 温度/degC；(b) 水分浓度云图，色标 C/(kg/kg)；(c) 干壳薄层放大。
    # 与 (a) 对照可见温度场在径向已基本均匀，水分却只在表面附近形成陡峭的低浓度薄层，
    # 即所谓干壳。加花与 fig08 同一套：叠等值线 + 线上标数值；把"干壳边界在哪"直接摊开
    # ——这是问题 2 的核心结论，原先只能从颜色深浅猜。
    Cg2 = np.array([s2.interp_u(y2[:, i], uu)[0] for i in range(len(t2))])
    Tg2 = np.array([s2.interp_u(y2[:, i], uu)[1] for i in range(len(t2))])
    # 三幅并排：(c) 把原先压在 (b) 左上角的内嵌放大图拉出来独立成面板——小图不再遮挡
    # 主图数据，放大后的干壳薄层也看得清。画幅加宽到 13 in，按 \textwidth 排入时高度反而
    # 比原来的 9.6×3.4 更矮，不会挤占正文版面；字号相应放大以补偿缩印比例。
    fig, ax = plt.subplots(1, 3, figsize=(13.0, 3.5))
    # colorbar 刻度是四位（47.5/42.5…），比 fig08 更宽，故 wspace 留得比默认大。
    fig.subplots_adjust(wspace=0.42)
    t2h = t2 / 3600
    m0 = ax[0].pcolormesh(rr, t2h, Tg2, shading="auto", cmap="inferno")
    plt.colorbar(m0, ax=ax[0], label="温度 / $^\\circ$C")
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_ylabel("时间 / h")
    ax[0].set_title("(a) 问题2 温度时空云图")
    m1 = ax[1].pcolormesh(rr, t2h, Cg2, shading="auto", cmap="viridis")
    plt.colorbar(m1, ax=ax[1], label="C / (kg/kg)")
    ax[1].set_xlabel("$r$ / cm"); ax[1].set_ylabel("时间 / h")
    ax[1].set_title("(b) 问题2 水分时空云图")
    m2 = ax[2].pcolormesh(rr, t2h, Cg2, shading="auto", cmap="viridis")
    plt.colorbar(m2, ax=ax[2], label="C / (kg/kg)")
    ax[2].set_xlim(1.7, 2.0)
    ax[2].set_xlabel("$r$ / cm"); ax[2].set_ylabel("时间 / h")
    ax[2].set_title("(c) 干壳薄层放大")     # 放大区间写进图注，标题不再压到 y 轴刻度上
    ax[0].grid(False); ax[1].grid(False); ax[2].grid(False)
    lv_T2 = [32, 38, 44, 48]
    lv_C2 = [1.2, 1.6, 2.0, 2.3]
    # 等值线：(a) 温度场标数值；(b) 只留线、不标值——(b) 与 (c) 数据相同，同一个量标两遍
    # 既重复，又会在 (b) 右边界上挤出框外（实测最外侧那个标值压框 0.063 in）；放大后的
    # (c) 里标值看得最清楚，标值统一放那里。
    for a, Z, lv, do_label in ((ax[0], Tg2, lv_T2, True),
                               (ax[1], Cg2, lv_C2, False),
                               (ax[2], Cg2, lv_C2, True)):
        cs = a.contour(rr, t2h, Z, levels=lv, colors="white", linewidths=0.9)
        if not do_label:
            continue
        for t in a.clabel(cs, fmt="%g", fontsize=11.5, inline=True, inline_spacing=2):
            t.set_color("#1a1a1a")
            t.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white")])
    for a in ax:
        a.tick_params(labelsize=11.5)
        a.xaxis.label.set_size(12); a.yaxis.label.set_size(12)
        a.title.set_size(12.5)
    save(fig, "fig13_p2_spacetime.png")

    print("[fig14] 表面传质通量（降速特征）")
    # ---------------- fig14：表面传质通量的降速特征与表面/中心水分浓度 ----------------
    # 数据不由缓存提供：这里另解一次问题 3（prob=3，控制方程、边界与附录 3 物性都与问题 2
    # 相同）并积分到 72 h，再从解的 y 中取 u=1 的表面 C_s 与 u=0 的中心 C_0 两条序列。
    # 横轴均为 时间/h。(a) 纵轴 表面传质通量 |h_m*(C_s-C_air)|/(kg/(m2 s))，半对数，用来看
    # 干壳形成后通量如何持续衰减；(b) 纵轴 C/(kg/kg)，两条曲线分别标为"中心"与"表面"，
    # 中心浓度的长平台正是 (a) 里通量衰减那么快的原因。

    # 问题 2 只要求输出 0~3 h，但正文用 6~72 h 的通量衰减来说明降速特征。
    # 问题 2 与问题 3 的控制方程、边界和附录 3 物性完全相同，故这里直接用
    # 问题 3 的解（它本来就积分到 t_dry 之后）把通量曲线画满 72 h，
    # 避免"正文讲 6~72 h、图只画 0~3 h"的口径不一致。
    FLUX_T_END = 72 * 3600.0
    s_flux = HerbSolver(64, prob=3)
    at_f = np.concatenate([np.full(s_flux.Np, 1e-13), np.full(s_flux.Np, 1e-11)])
    t_flux = np.linspace(0.0, FLUX_T_END, 433)
    sol_flux = s_flux.solve(FLUX_T_END, t_eval=list(t_flux), rtol=1e-11, atol=at_f)
    Cs_flux = sol_flux.y[s_flux.Np - 1, :]
    C0_flux = sol_flux.y[0, :]
    flux = 8e-7 * (Cs_flux - cm.C_air(t_flux))          # kg/(m^2 s)
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].semilogy(t_flux / 3600, np.abs(flux), lw=1.8, color="C3")
    ax[0].set_xlabel("时间 / h")
    ax[0].set_ylabel("表面传质通量 $|h_m(C_s-C_{\\mathrm{air}})|$ / (kg/(m$^2\\cdot$s))")
    ax[0].set_title("(a) 表面传质通量：0~72 h 降速衰减")
    ax[1].plot(t_flux / 3600, C0_flux, lw=1.8, label="中心")
    ax[1].plot(t_flux / 3600, Cs_flux, lw=1.8, label="表面")
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("C / (kg/kg)")
    ax[1].set_title("(b) 中心与表面水分浓度（同期）"); ax[1].legend(fontsize=8)
    save(fig, "fig14_p2_flux.png")

    # ---------------- 第 5 组：问题 3 的烘干时长与降速特征（fig15~fig17） ----------------
    C3, T3 = _samples(S3, S3)
    # 截面平均按论文 S11.6 的定义 Cbar=int_0^1 C(u,t)du 计算（面积平均）。
    # 不能用 5 个物理半径采样点的算术平均代替，达标时刻也不能写死。
    Cbar3 = y3[:s3.Np].T @ s3.area_weights()
    t_bar3 = _cross_time(t3, Cbar3, CRIT)
    lag_h = (td3 - t_bar3) / 3600.0
    under_pct = (td3 - t_bar3) / td3 * 100.0
    print("  面积平均达标 %.4f h（中心 %.4f h）；滞后 %.2f h，低估 %.2f%%"
          % (t_bar3 / 3600, td3 / 3600, lag_h, under_pct))
    print("[fig15-16] 问题3")
    # ---------------- fig15：问题 3 的半对数干燥曲线与"平均判据"的迟滞 ----------------
    # 数据取自 solutions.npz 的 p3_t [s]、p3_y，经 _samples() 插值到 RCOLS 得 C3 [kg/kg]。
    # (a) 横轴 时间/h（限定 0~62 h）、纵轴 C/(kg/kg) 半对数：画 r=0、1.0、2.0 cm 三条
    # 干燥曲线，红虚线是判据 0.15，黑点线是中心达标时刻 t_dry；
    # (b) 横轴 时间/h、纵轴 C/(kg/kg)，把中心 C_0（实线）与截面面积平均 Cbar（虚线）画在
    # 一起：Cbar 比 C_0 早达标，即按平均判据算会把烘干时长低估 under_pct%（写在子图标题里）。
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    # 加花：把"预热平衡 / 恒温干燥"两阶段直接铺成底纹。本题的工艺结构就写在这条
    # 4 h 分界上，而原先 0~57 h 的横轴上完全看不出来，读者只能回正文找。底纹 alpha
    # 压到 0.075，低于判据线与曲线的对比度，不会盖住任何数据。
    ax[0].axvspan(0, 4, color=COL_PHASE, alpha=0.075)
    ax[0].axvspan(4, 62, color=COL_MOIST, alpha=0.055)
    # 阶段名原写在底纹带里（y=2.05）：4 h 宽的带放不下 4 个中文字，文字必然横向伸出带外，
    # 纵向又正好落在 r=0 与 r=1.0 cm 两条曲线上（实测压住 55 个数据点，左端还越出坐标框
    # 0.15 in）。改写到曲线上方的空白带（y 约 2.9，曲线最高只到 2.55），各接一段短引线指回
    # 自己的阶段带；引线只经过 y=2.68~2.80 这段无数据的区间。
    ax[0].annotate("预热平衡", xy=(2.0, 2.68), xytext=(1.0, 2.92), ha="left", va="center",
                   fontsize=8.5, color=COL_PHASE,
                   arrowprops=dict(arrowstyle="-", color=COL_PHASE, lw=0.8,
                                   alpha=0.75, shrinkA=1, shrinkB=1))
    ax[0].annotate("恒温干燥", xy=(20.0, 2.68), xytext=(19.0, 2.92), ha="left", va="center",
                   fontsize=8.5, color=COL_MOIST,
                   arrowprops=dict(arrowstyle="-", color=COL_MOIST, lw=0.8,
                                   alpha=0.75, shrinkA=1, shrinkB=1))
    for j, lbl, is_main in [(0, "r=0 中心", True), (2, "r=1.0 cm", False), (4, "r=2.0 cm 表面", False)]:
        # 中心干燥最慢、决定烘干时长，故 r=0 这条是主结论线，发光并加粗；其余两条灰细
        if is_main:
            _glow(ax[0], t3 / 3600, C3[:, j], COL_MOIST, lw=1.9, glow=3.6, label=lbl)
        else:
            ax[0].semilogy(t3 / 3600, C3[:, j], lw=1.4, color=COL_SIDE, label=lbl)
    ax[0].axhline(CRIT, color=COL_CRIT, ls="--", lw=1.3, label="判据 $C\\leq0.15$")
    ax[0].axvline(td3 / 3600, color="k", ls=":", lw=1.2)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_ylim(0.04, 3.5); ax[0].set_xlim(0, 62)
    ax[0].set_yscale("log")     # _glow 用的是 ax.plot，纵轴对数必须显式设定
    # t_dry 的说明原贴在左下角（x=36、y=0.055）：横跨 x=36~56 h 的那一格里，表面干燥
    # 曲线恰好也在 0.054~0.058（实测压住 556 个数据点）。改到两条干燥曲线之间的空带
    # （y 约 0.10：中心线与 r=1.0 cm 线在该区间都在 0.15 以上、表面线在 0.06 以下），
    # 再看箭头指回终点。
    ax[0].annotate("$t_{\\mathrm{dry}}=%.4f$ h" % (td3 / 3600),
                   xy=(td3 / 3600, CRIT), xytext=(38.0, 0.10), ha="left", va="center",
                   fontsize=9.5, color="k",
                   arrowprops=dict(arrowstyle="->", color="k", lw=1.0,
                                   shrinkA=2, shrinkB=2))
    ax[0].set_title("(a) 半对数干燥曲线")
    # 判据的说明原先单独贴在虚线左下方，横向正好横穿表面干燥曲线（压住 10 个数据点）；
    # 收进图例——图例本就是回答"这条线是什么"的地方，且框内没有数据。四条目使图例底边
    # 落到 y 约 0.9，与上面那行 t_dry 说明（顶边 0.5）仍留有余量。
    ax[0].legend(fontsize=8, loc="upper right")
    _glow(ax[1], t3 / 3600, C3[:, 0], COL_MOIST, lw=1.8, label="中心 $C_0$")
    ax[1].semilogy(t3 / 3600, Cbar3, lw=1.6, ls="--", color=COL_SIDE,
                   label="截面平均 $\\bar C=\\int_0^1 C\\,\\mathrm{d}u$")
    ax[1].axhline(CRIT, color=COL_CRIT, ls="--", lw=1.3)
    ax[1].axvline(t_bar3 / 3600, color=COL_PHASE, ls=":", lw=1.2)
    ax[1].axvline(td3 / 3600, color="k", ls=":", lw=1.2)
    ax[1].annotate("$\\bar C$ 达标 %.1f h" % (t_bar3 / 3600),
                   xy=(t_bar3 / 3600 + 1.0, 0.5), fontsize=9.5, color=COL_PHASE)
    ax[1].annotate("$C_0$ 达标 %.1f h" % (td3 / 3600),
                   xy=(td3 / 3600 - 21, 0.9), fontsize=9.5)
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("C / (kg/kg)")
    ax[1].set_yscale("log")     # 同上：_glow 走的是线性 plot，对数轴要显式开
    ax[1].set_title("(b) 判据迟滞：平均判据低估 %.1f%%" % under_pct); ax[1].legend(fontsize=8)
    save(fig, "fig15_p3_drying.png")

    # ---------------- fig16：烘干终点附近的放大与中心干燥速率 ----------------
    # 数据同 fig15（C3 取自 p3_y），横轴均为 时间/h。由相邻两点差分得中心干燥速率
    # rate = -dC_0/dt [kg/(kg/h)]。
    # (a) 只保留 52~60 h、纵轴 C/(kg/kg)：画中心与表面两条曲线，红虚线是判据 0.15，红星
    # 标出 C_0(t_dry)=0.1500 的终点；(b) 纵轴 |rate|/(kg/kg/h) 半对数，用来看中心干燥
    # 速率在降速期是否按指数衰减。

    # 加花：这是全文最关键的结论图——烘干时长本身原先只能回正文找。这里把 t_dry 直接
    # 标注在图上、判据线标出阈值、中心曲线发光（它是判据取的那条线），表面曲线降为灰细。
    # 注意 t_dry 必须由 td3 变量格式化：论文摘要用的是 N=128 连续极限 57.1043 h，本图画
    # 的是交付口径 N=64 的生产解，写死常量就会与图里的曲线对不上。
    # 画布与布局：默认间距下 (b) 的纵轴标签会压进 (a) 的绘图区（实测 0.015 in），故用
    # constrained 布局；画布宽 8.2 in 按"存出 PNG 宽 ÷ 版心宽"反推，印后轴标签约 8 pt。
    fig = plt.figure(figsize=(8.2, 3.05), layout="constrained")
    ax = fig.subplots(1, 2)
    m = (t3 / 3600 > 52) & (t3 / 3600 < 60)
    _glow(ax[0], t3[m] / 3600, C3[m, 0], COL_MOIST, lw=2.2, glow=4.0, label="中心 $C_0$")
    ax[0].plot(t3[m] / 3600, C3[m, -1], lw=1.4, color=COL_SIDE, label="表面 $C_s$")
    ax[0].axhline(CRIT, color=COL_CRIT, ls="--", lw=1.4)
    ax[0].axvline(td3 / 3600, color="k", ls=":", lw=1.4)
    ax[0].plot([td3 / 3600], [CRIT], "*", color=COL_TEMP, ms=14, zorder=6)
    # 纵轴上下各留一条空带：下带给图例、上带给"判据"标签，中间那条才是标注框。
    # 不留的话，标注框的引线会从图例上穿过去。
    ax[0].set_ylim(0.02, 0.185)
    # 标注框落在这两行之间的空白带（中心线在 0.15 以上、表面线在 0.054 以下，中间无数据）
    ax[0].annotate("$t_{\\mathrm{dry}}=%.4f$ h\n$C_0(t_{\\mathrm{dry}})=0.1500$" % (td3 / 3600),
                   xy=(td3 / 3600, CRIT), xytext=(52.25, 0.088), fontsize=8.5,
                   color=COL_TEMP, ha="left", va="center",
                   bbox=dict(boxstyle="round,pad=0.28", fc="white", ec=COL_TEMP, lw=0.8),
                   arrowprops=dict(arrowstyle="->", color=COL_TEMP, lw=1.1))
    ax[0].annotate("判据 $C\\leq0.15$", xy=(56.0, CRIT), xytext=(52.05, 0.181),
                   fontsize=9.5, color=COL_CRIT, va="top")
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_title("(a) 干燥终点放大：$C_0(t_{\\mathrm{dry}})=0.1500$")
    ax[0].legend(fontsize=8, loc="lower left")
    dt_ = np.diff(t3 / 3600); dC = np.diff(C3[:, 0])
    rate = -dC / dt_
    _glow(ax[1], (t3[:-1]) / 3600, np.abs(rate), COL_MOIST, lw=1.6)
    ax[1].set_yscale("log")     # 同上：_glow 走线性 plot，半对数图的纵轴要显式开
    # 标出拟合用的降速段 30~54 h：这段的斜率才是图 15 里那个特征时间
    ax[1].axvspan(30, 54, color=COL_PHASE, alpha=0.075)
    ax[1].text(42, 0.06, "降速段 30~54 h", fontsize=9.5, color=COL_PHASE, ha="center")
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("$|-\\dot C_0|$ / (kg/kg/h)")
    ax[1].set_title("(b) 中心干燥速率指数衰减（降速期）")
    save(fig, "fig16_p3_endpoint.png")

    print("[fig17] 干燥速率对数的线性拟合")
    # ---------------- fig17：干燥速率对数与浓度对数两条斜率的对比 ----------------
    # 把 fig16(b) 的速率取自然对数，再对时间做一次最小二乘（窗口取 30~54 h）。
    # (a) 横轴 时间/h、纵轴 ln(-dC_0/dt)：散点是数据、红实线是拟合直线，斜率单位为
    # h^-1，其负倒数即图上标的速率衰减特征时间；
    # (b) 横轴 时间/h、纵轴 ln C_0：浓度自身对数的斜率，由 30 h 之后的数据拟合得到。
    # 两幅斜率相差很多，说明"速率按指数衰减"与"浓度按指数衰减"不是同一件事。
    mm = (t3[:-1] / 3600 > 30) & (t3[:-1] / 3600 < 54) & (np.abs(rate) > 0)
    xs = t3[:-1][mm] / 3600
    ys = np.log(np.abs(rate[mm]))
    co = np.polyfit(xs, ys, 1)
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(xs, ys, ".", ms=3, alpha=.5, label="数据")
    ax[0].plot(xs, np.polyval(co, xs), "r-", lw=1.8,
               label="拟合斜率 %.4f h$^{-1}$" % co[0])
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("$\\ln(-\\dot C_0)$")
    ax[0].set_title("(a) 30~54 h：速率衰减的特征时间 %.1f h" % (-1 / co[0]))
    ax[0].legend(fontsize=8)
    ax[1].plot(t3 / 3600, np.log(C3[:, 0]), lw=1.8, color="C0")
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("$\\ln C_0$")
    ax[1].set_title("(b) 浓度自身的对数斜率仅 %.4f h$^{-1}$" % np.polyfit(
        t3[t3 / 3600 > 30] / 3600, np.log(C3[t3 / 3600 > 30, 0]), 1)[0])
    save(fig, "fig17_p3_rate.png")

    # ---------------- 第 6 组：问题 4 的收缩效应（fig18~fig20） ----------------
    # 问题 4 物性取附录 4，并叠加附件 2 的实测收缩 R(t)。u=xi^2 是物质坐标：收缩过程
    # 中求解域在该坐标下固定为 0~1，换算回物理距离要乘当前半径 R(t)。
    C4, T4 = _samples(S4, S4)
    Rt4 = np.array([s4.R_of(x) for x in t4])
    print("[fig18-19] 问题4")
    # ---------------- fig18：问题 4 在物理坐标下的剖面与收缩中的表面半径 ----------------
    # 数据取自 p4_t [s]、p4_y；先用 Rt4[i] 取出该时刻的当前半径，再把 u 节点换算成物理
    # 距离 r=sqrt(u)*R(t)，浓度由 interp_u 取回。
    # (a) 横轴 当前物理距离 r/cm、纵轴 C/(kg/kg)：画 6、12、18、24、48 h 五条剖面，图例
    # 括号里给的是该时刻的 R（cm），可见剖面整体向左移动；
    # (b) 横轴 时间/h、纵轴 R(t)/cm：收缩中的表面半径本身，输入来自附件 2。

    # 加花：这张图就是"模型到底有没有真的做收缩"的直接证据，而原先每条剖面只是默默
    # 左移了 13%，视觉上几乎注意不到。改进三处：(1)每条剖面的右端点（= 当时的药材表面）
    # 加实心圆点，并用由浅到深的水分蓝区分早晚时刻；(2)在 48 h 的端点直标当时的半径；
    # (3)用一条弧形箭头把 6 h 与 48 h 的端点连起来，把"表面左移了多少"写成数字。
    # 左移量由 Rt4 现算，不写死——图上任何数字都必须对应代码里的一次真实计算。
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    shades = ["#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08306b"]
    ends = {}
    for k, tt in enumerate([6, 12, 18, 24, 48]):
        i = int(np.argmin(np.abs(t4 - tt * 3600)))
        Rn = np.sqrt(s4.u) * Rt4[i]
        Cn, _ = s4.interp_u(y4[:, i], s4.u)
        ax[0].plot(Rn * 100, Cn, lw=1.6, color=shades[k],
                   label="%.0f h (R=%.3f)" % (tt, Rt4[i] * 100))
        ax[0].plot([Rn[-1] * 100], [Cn[-1]], "o", ms=5.5, color=shades[k],
                   mec="white", mew=0.7, zorder=5)
        ends[tt] = float(Rn[-1] * 100)          # 该时刻的表面半径 [cm]
    r6, r48 = ends[6], ends[48]
    # 收缩的视觉证据：把 6 h 与 48 h 的两个表面半径用一条双向箭头连起来。
    # 版面：纵轴上方专门留一条空带放结论文字与图例，双向箭头放在 6 h 与 12 h 两条
    # 剖面之间的空白带（实测该带内无数据），三者的位置都按"不与任何曲线相交"挑的。
    ax[0].set_ylim(-0.05, 2.45)
    # 结论文字分三行并左对齐：单行时它有 30 多个字符、宽度超过绘图区，左端越出坐标框
    # 0.36 in（印后 26 pt）、压到 y 轴刻度上；右对齐拆两行后，第一行宽 2.10 in 仍伸进
    # 右上角的图例框（实测重叠 0.077 × 0.133 in）。三行后最长一行只有约 1.1 in。
    ax[0].text(0.02, 2.38,
               "收缩：表面半径\n$%.3f\\to%.3f$ cm\n（左移 $%.3f$ cm，$%.1f\\%%$）"
               % (r6, r48, r6 - r48, (r6 - r48) / r6 * 100),
               fontsize=9.5, color="#08306b", ha="left", va="top", linespacing=1.35)
    ybar = 0.95
    ax[0].annotate("", xy=(r48, ybar), xytext=(r6, ybar),
                   arrowprops=dict(arrowstyle="<|-|>", color="#08306b", lw=1.5,
                                   mutation_scale=11, shrinkA=0, shrinkB=0))
    ax[0].text((r6 + r48) / 2, ybar + 0.05, "$\\Delta R$", fontsize=9.5,
               color="#08306b", ha="center", va="bottom")
    ax[0].set_xlabel("当前物理距离 $r$ / cm"); ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_title("(a) 物理坐标下的水分剖面（随收缩移动）")
    ax[0].legend(fontsize=7.5, loc="upper right", bbox_to_anchor=(1.0, 0.92),
                 title="干燥时间 (当时半径)", title_fontsize=7.5)
    _glow(ax[1], t4 / 3600, Rt4 * 100, COL_MOIST, lw=1.9)
    # 实测半径按 0.001 cm 步进下降，20 h 之后基本落定；这一段是问题 4 的几何前提
    ax[1].axvspan(20, t4[-1] / 3600, color=COL_PHASE, alpha=0.075)
    ax[1].text((20 + t4[-1] / 3600) / 2, Rt4.min() * 100 + 0.035,
               "20 h 后收缩基本停止", fontsize=9.5, color=COL_PHASE, ha="center")
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("$R(t)$ / cm")
    ax[1].set_title("(b) 收缩中的表面半径")
    save(fig, "fig18_p4_physics.png")

    # ---------------- fig19：物质坐标下的剖面与论文表 6 的采样点 ----------------
    # (a) 横轴 物质坐标 u=xi^2（域固定为 0~1）、纵轴 C/(kg/kg)：画 6、12、24、48 h 四条
    # 剖面，与 fig18(a) 是同一批数据的两种坐标画法；
    # (b) 横轴 时间/h、纵轴 C/(kg/kg)：取 r=0、0.5、1.0、1.5 cm 四个物理距离在四个时刻的
    # 浓度值，即论文表 6 的采样；r > R(t) 处该位置不存在，不存点故曲线留空。
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    for tt in [6, 12, 24, 48]:
        i = int(np.argmin(np.abs(t4 - tt * 3600)))
        Rn = np.sqrt(s4.u) * Rt4[i]
        Cn, _ = s4.interp_u(y4[:, i], s4.u)
        ax[0].plot(s4.u, Cn, lw=1.5, label="%.0f h" % tt)
    ax[0].set_xlabel("物质坐标 $u=\\xi^2$"); ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_title("(a) 物质坐标下的剖面（域固定）"); ax[0].legend(fontsize=8)
    # 表6 采样
    rcs = np.array([0.0, 0.005, 0.010, 0.015])
    vals = np.full((len(rcs), 4), np.nan)
    for tt in [6, 12, 24, 48]:
        i = int(np.argmin(np.abs(t4 - tt * 3600)))
        for j, r in enumerate(rcs):
            if r <= Rt4[i]:
                v, _ = s4.interp_u(y4[:, i], (r / Rt4[i]) ** 2)
                vals[j, [6, 12, 24, 48].index(tt)] = v[0]
    for j, r in enumerate(rcs):
        ax[1].plot([6, 12, 24, 48], vals[j], "o-", ms=4, lw=1.4,
                   label="r=%.1f cm" % (r * 100))
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("C / (kg/kg)")
    ax[1].set_title("(b) 表6 各物理距离处的干燥曲线"); ax[1].legend(fontsize=8)
    save(fig, "fig19_p4_profiles.png")

    print("[fig20] 问题3 与问题4 对比（备用图，论文未采用）")
    # ---------------- fig20（备用）：问题 3 与问题 4 的净差异 ----------------
    # (a) 横轴 时间/h、纵轴 中心水分浓度 C/(kg/kg) 半对数：实线是问题 3，虚线是问题 4，
    # 红点线为判据 0.15，两条竖点线分别是 td3（蓝）与 td4（橙）；
    # (b) 纵轴 t_dry/h，两根柱对比两问的烘干时长，柱顶标出小时数。

    # 该图把问题 3 与问题 4 直接相减，而两者的物性关系与几何同时改变，
    # 其差值不是"纯收缩效应"；论文正式使用 fig20_p4_controlled.png 的受控对照。
    # 此处标题只陈述"净差异"，不做因果表述。
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].semilogy(t3 / 3600, C3[:, 0], lw=1.8, label="问题3（附录3物性、不收缩）")
    ax[0].semilogy(t4 / 3600, C4[:, 0], lw=1.8, ls="--", label="问题4（附录4物性、含收缩）")
    ax[0].axhline(CRIT, color="r", ls=":", lw=1.2)
    ax[0].axvline(td3 / 3600, color="C0", ls=":", lw=1.2)
    ax[0].axvline(td4 / 3600, color="C1", ls=":", lw=1.2)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("中心水分浓度 / (kg/kg)")
    ax[0].set_title("(a) 物性与几何同时改变后的净差异 %.1f%%"
                    % ((td3 - td4) / td3 * 100)); ax[0].legend(fontsize=8)
    ax[1].bar(["问题3\n$R\\equiv2$ cm", "问题4\n含收缩"],
              [td3 / 3600, td4 / 3600], color=["C0", "C1"])
    for i, v in enumerate([td3 / 3600, td4 / 3600]):
        ax[1].text(i, v + 0.6, "%.2f h" % v, ha="center", fontsize=10)
    ax[1].set_ylabel("$t_{\\mathrm{dry}}$ / h"); ax[1].set_ylim(0, 65)
    ax[1].set_title("(b) 烘干时长对比（净差异，非纯收缩）")
    save(fig, "fig20_p4_vs_p3.png")

    # ---------------- 第 7 组：三层验证、收敛性与误差分析（fig21~fig25） ----------------
    print("[fig21] 三族离散互证")
    # ---------------- fig21：谱配置法 / 有限体积 / 闭式解三族离散互证 ----------------
    # 前两族都在 t=1800 s 报出表面水分浓度 C(R_0,1800) [kg/kg]，闭式解仍是基准。
    # (a) 横轴 四组均匀网格 FV 的配置（N=200/200/400/400，dt=1.0/0.5/0.5/0.25 s）、纵轴
    # C(R_0,1800)/(kg/kg)：四种网格都落在红虚线（谱配置法的 1.510918）上；
    # (b) 横轴 与闭式解析解的最大偏差/K（对数）、纵轴 三条方法：闭式解按定义为基准，
    # 谱配置法 1.1e-10、均匀网格 FV 2.6e-5。
    from fv_uniform import solve_fv, sample_xi
    lab, vals_spec, vals_fv = [], [], []
    xis = []; fvv = []
    for Nn, dtv in [(200, 1.0), (200, 0.5), (400, 0.5), (400, 0.25)]:
        xi, Ctab, Ttab = solve_fv(Nn, dtv, 1800.0, prob=1, t_report=[1800.0])
        fvv.append(float(np.atleast_1d(sample_xi(xi, Ctab[0], 1.0))[0]))
        lab.append("$N$=%d\n$\\Delta t$=%.2f" % (Nn, dtv))
    # 画布与布局：(b) 的纵轴刻度是两行中文方法名（"均匀网格/守恒型 FV"、"闭式解析解/
    # (问题1 温度)"），默认间距下会压进 (a) 的绘图区（实测 0.153 in），故用 constrained
    # 布局自动按最宽的刻度留位；画布宽 8.2 in 按"存出 PNG 宽 ÷ 版心宽"反推，印后约 8 pt。
    fig = plt.figure(figsize=(8.2, 3.05), layout="constrained")
    ax = fig.subplots(1, 2)
    SPEC = 1.5109181                      # 本文谱配置法在该点的值（与表 7 同口径）
    ax[0].plot(range(len(fvv)), fvv, "o-", lw=1.5, label="均匀网格守恒型 FV")
    ax[0].axhline(SPEC, color="r", ls="--", lw=1.6, label="本文谱配置法 1.510918")
    # 刻度原来是 45° 倾斜的长串（"N=200,dt=1.00"），是全篇最挤的一处；拆成两行短标签
    ax[0].set_xticks(range(len(lab))); ax[0].set_xticklabels(lab, rotation=0, fontsize=8.5)
    ax[0].set_ylabel("$C(R_0,1800)$ / (kg/kg)")
    ax[0].set_title("(a) 两族离散收敛到同一值"); ax[0].legend(fontsize=8)
    ax[0].set_ylim(1.51085, 1.51115)
    # 把结论量直接写在图上：四组网格与谱解的最大相对差，用当前 fvv 现算，不写死。
    # 用"相对差"而不是绝对极差：绝对极差 7.9e-5 kg/kg 单看数字会比实际收敛质量显得差。
    rel = max(abs(np.array(fvv) - SPEC)) / SPEC
    ax[0].annotate("四组网格与谱解\n最大相对差 $%.1f\\times10^{-5}$" % (rel * 1e5),
                   xy=(0.02, 0.97), xycoords="axes fraction", fontsize=9.5,
                   ha="left", va="top", color="#08306b",
                   bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#08306b", lw=0.8))
    # 闭式解析解是"偏差"的基准，其偏差恰为 0；对数轴上 0 高的柱不可见，
    # 故改用水平柱 + 一个下限刻度，并把"0（基准）"直接写在图上。
    FLOOR = 1e-12
    names = ["闭式解析解\n(问题1 温度)", "谱配置法", "均匀网格\n守恒型 FV"]
    errs = [FLOOR, 1.1e-10, 2.6e-5]
    ax[1].barh(names, errs, color=["0.7", "C0", "C1"], height=0.55)
    ax[1].set_xscale("log"); ax[1].set_xlim(FLOOR / 2, 1e-3)
    ax[1].set_xlabel("与闭式解的最大偏差 / K")
    ax[1].set_title("(b) 三族方法的精度量级")
    # 纵轴上限抬到 2.7：数值标签改为"条右端上方 6 pt"，需要这 0.15 的顶部余量。
    ax[1].set_ylim(-0.6, 2.7)
    ax[1].tick_params(axis="y", pad=5.0)        # 方法名与色条左端之间留出约 5 pt
    # 数值标签原来贴在色条右端外侧：对数轴上 1.3 倍的位置偏移折算下来离条端只有约
    # 1.6 pt（看上去与色条粘在一起），而最长的 $2.6\times10^{-5}$ 还越出坐标框 0.18 in。
    # 改成"居中于条右端上方 7 pt"：标签不再出框（窄的谱配置法那条左端距框也还有余量），
    # 文字中心正对条端，数值与条长的对应关系比外侧更直接。
    ax[1].text(FLOOR * 1.6, 0, "0（定义基准）", va="center", fontsize=9.5, color="0.35")
    for yy, v, s in [(1, 1.1e-10, "$1.1\\times10^{-10}$"), (2, 2.6e-5, "$2.6\\times10^{-5}$")]:
        ax[1].annotate(s, xy=(v, yy), xytext=(0, 7), textcoords="offset points",
                       ha="center", va="bottom", fontsize=9.5)
    save(fig, "fig21_three_methods.png")

    print("[fig22] 收敛性（数据来自 中间数据/verification.json，不写死）")
    # ---------------- fig22：空间与时间收敛性（数据取自 中间数据/verification.json） ----------------
    # L 取最细网格的 t_dry 当作连续极限；经典 FV 的 O(dt) 时间误差由 VD.time_slope() 的
    # 斜率乘 dt 扣掉，不扣的话空间阶会被时间误差盖住。
    # (a) 横轴 节点数 N、纵轴 L-t_dry/s（双对数）：三角是谱解、方块是经典 FV+theta，另有
    # 斜率 -2 的参考线；(b) 横轴 时间步 dt/s（坐标轴反向）、纵轴 L-t_dry/s：一阶外推极限
    # L_fv 由 t(dt)=L-c*dt 对最细三档做最小二乘得到，黑虚线是斜率 1 的参考线。
    vd = VD.load()
    VD.check_consistency(vd, tdry_cached=td3)
    Ns = list(vd["space"]["N"])
    tds = np.asarray(vd["space"]["tdry"], dtype=float)
    Nf = list(vd["fv_space"]["N"])
    tdf = np.asarray(vd["fv_space"]["tdry"], dtype=float)
    c_fv = VD.time_slope(vd)                 # 该格式的 O(dt) 系数，用于扣时间误差
    dt_fv = vd["fv_space"]["dt"]
    LIM = float(tds.max())                   # 谱解的连续极限（取最细网格值）
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].loglog(Ns, np.maximum(LIM - tds, 0.01), "^-", lw=1.7, ms=6,
                 label="本文谱配置法")
    ax[0].loglog(Nf, LIM - (tdf + c_fv * dt_fv), "s-", lw=1.7, ms=6,
                 label="经典 FV+$\\theta$ 格式（扣去时间误差 %.2f s）" % (c_fv * dt_fv))
    ax[0].loglog([Nf[0], Nf[-1]],
                 (LIM - (tdf[-1] + c_fv * dt_fv)) * (np.array([Nf[0], Nf[-1]],
                                                             dtype=float) / Nf[-1]) ** -2,
                 "k--", lw=1.1, label="斜率 $-2$")
    ax[0].set_xlabel("节点数 $N$"); ax[0].set_ylabel("$L-t_{\\mathrm{dry}}$ / s")
    ax[0].set_title("(a) 空间收敛：谱解远快于二阶")
    # 谱解的最后一点已压到时间误差平台（纵轴下端），默认下半轴没有空位，
    # 图例放左下角会盖住该点；把纵轴下限再放低约一个半数量级腾出空白。
    ax[0].set_ylim(bottom=3e-4)
    ax[0].legend(fontsize=8, loc="lower left")
    # 把结论量直接写在图上：谱解最后一点已落进时间误差平台，而经典 FV 还差三个数量级。
    # 数值取 (a) 里谱解曲线真正画出的末点值——注意不能直接用 LIM - tds[-1]：LIM 就是取
    # tds 的最大值，那一项恒为 0，而 y=0 在对数轴上不可表示，标注会整条画不出来（实测）。
    plateau = float(np.maximum(LIM - tds, 0.01)[-1])
    ax[0].annotate("谱解误差已落到\n$%.1f\\times10^{-2}$ s 的平台" % (plateau * 1e2),
                   xy=(0.97, 0.96), xycoords="axes fraction", fontsize=9.5,
                   ha="right", va="top", color="#08306b",
                   bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#08306b", lw=0.8))
    dtl = np.asarray(vd["time"]["dt"], dtype=float)
    tdl = np.asarray(vd["time"]["tdry"], dtype=float)
    # 一阶外推极限：由 t(dt)=L-c*dt 的最小二乘截距给出，而不是写死的常数
    sl, L_fv = np.polyfit(dtl[-3:], tdl[-3:], 1)
    print("  经典 FV+theta 外推极限 L=%.3f s，c=%.4f s/s" % (L_fv, -sl))
    ax[1].loglog(dtl, L_fv - tdl, "o-", lw=1.7, ms=6, label="经典 FV+$\\theta$ 格式")
    ax[1].loglog(dtl, (L_fv - tdl[-1]) * (dtl / dtl[-1]), "k--", lw=1.1, label="斜率 1")
    ax[1].set_xlabel("时间步 $\\Delta t$ / s"); ax[1].set_ylabel("$L-t_{\\mathrm{dry}}$ / s")
    ax[1].set_title("(b) 时间收敛：严格一阶")
    ax[1].invert_xaxis()
    # 实测斜率直接标在图上（与 -2 参考线对照）。单位必须写成 s/s：这个 1.42 是
    # t(dt)=L-c*dt 的斜率 c，而正文另有一个"相邻增量-步长双对数斜率 1.01"（收敛阶），
    # 两个数不是一回事，标上单位才不会让人以为图上写错了。
    ax[1].annotate("最小二乘斜率 $c=%.4f$ s/s" % (-sl), xy=(0.03, 0.96),
                   xycoords="axes fraction", fontsize=9.5, ha="left", va="top",
                   color="#08306b",
                   bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#08306b", lw=0.8))
    ax[1].legend(fontsize=8, loc="lower left")
    save(fig, "fig22_convergence.png")

    print("[fig23] 积分器与容差稳定性")
    # ---------------- fig23：求解器与容差配置的稳定性 ----------------
    # 数据取自 verification.json 的 integ 组，共 6 种配置（BDF 1e-9/1e-10/1e-12、Radau
    # 1e-11、LSODA 1e-10 以及 1 s 密采样）各自报出的 t_dry [s]。
    # (a) 横轴 六种配置、纵轴 t_dry/s，数值直接标在点上；(b) 横轴 t_dry 相对均值的偏差/s、
    # 纵轴 频数。两幅标题都只报极差。
    vvv = [r["tdry"] for r in vd["integ"]]
    cfgs = ["BDF\n1e-9", "BDF\n1e-10", "BDF\n1e-12", "Radau\n1e-11", "LSODA\n1e-10", "1s密采样"]
    spread = max(vvv) - min(vvv)
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(range(len(vvv)), vvv, "o-", lw=1.5, ms=7)
    ax[0].set_xticks(range(len(cfgs))); ax[0].set_xticklabels(cfgs, fontsize=8)
    ax[0].set_ylabel("$t_{\\mathrm{dry}}$ / s")
    ax[0].set_title("(a) 六种积分器/容差配置：极差 %.3f s" % spread)
    ax[0].set_ylim(min(vvv) - 0.02, max(vvv) + 0.02)
    for i, v in enumerate(vvv):
        # 首末标签靠内对齐，否则会压到 y 轴刻度数字上；
        # 相邻两点数值只差 0.0001 s，标签上下错开以免互相叠字。
        ha = "left" if i == 0 else ("right" if i == len(vvv) - 1 else "center")
        ax[0].text(i, v + (0.004 if i % 2 == 0 else 0.009), "%.3f" % v,
                   ha=ha, fontsize=6.0)
    ax[0].margins(x=0.07)
    ax[1].hist(np.array(vvv) - np.mean(vvv), bins=8, color="C0", alpha=.8)
    ax[1].set_xlabel("$t_{\\mathrm{dry}}$ 相对均值的偏差 / s"); ax[1].set_ylabel("频数")
    # 注意：极差 0.007 s 不等于"各点落在均值 ±0.004 s 内"——6 个配置里 5 个聚在
    # 一侧、1 个偏离，相对均值的最大偏差可达 0.0056 s。此处只陈述极差。
    ax[1].set_title("(b) 相对均值的偏差：极差 %.3f s" % spread)
    save(fig, "fig23_integrator.png")

    print("[fig24] 谱解-闭式解残差统计")
    # ---------------- fig24：谱解与闭式解偏差的统计特征 ----------------
    # 残差 res 沿用 fig11 里算出的"谱解 − 闭式解"（单位 K），这里放大 1e10 后统计。
    # (a) 横轴 残差/(1e-10 K)、纵轴 概率密度：柱状为直方图，红实线为同均值同标准差的
    # 正态密度；(b) Q-Q 图，横轴 理论分位、纵轴 样本分位。
    res2 = res * 1e10
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].hist(res2, bins=26, color="C0", alpha=.8, density=True)
    xs = np.linspace(res2.min(), res2.max(), 200)
    from scipy import stats as st
    ax[0].plot(xs, st.norm.pdf(xs, res2.mean(), res2.std()), "r-", lw=1.8)
    ax[0].set_xlabel("残差 / ($10^{-10}$ K)"); ax[0].set_ylabel("概率密度")
    ax[0].set_title("(a) 谱解与闭式解偏差的分布")
    st.probplot(res2, dist="norm", plot=ax[1])
    ax[1].set_xlabel("理论分位"); ax[1].set_ylabel("样本分位")
    # 残差由时间积分的系统性偏差主导（201 个点同号、Shapiro-Wilk p~1e-12），
    # 不应表述为"与正态无偏离"；此处只描述其确定性特征。
    ax[1].set_title("(b) Q-Q 图：确定性离散误差，非随机噪声")
    ax[1].get_lines()[0].set_marker("."); ax[1].get_lines()[0].set_markersize(5)
    ax[1].get_lines()[1].set_color("r")
    save(fig, "fig24_residual_qc.png")

    print("[fig25] 方法流程图")
    # ---------------- fig25：技术路线流程图 ----------------
    # 本图没有数据来源，纯手工排版的方框 + 箭头，坐标都用 axes fraction 给出，与任何
    # 物理量无关。方框自上而下是：题目条件 + 附件 1/2 数据 → 几何简化 → 控制方程 →
    # 环境标定 / 坐标正则化 / 物性经验式 → 问题 1~4 的求解路线 → 三层验证 →
    # 灵敏度与结论。
    # 底色承载语义（三类含义固定的浅底，与全篇语义色一致；同一类节点在任何一版里同色）：

    #   #cfe2f3 蓝   = 数据、模型与求解流程（主干）
    #   #9dc3e6 深蓝 = 两个关键节点（控制方程、三层验证），边框加粗以形成视觉层级
    #   #fdf1cf 米黄 = 三项建模准备（中间产物，对应参考图里的黄色强调框）
    # 三处必须守住的版面底线：
    #   (1)"最后粉色的框框下面的黑线没有了"——原代码最后一个框中心 y=0.04、高 0.11，
    #     底边落在 y=-0.015 已被坐标轴裁掉（Rectangle 默认 clip_on=True）。新版把整列

    #     重排成 7 行（0.940 / 0.812 / 0.684 / 0.535 / 0.370 / 0.212 / 0.062），
    #     最低一行的底边在 y=0.018，四周都留出余量，不会再被裁。
    #   (2)"这个箭头必须有 - 再箭头"——原相邻框间距只有 0.03 axes fraction（约 0.16 in），
    #     箭头只看得见一个三角。新版行距放到 0.04~0.077，箭头 lw=1.6、mutation_scale=13，
    #     杆与头都清楚可见。
    #   (3)"不要怕占地方、要圆角"——直角矩形换 FancyBboxPatch 圆角框，方框文字 9.5→10.5 pt。

    # 画布宽仍是 9.6 in（不加大）：插图按 0.92\textwidth 排版，画布一加宽缩放系数就变小，
    # 印后字号反而更小；要"字更大"只能把字号做大于画布。画布只加高到 5.9 in，给行距让位。
    fig, ax = plt.subplots(figsize=(9.6, 5.9))
    ax.axis("off")
    BLUE, BLUE_D, YELLOW = "#cfe2f3", "#9dc3e6", "#fdf1cf"
    EDGE, ACCENT = "#2e75b6", "#bf9000"
    ROUND = "round,pad=0,rounding_size=0.012"
    rows = [
        (0.5, 0.940, 0.76, "题目条件 + 附件 1/2 实测数据", BLUE, EDGE, 1.6),
        (0.5, 0.812, 0.76, "几何简化：$L/R_0=12.5\\ \\Rightarrow$ 一维轴对称", BLUE, EDGE, 1.6),
        (0.5, 0.684, 0.76, "控制方程：热传导 + Fick 扩散 + Robin 边界", BLUE_D, EDGE, 2.2),
        (0.185, 0.535, 0.26, "环境标定\n一阶惯性拟合", YELLOW, ACCENT, 1.4),
        (0.500, 0.535, 0.26, "坐标正则化\n$u=(r/R)^2$", YELLOW, ACCENT, 1.4),
        (0.815, 0.535, 0.26, "物性经验式\n附录 2/3/4", YELLOW, ACCENT, 1.4),
        (0.185, 0.370, 0.26, "问题 1\n解耦 + 闭式解", BLUE, EDGE, 1.6),
        (0.500, 0.370, 0.26, "问题 2\n变物性耦合", BLUE, EDGE, 1.6),
        (0.815, 0.370, 0.26, "问题 3/4\n事件求交 + 收缩", BLUE, EDGE, 1.6),
        (0.5, 0.212, 0.76, "三层验证：闭式解析解 / 谱配置 / 两族有限体积",
         BLUE_D, EDGE, 2.2),
        (0.5, 0.062, 0.76, "灵敏度与误差分析 $\\Rightarrow$ 结论与工艺建议", BLUE, EDGE, 1.6),
    ]
    BH = 0.088                                    # 方框高度（axes fraction），全列统一
    for (x, yv, w, txt, fc, ec, lw) in rows:
        ax.add_patch(FancyBboxPatch((x - w / 2, yv - BH / 2), w, BH,
                                    boxstyle=ROUND, facecolor=fc, edgecolor=ec,
                                    lw=lw, transform=ax.transAxes, zorder=2))
        ax.text(x, yv, txt, ha="center", va="center", fontsize=10.5, linespacing=1.35,
                transform=ax.transAxes, zorder=3)
    # 三个"问题"框用虚线分组框圈起来，与参考图的分组画法一致（纯装饰，不含数值）；
    # 标签放在分组框左下角外，避开所有箭头（最左的收敛箭头从 x=0.185 才起步）
    ax.add_patch(FancyBboxPatch((0.045, 0.370 - BH / 2 - 0.030), 0.915, BH + 0.060,
                                boxstyle=ROUND, facecolor="none", edgecolor="0.55",
                                lw=1.1, ls=(0, (5, 3)), transform=ax.transAxes, zorder=1))
    ax.text(0.050, 0.370 - BH / 2 - 0.038, "四问求解路线", fontsize=9.5, color="0.40",
            ha="left", va="top", transform=ax.transAxes, zorder=1)
    arr = dict(arrowstyle="-|>", color=ACCENT, lw=1.6, shrinkA=0, shrinkB=0,
               mutation_scale=13)
    # 箭头一律从上一行的底边画到下一行的顶边，故杆长就是行距，肉眼可见
    A = [(0.5, 0.896, 0.5, 0.856), (0.5, 0.768, 0.5, 0.728),
         (0.5, 0.640, 0.185, 0.579), (0.5, 0.640, 0.5, 0.579), (0.5, 0.640, 0.815, 0.579),
         (0.185, 0.491, 0.185, 0.414), (0.5, 0.491, 0.5, 0.414), (0.815, 0.491, 0.815, 0.414),
         (0.185, 0.326, 0.5, 0.256), (0.5, 0.326, 0.5, 0.256), (0.815, 0.326, 0.5, 0.256),
         (0.5, 0.168, 0.5, 0.106)]
    for (x1, y1_, x2, y2_) in A:
        ax.annotate("", xy=(x2, y2_), xytext=(x1, y1_),
                    arrowprops=arr, xycoords="axes fraction")
    ax.set_title("本文技术路线", fontsize=13.5, pad=6)
    save(fig, "fig25_flow.png")

    print("[fig26] 物理模型与坐标系示意图")
    # ---------------- fig26：物理模型与坐标系示意图（无数据，纯示意） ----------------
    # 全篇此前没有一张物理示意图——翻到第五章会想知道"这个问题被想象成什么形状、
    # 边界加在哪里、收缩体现在哪个量上"，原先只能靠读文字。本图把这条建模主线画出来：
    # (1)药材横截面（热风对流 + 表面蒸发 + 收缩后的 R(t)）；(2)取 1/4 截面对称化，得到
    # 一维径向坐标 r 属于 [0, R(t)]；(3)控制方程与 Robin 边界写在右侧。
    # 本图不含任何数值计算，画的都是模型假设本身；图上虚线圆用的 0.60 是附件 2 实测

    # 收缩的终值比例（R 由 2.00 cm 降到 1.20 cm），不是拟合或外推出来的。
    # 坐标一律用 data 单位并 set_aspect("equal")：axes fraction 下圆会被拉成椭圆。
    # 版面纪律：三块内容各自留出"文字带"，说明文字一律放在图形外，不与图元交叠。
    fig, ax = plt.subplots(figsize=(8.2, 3.45))
    ax.set_xlim(0, 12.6); ax.set_ylim(0, 5.30)
    ax.set_aspect("equal"); ax.axis("off")

    # ---- (1) 药材横截面：热风、蒸发、收缩 ----
    # 外圆由 1.05 放大到 1.15：内圆半径固定为 0.60 R（附件 2 实测收缩终值比例，不能改），
    # 原尺寸下内圆直径只有 47.8 pt，而圆内"（湿物料）"一行在 9.5 pt 下就要 47.5 pt，
    # 两端正好压在红虚线圆上。放大外圆并把该行降到 9 pt（印后 8.5 pt，仍过 7 pt 线）后，
    # 文字两端各留约 3.6 pt；再加一层白色描边，即使与虚线相切也读得清。
    cx, cy, Ro = 1.85, 3.20, 1.15
    Ri = Ro * 0.60                       # R(t_end)/R_0 = 1.20/2.00，取自附件 2 实测收缩
    ax.add_patch(plt.Circle((cx, cy), Ro, facecolor="#eaf3fb", edgecolor="#2e75b6",
                            lw=1.8, zorder=2))
    ax.add_patch(plt.Circle((cx, cy), Ri, facecolor="none", edgecolor="#c0392b",
                            lw=1.6, ls=(0, (4, 3)), zorder=3))
    ax.text(cx, cy, "药材\n（湿物料）", ha="center", va="center", fontsize=9.0, zorder=4,
            path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])
    # 热风：从左侧斜向吹来，箭头止于表面
    for (dx, dy) in [(-1.60, 0.75), (-1.75, -0.10), (-1.30, -0.90)]:
        s = np.hypot(dx, dy)
        ax.annotate("", xy=(cx + dx / s * Ro, cy + dy / s * Ro),
                    xytext=(cx + dx, cy + dy),
                    arrowprops=dict(arrowstyle="-|>", color=COL_TEMP, lw=1.8,
                                    mutation_scale=13, shrinkA=0, shrinkB=0))
    # 蒸发：从表面向外
    for ang in (-30, 5, 45):
        a = np.deg2rad(ang)
        ax.annotate("", xy=(cx + np.cos(a) * (Ro + 0.55), cy + np.sin(a) * (Ro + 0.55)),
                    xytext=(cx + np.cos(a) * Ro, cy + np.sin(a) * Ro),
                    arrowprops=dict(arrowstyle="-|>", color=COL_MOIST, lw=1.6,
                                    mutation_scale=12, shrinkA=0, shrinkB=0))
    ax.text(0.10, 5.15, "热风：$T_{\\mathrm{air}}(t),\\,C_{\\mathrm{air}}(t)$",
            ha="left", va="top", fontsize=9.5, color=COL_TEMP)
    ax.text(4.60, 5.15, "表面蒸发\n（Robin 边界）", ha="right", va="top",
            fontsize=9.5, color=COL_MOIST)
    # 收缩：外圆 → 内虚线圆，箭头画在右下 45°（该方向既无热风入射也无蒸发射出）
    ax.annotate("", xy=(cx + Ri * 0.70, cy - Ri * 0.70),
                xytext=(cx + Ro * 0.70, cy - Ro * 0.70),
                arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.6,
                                mutation_scale=12, shrinkA=0, shrinkB=0))
    ax.text(cx, 1.82, "$R(t)$：表面半径随含水率收缩", ha="center", va="top",
            fontsize=9.5, color="#c0392b")
    ax.text(cx, 1.34, "(a) 圆柱药材横截面（长 $L=25$ cm）", ha="center", va="top",
            fontsize=9.5, color="0.25")

    # ---- (2) 取 1/4 截面对称化 → 一维径向 ----
    # 圆心高度与 (a) 对齐（cy=3.20），两块示意图才在同一水平线上
    mx, my, Rm = 5.95, 3.20, 1.15
    ax.add_patch(plt.Circle((mx, my), Rm, facecolor="none", edgecolor="0.55",
                            lw=1.3, ls=(0, (4, 3)), zorder=2))
    # 两条对称轴用点线画出：1/4 截面就是靠它们切出来的
    ax.plot([mx, mx], [my, my + Rm], ls=":", lw=1.0, color="0.6")
    ax.plot([mx, mx + Rm], [my, my], ls=":", lw=1.0, color="0.6")
    ax.annotate("", xy=(mx + Rm, my), xytext=(mx, my),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.8,
                                mutation_scale=13, shrinkA=0, shrinkB=0))
    ax.plot([mx], [my], "o", ms=5, color="k", zorder=4)
    ax.text(mx + Rm + 0.10, my + 0.12, "$r=R(t)$", fontsize=9.5, va="bottom")
    ax.text(mx + Rm * 0.45, my - 0.16, "$r$", fontsize=11, ha="center", va="top")
    ax.text(mx, 1.82, "$r=0$：对称（$\\partial_rT=\\partial_rC=0$）", ha="center",
            va="top", fontsize=9.5, color="0.35")
    ax.text(mx, 1.34, "(b) 化为轴对称一维径向问题", ha="center", va="top",
            fontsize=9.5, color="0.25")

    # ---- (3) 控制方程与边界条件 ----
    ax.add_patch(FancyBboxPatch((7.75, 0.30), 4.70, 4.70,
                                boxstyle="round,pad=0,rounding_size=0.16",
                                facecolor="#f4f8fc", edgecolor="#2e75b6", lw=1.5))
    ex = 10.10
    ax.text(ex, 4.74, "(c) 控制方程与边界条件", ha="center", va="center",
            fontsize=9.5, color="0.25")
    # 每一节"小标题在上、公式在下"，左右都不加标签，避免长公式越过小标题的左边界
    eqs = [
        ("表观体积热容模型",
         r"$C_v(C)\,\partial_tT=r^{-1}\partial_r\!\left(rk(C)\,\partial_rT\right)$"),
        ("干物质基准守恒",
         r"$\rho_d\,\partial_tC=r^{-1}\partial_r\!\left(r\rho_dD\,\partial_rC\right)$"),
        ("边界条件",
         r"$r=0$：$\partial_rT=\partial_rC=0$"),
        ("",
         r"$r=R(t)$：$-k\,\partial_rT=h\,(T_{\mathrm{air}}-T_s)$"),
        ("",
         r"$-D\,\partial_rC=h_m(C_{\mathrm{air}}-C_s)$"),
    ]
    yy = 4.40
    # 公式块整体上移 0.12：原起点 4.28 使最后一行 $-D\partial_rC=\dots$ 距框底边只剩
    # 2 pt（实测），紧得像是被裁掉；上移后首行距 (c) 标题 12 px（约 3.6 pt）、末行距框底
    # 净留白约 25 px，上下都不再贴边。
    for (tag, eq) in eqs:
        if tag:
            ax.text(ex, yy, tag, fontsize=9.5, ha="center", va="center", color="0.30")
            yy -= 0.46
        ax.text(ex, yy, eq, fontsize=9.5, ha="center", va="center")
        yy -= 0.60
    save(fig, "fig26_model_schematic.png")

    print("全部图完成，输出目录：", FIGDIR)
    return td3, td4


if __name__ == "__main__":
    main()




