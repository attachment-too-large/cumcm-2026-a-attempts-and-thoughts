# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------
# make_figures.py —— 绘制论文全部 21 幅插图（200 dpi PNG，输出到 figures/）
#
# 本脚本不做任何物理计算，只把算好的结果读进来画图。数据来源：
#   data/p1_fields.npz   问题 1（预热阶段）：时间轴 t_fine[s]、半径 radii_cm[cm]、
#                        温度场 T_fine、水分浓度场 C_fine 及抽样表 *_tab
#   data/p23_fields.npz  问题 2、3（整段烘干、不收缩）：t[s]、radii_cm[cm]、
#                        T、C、烘干时刻 t_dry；T/C 的行是时间、列是半径
#   data/p4_fields.npz   问题 4（收缩）：物质坐标场 C/T（列按 xi_mat 排列）、
#                        物理坐标场 C_phys、半径历史 R[m]、t_dry
#   data/sensitivity.csv 灵敏度分析结果，由 sensitivity.py 生成（图 18、19 用）
#   附件 1、附件 2       烘房温湿度实测序列与药材半径实测序列，经 common 读取
#
# 每个 figNN_xxx() 函数画一幅 figures/figNN_xxx.png，函数名与文件名一一对应；
# 每幅图画的是什么、取自哪个键，写在各自函数的 docstring 里。
# 运行方式：先跑 p1_preheat.py、p23_process.py、p4_shrinkage.py、sensitivity.py
# 生成上述数据文件，再执行 python make_figures.py；main() 按顺序调用各绘图函数。
# 中文字体、负号与出图分辨率由 common.setup_matplotlib() 统一设置。
# --------------------------------------------------------------------------
"""
论文全部 21 幅插图的生成入口：数据来源、输出位置与运行方式见上方文件头注释。
"""

import os
import sys
import io
import numpy as np
import pandas as pd

if sys.stdout.encoding is None or "utf" not in sys.stdout.encoding.lower():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm

plt = cm.setup_matplotlib()
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle
import matplotlib.gridspec as gridspec
from scipy import stats

C_TARGET = 0.15
D = cm.DATA_DIR
F = cm.FIG_DIR

# 画布整体缩放。插图在论文里按接近版心宽（约 14.5 cm）排版，而这里的画布宽是
# 10~13 in（约 25~33 cm），直接缩放排版会把坐标轴文字压到设计字号的 45% 以下，
# 印刷后只剩 4~5 pt、读不清。统一把画布缩到 0.8 倍、字号（10.5 pt）保持不变，
# 印刷后的实际字号约 6~9 pt，与期刊插图的常见量级一致。改动本常数后需要重新
# 运行本脚本，并逐图检查标签是否重叠。
FIG_SCALE = 0.80


def load(name):
    """读入 data/ 目录下的一个 npz 结果文件，各键的含义见使用它的绘图函数。"""
    return np.load(os.path.join(D, name))


# --------------------------------------------------------------------------
# fig01、fig02 —— 烘房环境这条边界条件本身，以及拟合式相对实测序列的残差诊断
# --------------------------------------------------------------------------
def fig01_room():
    """fig01：烘房环境边界条件，取附件 1 的实测序列（0~4 h，每 60 s 一点），与
    common.RoomConditions 的一阶惯性拟合式对照。左 (a) 烘房温度 T_air、右 (b) 烘房
    水分浓度 C_air：横轴均为时间/h，灰散点为附件 1 实测，实线为拟合式，横虚线是
    恒温干燥阶段的设定值（50.21 degC 与 0.0509 kg/kg），可见 4 h 后两者重合。
    """
    df = cm.load_attachment1()
    t = df["t"].to_numpy(float) / 3600.0
    Ta = df["T_air"].to_numpy(float)
    Ca = df["C_air"].to_numpy(float)
    room = cm.RoomConditions()
    tt = np.linspace(0, 4, 400)

    fig, ax = plt.subplots(1, 2, figsize=(9.6 * FIG_SCALE, 3.4 * FIG_SCALE))
    ax[0].plot(t, Ta, "o", ms=1.6, color="#888888", label="附件1 实测")
    ax[0].plot(tt, room.T_air(tt * 3600.0), "-", lw=1.8, color="#c0392b",
               label=r"拟合 $T_{air}=T_{set}-(T_{set}-28)e^{-t/\tau_T}$")
    ax[0].axhline(cm.T_SET_C, ls=":", lw=1.2, color="#2c3e50")
    ax[0].annotate(r"$T_{set}=50.21\,^\circ$C", xy=(2.6, 50.4), fontsize=9,
                   color="#2c3e50")
    ax[0].set_xlabel("时间 / h")
    ax[0].set_ylabel(r"烘房温度 $T_{air}$ / $^\circ$C")
    ax[0].set_title("(a) 烘房温度")
    ax[0].legend(loc="lower right", fontsize=8.5)

    ax[1].plot(t, Ca, "o", ms=1.6, color="#888888", label="附件1 实测")
    ax[1].plot(tt, room.C_air(tt * 3600.0), "-", lw=1.8, color="#1f6fb2",
               label=r"拟合 $C_{air}=C_{set}-(C_{set}-0.01963)e^{-t/\tau_C}$")
    ax[1].axhline(cm.C_SET, ls=":", lw=1.2, color="#2c3e50")
    ax[1].set_ylim(0.018, 0.0555)
    ax[1].annotate(r"$C_{set}=0.0509$ kg/kg", xy=(1.5, 0.0522), fontsize=9,
                   color="#2c3e50")
    ax[1].set_xlabel("时间 / h")
    ax[1].set_ylabel(r"烘房水分浓度 $C_{air}$ / (kg/kg)")
    ax[1].set_title("(b) 烘房水分浓度")
    ax[1].legend(loc="lower right", fontsize=8.5)
    return cm.savefig(fig, "fig01_room_conditions.png", plt)


def fig02_residuals():
    """fig02：拟合式相对附件 1 实测序列的残差诊断，2x3 共 6 个子图（残差 = 实测 - 拟合）。
    上行 (a)(b)(c) 是温度残差：时序（标出 RMSE）、直方图与正态密度、Q-Q 图；下行
    (d)(e)(f) 是水分浓度残差（放大 1e3 倍成 10^-3 kg/kg）的同样三张图，其中 (f) 还叠了
    温度与湿度残差的自相关函数（滞后 1~25 阶）与 ±1.96/sqrt(N) 置信带，检查白噪声性。
    """
    df = cm.load_attachment1()
    t = df["t"].to_numpy(float)
    Ta = df["T_air"].to_numpy(float)
    Ca = df["C_air"].to_numpy(float)
    room = cm.RoomConditions()
    rT = Ta - room.T_air(t)
    rC = Ca - room.C_air(t)

    fig = plt.figure(figsize=(9.8 * FIG_SCALE, 6.2 * FIG_SCALE))
    gs = gridspec.GridSpec(2, 3, hspace=0.42, wspace=0.32)

    # 子图位置：gs[0,0]=(a) 温度残差时序，gs[0,1]=(b) 直方图，gs[0,2]=(c) Q-Q 图；
    # 下一行 gs[1,0..2] 依次对应 (d)(e)(f)，画的是水分浓度残差的同样三张图
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(t / 3600.0, rT, lw=0.8, color="#c0392b")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title(r"(a) $T_{air}$ 残差 (RMSE=%.3f$^\circ$C)" %
                 np.sqrt((rT ** 2).mean()))
    ax.set_xlabel("时间 / h"); ax.set_ylabel("残差 / $^\\circ$C")

    ax = fig.add_subplot(gs[0, 1])
    ax.hist(rT, bins=30, color="#c0392b", alpha=0.75, density=True,
            edgecolor="white", linewidth=0.4)
    xs = np.linspace(rT.min(), rT.max(), 200)
    ax.plot(xs, stats.norm.pdf(xs, rT.mean(), rT.std()), "k--", lw=1.2)
    ax.set_title("(b) 温度残差直方图与正态密度")
    ax.set_xlabel("残差 / $^\\circ$C"); ax.set_ylabel("概率密度")

    ax = fig.add_subplot(gs[0, 2])
    stats.probplot(rT, dist="norm", plot=ax)
    ax.get_lines()[0].set_markerfacecolor("#c0392b")
    ax.get_lines()[0].set_markersize(2.5)
    ax.get_lines()[1].set_color("k")
    ax.set_title("(c) 温度残差 Q-Q 图")
    ax.set_xlabel("理论分位数"); ax.set_ylabel("样本分位数")

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(t / 3600.0, rC * 1e3, lw=0.8, color="#1f6fb2")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title(r"(d) $C_{air}$ 残差 (RMSE=%.2e kg/kg)" %
                 np.sqrt((rC ** 2).mean()))
    ax.set_xlabel("时间 / h"); ax.set_ylabel("残差 / $10^{-3}$ kg/kg")

    ax = fig.add_subplot(gs[1, 1])
    ax.hist(rC * 1e3, bins=30, color="#1f6fb2", alpha=0.75, density=True,
            edgecolor="white", linewidth=0.4)
    xs = np.linspace(rC.min(), rC.max(), 200) * 1e3
    ax.plot(xs, stats.norm.pdf(xs, rC.mean() * 1e3, rC.std() * 1e3), "k--",
            lw=1.2)
    ax.set_title("(e) 湿度残差直方图与正态密度")
    ax.set_xlabel("残差 / $10^{-3}$ kg/kg"); ax.set_ylabel("概率密度")

    ax = fig.add_subplot(gs[1, 2])
    for r, c, lb in ((rT, "#c0392b", "温度"), (rC / rC.std() * rT.std(),
                                               "#1f6fb2", "湿度(归一化)")):
        ac = [np.corrcoef(r[:-k], r[k:])[0, 1] for k in range(1, 26)]
        ax.plot(range(1, 26), ac, "o-", ms=3, lw=1.1, color=c, label=lb)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(1.96 / np.sqrt(len(rT)), color="gray", ls=":", lw=1)
    ax.axhline(-1.96 / np.sqrt(len(rT)), color="gray", ls=":", lw=1)
    ax.set_title("(f) 残差自相关函数")
    ax.set_xlabel("滞后阶数"); ax.set_ylabel("自相关系数")
    ax.legend(fontsize=8)

    return cm.savefig(fig, "fig02_room_residuals.png", plt)


# --------------------------------------------------------------------------
# fig03、fig04 —— 问题 1：预热平衡阶段（0~1800 s）的场与径向剖面
# --------------------------------------------------------------------------
def fig03_p1_fields():
    """fig03：问题 1 的温度场与水分浓度场，取自 data/p1_fields.npz 的 t_fine[s]、
    radii_cm[cm]、T_fine、C_fine。2x2 布局：上行 (a) 温度场、(b) 水分浓度场的时空云图
    （横轴时间/min，纵轴到中心距离/cm，色标分别为 degC 与 kg/kg）；下行 (c)(d) 是
    r=0、0.5、1、1.5、2 cm 五个位置的时间历程。场数组的行是时间、列是半径，故画云图时转置。
    """
    d = load("p1_fields.npz")
    t = d["t_fine"] / 60.0
    r = d["radii_cm"]
    fig, ax = plt.subplots(2, 2, figsize=(9.6 * FIG_SCALE, 6.4 * FIG_SCALE))
    for k, (dat, lab, cmap) in enumerate(
            [(d["T_fine"], "温度 / $^\\circ$C", "inferno"),
             (d["C_fine"], "水分浓度 / (kg/kg)", "viridis")]):
        mesh = ax[0, k].pcolormesh(t, r, dat.T, shading="auto", cmap=cmap)
        fig.colorbar(mesh, ax=ax[0, k], label=lab)
        ax[0, k].set_xlabel("时间 / min"); ax[0, k].set_ylabel("到中心距离 / cm")
        ax[0, k].set_title("(%s) %s 场" % ("ab"[k], lab.split(" / ")[0]))
    for k, (dat, lab) in enumerate([(d["T_fine"], "温度 / $^\\circ$C"),
                                    (d["C_fine"], "水分浓度 / (kg/kg)")]):
        for r_cm in (0.0, 0.5, 1.0, 1.5, 2.0):
            j = int(round(r_cm / 0.1))
            ax[1, k].plot(t, dat[:, j], lw=1.4,
                          label="$r=%.1f$ cm" % r_cm)
        ax[1, k].set_xlabel("时间 / min"); ax[1, k].set_ylabel(lab)
        ax[1, k].legend(fontsize=8, ncol=2)
        ax[1, k].set_title("(%s) 不同半径处的时间历程" % "cd"[k])
    fig.tight_layout()
    return cm.savefig(fig, "fig03_p1_fields.png", plt)


def fig04_p1_profiles():
    """fig04：问题 1 的径向剖面演化，取自 data/p1_fields.npz 的 t_fine、radii_cm、
    T_fine、C_fine。(a) 温度、(b) 水分浓度在 t=100、300、600、900、1200、1500、1800 s
    共 7 个时刻的剖面（横轴到中心距离/cm，颜色由浅到深表示时间先后）；(c) 只画 1800 s
    的水分剖面并标出表面值 C_s=1.5109 与尚未变化的中心值 C_0=2.55。
    """
    d = load("p1_fields.npz")
    r = d["radii_cm"]
    idx = [100, 300, 600, 900, 1200, 1500, 1800]
    times = d["t_fine"][np.array(idx) - 1]
    fig, ax = plt.subplots(1, 3, figsize=(11.4 * FIG_SCALE, 3.5 * FIG_SCALE))
    cols = plt.cm.plasma(np.linspace(0.05, 0.9, len(idx)))
    for c, i, tt in zip(cols, idx, times):
        ax[0].plot(r, d["T_fine"][i - 1], "-o", ms=3, lw=1.3, color=c,
                   label="%d s" % tt)
        ax[1].plot(r, d["C_fine"][i - 1], "-o", ms=3, lw=1.3, color=c,
                   label="%d s" % tt)
    ax[0].set_xlabel("到中心距离 / cm"); ax[0].set_ylabel("温度 / $^\\circ$C")
    ax[0].set_title("(a) 温度径向分布"); ax[0].legend(fontsize=7.5, ncol=2)
    ax[1].set_xlabel("到中心距离 / cm")
    ax[1].set_ylabel("水分浓度 / (kg/kg)")
    ax[1].set_title("(b) 水分浓度径向分布"); ax[1].legend(fontsize=7.5, ncol=2)

    # (c) 单独放大近表面那一层：D=7e-9 exp(-0.89/C) 在 C 变小时急剧变小，表面一旦
    # 降到 1.5 附近，紧贴表面的几个输出点就被拉开很陡的梯度，这是全章数值精度的关键处
    ax[2].plot(r, d["C_fine"][-1], "-o", ms=3, lw=1.4, color="#c0392b",
               label="模型 (表观界面)")
    ax[2].set_xlabel("到中心距离 / cm")
    ax[2].set_ylabel("水分浓度 / (kg/kg)")
    ax[2].set_title("(c) 1800 s 时近表面水分梯度")
    ax[2].annotate(r"$C_s=1.5109$", xy=(2.0, 1.5109), xytext=(1.0, 1.75),
                   arrowprops=dict(arrowstyle="->", lw=1),
                   fontsize=9, color="#c0392b")
    ax[2].annotate(r"$C_0=2.5500$ (未变化)", xy=(0.0, 2.55), xytext=(0.15, 2.30),
                   arrowprops=dict(arrowstyle="->", lw=1), fontsize=9)
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    return cm.savefig(fig, "fig04_p1_profiles.png", plt)


def fig05_p1_verification():
    """fig05：数值格式的解析验证。算例取常物性（附录 2）、恒定环境 T_air=60 degC、
    C_air=0.02 kg/kg 的圆柱导热问题，与 verify_solver.analytical_cylinder 的 Bessel 级数
    解对照。(a) 数值解与解析解之差在 (时间, 归一化半径 xi) 平面的分布（对数色标）；(b) 空间
    收敛：网格数 100~1600、dt 取 1 s 与 0.25 s；(c) 时间收敛：N=800、dt 4~0.25 s，黑虚线为二阶参考线。
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from verify_solver import analytical_cylinder

    class ConstRoom:
        def __init__(self, ta, ca):
            self.ta, self.ca = ta, ca

        def T_air(self, t):
            return np.full_like(np.asarray(t, float), self.ta)

        def C_air(self, t):
            return np.full_like(np.asarray(t, float), self.ca)

    alpha = 0.36 / (820.0 * 2600.0)
    bi = cm.H_CONV * cm.R0 / 0.36
    room = ConstRoom(60.0, 0.02)
    times = np.array([60.0, 300.0, 900.0, 1800.0, 3600.0])
    xi_out = np.linspace(0.0, 1.0, 41)
    res = cm.simulate(800, cm.props_problem1, 3600.0, 1.0, room,
                      t_report=times, xi_report=xi_out)
    exact = analytical_cylinder(xi_out * cm.R0, times, cm.T_INIT_C, 60.0,
                                alpha, bi)
    err = np.abs(res["T"] - exact)

    fig, ax = plt.subplots(1, 3, figsize=(11.4 * FIG_SCALE, 3.4 * FIG_SCALE))
    lo = max(err[err > 0].min(), 1e-8) if np.any(err > 0) else 1e-8
    m = ax[0].pcolormesh(times / 60.0, xi_out, err.T, shading="auto",
                         cmap="magma",
                         norm=plt.matplotlib.colors.LogNorm(vmin=lo,
                                                            vmax=err.max()))
    fig.colorbar(m, ax=ax[0], label="绝对误差 / K", extend="both")
    ax[0].set_xlabel("时间 / min"); ax[0].set_ylabel(r"归一化半径 $\xi=r/R$")
    ax[0].set_title(r"(a) 数值解与 Bessel 级数解之差 (Bi=%.2f)" % bi)

    ns = [100, 200, 400, 800, 1600]
    series = {}
    for dt_sp in (1.0, 0.25):
        e = []
        for n in ns:
            rr = cm.simulate(n, cm.props_problem1, 3600.0, dt_sp, room,
                             t_report=np.array([3600.0]), xi_report=xi_out)
            ee = analytical_cylinder(xi_out * cm.R0, np.array([3600.0]),
                                     cm.T_INIT_C, 60.0, alpha, bi)
            e.append(np.abs(rr["T"] - ee).max())
        series[dt_sp] = np.array(e)
    for dt_sp, c, mk in ((1.0, "#c0392b", "o"), (0.25, "#1f6fb2", "s")):
        ax[1].loglog(ns, series[dt_sp], mk + "-", color=c, lw=1.4,
                     label=r"$\Delta t=%.2f$ s" % dt_sp)
    ref = series[0.25][0] * (np.array(ns, float) / ns[0]) ** -2.0
    ax[1].loglog(ns, ref, "k--", lw=1.1, label=r"$\mathcal{O}(h^2)$ 参考线")
    ax[1].set_xlabel("径向网格数 $N$"); ax[1].set_ylabel("最大绝对误差 / K")
    ax[1].set_title("(b) 空间收敛性（$t=3600$ s）")
    ax[1].legend(fontsize=7.5)

    dts = [4.0, 2.0, 1.0, 0.5, 0.25]
    e_dt = []
    for dt in dts:
        rr = cm.simulate(800, cm.props_problem1, 3600.0, dt, room,
                         t_report=np.array([3600.0]), xi_report=xi_out)
        ee = analytical_cylinder(xi_out * cm.R0, np.array([3600.0]),
                                 cm.T_INIT_C, 60.0, alpha, bi)
        e_dt.append(np.abs(rr["T"] - ee).max())
    ax[2].loglog(dts, e_dt, "s-", color="#1f6fb2", lw=1.4)
    ax[2].loglog(dts, e_dt[0] * (np.array(dts) / dts[0]) ** 2.0, "k--", lw=1.2,
                 label=r"$\mathcal{O}(\Delta t^2)$ 参考线")
    ax[2].set_xlabel(r"时间步长 $\Delta t$ / s")
    ax[2].set_ylabel("最大绝对误差 / K")
    ax[2].set_title("(c) 时间收敛性（$N=800$）")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    return cm.savefig(fig, "fig05_verification.png", plt)


# --------------------------------------------------------------------------
# fig06 —— 附录给出的物性经验关系
# --------------------------------------------------------------------------
def fig06_properties():
    """fig06：附录 2 / 3 / 4 三条物性经验关系的对比，由 common.props_problem1、
    props_problem23、props_problem4 在 T=323.15 K（50 degC）下现算，不读结果文件。
    四个子图横轴都是水分浓度 C：(a) 密度 rho、(b) 比热容 cp、(c) 热传导系数 k、
    (d) 扩散系数 D（对数纵轴）。D 随 C 减小而迅速下降，是干燥速率被内部扩散控制的根源。
    """
    C = np.linspace(0.02, 2.55, 400)
    Tk = np.full_like(C, 323.15)
    r1, cp1, k1, D1 = cm.props_problem1(C)
    r3, cp3, k3, D3 = cm.props_problem23(C, Tk)
    r4, cp4, k4, D4 = cm.props_problem4(C, Tk)
    fig, ax = plt.subplots(1, 4, figsize=(13.2 * FIG_SCALE, 3.2 * FIG_SCALE))
    for a, y1, y3, y4, lab, ttl in [
            (ax[0], r1, r3, r4, r"$\rho$ / (kg/m$^3$)", "(a) 密度"),
            (ax[1], cp1, cp3, cp4, r"$c_p$ / (J kg$^{-1}$K$^{-1}$)", "(b) 比热容"),
            (ax[2], k1, k3, k4, r"$k$ / (W m$^{-1}$K$^{-1}$)", "(c) 热传导系数")]:
        a.plot(C, y1, lw=1.6, label="附录2 (问题1)")
        a.plot(C, y3, lw=1.6, label="附录3 (问题2,3)")
        a.plot(C, y4, lw=1.6, label="附录4 (问题4)")
        a.set_xlabel(r"水分浓度 $C$ / (kg/kg)"); a.set_ylabel(lab)
        a.set_title(ttl); a.legend(fontsize=7.5)
    ax[3].semilogy(C, D1, lw=1.6, label="附录2")
    ax[3].semilogy(C, D3, lw=1.6, label="附录3")
    ax[3].semilogy(C, D4, lw=1.6, label="附录4")
    ax[3].set_xlabel(r"水分浓度 $C$ / (kg/kg)")
    ax[3].set_ylabel(r"$D$ / (m$^2$/s)")
    ax[3].set_title(r"(d) 扩散系数 $D(C,T=50\,^\circ$C)")
    ax[3].legend(fontsize=8)
    fig.tight_layout()
    return cm.savefig(fig, "fig06_properties.png", plt)


# --------------------------------------------------------------------------
# fig07~fig10 —— 问题 2、3：整段烘干过程的场、干燥曲线、通量与终点附近的细节
# --------------------------------------------------------------------------
def fig07_p23_fields():
    """fig07：问题 2 整段 0~72 h 烘干的温度场与水分浓度场，取自 data/p23_fields.npz 的
    t[s]、radii_cm[cm]、T、C、t_dry。2x2 布局：上行 (a) 温度场、(b) 水分浓度场的时空云图
    （横轴时间/h）；下行 (c)(d) 是 r=0、0.5、1、1.5、2 cm 处的时间历程，(d) 中标出判据
    C=0.15（红虚线）与烘干时刻 t_dry（竖点线），是展示问题 2、3 全过程的主图。
    """
    d = load("p23_fields.npz")
    t = d["t"] / 3600.0
    r = d["radii_cm"]
    fig, ax = plt.subplots(2, 2, figsize=(10.2 * FIG_SCALE, 6.6 * FIG_SCALE))
    for k, (dat, lab, cmap) in enumerate(
            [(d["T"], "温度 / $^\\circ$C", "inferno"),
             (d["C"], "水分浓度 / (kg/kg)", "viridis")]):
        mesh = ax[0, k].pcolormesh(t, r, dat.T, shading="auto", cmap=cmap)
        fig.colorbar(mesh, ax=ax[0, k], label=lab)
        ax[0, k].set_xlabel("时间 / h"); ax[0, k].set_ylabel("到中心距离 / cm")
        ax[0, k].set_title("(%s) %s 场（整段烘干过程）" %
                           ("ab"[k], lab.split(" / ")[0]))
    for k, (dat, lab) in enumerate([(d["T"], "温度 / $^\\circ$C"),
                                    (d["C"], "水分浓度 / (kg/kg)")]):
        for r_cm in (0.0, 0.5, 1.0, 1.5, 2.0):
            j = int(round(r_cm / 0.1))
            ax[1, k].plot(t, dat[:, j], lw=1.4, label="$r=%.1f$ cm" % r_cm)
        ax[1, k].set_xlabel("时间 / h"); ax[1, k].set_ylabel(lab)
        ax[1, k].legend(fontsize=8, ncol=2)
        ax[1, k].set_title("(%s) 不同半径处的时间历程" % "cd"[k])
    ax[1, 1].axhline(C_TARGET, color="r", ls="--", lw=1.2)
    ax[1, 1].annotate(r"判据 $C=0.15$", xy=(40, 0.17), color="r", fontsize=9)
    tdry = float(d["t_dry"][0]) / 3600.0
    ax[1, 1].axvline(tdry, color="k", ls=":", lw=1.2)
    ax[1, 1].annotate("烘干结束 %.2f h" % tdry, xy=(tdry * 0.55, 1.6),
                      fontsize=9)
    fig.tight_layout()
    return cm.savefig(fig, "fig07_p23_fields.png", plt)


def fig08_p23_curves():
    """fig08：问题 2、3 的干燥曲线与径向剖面，取自 data/p23_fields.npz 的 t、radii_cm、
    C、t_dry。(a) 半对数干燥曲线：5 个半径处的 C(t)，红虚线与竖点线为判据 0.15 和 t_dry；
    (b) 6、12、18、24、36、48、57 h 共 7 条剖面 C(r)；(c) 只取 6、24、48 h，显示剖面
    形状由平缓转为陡峭，即水分蒸发的“前沿”逐步向中心推进。
    """
    d = load("p23_fields.npz")
    t = d["t"] / 3600.0
    tdry = float(d["t_dry"][0]) / 3600.0
    fig, ax = plt.subplots(1, 3, figsize=(12.2 * FIG_SCALE, 3.6 * FIG_SCALE))
    for r_cm in (0.0, 0.5, 1.0, 1.5, 2.0):
        j = int(round(r_cm / 0.1))
        ax[0].semilogy(t, d["C"][:, j], lw=1.4, label="$r=%.1f$ cm" % r_cm)
    ax[0].axhline(C_TARGET, color="r", ls="--", lw=1.2)
    ax[0].axvline(tdry, color="k", ls=":", lw=1.2)
    ax[0].set_xlabel("时间 / h")
    ax[0].set_ylabel("水分浓度 $C$ / (kg/kg)")
    ax[0].set_title("(a) 半对数干燥曲线")
    ax[0].legend(fontsize=8)
    ax[0].annotate("$t_{dry}=%.2f$ h" % tdry, xy=(tdry + 1, 0.9), fontsize=9)

    for h in (6, 12, 18, 24, 36, 48, 57):
        j = int(np.argmin(np.abs(t - h)))
        ax[1].plot(d["radii_cm"], d["C"][j], "-o", ms=3, lw=1.3,
                   label="%d h" % h)
    ax[1].axhline(C_TARGET, color="r", ls="--", lw=1.1)
    ax[1].set_xlabel("到中心距离 / cm")
    ax[1].set_ylabel("水分浓度 $C$ / (kg/kg)")
    ax[1].set_title("(b) 径向剖面演化")
    ax[1].legend(fontsize=7.5, ncol=2)

    for h, c in ((6, "#c0392b"), (24, "#e67e22"), (48, "#1f6fb2")):
        j = int(np.argmin(np.abs(t - h)))
        ax[2].plot(d["radii_cm"], d["C"][j], "-o", ms=3, lw=1.4, color=c,
                   label="%d h" % h)
    ax[2].set_xlabel("到中心距离 / cm")
    ax[2].set_ylabel("水分浓度 $C$ / (kg/kg)")
    ax[2].set_title("(c) 剖面形状：由平缓转为陡峭")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    return cm.savefig(fig, "fig08_p23_curves.png", plt)


def fig09_p23_flux():
    """fig09：问题 2、3 的传质与传热过程，取自 data/p23_fields.npz 的 t、T、C。
    (a) 截面平均水分 Cbar 与中心水分 C(0,t) 的对比：Cbar 按 2*xi*dxi 加权积分（xi=radii_cm/2），
    它先降到 0.15，故用截面平均作判据会低估烘干时间；(b) 表面传质通量 h_m(C_s-C_air)
    （对数纵轴，反映先恒速后降速的过程）；(c) T_air、表面温度 T_s 与中心温度 T(0,t)。
    """
    d = load("p23_fields.npz")
    t = d["t"]
    th = t / 3600.0
    xi = d["radii_cm"] / 2.0
    # 截面平均的梯形权重：Cbar = 2*int_0^1 C xi dxi，xi 间距为 0.05，
    # 轴心与表面各只有半个控制体厚，故两端取半权
    w = 2.0 * xi * 0.05
    w[0] *= 0.5
    w[-1] *= 0.5
    Cbar = d["C"] @ w
    Ts = d["T"][:, -1]
    Tair = cm.RoomConditions().T_air(t)
    Csurf = d["C"][:, -1]
    Cair = cm.RoomConditions().C_air(t)
    flux = cm.HM_CONV * (Csurf - Cair)              # model-unit surface flux

    fig, ax = plt.subplots(1, 3, figsize=(12.2 * FIG_SCALE, 3.5 * FIG_SCALE))
    ax[0].plot(th, Cbar, lw=1.6, color="#1f6fb2", label=r"截面平均 $\bar C$")
    ax[0].plot(th, d["C"][:, 0], lw=1.6, color="#c0392b", label="$C(r=0)$")
    ax[0].axhline(C_TARGET, color="k", ls="--", lw=1.1)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("水分浓度 / (kg/kg)")
    ax[0].set_title(r"(a) 平均水分与中心水分（判据 $\bar C$ 提前满足）")
    ax[0].legend(fontsize=8)
    ax[0].axvline(35.17, color="gray", ls=":", lw=1.1)
    ax[0].annotate(r"$\bar C=0.15$ @ 35.2 h", xy=(20, 0.35), fontsize=8.5,
                   color="gray")

    ax[1].plot(th, flux * 1e9, lw=1.6, color="#16a085")
    ax[1].set_xlabel("时间 / h")
    ax[1].set_ylabel(r"表面水分通量 $h_m(C_s-C_{air})$ / ($10^{-9}$ kg kg$^{-1}$m s$^{-1}$)")
    ax[1].set_title("(b) 表面传质通量：先快后慢的降速期")
    ax[1].set_yscale("log")

    ax[2].plot(th, Tair, lw=1.6, color="#7f8c8d", label=r"$T_{air}$")
    ax[2].plot(th, Ts, lw=1.6, color="#c0392b", label=r"$T_s$ (表面)")
    ax[2].plot(th, d["T"][:, 0], lw=1.6, color="#2c3e50", label=r"$T(r=0)$")
    ax[2].set_xlabel("时间 / h"); ax[2].set_ylabel("温度 / $^\\circ$C")
    ax[2].set_title("(c) 传热：药温约 1.5 h 内与环境平衡")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    return cm.savefig(fig, "fig09_p23_flux.png", plt)


def fig10_p23_zoom():
    """fig10：烘干终点附近的细节与降速期特征，取自 data/p23_fields.npz 的 t、C、t_dry。
    (a) 放大 t_dry-12 h ~ t_dry+6 h 窗口的中心水分 C_0(t)，标出 t_dry 处恰好为 0.1500，
    说明判据是由中心点最后满足的；(b) 中心干燥速率 -dC_0/dt（由 np.gradient 数值微分）
    的双纵轴图：左轴对时间、右轴对 C_0，展示降速期速率随 C 近似指数衰减。
    """
    d = load("p23_fields.npz")
    t = d["t"] / 3600.0
    tdry = float(d["t_dry"][0]) / 3600.0
    cc = d["C"][:, 0]
    m = (t > tdry - 12) & (t < tdry + 6)
    fig, ax = plt.subplots(1, 2, figsize=(10.0 * FIG_SCALE, 3.6 * FIG_SCALE))
    ax[0].plot(t[m], cc[m], "-o", ms=2.5, lw=1.5, color="#c0392b")
    ax[0].axhline(C_TARGET, color="k", ls="--", lw=1.2)
    ax[0].axvline(tdry, color="k", ls=":", lw=1.2)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("中心水分浓度 / (kg/kg)")
    ax[0].set_title("(a) 判据附近的放大：$t_{dry}=%.3f$ h" % tdry)
    ax[0].annotate(r"$C_0(t_{dry})=0.1500$", xy=(tdry - 11.5, 0.1504),
                   fontsize=9)

    dt = np.gradient(t) * 3600.0
    rate = np.gradient(cc) / dt
    ax[1].semilogy(t[1:], -rate[1:], lw=1.6, color="#1f6fb2")
    ax[1].set_xlabel("时间 / h")
    ax[1].set_ylabel(r"中心干燥速率 $-\mathrm{d}C_0/\mathrm{d}t$ / s$^{-1}$")
    ax[1].set_title(r"(b) 降速期：速率随 $C$ 下降呈指数衰减")
    ax2 = ax[1].twinx()
    ax2.semilogy(cc[1:], -rate[1:], lw=1.0, ls="--", color="#e67e22")
    ax2.set_ylabel(r"对 $C_0$ 的依赖", color="#e67e22")
    ax2.tick_params(axis="y", colors="#e67e22")
    ax2.grid(False)
    fig.tight_layout()
    return cm.savefig(fig, "fig10_p23_zoom.png", plt)


# --------------------------------------------------------------------------
# fig11~fig14 —— 问题 4：收缩半径、收缩工况下的场，以及与问题 3 的对比和外部检验
# --------------------------------------------------------------------------
def fig11_p4_radius():
    """fig11：问题 4 的半径收缩数据，来自附件 2（common.load_attachment2）与
    common.RadiusHistory。(a) 实测半径、单指数拟合 R_inf+(R_0-R_inf)e^{-t/tau} 与本文实际
    采用的分段线性插值三条曲线（横轴时间/h，纵轴半径/cm）；(b) 收缩速率 dR/dt（μm/s，
    对实测表差分），90% 的收缩在前 10 h 内完成；(c) 三个无量纲比值：半径比 R/R_0、
    横截面积比 (R/R_0)^2（本文无限长圆柱口径）与各向同性体积比 (R/R_0)^3。
    """
    df = cm.load_attachment2()
    t = df["t"].to_numpy(float) / 3600.0
    R = df["R_cm"].to_numpy(float)
    rad = cm.RadiusHistory(mode="table")
    rs = cm.RadiusHistory(mode="smooth")
    tt = np.linspace(0, 72, 600)
    fig, ax = plt.subplots(1, 3, figsize=(12.2 * FIG_SCALE, 3.5 * FIG_SCALE))
    ax[0].plot(t, R, "o", ms=3, color="#c0392b", label="附件2 实测")
    ax[0].plot(tt, rs.R(tt * 3600.0) * 100, "--", lw=1.4, color="#2c3e50",
               label=r"指数拟合 $R_\infty+(R_0-R_\infty)e^{-t/\tau_R}$")
    ax[0].plot(tt, rad.R(tt * 3600.0) * 100, "-", lw=1.2, color="#1f6fb2",
               label="分段线性插值（本文采用）")
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("半径 $R$ / cm")
    ax[0].set_title("(a) 药材半径的收缩过程")
    ax[0].legend(fontsize=8)

    dR = np.gradient(R * 1e-2, df["t"].to_numpy(float))
    ax[1].plot(t, dR * 1e6, "-", lw=1.5, color="#16a085")
    ax[1].set_xlabel("时间 / h")
    ax[1].set_ylabel(r"$\mathrm{d}R/\mathrm{d}t$ / ($\mu$m/s)")
    ax[1].set_title("(b) 收缩速率：前 10 h 完成 90%")

    ax[2].plot(t, (R / R[0]) ** 3, "--", lw=1.5, color="#8e44ad",
               label=r"各向同性体积比 $(R/R_0)^3$")
    ax[2].plot(t, (R / R[0]) ** 2, "-", lw=1.6, color="#1f6fb2",
               label=r"横截面积比 $(R/R_0)^2$（本文口径）")
    ax[2].plot(t, R / R[0], "-", lw=1.6, color="#e67e22",
               label=r"半径比 $R/R_0$")
    ax[2].set_xlabel("时间 / h"); ax[2].set_ylabel("无量纲比值")
    ax[2].set_title("(c) 半径减至 59.9%，横截面积减至 35.9%")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    return cm.savefig(fig, "fig11_p4_radius.png", plt)


def fig12_p4_fields():
    """fig12：问题 4 收缩工况下的两个场，取自 data/p4_fields.npz：t[s]、xi_mat（物质坐标
    xi=r/R(t)，0~1 等距）、radii_cm（物理半径 0~2 cm）、C/T（物质坐标场）、C_phys（物理坐标场）、
    R[m]、t_dry。2x2 布局：上行 (a) 水分、(b) 温度的物质坐标云图 C(xi,t)、T(xi,t)；下行
    (c) 物理坐标下的 C 场（已收缩出界的点不画，故域随时间变窄）、(d) 物理坐标下的径向剖面。
    """
    d = load("p4_fields.npz")
    t = d["t"] / 3600.0
    xi = d["xi_mat"]                       # 物质坐标轴 xi=r/R(t)，剖面按它排列
    Rcm = d["R"] * 100.0
    tdry = float(d["t_dry"][0]) / 3600.0
    fig, ax = plt.subplots(2, 2, figsize=(10.2 * FIG_SCALE, 6.6 * FIG_SCALE))
    mesh = ax[0, 0].pcolormesh(t, xi, d["C"].T, shading="auto", cmap="viridis")
    fig.colorbar(mesh, ax=ax[0, 0], label="水分浓度 / (kg/kg)")
    ax[0, 0].set_xlabel("时间 / h")
    ax[0, 0].set_ylabel(r"物质坐标 $\xi=r/R(t)$")
    ax[0, 0].set_title(r"(a) 物质坐标下的水分浓度场 $C(\xi,t)$")
    ax[0, 0].plot(tdry, 0, "r*", ms=10)
    ax[0, 0].plot(tdry, 1, "r*", ms=10)

    mesh = ax[0, 1].pcolormesh(t, xi, d["T"].T, shading="auto", cmap="inferno")
    fig.colorbar(mesh, ax=ax[0, 1], label="温度 / $^\\circ$C")
    ax[0, 1].set_xlabel("时间 / h")
    ax[0, 1].set_ylabel(r"物质坐标 $\xi=r/R(t)$")
    ax[0, 1].set_title(r"(b) 物质坐标下的温度场 $T(\xi,t)$")

    # 坐标换算：求解器工作在物质坐标 xi 上，这里要把每个时刻的剖面画到物理半径上。
    # C_phys 的第 j 列已经是按 xi = r_j/R(t) 采样得到的，所以列标签 r_j 就是该时刻的
    # 真实距离；反过来 (d) 子图里的 C 是物质坐标下的场，需要把 xi 乘上该时刻的 R(t)
    # 才变回物理半径。r_j > R(t) 的点已经随收缩离开药材，用 nan 掩掉即为收缩域。
    r_axis = d["radii_cm"]
    Cphys = np.where(r_axis[None, :] <= Rcm[:, None], d["C_phys"], np.nan)
    # 色标固定为 0~2.55（初始含水率），与 (a) 子图同一含义，便于比较两个坐标下的场
    mesh = ax[1, 0].pcolormesh(t, r_axis, Cphys.T, shading="auto",
                               cmap="viridis", vmin=0, vmax=2.55)
    fig.colorbar(mesh, ax=ax[1, 0], label="水分浓度 / (kg/kg)")
    ax[1, 0].plot(t, Rcm, "r-", lw=2, label="药材表面 $R(t)$")
    ax[1, 0].set_xlabel("时间 / h"); ax[1, 0].set_ylabel("到中心距离 / cm")
    ax[1, 0].set_title("(c) 物理坐标下的水分浓度场（域随时间收缩）")
    ax[1, 0].legend(fontsize=8)

    for h, c in ((6, "#c0392b"), (12, "#e67e22"), (24, "#f1c40f"),
                 (48, "#1f6fb2")):
        j = int(np.argmin(np.abs(t - h)))
        rr = d["radii_cm"] * (d["R"][j] * 100.0 / 2.0)
        ax[1, 1].plot(rr, d["C"][j], "-o", ms=3, lw=1.5, color=c,
                      label="%d h (R=%.2f cm)" % (h, d["R"][j] * 100))
    ax[1, 1].axhline(C_TARGET, color="r", ls="--", lw=1.1)
    ax[1, 1].set_xlabel("到中心距离 / cm")
    ax[1, 1].set_ylabel("水分浓度 / (kg/kg)")
    ax[1, 1].set_title("(d) 物理坐标下的径向剖面")
    ax[1, 1].legend(fontsize=7.5)
    fig.tight_layout()
    return cm.savefig(fig, "fig12_p4_fields.png", plt)


def fig13_p3_vs_p4():
    """fig13：收缩对干燥过程的影响，对比问题 3（不收缩，data/p23_fields.npz）与问题 4
    （收缩，data/p4_fields.npz）。(a) 中心水分 C(0,t) 的半对数曲线，两条竖点线是各自的
    t_dry；(b) 表面水分与内部一点的水分：问题 3 取 r=1.0 cm、问题 4 取 xi=0.5，在各自
    的坐标系里都是第 11 列但物理含义不同，红虚线为判据 0.15。
    """
    d3 = load("p23_fields.npz")
    d4 = load("p4_fields.npz")
    t3 = d3["t"] / 3600.0
    t4 = d4["t"] / 3600.0
    td3 = float(d3["t_dry"][0]) / 3600.0
    td4 = float(d4["t_dry"][0]) / 3600.0
    fig, ax = plt.subplots(1, 2, figsize=(10.2 * FIG_SCALE, 3.8 * FIG_SCALE))
    ax[0].semilogy(t3, d3["C"][:, 0], lw=1.7, color="#c0392b",
                   label="问题3（不收缩）$t_{dry}=%.2f$ h" % td3)
    ax[0].semilogy(t4, d4["C"][:, 0], lw=1.7, color="#1f6fb2",
                   label="问题4（收缩）$t_{dry}=%.2f$ h" % td4)
    ax[0].axhline(C_TARGET, color="k", ls="--", lw=1.2)
    ax[0].axvline(td3, color="#c0392b", ls=":", lw=1.1)
    ax[0].axvline(td4, color="#1f6fb2", ls=":", lw=1.1)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("中心水分浓度 / (kg/kg)")
    ax[0].set_title("(a) 中心水分浓度对比")
    ax[0].legend(fontsize=8.5)

    ax[1].plot(t3, d3["C"][:, -1], lw=1.7, color="#c0392b", label="问题3 表面")
    ax[1].plot(t4, d4["C"][:, -1], lw=1.7, color="#1f6fb2", label="问题4 表面")
    ax[1].plot(t3, d3["C"][:, 10], lw=1.3, ls="--", color="#c0392b",
               label="问题3 $r=1.0$ cm")
    ax[1].plot(t4, d4["C"][:, 10], lw=1.3, ls="--", color="#1f6fb2",
               label="问题4 $\\xi=0.5$")
    ax[1].axhline(C_TARGET, color="k", ls="--", lw=1.1)
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("水分浓度 / (kg/kg)")
    ax[1].set_title("(b) 表面与内部点对比")
    ax[1].legend(fontsize=7.5)
    fig.tight_layout()
    return cm.savefig(fig, "fig13_p3_vs_p4.png", plt)


def fig14_shrink_consistency():
    """fig14：外部数据交叉检验——附件 2 的实测收缩与附录 4 的密度式相容吗？取
    data/p4_fields.npz 的 t、C（物质坐标下的剖面）、R（实测半径[m]）、radii_cm。用
    rho=760+90C 把各时刻的干基质量沿截面 2*int rho_d xi dxi 积分，再由干物质守恒反推
    半径 R_pred。(a) R_pred 与实测 R(t)；(b) 两者的相对偏差；(c) 两种体积收缩比。
    """
    d = load("p4_fields.npz")
    t = d["t"] / 3600.0
    C = d["C"]
    xi = d["radii_cm"] / 2.0
    # 梯形权重：w 积的是 2*int_0^1 f xi dxi，xi 间距为 0.05；轴心与表面处只有半个控制体厚，
    # 故两端的权取一半。基函数是 xi，所以 x=0 处的权本来就是 0
    w = 2.0 * xi * 0.05
    w[0] *= 0.5
    w[-1] *= 0.5
    rho = 760.0 + 90.0 * C
    rho_d = rho / (1.0 + C)                       # 单位湿体积内的干物质质量
    # 截面平均的扇形积分：把剖面按半径加权求和即得单位长度上的干物质质量 /(2 pi)
    integral = 0.5 * (rho_d @ w)
    # 同一常数由初始均匀状态算出：R_0^2 * rho_d(C_INIT)/2；半径不收缩时两者应相等
    M0 = cm.R0 ** 2 * ((760.0 + 90.0 * cm.C_INIT) / (1.0 + cm.C_INIT)) / 2.0
    R_pred = np.sqrt(M0 / integral) * 100.0
    R_meas = d["R"] * 100.0

    fig, ax = plt.subplots(1, 3, figsize=(12.2 * FIG_SCALE, 3.5 * FIG_SCALE))
    ax[0].plot(t, R_meas, lw=1.8, color="#c0392b", label="附件2 实测 $R(t)$")
    ax[0].plot(t, R_pred, lw=1.8, ls="--", color="#1f6fb2",
               label=r"由 $\rho(C)$ 与干基质量守恒反推")
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("半径 / cm")
    ax[0].set_title("(a) 收缩数据与物性公式的一致性检验")
    ax[0].legend(fontsize=8.5)

    ax[1].plot(t, 100.0 * (R_pred - R_meas) / R_meas, lw=1.6, color="#8e44ad")
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("相对偏差 / %")
    ax[1].set_title("(b) 相对偏差：物性式与实测收缩不自洽")
    ax[1].annotate("最大偏差 %.1f%%" % np.max(np.abs(100 * (R_pred - R_meas) / R_meas)),
                   xy=(20, np.max(100 * (R_pred - R_meas) / R_meas) * 0.7),
                   fontsize=9)

    ax[2].plot(t, (R_meas / R_meas[0]) ** 3, lw=1.7, color="#c0392b",
               label="实测体积比")
    ax[2].plot(t, (R_pred / R_pred[0]) ** 3, lw=1.7, ls="--", color="#1f6fb2",
               label="由 $\\rho(C)$ 反推的体积比")
    ax[2].set_xlabel("时间 / h"); ax[2].set_ylabel(r"$V/V_0$")
    ax[2].set_title("(c) 体积收缩比对比")
    ax[2].legend(fontsize=8.5)
    fig.tight_layout()
    return cm.savefig(fig, "fig14_shrink_consistency.png", plt)


# --------------------------------------------------------------------------
# fig15、fig16 —— 独立方法的互证，以及全局水分守恒检验
# --------------------------------------------------------------------------
def fig15_crosscheck():
    """fig15：用三条彼此独立的求解路径互证问题 2、3 的结果（由 common.simulate 与
    crosscheck.mol_solve 现算，不读 npz）。(a) 中心水分 C(0,t) 的干燥曲线：有限体积节点
    格式、有限体积单元中心格式、直线法（MOL + BDF）；(b) 后两者与节点格式之差（对数纵轴，
    全时段 <1e-4 kg/kg）；(c) 两种格式各自的网格收敛，纵轴为 t=2e5 s 时的中心值之差。
    """
    import crosscheck as cc
    room = cm.RoomConditions()
    rad = cm.RadiusHistory(mode="table")
    t_ev = np.arange(0.0, 2.0e5 + 1.0, 1200.0)
    xi = np.array([0.0, 0.25, 0.5, 0.75, 0.9, 1.0])

    res = {}
    for tag, prop, radius, grid, n in [("FV (node)", cm.props_problem23, None, "node", 800),
                                       ("FV (cell)", cm.props_problem23, None, "cell", 800),
                                       ("MOL (BDF)", None, None, None, 300)]:
        if tag.startswith("MOL"):
            # 直线法这里用 300 个单元、共 301 个节点，取出 xi=0、0.25、0.5、0.75、
            # 0.9、1 对应的 6 列，与上面 xi 数组给出的采样位置一一对齐
            r = cc.mol_solve(prop=cm.props_problem23, t_end=2.0e5, n_node=n,
                             room=room, t_eval=t_ev)
            C = r["C"][:, [0, 75, 150, 225, 270, 300]]
        else:
            r = cm.simulate(n, prop, 2.0e5, 1.0, room, grid=grid,
                            dt_schedule=[(7200.0, 1.0), (2.0e5, 10.0)],
                            t_report=t_ev, xi_report=xi)
            C = r["C"]
        res[tag] = C

    fig, ax = plt.subplots(1, 3, figsize=(12.2 * FIG_SCALE, 3.5 * FIG_SCALE))
    for tag, c in res.items():
        ax[0].semilogy(t_ev / 3600.0, c[:, 0], lw=1.6, label=tag)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("中心水分浓度 / (kg/kg)")
    ax[0].set_title("(a) 三种独立计算的干燥曲线")
    ax[0].legend(fontsize=8.5)

    ref = res["FV (node)"]
    for tag in ("FV (cell)", "MOL (BDF)"):
        ax[1].plot(t_ev / 3600.0, np.abs(res[tag][:, 0] - ref[:, 0]), lw=1.5,
                   label=tag + " vs FV(node)")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("时间 / h")
    ax[1].set_ylabel(r"$|C_0|$ 之差 / (kg/kg)")
    ax[1].set_title(r"(b) 中心点差异（全时段 $<10^{-4}$ kg/kg）")
    ax[1].legend(fontsize=8)

    ns = np.array([100, 200, 400, 800, 1600])
    e1 = []
    for n in ns:
        r = cm.simulate(n, cm.props_problem23, 2.0e5, 1.0, room, grid="node",
                        dt_schedule=[(7200.0, 1.0), (2.0e5, 10.0)],
                        t_report=np.array([2.0e5]), xi_report=np.array([0.0]))
        e1.append(r["C"][0, 0])
    ns2 = np.array([75, 150, 300, 600, 1200])
    e2 = []
    for n in ns2:
        r = cc.mol_solve(cm.props_problem23, 2.0e5, n_node=n, room=room,
                         t_eval=np.array([2.0e5]))
        e2.append(r["C"][0, 0])
    ax[2].semilogx(ns, (np.array(e1) - e1[-1]), "o-", lw=1.5,
                   label="FV (node)")
    ax[2].semilogx(ns2, (np.array(e2) - e2[-1]), "s-", lw=1.5, label="MOL (BDF)")
    ax[2].set_xlabel("网格数")
    ax[2].set_ylabel(r"$C_0(2.0\times10^5 \mathrm{s})$ 相对最细网格之差")
    ax[2].set_title("(c) 两种格式的网格收敛")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    return cm.savefig(fig, "fig15_crosscheck.png", plt)


def fig16_water_balance():
    """fig16：问题 2、3 的全局水分守恒检验，取自 data/p23_fields.npz 的 t、C，环境浓度由
    common.RoomConditions 给出。(a) 截面含水量 W=pi*R_0^2*Cbar 的衰减（Cbar 为按 2*xi*dxi
    加权的截面平均，故 W 只差一个单位长度因子，为模型单位）；(b) 两条累计失水曲线：表面通量
    的时间积分（乘 2*pi*R_0）与内部储量变化 W(0)-W(t)，两者重合即说明离散格式整体守恒。
    """
    d = load("p23_fields.npz")
    t = d["t"]
    xi = d["radii_cm"] / 2.0
    # 与 fig09 相同的截面平均权重：w 积 2*int f xi dxi，两端各取半权
    w = 2.0 * xi * 0.05
    w[0] *= 0.5
    w[-1] *= 0.5
    Cbar = d["C"] @ w
    room = cm.RoomConditions()
    Csurf = d["C"][:, -1]
    flux = cm.HM_CONV * (Csurf - room.C_air(t))
    W = np.pi * cm.R0 ** 2 * Cbar                       # model-unit water
    loss = W[0] - W
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (flux[1:] + flux[:-1])
                                           * np.diff(t))]) * 2 * np.pi * cm.R0
    fig, ax = plt.subplots(1, 2, figsize=(9.6 * FIG_SCALE, 3.5 * FIG_SCALE))
    ax[0].plot(t / 3600.0, W * 1e3, lw=1.7, color="#1f6fb2")
    ax[0].set_xlabel("时间 / h")
    ax[0].set_ylabel(r"截面含水（模型单位）/ $10^{-3}$")
    ax[0].set_title("(a) 药材内部水量的衰减")
    ax[1].plot(t / 3600.0, cum * 1e3, lw=1.7, color="#c0392b",
               label="表面通量累积（式 (24)）")
    ax[1].plot(t / 3600.0, loss * 1e3, "--", lw=1.5, color="#2c3e50",
               label="内部储量变化")
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel(r"累计失水 / $10^{-3}$")
    ax[1].set_title("(b) 全局水分守恒检验")
    ax[1].legend(fontsize=8.5)
    fig.tight_layout()
    return cm.savefig(fig, "fig16_water_balance.png", plt)


# --------------------------------------------------------------------------
# fig17~fig19 —— 模型形式变体与烘房环境条件的灵敏度
# --------------------------------------------------------------------------
def fig17_variants():
    """fig17：模型形式变体的影响，各算例都由 common.simulate 现算（400 个控制体）。
    (a) 中心温度 T(0,t)：基准（问题 3）与计入蒸发潜热冷却的变体，后者表面散热超过对流供热、
    出现不物理的过冷，且该变体在网格加密下发散，故只能说明"潜热无法由题面数据定量"，
    不能当作敏感性上界；(b) 两者的中心水分：潜热变体使干燥变慢约 9.4%（数值不稳定，
    不作为结论）；(c) 问题 4 的基准收缩、强制无收缩、能量方程含体积项三个变体，说明收缩
    才是主导效应（不收缩时约 129 h）。"不收缩"算例要跑到 139 h 才会达到判据，故单独用
    更长的时限，避免曲线被截断。
    """
    room = cm.RoomConditions()
    rad = cm.RadiusHistory(mode="table")
    T_END = 2.6e5
    T_END_LONG = 5.0e5
    sched = [(7200.0, 1.0), (T_END, 10.0)]
    sched_long = [(7200.0, 1.0), (T_END_LONG, 10.0)]
    t_rep = np.arange(1.0, T_END + 0.5, 120.0)
    t_rep_long = np.arange(1.0, T_END_LONG + 0.5, 120.0)
    xi = np.array([0.0, 1.0])              # 只输出中心与表面两点，够画对比曲线
    # 各算例：前两个用问题 3 物性（基准 / 计入蒸发潜热冷却），后三个用问题 4 物性
    # （基准收缩 / 强制不收缩 / 能量方程额外计入体积变化项）
    runs = {
        "基准（问题3）": dict(prop=cm.props_problem23),
        "含蒸发潜热冷却": dict(prop=cm.props_problem23, evap_cooling=True),
        "基准（问题4）": dict(prop=cm.props_problem4, radius=rad),
        "问题4 无收缩": dict(prop=cm.props_problem4, long=True),
        "问题4 能量方程体积项": dict(prop=cm.props_problem4, radius=rad,
                                     volume_term=True),
    }
    out = {}
    for k, kw in runs.items():
        long_run = kw.get("long", False)
        r = cm.simulate(400, kw["prop"],
                        T_END_LONG if long_run else T_END, 1.0, room,
                        radius=kw.get("radius"), grid="node",
                        volume_term=kw.get("volume_term", False),
                        evap_cooling=kw.get("evap_cooling", False),
                        dt_schedule=sched_long if long_run else sched,
                        t_report=t_rep_long if long_run else t_rep, xi_report=xi)
        out[k] = r

    fig, ax = plt.subplots(1, 3, figsize=(12.2 * FIG_SCALE, 3.6 * FIG_SCALE))
    for k in ("基准（问题3）", "含蒸发潜热冷却"):
        ax[0].plot(out[k]["t"] / 3600.0, out[k]["T"][:, 0], lw=1.6, label=k)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("中心温度 / $^\\circ$C")
    ax[0].set_title("(a) 潜热变体出现不物理的表面过冷")
    ax[0].legend(fontsize=8, loc="lower right")
    ax[0].annotate("题给 $(h,h_m)$ 与 Lewis 关系不自洽：\n"
                   "计入潜热后表面散热远超对流供热，\n"
                   "中心最低降至 $8.9\\,^\\circ$C（不物理）；\n"
                   "该变体在网格加密下发散，不作定量结论",
                   xy=(0.7, 11.0), xytext=(8.5, 17.5), fontsize=7.4,
                   color="#b9770e",
                   arrowprops=dict(arrowstyle="->", color="#b9770e", lw=1.0))

    for k in ("基准（问题3）", "含蒸发潜热冷却"):
        ax[1].semilogy(out[k]["t"] / 3600.0, out[k]["C"][:, 0], lw=1.6, label=k)
    ax[1].axhline(C_TARGET, color="k", ls="--", lw=1.1)
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("中心水分浓度 / (kg/kg)")
    ax[1].set_title("(b) 潜热变体使干燥变慢 $+9.4\\%$（数值不稳定，不作上界）")
    ax[1].legend(fontsize=8.5)

    for k in ("基准（问题4）", "问题4 无收缩", "问题4 能量方程体积项"):
        ax[2].semilogy(out[k]["t"] / 3600.0, out[k]["C"][:, 0], lw=1.6, label=k)
    ax[2].axhline(C_TARGET, color="k", ls="--", lw=1.1)
    ax[2].set_xlabel("时间 / h"); ax[2].set_ylabel("中心水分浓度 / (kg/kg)")
    ax[2].set_title("(c) 收缩项是主导效应（不收缩 $\\to 129$ h）")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    return cm.savefig(fig, "fig17_variants.png", plt)


def fig18_sensitivity():
    """fig18：烘干时间对关键参数的灵敏度（龙卷风图），读 data/sensitivity.csv（由
    sensitivity.py 生成，缺该文件时本图直接跳过）。左 (a) 问题 3、右 (b) 问题 4：取每个
    参数 +2% 扰动那一行，横轴是烘干时间的变化量 Δt_dry/min（正=变慢，红色为延长、蓝色为
    缩短），纵轴是参数名并按 |Δt_dry| 从大到小排列，用来指出主导不确定性的参数。
    """
    path = os.path.join(D, "sensitivity.csv")
    if not os.path.exists(path):
        print("   sensitivity.csv missing - run sensitivity.py first")
        return None
    df = pd.read_csv(path)
    # sensitivity.csv 的列：case（P3/P4 对应问题 3、4）、quantity（被扰动量的英文标识）、
    # change（扰动量，如 "+2%"）、t_dry_h（该扰动下的烘干时间/h）、delta_h（相对基准的变化量/h）
    fig, ax = plt.subplots(1, 2, figsize=(11.6 * FIG_SCALE, 5.2 * FIG_SCALE))
    for k, case in enumerate(("P3", "P4")):
        sub = df[(df.case == case) & (df.change == "+2%")].copy()
        sub["absd"] = sub.delta_h.abs()
        sub = sub.sort_values("absd")
        y = np.arange(len(sub))
        ax[k].barh(y, sub.delta_h * 60.0, color=np.where(sub.delta_h > 0,
                                                         "#c0392b", "#1f6fb2"))
        ax[k].set_yticks(y)
        ax[k].set_yticklabels(sub.quantity, fontsize=8)
        ax[k].axvline(0, color="k", lw=0.8)
        ax[k].set_xlabel(r"$\Delta t_{dry}$ / min  （参数 +2%）")
        ax[k].set_title("(%s) %s 的主导不确定性" % ("ab"[k], case))
    fig.tight_layout()
    return cm.savefig(fig, "fig18_sensitivity.png", plt)


def fig19_room_sensitivity():
    """fig19：烘房环境条件对烘干时间的影响，同样读 data/sensitivity.csv（缺文件时跳过）。
    (a) 设定温度 T_set、(b) 设定含湿量 C_set、(c) 预热/恒温切换时刻三种扰动，纵轴都是
    烘干时间 t_dry/h，红色与蓝色曲线分别是问题 3 与问题 4；横轴直接取 CSV 里 change 列的
    字符串，所以刻度是离散的（如 "-2%"、"0%"、"+2%"）。
    """
    path = os.path.join(D, "sensitivity.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    fig, ax = plt.subplots(1, 3, figsize=(12.2 * FIG_SCALE, 3.6 * FIG_SCALE))
    for case, c in (("P3", "#c0392b"), ("P4", "#1f6fb2")):
        s = df[(df.case == case) & (df.quantity == "set-point temperature")]
        s = s.sort_values("change")
        ax[0].plot(s.change.astype(str), s.t_dry_h, "o-", lw=1.6, color=c,
                   label=case)
        s = df[(df.case == case) & (df.quantity == "set-point moisture")]
        s = s.sort_values("change")
        ax[1].plot(s.change.astype(str), s.t_dry_h, "o-", lw=1.6, color=c,
                   label=case)
        s = df[(df.case == case) & (df.quantity == "preheat->constant switch")]
        s = s.sort_values("change")
        ax[2].plot(s.change.astype(str), s.t_dry_h, "o-", lw=1.6, color=c,
                   label=case)
    for a, t in zip(ax, ["(a) 设定温度 $T_{set}$ 敏感度",
                         "(b) 设定湿度 $C_{set}$ 敏感度",
                         "(c) 预热/恒温切换时刻敏感度"]):
        a.set_ylabel("烘干时间 / h"); a.set_title(t); a.legend(fontsize=8.5)
    ax[0].set_xlabel(r"$\Delta T_{set}$ / K")
    ax[1].set_xlabel("相对变化")
    ax[2].set_xlabel("切换时刻")
    fig.tight_layout()
    return cm.savefig(fig, "fig19_room_sensitivity.png", plt)


# --------------------------------------------------------------------------
# fig20 —— 建模示意图（不含计算结果）
# --------------------------------------------------------------------------
def fig20_schematic():
    """fig20：三个子图说明建模对象与工艺设定。(a) 物理对象与环境：长 25 cm、半径 2 cm 的
    圆柱药材置于热风烘房中，L≫R 是无穷长圆柱假设的依据；(b) 横截面：轴对称且无轴向梯度，
    问题化简为一维径向模型（r、theta 的含义与对称轴）；(c) 两阶段工艺：纵轴为按附件 2
    实测半径画的 R(t)/R_0，横轴是真实时间 0~72 h，灰虚线为 4 h 的阶段分界。
    """
    fig = plt.figure(figsize=(11.0 * FIG_SCALE, 4.2 * FIG_SCALE))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.0, 1.15, 1.25], wspace=0.28)

    ax = fig.add_subplot(gs[0, 0])
    ax.add_patch(Rectangle((0.06, 0.10), 0.72, 0.80, fc="#f5f5f5",
                           ec="#7f8c8d", lw=1.2))
    ax.add_patch(Rectangle((0.20, 0.18), 0.44, 0.62, fc="#fdebd0",
                           ec="#b9770e", lw=1.4))
    ax.add_patch(Rectangle((0.30, 0.18), 0.24, 0.62, fc="#fad7a0",
                           ec="none"))
    ax.add_patch(Circle((0.42, 0.49), 0.028, fc="#b9770e", ec="none"))
    for xx in (0.28, 0.34, 0.50, 0.56):
        ax.annotate("", xy=(xx, 0.84), xytext=(xx, 0.92),
                    arrowprops=dict(arrowstyle="-|>", color="#2980b9", lw=1.5))
    ax.text(0.51, 0.955, "热风", ha="center", fontsize=9, color="#2980b9")
    ax.text(0.42, 0.49, "药材\n圆柱", ha="center", va="center", fontsize=9)
    ax.annotate("", xy=(0.86, 0.49), xytext=(0.66, 0.49),
                arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.2))
    ax.text(0.86, 0.545, "$R=2$ cm", ha="center", fontsize=8.5,
            color="#c0392b")
    ax.annotate("", xy=(0.42, 0.80), xytext=(0.42, 0.18),
                arrowprops=dict(arrowstyle="<->", color="#16a085", lw=1.2))
    ax.text(0.44, 0.30, "$L=25$ cm", fontsize=8.5, color="#16a085",
            ha="left", va="center")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_title("(a) 物理对象：长圆柱，$L\\gg R$")

    ax = fig.add_subplot(gs[0, 1])
    ax.set_aspect("equal")
    th = np.linspace(0, 2 * np.pi, 200)
    for rr, al in ((1.0, 0.18), (0.66, 0.30), (0.33, 0.45)):
        ax.fill(rr * np.cos(th), rr * np.sin(th), color="#b9770e", alpha=al)
    ax.plot(np.cos(th), np.sin(th), color="#7e5109", lw=1.6)
    ax.plot([0, 0], [0, 0], "k+")
    ax.annotate("", xy=(1.0, 0.0), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.3))
    ax.text(0.55, -0.13, "$r$", color="#c0392b", fontsize=11)
    ax.annotate("", xy=(np.cos(0.55), np.sin(0.55)), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color="#16a085", lw=1.1))
    ax.text(0.30, 0.52, r"$\theta$", color="#16a085", fontsize=10)
    ax.text(0, 0, "对称轴", fontsize=8, ha="center", va="bottom")
    ax.text(1.18, 0.10, r"$R$", fontsize=10, color="#7e5109")
    ax.annotate("", xy=(1.18, 0.0), xytext=(0.0, 0.0),
                arrowprops=dict(arrowstyle="-", color="#7e5109", lw=1.0))
    ax.set_xlim(-1.45, 1.5); ax.set_ylim(-1.2, 1.2); ax.axis("off")
    ax.set_title(r"(b) 横截面与一维径向简化（轴对称、无轴向梯度）")

    ax = fig.add_subplot(gs[0, 2])
    # (c) 子图的横轴用真实小时，纵轴是归一化半径：把附件 2 实测的 R(t) 除以 R_0 得到，
    # 4 h 处的灰虚线即预热平衡阶段与恒温干燥阶段的分界
    _rh = cm.RadiusHistory(mode="table")
    th_h = np.linspace(0.0, 72.0, 400)
    rr = np.array([_rh.R(x * 3600.0) for x in th_h]) / cm.R0
    ax.plot(th_h, rr, lw=2.0, color="#c0392b")
    ax.fill_between(th_h, 0, rr, color="#c0392b", alpha=0.12)
    ax.axvline(4.0, color="gray", ls="--", lw=1.2)
    ax.text(14.0, 0.22, "预热平衡阶段\n$0\\leq t\\leq 4$ h\n环境温湿度按附件1上升",
            fontsize=8.5, ha="center")
    ax.text(22.0, 0.62, "恒温干燥阶段\n$t>4$ h\n$T_{air}=50.21\\,^\\circ$C\n"
            "$C_{air}=0.0509$ kg/kg", fontsize=8.5, ha="left")
    ax.set_xlabel("时间 / h"); ax.set_ylabel("半径 $R(t)$ / $R_0$")
    ax.set_title("(c) 两阶段工艺与收缩（半径取附件2 实测）")
    ax.set_xlim(0, 72); ax.set_ylim(0, 1.12)
    return cm.savefig(fig, "fig20_schematic.png", plt)


# --------------------------------------------------------------------------
# fig21 —— 无量纲数：模型所处的工作区间
# --------------------------------------------------------------------------
def fig21_regime():
    """fig21：无量纲分析，物性由 common 的物性函数在 T=323.15 K（50 degC）下现算，横轴均为 C。
    (a) 传质 Biot 数 Bi_m=h_m*R_0/D（问题 1、3、4 各一条；问题 3、4 取 T=323.15 K，
    问题 1 的 D 只随 C 变化），Bi_m 不小于 1 说明过程由内部扩散控制；(b) 特征时间 tau_D=R_0^2/D 与 tau_T=R_0^2/alpha（单位 h），tau_T≪tau_D 说明
    药温先与环境平衡、干燥后完成；(c) 热 Biot 数 Bi_T=h*R_0/k 与 Lewis 数 Le=D/alpha。
    """
    C = np.linspace(0.05, 2.55, 300)
    Tk = np.full_like(C, 323.15)
    D3 = cm.props_problem23(C, Tk)[3]
    D4 = cm.props_problem4(C, Tk)[3]
    bi3 = cm.HM_CONV * cm.R0 / D3
    bi4 = cm.HM_CONV * cm.R0 / D4
    alpha = 0.36 / (820.0 * 2600.0)
    k3 = cm.props_problem23(C, Tk)[2]
    cp3 = cm.props_problem23(C, Tk)[1]
    rho3 = cm.props_problem23(C, Tk)[0]
    a3 = k3 / (rho3 * cp3)

    fig, ax = plt.subplots(1, 3, figsize=(12.2 * FIG_SCALE, 3.5 * FIG_SCALE))
    ax[0].semilogy(C, bi3, lw=1.8, label="问题3 ($Bi_m$)")
    ax[0].semilogy(C, bi4, lw=1.8, label="问题4 ($Bi_m$)")
    ax[0].axhline(1, color="k", ls=":", lw=1.1)          # Bi=1：内外阻力相当的分界
    # 问题 1 的传质 Biot 数：其 D 只随含水率变化，故画成曲线而不是水平线。
    # 不能用 h*R_0/k=1.39 代替——那是热 Biot 数，与纵轴定义的 Bi_m 不是同一个量。
    D1 = cm.props_problem1(C)[3]
    ax[0].semilogy(C, cm.HM_CONV * cm.R0 / D1, color="#c0392b", ls="--", lw=1.4,
                   label="问题1 ($Bi_m$)")
    ax[0].set_xlabel(r"水分浓度 $C$ / (kg/kg)")
    ax[0].set_ylabel(r"$Bi_m=h_mR/D$")
    ax[0].set_title("(a) 传质 Biot 数：内部扩散控制")
    ax[0].legend(fontsize=8)

    tau3 = cm.R0 ** 2 / D3 / 3600.0
    tau4 = cm.R0 ** 2 / D4 / 3600.0
    ax[1].semilogy(C, tau3, lw=1.8, label=r"$\tau_D=R^2/D$, 问题3")
    ax[1].semilogy(C, tau4, lw=1.8, label=r"$\tau_D$, 问题4")
    ax[1].semilogy(C, cm.R0 ** 2 / a3 / 3600.0, lw=1.8, ls="--",
                   label=r"$\tau_T=R^2/\alpha$")
    ax[1].set_xlabel(r"水分浓度 $C$ / (kg/kg)")
    ax[1].set_ylabel(r"特征时间 / h")
    ax[1].set_title(r"(b) 时间尺度：$\tau_T\ll\tau_D$（热平衡远快于干燥）")
    ax[1].legend(fontsize=8)

    ax[2].semilogy(C, cm.H_CONV * cm.R0 / k3, lw=1.8, label=r"$Bi_T=hR/k$")
    ax[2].semilogy(C, D3 / a3, lw=1.8, label=r"$Le=D/\alpha$")
    ax[2].set_xlabel(r"水分浓度 $C$ / (kg/kg)")
    ax[2].set_ylabel("无量纲数")
    ax[2].set_title("(c) 热 Biot 数与 Lewis 数")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    return cm.savefig(fig, "fig21_regime.png", plt)


def main():
    """按固定顺序调用 21 个绘图函数，每个函数把一幅 200 dpi 的 PNG 写进 figures/。"""
    print("=" * 70)
    print("GENERATING PAPER FIGURES")
    print("=" * 70)
    jobs = [fig01_room, fig02_residuals, fig03_p1_fields, fig04_p1_profiles,
            fig05_p1_verification, fig06_properties, fig07_p23_fields,
            fig08_p23_curves, fig09_p23_flux, fig10_p23_zoom,
            fig11_p4_radius, fig12_p4_fields, fig13_p3_vs_p4,
            fig14_shrink_consistency, fig15_crosscheck, fig16_water_balance,
            fig17_variants, fig18_sensitivity, fig19_room_sensitivity,
            fig20_schematic, fig21_regime]
    for f in jobs:
        print("->", f.__name__, flush=True)
        f()
    print("\nALL FIGURES WRITTEN TO", cm.FIG_DIR)


if __name__ == "__main__":
    main()
