# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s07_figures.py
# 作用  : 生成论文全部插图(30 张)。所有图均由真实求解结果绘制, 不使用示意数据;
#         示意图(fig06)除外, 其中标注的几何与边界条件与题面/附录完全一致。
#         绘图统一使用 plotstyle 的中文字体设置, 输出 220 dpi PNG。
# 产出  : figures/fig05_*.png ... fig30_*.png
# 运行  : python s07_figures.py
# =============================================================================
"""Generate every figure used by the paper (真实结果, 非示意数据)."""

from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drying_core as dc      # noqa: E402
import plotstyle as ps        # noqa: E402
import scenarios as sc        # noqa: E402
from s05_verification import (analytic_cylinder_robin, constant_env,  # noqa: E402
                              frozen_props)

P = ps.PALETTE
FIG = sc.DIR_FIG


def run_p1_snapshots():
    """问题 1: 短时全场算例, 含 T/C 快照。"""
    env = sc.load_env()
    ts = list(np.arange(60.0, 1801.0, 60.0))
    res = dc.simulate(400, 1800.0, 0.5, dc.props_appendix2, env, sc.fixed_radius(),
                      theta=1.0, record_dt=60.0, snapshot_times=ts)
    return res, np.array(ts)


def run_p2_snapshots():
    """问题 2: 3 h 全场算例。"""
    env = sc.load_env()
    ts = list(np.arange(900.0, 10801.0, 900.0))
    res = dc.simulate(400, 10800.0, 1.0, dc.props_appendix3, env, sc.fixed_radius(),
                      theta=1.0, record_dt=60.0, snapshot_times=ts)
    return res, np.array(ts)


def run_p3_full():
    """问题 3: 全过程算例(用于全场图与曲线)。"""
    env = sc.load_env()
    ts = list(np.arange(3600.0, 220000.0, 3600.0))
    res = dc.simulate(400, 400000.0, 5.0, dc.props_appendix3, env, sc.fixed_radius(),
                      theta=1.0, c_stop=dc.C_DRY, record_dt=60.0, snapshot_times=ts)
    return res


def run_p4_full():
    """问题 4: 全过程算例(物质坐标)。"""
    env = sc.load_env()
    ts = list(np.arange(1800.0, 220000.0, 1800.0))
    res = dc.simulate(400, 400000.0, 5.0, dc.props_appendix4, env,
                      sc.load_radius("linear"), theta=1.0, c_stop=dc.C_DRY,
                      record_dt=60.0, snapshot_times=ts)
    return res


# ---------------------------------------------------------------------------
def fig05_plateau():
    """平台边界口径: 实测末段均值 vs 一阶指数拟合渐近值。"""
    from scipy.optimize import curve_fit
    df = sc.load_env_table()
    t = df["t"].to_numpy(float)
    T = df["T"].to_numpy(float)
    C = df["C"].to_numpy(float)
    f = lambda x, a, b, c: a - (a - b) * np.exp(-x / c)   # noqa: E731
    pT, _ = curve_fit(f, t, T, p0=[50.2, 28.0, 1800.0], maxfev=20000)
    pC, _ = curve_fit(f, t, C, p0=[0.0509, 0.0196, 2700.0], maxfev=20000)
    tt = np.linspace(0, 14400, 400)
    fig, ax = plt.subplots(1, 2, figsize=(10.4, 3.4))
    ax[0].plot(t / 3600, T, color=P["gray"], lw=1.4, label="附件 1 实测")
    ax[0].plot(tt / 3600, f(tt, *pT), color=P["red"], lw=1.6, ls="--",
               label="一阶惯性拟合 $T_\\infty$=%.3f $^\\circ$C" % pT[0])
    ax[0].axhline(np.mean(T[-60:]), color=P["blue"], lw=1.4,
                  label="末段 3600 s 均值 %.3f $^\\circ$C" % np.mean(T[-60:]))
    ax[0].set_xlabel("时间 / h")
    ax[0].set_ylabel("烘房温度 / $^\\circ$C")
    ax[0].set_title("温度边界平台口径")
    ax[0].legend(loc="lower right", fontsize=8.5)
    ax[1].plot(t / 3600, C, color=P["gray"], lw=1.4, label="附件 1 实测")
    ax[1].plot(tt / 3600, f(tt, *pC), color=P["red"], lw=1.6, ls="--",
               label="一阶惯性拟合 $C_\\infty$=%.5f" % pC[0])
    ax[1].axhline(np.mean(C[-60:]), color=P["blue"], lw=1.4,
                  label="末段 3600 s 均值 %.5f" % np.mean(C[-60:]))
    ax[1].set_xlabel("时间 / h")
    ax[1].set_ylabel("烘房水分浓度 / (kg/kg)")
    ax[1].set_title("水分浓度边界平台口径")
    ax[1].legend(loc="lower right", fontsize=8.5)
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig05_plateau_estimators.png"))


def fig06_schematic():
    """模型示意图: 几何、径向坐标、边界条件与两阶段工艺。"""
    fig, ax = plt.subplots(1, 2, figsize=(10.6, 3.6))
    # 左: 圆柱截面
    ax[0].add_patch(Rectangle((0.32, 0.18), 0.36, 0.64, facecolor="#dbe5f1",
                              edgecolor=P["blue"], lw=1.5))
    ax[0].plot([0.5, 0.5], [0.18, 0.82], color=P["red"], lw=1.2, ls="--")
    ax[0].annotate("", xy=(0.5, 0.5), xytext=(0.32, 0.5),
                   arrowprops=dict(arrowstyle="->", color="k", lw=1.2))
    ax[0].text(0.40, 0.52, "$r$", fontsize=11)
    ax[0].text(0.50, 0.135, "轴线 $r=0$\\n$\\partial T/\\partial r=\\partial C/\\partial r=0$",
               ha="center", fontsize=9)
    ax[0].annotate("表面 $r=R(t)$\\n$-k\\,\\partial T/\\partial r=h(T_{air}-T_s)$\\n"
                   "$-D\\,\\partial C/\\partial r=h_m(C_{air}-C_s)$",
                   xy=(0.68, 0.5), xytext=(0.72, 0.60), fontsize=9,
                   arrowprops=dict(arrowstyle="->", color=P["red"], lw=1.2))
    for y in np.linspace(0.24, 0.76, 7):
        ax[0].annotate("", xy=(0.24, y), xytext=(0.10, y),
                       arrowprops=dict(arrowstyle="->", color=P["cyan"], lw=1.1))
    ax[0].text(0.06, 0.86, "热风", fontsize=9, color=P["cyan"])
    ax[0].text(0.30, 0.86, "长 25 cm $\\gg$ 半径 2 cm\\n→ 一维径向模型", fontsize=9)
    ax[0].set_xlim(0, 1)
    ax[0].set_ylim(0, 1)
    ax[0].axis("off")
    ax[0].set_title("(a) 几何与边界条件")
    # 右: 两阶段
    ax[1].axvspan(0, 4, color=P["orange"], alpha=0.13)
    ax[1].axvspan(4, 57.2, color=P["blue"], alpha=0.10)
    ax[1].plot([0, 1.2], [0, 0], color=P["gray"], lw=0)
    ax[1].annotate("", xy=(60, 0.2), xytext=(0, 0.2),
                   arrowprops=dict(arrowstyle="->", color="k", lw=1.2))
    ax[1].text(2, 0.28, "预热平衡阶段\\n0~4 h 环境由附件 1 给出", fontsize=9, ha="center")
    ax[1].text(30, 0.28, "恒温干燥阶段\\n4 h 以后环境取平台(模型假设)", fontsize=9, ha="center")
    ax[1].text(30, 0.05, "总时长 2~3 天(题面)", fontsize=9, ha="center", color=P["gray"])
    ax[1].set_xlim(-2, 62)
    ax[1].set_ylim(0, 0.5)
    ax[1].axis("off")
    ax[1].set_title("(b) 两阶段工艺与环境边界")
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig06_model_schematic.png"))


def fig07_properties():
    """物性经验关系曲线(附录 2/3/4)。"""
    Cv = np.linspace(0.05, 2.6, 400)
    T_C = 50.0
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.1))
    for tag, props, col in (("附录3", dc.props_appendix3, P["blue"]),
                            ("附录4", dc.props_appendix4, P["red"])):
        rho, cp, k, D = props(Cv, np.full_like(Cv, T_C))
        axes[0].plot(Cv, rho, color=col, lw=1.8, label=tag)
        axes[1].plot(Cv, cp, color=col, lw=1.8, label=tag)
        axes[2].plot(Cv, k, color=col, lw=1.8, label=tag)
        axes[3].plot(Cv, D, color=col, lw=1.8, label=f"{tag} (T=50$^\\circ$C)")
    r2, c2, k2, D2 = dc.props_appendix2(Cv, np.full_like(Cv, T_C))
    axes[3].plot(Cv, D2, color=P["green"], lw=1.8, ls="--", label="附录2 (T 无关)")
    axes[0].set_ylabel("密度 / (kg/m$^3$)")
    axes[1].set_ylabel("比热容 / (J/(kg·K))")
    axes[2].set_ylabel("导热系数 / (W/(m·K))")
    axes[3].set_ylabel("扩散系数 $D$ / (m$^2$/s)")
    axes[3].set_yscale("log")
    for ax, ttl in zip(axes, ["(a) 密度 $\\rho(C)$", "(b) 比热容 $c_p(C)$",
                              "(c) 导热系数 $k(C)$", "(d) 扩散系数 $D(C,T)$"]):
        ax.set_xlabel("含水率 $C$ / (kg/kg)")
        ax.set_title(ttl)
        ax.legend(fontsize=8.5)
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig07_property_curves.png"))


def fig08_dimensionless():
    """无量纲数与特征时间。"""
    Cv = np.linspace(0.15, 2.55, 200)
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.1))
    for tag, props, col in (("附录3", dc.props_appendix3, P["blue"]),
                            ("附录4", dc.props_appendix4, P["red"])):
        rows = [dc.dimensionless_numbers(props, C=float(c), T_C=50.0) for c in Cv]
        axes[0].plot(Cv, [r["Bi_mass"] for r in rows], color=col, lw=1.8, label=tag)
        axes[1].plot(Cv, [r["t_diff_mass"] / 3600 for r in rows], color=col, lw=1.8, label=tag)
        axes[2].plot(Cv, [r["D"] for r in rows], color=col, lw=1.8, label=tag)
    axes[0].set_ylabel("$Bi_m=h_mR/D$")
    axes[0].set_title("(a) 传质 Biot 数")
    axes[1].set_ylabel("$R^2/D$ / h")
    axes[1].set_title("(b) 传质特征时间")
    axes[2].set_ylabel("$D$ / (m$^2$/s)")
    axes[2].set_yscale("log")
    axes[2].set_title("(c) 扩散系数范围(与文献区间对照)")
    axes[2].axhspan(1e-11, 1e-8, color=P["green"], alpha=0.12)
    for ax in axes:
        ax.set_xlabel("含水率 $C$ / (kg/kg)")
        ax.legend(fontsize=8.5)
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig08_dimensionless.png"))


def fig09_10_p1(res, ts):
    """问题 1 剖面图与热图。"""
    R = res.R_final
    r_cm = res.xi_c * R * 100
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.5))
    cols = plt.cm.viridis(np.linspace(0.05, 0.9, 7))
    for k, tc in enumerate([60, 300, 600, 900, 1200, 1500, 1800]):
        i = int(np.argmin(np.abs(ts - tc)))
        axes[0].plot(r_cm, res.snap_T[i], color=cols[k], lw=1.6, label=f"{tc} s")
        axes[1].plot(r_cm, res.snap_C[i], color=cols[k], lw=1.6, label=f"{tc} s")
    axes[0].set_ylabel("温度 / $^\\circ$C")
    axes[1].set_ylabel("含水率 / (kg/kg)")
    for ax, ttl in zip(axes, ["(a) 温度径向剖面", "(b) 含水率径向剖面"]):
        ax.set_xlabel("到药材中心的距离 / cm")
        ax.set_title(ttl)
        ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    ps.save(fig, os.path.join(FIG, "fig09_p1_profiles.png"))

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.4))
    tt = ts
    for ax, key, lab, cmap in ((axes[0], "T", "温度 / $^\\circ$C", "inferno"),
                               (axes[1], "C", "含水率 / (kg/kg)", "viridis")):
        F = getattr(res, f"snap_{key}")
        im = ax.contourf(r_cm, tt / 60.0, F, levels=24, cmap=cmap)
        ax.set_xlabel("到药材中心的距离 / cm")
        ax.set_ylabel("时间 / min")
        fig.colorbar(im, ax=ax, label=lab)
        ax.grid(False)
    axes[0].set_title("(a) 温度场 $T(r,t)$")
    axes[1].set_title("(b) 含水率场 $C(r,t)$")
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig10_p1_heatmap.png"))


def fig11_surface_reconstruction(res):
    """表面重构: 最后单元中心值 vs 半单元重构表面值。"""
    xi_last = 1.0 - 0.5 / res.n_cells              # 最后一个控制体中心的归一化坐标
    idx, w = dc.lagrange_weights_grid(res.n_cells, [xi_last])
    last_cell = np.einsum("a,na->n", w[0], res.C[:, idx[0]])
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.3))
    t = res.times / 60.0
    axes[0].plot(t, last_cell, color=P["orange"], lw=1.7, label="最后一个单元中心值")
    axes[0].plot(t, res.c_surf, color=P["blue"], lw=1.7, label="半单元重构表面值 $C_s$")
    axes[0].set_xlabel("时间 / min")
    axes[0].set_ylabel("含水率 / (kg/kg)")
    axes[0].set_title("(a) 表面含水率两种取法的差别")
    axes[0].legend(fontsize=8.5)
    axes[1].plot(t, res.c_surf - last_cell, color=P["red"], lw=1.6)
    axes[1].set_xlabel("时间 / min")
    axes[1].set_ylabel("$C_s-C_{N-1}$ / (kg/kg)")
    axes[1].set_title("(b) 差值随时间变化(40 min 后趋于稳定)")
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig11_p1_surface_reconstruction.png"))


def fig12_p1_curves(res):
    """问题 1 中心/表面曲线。"""
    idx, w = dc.lagrange_weights_grid(res.n_cells, [0.0])
    c0 = dc.apply_weights(res.C, idx, w)
    t0 = dc.apply_weights(res.T, idx, w)
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.3))
    axes[0].plot(res.times / 60, t0, color=P["red"], lw=1.8, label="中心 $r=0$")
    axes[0].plot(res.times / 60, res.t_surf, color=P["blue"], lw=1.8, label="表面 $r=R$")
    axes[0].plot(res.times / 60, res.step_Tair[:res.times.size], color=P["gray"],
                 lw=1.3, ls="--", label="环境 $T_{air}$")
    axes[0].set_ylabel("温度 / $^\\circ$C")
    axes[1].plot(res.times / 60, c0, color=P["red"], lw=1.8, label="中心 $r=0$")
    axes[1].plot(res.times / 60, res.c_surf, color=P["blue"], lw=1.8, label="表面 $r=R$")
    axes[1].plot(res.times / 60, res.step_Cair[:res.times.size], color=P["gray"],
                 lw=1.3, ls="--", label="环境 $C_{air}$")
    axes[1].set_ylabel("含水率 / (kg/kg)")
    for ax, ttl in zip(axes, ["(a) 温度随时间变化", "(b) 含水率随时间变化"]):
        ax.set_xlabel("时间 / min")
        ax.set_title(ttl)
        ax.legend(fontsize=8.5)
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig12_p1_center_surface.png"))


def fig13_14_p2(res2, ts2):
    """问题 2: 剖面与热图。"""
    r_cm = res2.xi_c * res2.R_final * 100
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.5))
    for h in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
        i = int(np.argmin(np.abs(ts2 / 3600 - h)))
        axes[0].plot(r_cm, res2.snap_T[i], lw=1.6, label=f"{h:g} h")
        axes[1].plot(r_cm, res2.snap_C[i], lw=1.6, label=f"{h:g} h")
    axes[0].set_ylabel("温度 / $^\\circ$C")
    axes[1].set_ylabel("含水率 / (kg/kg)")
    for ax, ttl in zip(axes, ["(a) 温度径向剖面", "(b) 含水率径向剖面"]):
        ax.set_xlabel("到药材中心的距离 / cm")
        ax.set_title(ttl)
        ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    ps.save(fig, os.path.join(FIG, "fig13_p2_profiles.png"))

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.4))
    tt = ts2
    for ax, key, lab, cmap in ((axes[0], "T", "温度 / $^\\circ$C", "inferno"),
                               (axes[1], "C", "含水率 / (kg/kg)", "viridis")):
        F = getattr(res2, f"snap_{key}")
        im = ax.contourf(r_cm, tt / 3600.0, F, levels=24, cmap=cmap)
        ax.set_xlabel("到药材中心的距离 / cm")
        ax.set_ylabel("时间 / h")
        fig.colorbar(im, ax=ax, label=lab)
        ax.grid(False)
    axes[0].set_title("(a) 温度场 $T(r,t)$")
    axes[1].set_title("(b) 含水率场 $C(r,t)$")
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig14_p2_heatmap.png"))


def fig15_18_p3(res3):
    """问题 3: 全过程曲线、热图、通量、温度。"""
    idx, w = dc.lagrange_weights_grid(res3.n_cells, [0.0])
    c0 = dc.apply_weights(res3.C, idx, w)
    t0 = dc.apply_weights(res3.T, idx, w)
    th = res3.times / 3600.0
    td = sc.drying_time(res3)

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.plot(th, c0, color=P["red"], lw=2.0, label="中心 $r=0$")
    ax.plot(th, res3.c_surf, color=P["blue"], lw=2.0, label="表面 $r=R_0$")
    cbar = np.array([float(np.sum(res3.C[i] * res3.grid["vol"])) for i in range(res3.C.shape[0])])
    ax.plot(th, cbar, color=P["green"], lw=1.6, ls="-.", label="截面平均")
    ax.axhline(dc.C_DRY, color="k", ls=":", lw=1.4)
    ax.text(2, dc.C_DRY + 0.05, "干燥判据 $C=0.15$", fontsize=9)
    ax.axvline(td / 3600.0, color=P["orange"], ls="--", lw=1.5)
    ax.annotate(f"$t_3$ = {td/3600.0:.2f} h", xy=(td / 3600.0, 1.4),
                xytext=(td / 3600.0 - 22, 1.9), fontsize=10,
                arrowprops=dict(arrowstyle="->", color=P["orange"]))
    ax.set_xlabel("时间 / h")
    ax.set_ylabel("含水率 / (kg/kg)")
    ax.set_title("问题 3: 全过程含水率演化")
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    ps.save(fig, os.path.join(FIG, "fig15_p3_curves_full.png"))

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 3.5))
    r_cm = res3.xi_c * res3.R_final * 100
    tt = np.concatenate([[0.0], res3.snap_times])
    for ax, key, lab, cmap in ((axes[0], "T", "温度 / $^\\circ$C", "inferno"),
                               (axes[1], "C", "含水率 / (kg/kg)", "viridis")):
        F = np.vstack([getattr(res3, f"snap_{key}")[0], getattr(res3, f"snap_{key}")])
        im = ax.contourf(r_cm, tt / 3600.0, F, levels=24, cmap=cmap)
        ax.set_xlabel("到药材中心的距离 / cm")
        ax.set_ylabel("时间 / h")
        fig.colorbar(im, ax=ax, label=lab)
        ax.grid(False)
    axes[0].set_title("(a) 全场温度 $T(r,t)$")
    axes[1].set_title("(b) 全场含水率 $C(r,t)$(内部长期偏高 → 内扩散控制)")
    fig.tight_layout()
    ps.save(fig, os.path.join(FIG, "fig16_p3_heatmap.png"))

    t = res3.step_times
    rate = np.gradient(-res3.step_maxC, t)
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.3))
    axes[0].plot(t / 3600, -rate * 3600, color=P["blue"], lw=1.7)
    axes[0].set_xlabel("时间 / h")
    axes[0].set_ylabel("干燥速率 / (kg/kg/h)")
    axes[0].set_title("(a) 中心含水率下降速率(降速干燥特征)")
    axes[1].plot(t / 3600, res3.cum_flux, color=P["green"], lw=1.8,
                 label="累计表面逸出(归一化)")
    axes[1].plot(res3.times / 3600,
                 res3.water[0] - res3.water[:res3.times.size], color=P["red"],
                 lw=1.5, ls="--", label="总水量减少量")
    axes[1].set_xlabel("时间 / h")
    axes[1].set_ylabel("归一化水量")
    axes[1].set_title("(b) 水量收支闭合检验")
    axes[1].legend(fontsize=8.5)
    fig.tight_layout()
    ps.save(fig, os.path.join(FIG, "fig17_p3_flux.png"))

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.3))
    axes[0].plot(th, t0, color=P["red"], lw=1.8, label="中心")
    axes[0].plot(th, res3.t_surf, color=P["blue"], lw=1.8, label="表面")
    axes[0].plot(res3.step_times / 3600,
                 res3.step_Tair[:res3.step_times.size], color=P["gray"], lw=1.2,
                 ls="--", label="环境")
    axes[0].set_ylabel("温度 / $^\\circ$C")
    axes[0].set_title("(a) 全过程温度演化")
    axes[1].plot(th, res3.c_surf - res3.step_Cair[:th.size], color=P["purple"], lw=1.8)
    axes[1].set_ylabel("$C_s-C_{air}$ / (kg/kg)")
    axes[1].set_yscale("log")
    axes[1].set_title("(b) 传质推动力(对数坐标, 后期迅速衰减)")
    for ax in axes:
        ax.set_xlabel("时间 / h")
        ax.legend(fontsize=8.5)
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig18_p3_temperature.png"))


def fig19_21_p4(res4, res3):
    """问题 4: 曲线、物质坐标热图与 P3/P4 对比。"""
    idx, w = dc.lagrange_weights_grid(res4.n_cells, [0.0])
    c0 = dc.apply_weights(res4.C, idx, w)
    th = res4.times / 3600.0
    td4 = sc.drying_time(res4)
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 3.4))
    axes[0].plot(th, c0, color=P["red"], lw=2.0, label="中心")
    axes[0].plot(th, res4.c_surf, color=P["blue"], lw=2.0, label="表面")
    axes[0].axhline(dc.C_DRY, color="k", ls=":", lw=1.3)
    axes[0].axvline(td4 / 3600.0, color=P["orange"], ls="--", lw=1.4)
    axes[0].annotate(f"$t_4$ = {td4/3600.0:.2f} h", xy=(td4 / 3600.0, 1.4),
                     xytext=(td4 / 3600.0 - 20, 1.9), fontsize=10,
                     arrowprops=dict(arrowstyle="->", color=P["orange"]))
    axes[0].set_xlabel("时间 / h")
    axes[0].set_ylabel("含水率 / (kg/kg)")
    axes[0].set_title("(a) 问题 4 含水率演化(含收缩)")
    axes[0].legend(fontsize=8.5)
    ax2 = axes[1]
    ax2.plot(th, res4.R * 100, color=P["green"], lw=1.9, label="药材半径 $R(t)$")
    ax2.set_xlabel("时间 / h")
    ax2.set_ylabel("半径 / cm", color=P["green"])
    ax2.tick_params(axis="y", labelcolor=P["green"])
    ax2b = ax2.twinx()
    ax2b.plot(th, (res4.R / res4.R[0]) ** 2, color=P["purple"], lw=1.5, ls="--",
              label="相对体积 $V/V_0$")
    ax2b.set_ylabel("$V/V_0$", color=P["purple"])
    ax2b.tick_params(axis="y", labelcolor=P["purple"])
    ax2b.grid(False)
    ax2.set_title("(b) 收缩轨迹")
    fig.tight_layout()
    ps.save(fig, os.path.join(FIG, "fig19_p4_curves.png"))

    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    tt = np.concatenate([[0.0], res4.snap_times])
    F = np.vstack([res4.snap_C[0], res4.snap_C])
    im = ax.contourf(res4.xi_c, tt / 3600.0, F, levels=24, cmap="viridis")
    ax.set_xlabel("归一化半径 $\\xi=r/R(t)$")
    ax.set_ylabel("时间 / h")
    ax.set_title("问题 4: 物质坐标下的含水率场(全部材料点走完相同轨迹)")
    fig.colorbar(im, ax=ax, label="含水率 / (kg/kg)")
    ax.grid(False)
    fig.tight_layout()
    ps.save(fig, os.path.join(FIG, "fig20_p4_heatmap.png"))

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.3))
    i3, w3 = dc.lagrange_weights_grid(res3.n_cells, [0.0])
    c3 = dc.apply_weights(res3.C, i3, w3)
    axes[0].plot(res3.times / 3600, c3, color=P["blue"], lw=1.9, label="问题 3(不收缩)")
    axes[0].plot(th, c0, color=P["red"], lw=1.9, label="问题 4(收缩)")
    axes[0].axhline(dc.C_DRY, color="k", ls=":", lw=1.2)
    axes[0].set_xlabel("时间 / h")
    axes[0].set_ylabel("中心含水率 / (kg/kg)")
    axes[0].set_title("(a) 中心含水率对比")
    axes[0].legend(fontsize=8.5)
    axes[1].bar(["问题 3", "问题 4"], [sc.drying_time(res3) / 3600, td4],
                color=[P["blue"], P["red"]], width=0.45)
    for i, v in enumerate([sc.drying_time(res3) / 3600, td4]):
        axes[1].text(i, v + 0.6, f"{v:.2f} h", ha="center", fontsize=10)
    axes[1].set_ylabel("干燥时间 / h")
    axes[1].set_title("(b) 干燥时间对比(收缩使干燥时间缩短)")
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig21_p3_vs_p4.png"))


def fig22_analytic():
    """解析解对照图。"""
    R = dc.R0
    times = np.array([60.0, 300.0, 900.0, 1800.0])
    r = np.linspace(0, R, 60)
    props = frozen_props(5e-9)
    env = constant_env(50.0, 0.05)
    res = dc.simulate(800, 1800.0, 0.5, props, env, sc.fixed_radius(), theta=1.0,
                      record_dt=60.0, snapshot_times=list(times))
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.3))
    for k, tt in enumerate(times):
        i = int(np.argmin(np.abs(res.snap_times - tt)))
        ana = analytic_cylinder_robin(r, np.array([tt]), R, 5e-9, dc.HM_CONV,
                                      0.05, 2.55)[0]
        axes[0].plot(r * 100, ana, color=P["blue"], lw=2.4, alpha=0.45)
        axes[0].plot(res.xi_c * R * 100, res.snap_C[i], color=P["red"], lw=1.3, ls="--")
        axes[1].plot(r * 100, np.abs(np.interp(r, res.xi_c * R, res.snap_C[i]) - ana),
                     lw=1.6, label=f"{tt:.0f} s")
    axes[0].plot([], [], color=P["blue"], lw=2.4, alpha=0.6, label="解析解(Bessel 级数)")
    axes[0].plot([], [], color=P["red"], lw=1.3, ls="--", label="数值解(本文 FVM)")
    axes[0].set_xlabel("到药材中心的距离 / cm")
    axes[0].set_ylabel("含水率 / (kg/kg)")
    axes[0].set_title("(a) 常物性算例的解析/数值对照")
    axes[0].legend(fontsize=8.5)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("到药材中心的距离 / cm")
    axes[1].set_ylabel("绝对误差 / (kg/kg)")
    axes[1].set_title("(b) 逐点绝对误差")
    axes[1].legend(fontsize=8)
    tab = pd.read_csv(os.path.join(sc.DIR_DATA, "verification_analytic.csv"))
    for field, col in (("mass", P["blue"]), ("heat", P["red"])):
        sub = tab[tab["field"] == field]
        axes[2].loglog(sub["N"], sub["max_abs_err"], "o-", color=col, lw=1.7,
                       label=f"{'传质' if field=='mass' else '传热'} 最大误差")
    axes[2].loglog([200, 1600], [3e-2, 3e-2 / 64], "k--", lw=1.2,
                   label="参考斜率 $-2$(二阶)")
    axes[2].set_xlabel("控制体数 $N$")
    axes[2].set_ylabel("相对解析解的最大误差")
    axes[2].set_title("(c) 收敛精度(空间二阶)")
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig22_analytic_verification.png"))


def fig23_convergence():
    """网格与时间步收敛证据汇总。"""
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.3))
    m = pd.read_csv(os.path.join(sc.DIR_DATA, "p1_mesh_convergence.csv"))
    axes[0].loglog(m["N"], m["max|dC|"], "o-", color=P["blue"], lw=1.7, label="含水率")
    axes[0].loglog(m["N"], m["max|dT|"], "s-", color=P["red"], lw=1.7, label="温度")
    axes[0].loglog([200, 3200], [3e-4, 3e-4 / 256], "k--", lw=1.1, label="参考斜率 $-2$")
    axes[0].set_xlabel("控制体数 $N$")
    axes[0].set_ylabel("相对最细网格的最大偏差")
    axes[0].set_title("(a) 问题 1 网格收敛")
    axes[0].legend(fontsize=8.5)
    d = pd.read_csv(os.path.join(sc.DIR_DATA, "p1_dt_convergence.csv"))
    axes[1].loglog(d["dt/s"], np.abs(d["T_center"] - d["T_center"].iloc[-1]) + 1e-9,
                   "o-", color=P["red"], lw=1.7, label="$T$ 中心")
    axes[1].loglog(d["dt/s"], np.abs(d["C_surface"] - d["C_surface"].iloc[-1]) + 1e-12,
                   "s-", color=P["blue"], lw=1.7, label="$C$ 表面")
    axes[1].loglog([0.02, 10], [4e-6, 4e-6 * 500], "k--", lw=1.1, label="参考斜率 1")
    axes[1].set_xlabel("时间步 $\\Delta t$ / s")
    axes[1].set_ylabel("相对最细时间步的偏差")
    axes[1].set_title("(b) 问题 1 时间步收敛(隐式 Euler 一阶)")
    axes[1].legend(fontsize=8.5)
    p4 = pd.read_csv(os.path.join(sc.DIR_DATA, "p4_convergence.csv"))
    dd = p4[p4["kind"] == "dt"]
    mm = p4[p4["kind"] == "mesh"]
    axes[2].semilogx(dd["value"], dd["t_dry/h"], "o-", color=P["blue"], lw=1.7, label="时间步(N=1600)")
    axes[2].semilogx(mm["value"], mm["t_dry/h"], "s-", color=P["red"], lw=1.7, label="网格(dt=10 s)")
    axes[2].set_xlabel("$\\Delta t$ / s 或 $N$")
    axes[2].set_ylabel("问题 4 干燥时间 / h")
    axes[2].set_title("(c) 问题 4 干燥时间的双重收敛")
    axes[2].legend(fontsize=8.5)
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig23_convergence.png"))


def fig24_residuals():
    """合成数据反演残差的统计特性。"""
    from scipy import stats
    df = pd.read_csv(os.path.join(sc.DIR_DATA, "inverse_residuals.csv"))
    r = df["resid"].to_numpy(float)
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.2))
    axes[0].hist(r, bins=18, density=True, color=P["blue"], alpha=0.75)
    xs = np.linspace(r.min(), r.max(), 200)
    axes[0].plot(xs, stats.norm.pdf(xs, r.mean(), r.std(ddof=1)), color=P["red"], lw=1.8)
    axes[0].set_xlabel("反演残差 / (kg/kg)")
    axes[0].set_ylabel("概率密度")
    axes[0].set_title("(a) 残差直方图(对照正态)")
    stats.probplot(r, dist="norm", plot=axes[1])
    axes[1].get_lines()[0].set(color=P["blue"], ms=3.2)
    axes[1].get_lines()[1].set(color=P["red"], lw=1.6)
    axes[1].set_title("(b) 残差 Q-Q 图")
    axes[1].set_xlabel("理论分位数")
    axes[1].set_ylabel("样本分位数")
    ac = [np.corrcoef(r[:-k], r[k:])[0, 1] for k in range(1, 13)]
    axes[2].bar(range(1, 13), ac, color=P["green"], alpha=0.85)
    axes[2].axhline(1.96 / np.sqrt(r.size), color=P["gray"], ls="--", lw=1.1)
    axes[2].axhline(-1.96 / np.sqrt(r.size), color=P["gray"], ls="--", lw=1.1)
    axes[2].set_xlabel("滞后阶数")
    axes[2].set_ylabel("自相关")
    axes[2].set_title("(c) 残差自相关(虚线为 95% 白噪带)")
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig24_inversion_residuals.png"))


def fig25_inverse():
    """反演目标函数面与置信椭圆。"""
    df = pd.read_csv(os.path.join(sc.DIR_DATA, "inverse_objective.csv"), index_col=0)
    est = pd.read_csv(os.path.join(sc.DIR_DATA, "inverse_estimate.csv")).iloc[0]
    hm = df.columns.to_numpy(float)          # 1e-7 m/s
    ds = df.index.to_numpy(float)
    Z = df.to_numpy(float).T
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.6))
    cs = axes[0].contourf(hm, ds, Z, levels=22, cmap="viridis_r")
    axes[0].contour(hm, ds, Z, levels=22, colors="w", linewidths=0.4, alpha=0.5)
    axes[0].plot(hm[np.unravel_index(Z.argmin(), Z.shape)[0]],
                 ds[np.unravel_index(Z.argmin(), Z.shape)[1]], "r*", ms=13,
                 label="网格极小点")
    axes[0].plot(est["hm_true"] * 1e7, est["dscale_true"], "w^", ms=10,
                 markeredgecolor="k", label="真值")
    axes[0].set_xlabel("$h_m$ / $10^{-7}$ m/s")
    axes[0].set_ylabel("$D$ 乘性尺度")
    axes[0].set_title("(a) 残差平方和(SSE)目标函数面")
    axes[0].legend(fontsize=8.5)
    fig.colorbar(cs, ax=axes[0], label="SSE / (kg/kg)$^2$")
    axes[0].grid(False)
    # 置信椭圆(线性化)
    d_scale = est["dscale_fit"]
    sd = np.array([est["hm_sd"] / 1e-7, est["dscale_sd"]])
    s = est["corr"]
    th = np.linspace(0, 2 * np.pi, 200)
    for lev, col, lab in ((1.0, P["orange"], "1$\\sigma$"), (2.447, P["red"], "95%")):
        x = est["hm_fit"] / 1e-7 + lev * sd[0] * np.cos(th)
        y = d_scale + lev * sd[1] * (s * np.cos(th) + np.sqrt(1 - s ** 2) * np.sin(th))
        axes[1].plot(x, y, color=col, lw=1.8, label=lab)
    axes[1].plot(est["hm_true"] * 1e7, est["dscale_true"], "k^", ms=10, label="真值")
    axes[1].plot(est["hm_fit"] / 1e-7, est["dscale_fit"], "r*", ms=13, label="估计值")
    axes[1].set_xlabel("$h_m$ / $10^{-7}$ m/s")
    axes[1].set_ylabel("$D$ 乘性尺度")
    axes[1].set_title(f"(b) 参数置信椭圆(相关系数 {s:+.2f})")
    axes[1].legend(fontsize=8.5)
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig25_inverse_objective.png"))


def fig26_28_sensitivity():
    """灵敏度: 龙卷风图、平台口径、Monte Carlo。"""
    oat = pd.read_csv(os.path.join(sc.DIR_DATA, "sensitivity_oat.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.6))
    for ax, problem, col in ((axes[0], "P3", P["blue"]), (axes[1], "P4", P["red"])):
        sub = oat[(oat["problem"] == problem) & (oat["param"] != "base")]
        # 每个参数取 |Δ| 最大的一次(即 ±2% 那一档), 构成龙卷风图
        sub = (sub.loc[sub.groupby("param")["delta_h"].apply(
            lambda s: s.abs().idxmax())].sort_values("delta_h"))
        ax.barh(sub["param"], sub["delta_h"], color=col, alpha=0.85)
        ax.axvline(0, color="k", lw=1.0)
        ax.set_xlabel("干燥时间变化 $\\Delta t$ / h")
        ax.set_title(f"({'(a)' if problem=='P3' else '(b)'} ) 问题 {problem[-1]} 单因素扰动响应")
    fig.tight_layout()
    ps.save(fig, os.path.join(FIG, "fig26_sensitivity_tornado.png"))

    pl = pd.read_csv(os.path.join(sc.DIR_DATA, "sensitivity_plateau.csv"))
    ca = pd.read_csv(os.path.join(sc.DIR_DATA, "sensitivity_cair.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.4))
    for problem, col in (("P3", P["blue"]), ("P4", P["red"])):
        sub = pl[pl["problem"] == problem]
        axes[0].plot(range(len(sub)), sub["t_dry_h"], "o-", color=col, lw=1.7, label=problem)
    axes[0].set_xticks(range(len(pl[pl["problem"] == "P3"])))
    axes[0].set_xticklabels(pl[pl["problem"] == "P3"]["variant"], rotation=20,
                            ha="right", fontsize=8)
    axes[0].set_ylabel("干燥时间 / h")
    axes[0].set_title("(a) 平台边界口径的影响")
    axes[0].legend(fontsize=8.5)
    for problem, col in (("P3", P["blue"]), ("P4", P["red"])):
        sub = ca[ca["problem"] == problem]
        axes[1].plot(sub["C_air_scale"], sub["t_dry_h"], "o-", color=col, lw=1.7,
                     label=problem)
    axes[1].set_xlabel("$C_{air}$ 乘性系数")
    axes[1].set_ylabel("干燥时间 / h")
    axes[1].set_title("(b) 空气含湿量口径的影响")
    axes[1].legend(fontsize=8.5)
    fig.tight_layout()
    ps.save(fig, os.path.join(FIG, "fig27_sensitivity_assumptions.png"))

    mc_path = os.path.join(sc.DIR_DATA, "sensitivity_montecarlo.csv")
    ct_path = os.path.join(sc.DIR_DATA, "sensitivity_mc_contrib.csv")
    if not (os.path.exists(mc_path) and os.path.exists(ct_path)):
        print("[WARN] missing Monte Carlo results, skip fig28")
        return None
    mc = pd.read_csv(mc_path)
    ct = pd.read_csv(ct_path)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.2))
    for problem, col in (("P3", P["blue"]), ("P4", P["red"])):
        sub = mc[mc["problem"] == problem]
        axes[0].hist(sub["t_dry_h"], bins=12, alpha=0.6, color=col, density=True,
                     label=problem)
    axes[0].set_xlabel("干燥时间 / h")
    axes[0].set_ylabel("概率密度")
    axes[0].set_title("(a) Monte Carlo 输出分布")
    axes[0].legend(fontsize=8.5)
    ents = {"h": "$h$", "h_m": "$h_m$", "D_scale": "$D$ 尺度",
            "T_plateau": "$T_{air}$ 平台", "C_plateau": "$C_{air}$ 平台"}
    for ax, problem, col in ((axes[1], "P3", P["blue"]), (axes[2], "P4", P["red"])):
        sub = ct[ct["problem"] == problem].set_index("param").loc[list(ents)]
        ax.bar([ents[k] for k in sub.index], sub["share_pct"], color=col, alpha=0.85)
        ax.set_ylabel("方差贡献 / %")
        ax.set_title(f"({'bc'[problem=='P4']}) {problem} 方差分解")
        ax.tick_params(axis="x", labelrotation=20)
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig28_montecarlo.png"))


def fig29_drymass():
    """干物质守恒检验图。"""
    df = pd.read_csv(os.path.join(sc.DIR_DATA, "drymass_consistency.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.3))
    axes[0].plot(df["t_h"], df["R_measured_cm"], color=P["blue"], lw=1.9,
                 label="附件 2 实测半径")
    axes[0].plot(df["t_h"], df["R_implied_cm"], color=P["red"], lw=1.7, ls="--",
                 label="由 $\\rho(C)$ 与干物质守恒反推")
    axes[0].set_xlabel("时间 / h")
    axes[0].set_ylabel("半径 / cm")
    axes[0].set_title("(a) 实测收缩轨迹与密度式隐含轨迹")
    axes[0].legend(fontsize=8.5)
    axes[1].plot(df["t_h"], df["m_dry"] * 1000, color=P["green"], lw=1.9)
    axes[1].set_xlabel("时间 / h")
    axes[1].set_ylabel("单位长度干物质质量 / (g/m)")
    axes[1].set_title("(b) 干物质守恒残差(理想应为常数)")
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig29_drymass_consistency.png"))


def fig30_latent():
    """潜热扩展模型的影响。"""
    env = sc.load_env()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.3))
    for lat, col, lab in ((0.0, P["blue"], "基准(不含潜热)"),
                          (2.38e6, P["red"], "扩展(含汽化潜热 2.38 MJ/kg)")):
        r = dc.simulate(400, 200000.0, 10.0, dc.props_appendix3, env,
                        sc.fixed_radius(), theta=1.0, c_stop=dc.C_DRY,
                        record_dt=300.0, latent_heat=lat)
        idx, w = dc.lagrange_weights_grid(r.n_cells, [0.0])
        c0 = dc.apply_weights(r.C, idx, w)
        axes[0].plot(r.times / 3600, c0, color=col, lw=1.8, label=lab)
        axes[1].plot(r.times / 3600, dc.apply_weights(r.T, idx, w), color=col,
                     lw=1.8, label=lab)
        axes[0].axvline(sc.drying_time(r) / 3600, color=col, ls=":", lw=1.2)
    axes[0].axhline(dc.C_DRY, color="k", ls=":", lw=1.2)
    axes[0].set_xlabel("时间 / h")
    axes[0].set_ylabel("中心含水率 / (kg/kg)")
    axes[0].set_title("(a) 潜热对干燥进程的影响")
    axes[0].legend(fontsize=8.5)
    axes[1].set_xlabel("时间 / h")
    axes[1].set_ylabel("中心温度 / $^\\circ$C")
    axes[1].set_title("(b) 潜热对温度演化的影响")
    axes[1].legend(fontsize=8.5)
    fig.tight_layout()
    return ps.save(fig, os.path.join(FIG, "fig30_latent_heat.png"))


def _run_with_retry(fn, dts=(5.0, 2.0, 1.0)):
    """以逐步减小的时间步重试含移动边界的长时间算例。

    目的: 个别配置下隐式 Euler 在移动边界附近可能出现非有限值,
    此时用更小的时间步重新推进即可避开(收敛序列本身表明结果不受影响)。
    每一步的失败原因都会被打印, 不做静默吞掉。
    """
    last = None
    for dt in dts:
        try:
            return fn(dt)
        except (FloatingPointError, ValueError) as e:      # noqa: PERF203
            last = e
            print(f"[WARN] dt={dt} 失败: {type(e).__name__}: {e}; 改用更小步长重试")
    raise last


def _safe(name, fn):
    """执行单个绘图任务, 失败只记录不中断其余图。"""
    try:
        out = fn()
        print("OK", name, out)
        return out
    except Exception as e:                                  # noqa: BLE001
        print(f"[FAIL] {name}: {type(e).__name__}: {e}")
        return None


def main():
    ps.setup()
    _safe("fig05", fig05_plateau)
    _safe("fig06", fig06_schematic)
    _safe("fig07", fig07_properties)
    _safe("fig08", fig08_dimensionless)
    store = {}
    _safe("fig09/10", lambda: store.setdefault(
        "p1", run_p1_snapshots()) and fig09_10_p1(*store["p1"]))
    if "p1" in store:
        _safe("fig11", lambda: fig11_surface_reconstruction(store["p1"][0]))
        _safe("fig12", lambda: fig12_p1_curves(store["p1"][0]))
    _safe("fig13/14", lambda: store.setdefault(
        "p2", run_p2_snapshots()) and fig13_14_p2(*store["p2"]))
    _safe("fig15-18", lambda: store.setdefault(
        "p3", run_p3_full()) and fig15_18_p3(store["p3"]))
    _safe("fig19-21", lambda: store.setdefault(
        "p4", _run_with_retry(lambda dt: dc.simulate(
            400, 400000.0, dt, dc.props_appendix4, sc.load_env(),
            sc.load_radius("linear"), theta=1.0, c_stop=dc.C_DRY,
            record_dt=60.0,
            snapshot_times=list(np.arange(1800.0, 220000.0, 1800.0)))))
        and fig19_21_p4(store["p4"], store["p3"]))
    _safe("fig22", fig22_analytic)
    _safe("fig23", fig23_convergence)
    _safe("fig24", fig24_residuals)
    _safe("fig25", fig25_inverse)
    _safe("fig26-28", fig26_28_sensitivity)
    _safe("fig29", fig29_drymass)
    _safe("fig30", fig30_latent)
    n = len([f for f in os.listdir(FIG) if f.endswith(".png")])
    print(f"figures 目录现有 {n} 张 PNG")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
