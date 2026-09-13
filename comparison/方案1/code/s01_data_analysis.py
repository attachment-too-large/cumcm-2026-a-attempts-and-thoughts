# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s01_data_analysis.py
# 作用  : 第 2 步「数据分析」。读取附件 1(烘房温湿度)与附件 2(药材半径),
#         完成数据规模/缺失/异常/分段/平台/单调性检查, 并输出:
#           data/data_overview.txt            数据概况报告
#           data/attachment1_clean.csv        附件 1 原始表(带派生列)
#           data/phase_summary.csv            分段统计(升温段 / 平台段)
#           figures/fig01_room_environment.png
#           figures/fig02_radius_history.png
#           figures/fig03_room_increments.png
# 运行  : python s01_data_analysis.py
# 说明  : 本脚本只做"用数据反推模型结构"的诊断, 不参与求解; 所有结论都以
#         打印与落盘数字为准, 不使用目测值。
# =============================================================================
"""Step 2: exploratory data analysis of the two attachments."""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotstyle  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_IN = os.path.join(BASE, "data", "input")
DIR_DATA = os.path.join(BASE, "data")
DIR_FIG = os.path.join(BASE, "figures")

R0_CM = 2.0
KELVIN = 273.15
P_ATM = 101325.0            # 标准大气压 [Pa]


def load_attachment1():
    """读取附件 1: 时间 [s]、温度 [degC]、水分浓度 [kg/kg]。"""
    df = pd.read_excel(os.path.join(DIR_IN, "attachment1_room.xlsx"))
    df.columns = ["t", "T_air", "C_air"]
    return df


def load_attachment2():
    """读取附件 2: 时间 [s]、半径 [cm]。"""
    df = pd.read_excel(os.path.join(DIR_IN, "attachment2_radius.xlsx"))
    df.columns = ["t", "R_cm"]
    return df


def saturation_pressure(T_C):
    """饱和水蒸气压 [Pa], Magnus 公式(Alduchov & Eskridge, 1996 形式)。

    适用范围: -40 ~ 50 degC, 常压空气; 50 degC 处不确定度约 0.1%。
    """
    return 611.21 * np.exp((18.678 - T_C / 234.5) * T_C / (257.14 + T_C))


def humidity_diagnostics(t, T_air, C_air, label, out):
    """把附件 1 的"水分浓度"分别按湿度比与绝对湿度两种口径换成相对湿度。

    目的: 判断该列更可能是 kg/kg(干空气) 还是 kg/m3; 两种口径下平台段相对
    湿度都应落在 50%~70%, 若某一口径给出 >100% 或 <10% 则应排除。
    """
    p_sat = saturation_pressure(T_air)
    # 口径 A: 湿度比 w [kg/kg 干空气] -> 水蒸气分压 -> 相对湿度
    p_v_A = P_ATM * C_air / (0.622 + C_air)
    rh_A = 100.0 * p_v_A / p_sat
    # 口径 B: 绝对湿度 [kg/m3] -> 理想气体 -> 水蒸气分压
    R_v = 461.5
    p_v_B = C_air * R_v * (T_air + KELVIN)
    rh_B = 100.0 * p_v_B / p_sat
    out.append(f"[{label}] 首点 T={T_air[0]:.3f} degC, C_air={C_air[0]:.5f}")
    out.append(f"[{label}] 末点 T={T_air[-1]:.3f} degC, C_air={C_air[-1]:.5f}")
    out.append(f"[{label}] 口径 A(湿度比 kg/kg) 相对湿度: 首 {rh_A[0]:.1f}%, 末 {rh_A[-1]:.1f}%")
    out.append(f"[{label}] 口径 B(绝对湿度 kg/m3) 相对湿度: 首 {rh_B[0]:.1f}%, 末 {rh_B[-1]:.1f}%")
    return rh_A, rh_B


def main():
    plotstyle.setup()
    lines = []
    push = lines.append

    # ---------------------------------------------------------------- 附件 1
    df1 = load_attachment1()
    t = df1["t"].to_numpy(float)
    Ta = df1["T_air"].to_numpy(float)
    Ca = df1["C_air"].to_numpy(float)

    push("=" * 78)
    push("附件 1(烘房环境)数据概况")
    push("=" * 78)
    push(f"记录条数 = {len(df1)}; 字段 = {list(df1.columns)}")
    push(f"时间范围 = {t[0]:.0f} ~ {t[-1]:.0f} s ({(t[-1]-t[0])/3600:.3f} h)")
    dt = np.diff(t)
    push(f"时间步长: 唯一值 = {np.unique(dt)}, 等间隔 = {bool(np.all(dt == dt[0]))}")
    push(f"缺失值: T_air {int(df1['T_air'].isna().sum())} 个, C_air {int(df1['C_air'].isna().sum())} 个")
    push(f"温度范围   : {Ta.min():.4f} ~ {Ta.max():.4f} degC")
    push(f"水分浓度范围: {Ca.min():.6f} ~ {Ca.max():.6f} kg/kg")

    # 异常值: 用相邻二阶差分做尖峰检测(对 60 s 采样, 二阶差分应远小于量程)
    d2T = np.abs(np.diff(Ta, 2))
    d2C = np.abs(np.diff(Ca, 2))
    push(f"二阶差分最大值: T_air {d2T.max():.5f} degC, C_air {d2C.max():.6f} kg/kg")
    push(f"二阶差分 > 5x 中位数 的点数: T_air {int(np.sum(d2T > 5*np.median(d2T)))} 个, "
         f"C_air {int(np.sum(d2C > 5*np.median(d2C)))} 个")

    # 分段: 以 600 s 滑窗内的平均斜率为判据, 划分升温段与平台段
    win = 10
    slope_T = np.array([np.polyfit(t[i:i + win], Ta[i:i + win], 1)[0]
                        for i in range(len(t) - win + 1)])
    slope_C = np.array([np.polyfit(t[i:i + win], Ca[i:i + win], 1)[0]
                        for i in range(len(t) - win + 1)])
    thr_T = 0.20 * np.max(np.abs(slope_T))
    plateau_idx = np.where(np.abs(slope_T) < thr_T)[0] * 1 + win // 2
    t_plateau = float(t[plateau_idx[0]]) if plateau_idx.size else float(t[-1])
    push("-" * 78)
    push(f"分段判据: 600 s 滑窗线性斜率 |dT/dt| < {thr_T:.4f} degC/s 视为平台段")
    push(f"升温段(预热平衡阶段): 0 ~ {t_plateau:.0f} s")
    push(f"平台段(恒温干燥阶段): {t_plateau:.0f} ~ {t[-1]:.0f} s")

    # 平台取值的多种口径(平台假定是模型假设, 不是数据的唯一推论)
    plateau_variants = {
        "last_point": (Ta[-1], Ca[-1]),
        "tail_mean_1800s(30点)": (Ta[-30:].mean(), Ca[-30:].mean()),
        "tail_mean_3600s(60点)": (Ta[-60:].mean(), Ca[-60:].mean()),
        "tail_mean_7200s(120点)": (Ta[-120:].mean(), Ca[-120:].mean()),
    }
    push("-" * 78)
    push("平台段环境取值的多种口径(用于灵敏度分析):")
    for k, (tt, cc) in plateau_variants.items():
        push(f"  {k:26s}: T_air = {tt:.4f} degC, C_air = {cc:.6f} kg/kg")
    push(f"  平台段标准差          : T_air = {Ta[-60:].std(ddof=1):.4f} degC, "
         f"C_air = {Ca[-60:].std(ddof=1):.6f} kg/kg")

    rh_A, rh_B = humidity_diagnostics(t, Ta, Ca, "附件1", lines)

    # ---------------------------------------------------------------- 附件 2
    df2 = load_attachment2()
    tr = df2["t"].to_numpy(float)
    Rcm = df2["R_cm"].to_numpy(float)
    push("")
    push("=" * 78)
    push("附件 2(药材半径)数据概况")
    push("=" * 78)
    push(f"记录条数 = {len(df2)}; 时间范围 = {tr[0]:.0f} ~ {tr[-1]:.0f} s ({tr[-1]/3600:.3f} h)")
    dtr = np.diff(tr)
    push(f"时间步长: 唯一值 = {np.unique(dtr)}, 等间隔 = {bool(np.all(dtr == dtr[0]))}")
    push(f"半径范围 = {Rcm.min():.4f} ~ {Rcm.max():.4f} cm")
    push(f"R(0) = {Rcm[0]:.4f} cm (题面给定初始半径 2 cm, 相对偏差 "
         f"{100*(Rcm[0]-R0_CM)/R0_CM:+.2f}%)")
    push(f"R(末) = {Rcm[-1]:.4f} cm; 总收缩率 = {100*(1-Rcm[-1]/Rcm[0]):.2f}%")
    dR = np.diff(Rcm)
    push(f"单调不增 = {bool(np.all(dR <= 1e-12))}; 严格递减段数 = {int(np.sum(dR < -1e-12))}")
    push(f"|dR/dt| 最大 = {np.max(np.abs(dR/dtr)):.3e} cm/s, 中位 = "
         f"{np.median(np.abs(dR/dtr)):.3e} cm/s")
    push(f"半径二阶差分最大绝对值 = {np.max(np.abs(np.diff(Rcm, 2))):.5f} cm "
         "(若偏大说明半径历史存在折点, 需检查分段线性插值的适用性)")
    # 体积比与干物质守恒一致性检查
    vol_ratio = (Rcm / Rcm[0]) ** 2
    push(f"体积比 V/V0 = (R/R0)^2: 末端 {vol_ratio[-1]:.4f}, 即体积收缩 "
         f"{100*(1-vol_ratio[-1]):.2f}%")
    push("注: 半径历史为实测输入, 与附录 3/4 的 rho(C) 经验式是否自洽需另做检验"
         "(见 s05_verification.py 的干物质守恒检验)。")

    # ------------------------------------------------------------ 落盘与出图
    df1_out = df1.copy()
    df1_out["dT_dt"] = np.gradient(Ta, t)
    df1_out["dC_dt"] = np.gradient(Ca, t)
    df1_out["RH_percent(kg/kg dry air)"] = rh_A
    df1_out["RH_percent(kg/m3)"] = rh_B
    df1_out.to_csv(os.path.join(DIR_DATA, "attachment1_clean.csv"),
                   index=False, encoding="utf-8-sig")

    df2_out = df2.copy()
    df2_out["volume_ratio"] = vol_ratio
    df2_out.to_csv(os.path.join(DIR_DATA, "attachment2_clean.csv"),
                   index=False, encoding="utf-8-sig")

    pd.DataFrame([
        {"phase": "preheat", "t_start_s": 0.0, "t_end_s": t_plateau,
         "T_air_start": Ta[0], "T_air_end": float(np.interp(t_plateau, t, Ta)),
         "C_air_start": Ca[0], "C_air_end": float(np.interp(t_plateau, t, Ca))},
        {"phase": "constant_drying", "t_start_s": t_plateau, "t_end_s": t[-1],
         "T_air_start": float(np.interp(t_plateau, t, Ta)), "T_air_end": Ta[-1],
         "C_air_start": float(np.interp(t_plateau, t, Ca)), "C_air_end": Ca[-1]},
    ]).to_csv(os.path.join(DIR_DATA, "phase_summary.csv"),
              index=False, encoding="utf-8-sig")

    import matplotlib.pyplot as plt

    # 图 1: 烘房环境
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.4), sharex=True)
    axes[0].plot(t / 3600.0, Ta, color=plotstyle.PALETTE["red"], lw=1.6, label="烘房温度 $T_{air}$")
    axes[0].axvline(t_plateau / 3600.0, color=plotstyle.PALETTE["gray"], ls=":", lw=1.3)
    axes[0].axhline(Ta[-60:].mean(), color=plotstyle.PALETTE["blue"], ls="--", lw=1.1,
                    label="末段 3600 s 均值 %.3f $^\\circ$C" % Ta[-60:].mean())
    axes[0].set_ylabel("温度 / $^\\circ$C")
    axes[0].legend(loc="lower right")
    axes[0].set_title("附件 1: 烘房温湿度环境与工艺分段(虚线为预热/恒温分界)")
    axes[1].plot(t / 3600.0, Ca, color=plotstyle.PALETTE["blue"], lw=1.6,
                 label="烘房水分浓度 $C_{air}$")
    axes[1].axvline(t_plateau / 3600.0, color=plotstyle.PALETTE["gray"], ls=":", lw=1.3)
    axes[1].axhline(Ca[-60:].mean(), color=plotstyle.PALETTE["red"], ls="--", lw=1.1,
                    label="末段 3600 s 均值 %.5f" % Ca[-60:].mean())
    axes[1].set_xlabel("时间 / h")
    axes[1].set_ylabel("水分浓度 / (kg/kg)")
    axes[1].legend(loc="lower right")
    for ax in axes:
        ax.axvspan(0, t_plateau / 3600.0, color=plotstyle.PALETTE["orange"], alpha=0.07)
    fig.tight_layout()
    plotstyle.save(fig, os.path.join(DIR_FIG, "fig01_room_environment.png"))

    # 图 2: 半径历史
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.3))
    axes[0].plot(tr / 3600.0, Rcm, color=plotstyle.PALETTE["green"], lw=1.8, marker="o", ms=3)
    axes[0].set_xlabel("时间 / h")
    axes[0].set_ylabel("半径 / cm")
    axes[0].set_title("附件 2: 药材半径历史")
    axes[1].plot(tr / 3600.0, vol_ratio, color=plotstyle.PALETTE["purple"], lw=1.8)
    axes[1].set_xlabel("时间 / h")
    axes[1].set_ylabel("$V/V_0=(R/R_0)^2$")
    axes[1].set_title("相对体积收缩")
    axes[2].plot(tr[1:] / 3600.0, np.diff(Rcm) / np.diff(tr) * 3600.0,
                 color=plotstyle.PALETTE["orange"], lw=1.6)
    axes[2].set_xlabel("时间 / h")
    axes[2].set_ylabel("d$R$/d$t$ / (cm/h)")
    axes[2].set_title("收缩速率(递减说明收缩逐渐停止)")
    fig.tight_layout()
    plotstyle.save(fig, os.path.join(DIR_FIG, "fig02_radius_history.png"))

    # 图 3: 环境增量诊断(判断平台段是否真实平稳、有无噪声)
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.2))
    axes[0].plot(t[1:] / 3600.0, np.diff(Ta), color=plotstyle.PALETTE["red"], lw=1.2)
    axes[0].set_xlabel("时间 / h")
    axes[0].set_ylabel("$\\Delta T_{air}$ / $^\\circ$C")
    axes[0].set_title("相邻采样温度增量(60 s)")
    axes[1].plot(t[1:] / 3600.0, np.diff(Ca), color=plotstyle.PALETTE["blue"], lw=1.2)
    axes[1].set_xlabel("时间 / h")
    axes[1].set_ylabel("$\\Delta C_{air}$ / (kg/kg)")
    axes[1].set_title("相邻采样水分浓度增量(60 s)")
    axes[2].hist(Ta[-60:], bins=12, color=plotstyle.PALETTE["gray"], alpha=0.85)
    axes[2].set_xlabel("平台段 $T_{air}$ / $^\\circ$C")
    axes[2].set_ylabel("频数")
    axes[2].set_title("平台段温度分布(末 3600 s)")
    fig.tight_layout()
    plotstyle.save(fig, os.path.join(DIR_FIG, "fig03_room_increments.png"))

    # 图 4: 湿度口径检验
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    ax.plot(t / 3600.0, rh_A, color=plotstyle.PALETTE["blue"], lw=1.7,
            label="口径 A: 视为湿度比 $w$ (kg/kg 干空气)")
    ax.plot(t / 3600.0, rh_B, color=plotstyle.PALETTE["orange"], lw=1.7, ls="--",
            label="口径 B: 视为绝对湿度 (kg/m$^3$)")
    ax.axhline(100.0, color=plotstyle.PALETTE["red"], ls=":", lw=1.2)
    ax.set_xlabel("时间 / h")
    ax.set_ylabel("推算相对湿度 / %")
    ax.set_title("附件 1「水分浓度」两种物理口径下的相对湿度")
    ax.legend(loc="lower right")
    fig.tight_layout()
    plotstyle.save(fig, os.path.join(DIR_FIG, "fig04_humidity_interpretation.png"))

    text = "\n".join(lines)
    with open(os.path.join(DIR_DATA, "data_overview.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print(text)
    print("\n[OK] 数据分析完成, 报告: data/data_overview.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
