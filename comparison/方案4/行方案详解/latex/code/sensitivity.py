# -*- coding: utf-8 -*-
# ============================================================================
# sensitivity.py —— 干燥时间对参数与模型假设的灵敏度（单因子扰动）
#
#   (1) 物性与传递参数：D 的前置因子、活化温度 3850 K、浓度指数、rho、cp、k，
#       以及 h、h_m，做单因子扰动（±1%、±2%；传递系数取 ±5%、±10%）；
#   (2) 烘房条件：设定温度、设定含湿量、预热/恒温切换时刻、预热时间常数，以及
#       用附件 1 的原始含噪数据代替平滑拟合曲线；
#   (3) 模型形式：表面蒸发潜热冷却、问题 4 是否计入收缩、是否计入物质坐标下
#       能量方程的体积变化项、收缩曲线用实测表还是指数拟合；
#   (4) 汇总排序：按影响大小排出各因素，指出主导不确定性的参数。
#
# 所有结果都以「与同一网格基准的差值」给出：扰动算例与基准用完全相同的网格与
# 步长，相减后离散误差被抵消，剩下的差异才归因于被扰动的那个量。
#
# 输入数据：物性来自附录 3 / 4；环境历程与含噪原始数据来自附件 1，收缩半径
#           来自附件 2（均由 common.py 读取）。
# 输出：data/sensitivity.csv（每个算例一行）、data/sensitivity_summary.txt
#       （控制台报告全文）。
#
# 关键变量：
#   C_TARGET = 0.15 kg/kg 干燥终点判据；N_CELL = 400 基准网格数
#   SCHED 时间步安排（前 7200 s 用 1 s，之后 10 s）；base 基准干燥时间 [s]；
#   delta = 扰动算例 t_dry - 基准 t_dry [s]
# ============================================================================
"""
运行方式：在目录 A题解答 下执行 python code/sensitivity.py（约需数分钟）。

报告按四组逐行打印每个算例的 t_dry 与相对基准的差值，最后一段是按 |delta|
排序的主导因素表。
"""

import os
import sys
import io
import time
import numpy as np
import pandas as pd

if sys.stdout.encoding is None or "utf" not in sys.stdout.encoding.lower():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm

C_TARGET = 0.15
T_END = 2.6e5
N_CELL = 400
SCHED = [(7200.0, 1.0), (T_END, 10.0)]
T_REP = np.arange(1.0, T_END + 0.5, 30.0)
XI_REP = np.array([0.0, 0.25, 0.5, 0.75, 0.9, 1.0])


# --------------------------------------------------------------------------
# 把附录 3 / 4 的物性式改写成「常数全部暴露为关键字参数」的形式，
# 扰动哪一个就传哪一个，其余保持默认值，从而保证每次只动一个因子
# --------------------------------------------------------------------------
def props23_scaled(A=2.4e-3, Ea=3850.0, b=0.45, rho0=650.0, rho1=128.0,
                   cp0=1450.0, cp1=2736.0, k0=0.21, k1=0.38):
    """附录 3 的物性关系，常数暴露为关键字参数以便逐个扰动。

    A 为 D 的前置因子 [m2/s]，Ea 为活化温度 [K]，b 为浓度指数（无量纲），
    rho0/rho1、cp0/cp1、k0/k1 分别是 rho、cp、k 随含水率变化式的两个系数。
    返回物性函数 prop(C, T_K) -> (rho, cp, k, D)，签名与 common.py 中一致。
    """
    def prop(C, T_K):
        C = np.maximum(np.asarray(C, dtype=float), cm.C_FLOOR)
        f = C / (C + 1.0)
        return (rho0 + rho1 * C, cp0 + cp1 * f, k0 + k1 * f,
                A * np.exp(-b / C) * np.exp(-Ea / T_K))
    return prop


def props4_scaled(A=4.2e-4, Ea=3850.0, b=0.30, rho0=760.0, rho1=90.0,
                  cp0=1850.0, cp1=2150.0, k0=0.12, k1=0.20):
    """附录 4（收缩工况）的物性关系，常数暴露为关键字参数以便逐个扰动。

    参数含义与 props23_scaled 完全相同，只是默认取值换成附录 4 的系数。
    返回物性函数 prop(C, T_K) -> (rho, cp, k, D)。
    """
    def prop(C, T_K):
        C = np.maximum(np.asarray(C, dtype=float), cm.C_FLOOR)
        f = C / (C + 1.0)
        return (rho0 + rho1 * C, cp0 + cp1 * f, k0 + k1 * f,
                A * np.exp(-b / C) * np.exp(-Ea / T_K))
    return prop


def drying_time(prop, radius=None, room=None, grid="node", n_cell=N_CELL,
                h_conv=None, hm_conv=None, volume_term=False,
                evap_cooling=False):
    """用统一的网格与步长安排算一次干燥时间 t_dry [s]，供各扰动算例调用。

    room、radius 为环境与收缩历史（None 表示用默认值），h_conv/hm_conv 可覆盖
    对流换热与传质系数，volume_term/evap_cooling 打开模型形式的两个可选项。
    这里用 track_max=True 让 simulate 额外记录整场的最大含水率 Cmax（而不是只在
    少数采样点上取最大），再线性插值出 Cmax 首次降到 C_TARGET 的时刻。
    返回 (t_dry, 求解结果)；若 T_END 内没降到判据以下，t_dry 为 nan。
    """
    room = room if room is not None else cm.RoomConditions()
    res = cm.simulate(n_cell, prop, T_END, 1.0, room, radius=radius, grid=grid,
                      h_conv=h_conv, hm_conv=hm_conv, volume_term=volume_term,
                      evap_cooling=evap_cooling, dt_schedule=SCHED,
                      t_report=T_REP, xi_report=XI_REP, track_max=True)
    Cmax = res["Cmax"]                     # 整场最大值（含未输出的内部节点）
    idx = np.where(Cmax < C_TARGET)[0]
    if not len(idx):
        return np.nan, res
    i0 = idx[0]
    t0, t1 = res["t"][i0 - 1], res["t"][i0]
    v0, v1 = Cmax[i0 - 1], Cmax[i0]
    return t0 + (C_TARGET - v0) * (t1 - t0) / (v1 - v0), res


class MeasuredRoom:
    """直接用附件 1 的原始实测序列当环境，接口与 RoomConditions 相同。

    一阶惯性式拟合会抹掉测量噪声，本类保留噪声，用于量化「用平滑曲线近似实测
    烘房历程」对 t_dry 的影响。实测只记到 14400 s，之后保持最后一个值不变。
    """

    def __init__(self):
        df = cm.load_attachment1()
        self.t = df["t"].to_numpy(float)
        self.T = df["T_air"].to_numpy(float)
        self.C = df["C_air"].to_numpy(float)

    def T_air(self, t):
        t = np.asarray(t, dtype=float)
        v = np.interp(t, self.t, self.T)
        # 超过实测末时刻就维持末值（对应恒温干燥阶段）
        return np.where(t > self.t[-1], self.T[-1], v) if v.ndim else float(v)

    def C_air(self, t):
        t = np.asarray(t, dtype=float)
        v = np.interp(t, self.t, self.C)
        return np.where(t > self.t[-1], self.C[-1], v) if v.ndim else float(v)


def main():
    """跑完四组灵敏度算例：基准、参数扰动、环境与模型形式、主导因素排序。

    P(s) 既打印（带 flush，便于看后台进度）又把同一行收进 out，最后一次性写文件。
    """
    rad = cm.RadiusHistory(mode="table")
    out = []
    P = lambda s: (print(s, flush=True), out.append(s))
    rows = []

    P("=" * 84)
    P("SENSITIVITY OF THE DRYING TIME  (n=%d, dt = 1 s -> 10 s, t_end = %.0f s)"
      % (N_CELL, T_END))
    P("=" * 84)

    # 基准：用附录 3 / 4 系数原值各算一次，后面所有差值都以它为参照
    base = {}
    for tag, prop, radius in [("P3", props23_scaled(), None),
                              ("P4", props4_scaled(), rad)]:
        t0 = time.time()
        td, _ = drying_time(prop, radius=radius)
        base[tag] = td
        P("  baseline %s : t_dry = %.1f s = %.4f h   (%.1f s wall)"
          % (tag, td, td / 3600.0, time.time() - t0))
    P("")
    P("  %-6s %-28s %-10s %-14s %-12s" %
      ("case", "perturbed quantity", "change", "t_dry / h", "delta / h"))

    def run(tag, label, change, prop, radius, **kw):
        """跑一个扰动算例：打印 t_dry 与相对基准的差值，记入 rows 并返回该差值 [s]。"""
        td, _ = drying_time(prop, radius=radius, **kw)
        d = td - base[tag]
        P("  %-6s %-28s %-10s %-14.4f %+12.4f" %
          (tag, label, change, td / 3600.0, d / 3600.0))
        rows.append(dict(case=tag, quantity=label, change=change,
                         t_dry_h=td / 3600.0, delta_h=d / 3600.0))
        return d

    # ---- 1. 物性与传递参数：每次只扰动一个量（factory 只传被扰动的关键字） ----
    for tag, radius, factory in [("P3", None, props23_scaled),
                                 ("P4", rad, props4_scaled)]:
        for pct in (-0.02, -0.01, 0.01, 0.02):
            run(tag, "D0 pre-exponential", "%+.0f%%" % (100 * pct),
                factory(A=(2.4e-3 if tag == "P3" else 4.2e-4) * (1 + pct)),
                radius=radius)
        for pct in (-0.02, -0.01, 0.01, 0.02):
            run(tag, "activation T (3850 K)", "%+.0f%%" % (100 * pct),
                factory(Ea=3850.0 * (1 + pct)), radius=radius)
        for pct in (-0.02, -0.01, 0.01, 0.02):
            run(tag, "concentration exponent", "%+.0f%%" % (100 * pct),
                factory(b=(0.45 if tag == "P3" else 0.30) * (1 + pct)),
                radius=radius)
        # 密度、比热、导热系数只取 ±2% 两档，每个量各扰动一次
        for pct in (-0.02, 0.02):
            run(tag, "rho correlation", "%+.0f%%" % (100 * pct),
                factory(rho0=(650.0 if tag == "P3" else 760.0) * (1 + pct),
                        rho1=(128.0 if tag == "P3" else 90.0) * (1 + pct)),
                radius=radius)
            run(tag, "cp correlation", "%+.0f%%" % (100 * pct),
                factory(cp0=(1450.0 if tag == "P3" else 1850.0) * (1 + pct),
                        cp1=(2736.0 if tag == "P3" else 2150.0) * (1 + pct)),
                radius=radius)
            run(tag, "k correlation", "%+.0f%%" % (100 * pct),
                factory(k0=(0.21 if tag == "P3" else 0.12) * (1 + pct),
                        k1=(0.38 if tag == "P3" else 0.20) * (1 + pct)),
                radius=radius)
        # 传递系数按更大的幅度扰动：h_m 取 ±5%、±10% 四档
        for pct in (-0.10, -0.05, 0.05, 0.10):
            run(tag, "h_m (mass transfer)", "%+.0f%%" % (100 * pct),
                factory(), radius=radius, hm_conv=cm.HM_CONV * (1 + pct))
        # 换热系数只取 ±10% 两档
        for pct in (-0.10, 0.10):
            run(tag, "h (heat transfer)", "%+.0f%%" % (100 * pct),
                factory(), radius=radius, h_conv=cm.H_CONV * (1 + pct))

    # ---- 2. 烘房条件：设定值、切换时刻、预热时间常数、原始实测数据 ----
    P("")
    for tag, radius, factory in [("P3", None, props23_scaled),
                                 ("P4", rad, props4_scaled)]:
        for dT in (-0.5, -0.2, 0.2, 0.5):
            run(tag, "set-point temperature", "%+.1f K" % dT, factory(),
                radius=radius, room=cm.RoomConditions(t_set=cm.T_SET_C + dT))
        for dc in (-0.02, -0.01, 0.01, 0.02):
            run(tag, "set-point moisture", "%+.0f%%" % (100 * dc), factory(),
                radius=radius,
                room=cm.RoomConditions(c_set=cm.C_SET * (1 + dc)))
        # 预热段与恒温段的切换时刻：2 h / 6 h / 8 h（拟合值 4 h 即基准）
        for tsw in (7200.0, 21600.0, 28800.0):
            run(tag, "preheat->constant switch", "%.0f h" % (tsw / 3600.0),
                factory(), radius=radius, room=cm.RoomConditions(t_switch=tsw))
        # 温度与含湿量的一阶惯性时间常数同比例缩放 0.8 / 1.2 倍
        for tau in (0.8, 1.2):
            run(tag, "preheat time constant", "x%.1f" % tau, factory(),
                radius=radius,
                room=cm.RoomConditions(tau_t=cm.TAU_AIR_T * tau,
                                       tau_c=cm.TAU_AIR_C * tau))
        # 时间常数取 1e-6 s 相当于取消预热段：t>0 时环境立刻达到设定值，
        # 而初值仍是 28 degC，用于对比「预热阶段」本身对总干燥时间的贡献
        run(tag, "no preheat ramp (T_air=28)", "alt", factory(),
            radius=radius,
            room=cm.RoomConditions(t_set=cm.T_SET_C, c_set=cm.C_SET,
                                   tau_t=1e-6, tau_c=1e-6))
        run(tag, "raw measured chamber data", "noise", factory(),
            radius=radius, room=MeasuredRoom())

    # ---- 3. 模型形式：潜热冷却、收缩、能量方程体积项、R(t) 取法 ----
    P("")
    for tag, radius, factory in [("P3", None, props23_scaled),
                                 ("P4", rad, props4_scaled)]:
        run(tag, "evaporative cooling ON", "latent", factory(), radius=radius,
            evap_cooling=True)
    # radius=None 表示半径恒为 R0 = 2 cm，即不计收缩
    run("P4", "shrinkage OFF (R = 2 cm)", "alt", props4_scaled(),
        radius=None)
    # 物质坐标下对流项本来为零，这里额外保留 dR/dt 引起的体积变化项
    run("P4", "volume term in energy eq.", "alt", props4_scaled(), radius=rad,
        volume_term=True)
    # 收缩曲线改用单指数拟合的平滑 R(t) 代替分段线性插值的实测表
    rad_s = cm.RadiusHistory(mode="smooth")
    run("P4", "smoothed R(t) instead of table", "alt", props4_scaled(),
        radius=rad_s)

    # 逐算例结果落盘为 CSV，供后面的出图脚本复用
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(cm.DATA_DIR, "sensitivity.csv"), index=False)
    P("")
    P("  baseline t_dry : P3 = %.4f h, P4 = %.4f h"
      % (base["P3"] / 3600.0, base["P4"] / 3600.0))

    # ---- 4. 主导因素排序：只比同一扰动幅度（+2%）下的 |delta|，从大到小排 ----
    P("")
    P("=" * 84)
    P("RANKING OF THE DOMINANT UNCERTAINTIES (|delta t_dry| for a 2% change)")
    P("=" * 84)
    for tag in ("P3", "P4"):
        sub = df[(df.case == tag) & (df.change == "+2%")]
        if len(sub):
            sub = sub.reindex(sub.delta_h.abs().sort_values(ascending=False).index)
            for _, r in sub.iterrows():
                P("  %-4s %-28s  %+8.4f h  (%+6.2f %%)"
                  % (tag, r["quantity"], r["delta_h"],
                     100.0 * r["delta_h"] / (base[tag] / 3600.0)))

    with open(os.path.join(cm.DATA_DIR, "sensitivity_summary.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("\n   saved ->", os.path.join(cm.DATA_DIR, "sensitivity_summary.txt"))


if __name__ == "__main__":
    main()
