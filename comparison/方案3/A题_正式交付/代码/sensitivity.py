# -*- coding: utf-8 -*-
"""灵敏度分析（真实运行，非口头推断）。

扰动对象：活化温度 E=3850 K、活化浓度指数 a、扩散系数前因子 D0、
对流换热系数 h、对流传质系数 hm、烘房设定温度 T_set、以及环境口径（拟合 vs 原始插值）。
考察指标：问题 3 的烘干时长 t_dry，以及 t=24 h / 48 h 的中心含水率。

输出：<交付根>/中间数据/sensitivity.txt 与 sensitivity.npz
"""
import os
from console_utf8 import fix_console
fix_console()
import numpy as np
import common as cm
from spectral import HerbSolver
import fig_core as FC

CRIT = 0.15
C_FLOOR = 1e-3
OUT = os.path.join(FC.DATADIR, "sensitivity.txt")
NPZ = os.path.join(FC.DATADIR, "sensitivity.npz")
LOG = []


def log(s):
    print(s, flush=True)
    LOG.append(str(s))


def props_pert(prob, C, T, pert):
    """按扰动因子重算物性（与 common.props 结构一致）。"""
    C = np.asarray(C, dtype=float)
    T = np.asarray(T, dtype=float)
    Cf = np.maximum(C, C_FLOOR)
    TK = T + 273.15
    E = 3850.0 * pert.get("E", 1.0)
    fD = pert.get("D0", 1.0)
    fa = pert.get("a", 1.0)
    if prob == 1:
        rho = np.full_like(C, 820.0)
        cp = np.full_like(C, 2600.0)
        k = np.full_like(C, 0.36)
        D = 7e-9 * fD * np.exp(-0.89 * fa / Cf)
    elif prob in (2, 3):
        rho = 650.0 + 128.0 * C
        cp = 1450.0 + 2736.0 * C / (C + 1.0)
        k = 0.21 + 0.38 * C / (C + 1.0)
        D = 2.4e-3 * fD * np.exp(-0.45 * fa / Cf) * np.exp(-E / TK)
    else:
        rho = 760.0 + 90.0 * C
        cp = 1850.0 + 2150.0 * C / (C + 1.0)
        k = 0.12 + 0.20 * C / (C + 1.0)
        D = 4.2e-4 * fD * np.exp(-0.30 * fa / Cf) * np.exp(-E / TK)
    fh = pert.get("h", 1.0)
    fhm = pert.get("hm", 1.0)
    return rho, cp, k * (1.0 if prob == 1 else 1.0) * (1.0), D, fh, fhm


def run_case(tag, pert=None, Tair_mode="fit", N=48, prob=3):
    """在给定扰动下算 t_dry 与两个时刻的中心含水率。"""
    import spectral
    pert = pert or {}
    base_props = cm.props
    base_T, base_C = cm.T_air, cm.C_air
    base_h, base_hm = cm.h, cm.hm

    def _props(p, C, T):
        rho, cp, k, D, fh, fhm = props_pert(p, C, T, pert)
        return rho, cp, k, D

    cm.props = _props
    cm.h = base_h * pert.get("h", 1.0)
    cm.hm = base_hm * pert.get("hm", 1.0)
    if "Tset" in pert:
        dT = pert["Tset"]
        cm.T_air = lambda t: np.where(np.asarray(t) <= cm.T_END,
                                      base_T(t) + dT, cm.Tset + dT)

    try:
        s = spectral.HerbSolver(N, prob=prob)
        at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])

        def ev(t, y):
            return np.max(y[:s.Np]) - CRIT
        ev.terminal = True
        ev.direction = -1
        e = s.solve(320000.0, rtol=1e-10, atol=at, events=ev)
        td = float(e.t_events[0][0])
        # 固定时刻的中心含水率（只取落在积分区间内的时刻）
        targets = [24 * 3600.0, 48 * 3600.0]
        tmax = min(td, 48 * 3600.0)
        tv = [x for x in targets if x <= tmax + 1e-9]
        c24 = c48 = np.nan
        if tv:
            so = s.solve(tmax, t_eval=tv, rtol=1e-10, atol=at)
            for k, x in enumerate(tv):
                v = float(s.interp_u(so.y[:, k], np.array([0.0]))[0][0])
                if abs(x - targets[0]) < 1.0:
                    c24 = v
                if abs(x - targets[1]) < 1.0:
                    c48 = v
    finally:
        cm.props = base_props
        cm.T_air, cm.C_air = base_T, base_C
        cm.h, cm.hm = base_h, base_hm
    return td, c24, c48


def main():
    log("=" * 74)
    log("灵敏度分析（问题 3，烘干时长 t_dry；基准 t_dry 由 BDF 求解给出）")
    log("=" * 74)

    base = run_case("base")
    log("基准（无扰动）：t_dry = %.1f s = %.4f h ; C0(24h)=%.4f ; C0(48h)=%.4f"
        % (base[0], base[0] / 3600, base[1], base[2]))
    t0 = base[0]

    cases = [
        ("E   +2%  (3850 -> 3927 K)",      {"E": 1.02}),
        ("E   -2%  (3850 -> 3773 K)",      {"E": 0.98}),
        ("a   +2%  (0.45 -> 0.459)",       {"a": 1.02}),
        ("a   -2%  (0.45 -> 0.441)",       {"a": 0.98}),
        ("D0  +2%",                        {"D0": 1.02}),
        ("D0  -2%",                        {"D0": 0.98}),
        ("h   +2%  (25 -> 25.5 W/m2K)",    {"h": 1.02}),
        ("h   -2%",                        {"h": 0.98}),
        ("hm  +2%  (8e-7 -> 8.16e-7)",     {"hm": 1.02}),
        ("hm  -2%",                        {"hm": 0.98}),
        ("Tset +0.2 K (50.212 -> 50.412)", {"Tset": 0.2}),
        ("Tset -0.2 K",                    {"Tset": -0.2}),
    ]
    rows = []
    log("")
    log("%-34s %12s %10s %10s" % ("扰动", "t_dry / s", "相对变化", "绝对变化/s"))
    log("-" * 74)
    for name, pert in cases:
        td, c24, c48 = run_case(name, pert)
        rel = (td - t0) / t0 * 100
        rows.append((name, td, rel, td - t0, c24, c48))
        log("%-34s %12.1f %9.3f%% %10.1f" % (name, td, rel, td - t0))
        print("   C0(24h)=%.4f  C0(48h)=%.4f" % (c24, c48))

    # ---- 环境口径分叉 ----
    log("")
    log("环境口径分叉（原始序列插值 vs 一阶惯性拟合）：")
    base_Tair = cm.T_air

    def T_raw(t):
        t = np.asarray(t, dtype=float)
        return np.where(t <= cm.T_END, np.interp(t, cm._t_env, cm._Ta_env), cm.Tset)

    try:
        cm.T_air = T_raw
        td_raw, c24r, c48r = run_case("raw")
    finally:
        cm.T_air = base_Tair
    log("  拟合口径 t_dry = %.1f s ; 原始插值口径 t_dry = %.1f s ; 差 %+.1f s (%+.3f%%)"
        % (t0, td_raw, td_raw - t0, (td_raw - t0) / t0 * 100))
    log("  C0(24h): %.4f -> %.4f ; C0(48h): %.4f -> %.4f" % (base[1], c24r, base[2], c48r))

    # ---- 排序 ----
    log("")
    log("按 |相对变化| 排序（主导不确定性的参数）：")
    for nm, td, rel, ab, _, _ in sorted(rows, key=lambda r: -abs(r[2])):
        log("   %-34s %+8.3f%%" % (nm, rel))

    # ---- 解析局部灵敏度对照 ----
    log("")
    log("解析局部灵敏度对照（D = D0 exp(-a/C) exp(-E/T)）：")
    for T_C in [28.0, 40.0, 50.0]:
        TK = T_C + 273.15
        log("   T=%.0f degC: dlnD/dT = E/T^2 = %.5f /K  (即 1 K 偏差 -> D 变 %.2f%%)"
            % (T_C, 3850.0 / TK ** 2, 3850.0 / TK ** 2 * 100))
    for Cv in [2.55, 1.0, 0.15]:
        log("   C=%.2f kg/kg: dlnD/da = -1/C = %.4f  (a 相对变 2%% -> lnD 变 %.3f%%)"
            % (Cv, -1.0 / Cv, -1.0 / Cv * 0.45 * 0.02 * 100))

    np.savez(os.path.join(FC.DATADIR, "sensitivity.npz"),
             rows=np.array([[r[1], r[2], r[3]] for r in rows]),
             names=np.array([r[0] for r in rows]),
             base=np.array([t0, base[1], base[2]]),
             raw=np.array([td_raw, c24r, c48r]))
    open(OUT, "w", encoding="utf-8").write("\n".join(LOG))
    log("")
    log("已写出：" + OUT)


if __name__ == "__main__":
    main()

