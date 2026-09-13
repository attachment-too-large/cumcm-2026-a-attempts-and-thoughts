# -*- coding: utf-8 -*-
# ==============================================================================
# sensitivity.py —— 关键参数扰动的灵敏度分析（真实重算，不是口头推断）
# ==============================================================================
# 做什么：逐个扰动物性经验式与传递系数，把问题 3 的整个烘干过程重新解一遍，看烘干时长
#         t_dry 与两个时刻的中心含水率随之变化多少，据此排定"哪个参数最要紧"。
# 输入：common.py（物性、一阶惯性环境、初值）、spectral.HerbSolver（主求解器）、
#       fig_core.DATADIR（输出目录）；基准解与 result3.xlsx 用同一网格 BASE_N。
# 输出：<交付根>/中间数据/sensitivity.txt（可读表格）与 sensitivity.npz（供绘图与核对）。
# 运行方式与依赖顺序：python sensitivity.py，可单独运行，与 run_all.py 互不依赖；
#       每个扰动都要从头重解一次（共 14 次），耗时主要在求解上。
# 关键变量：
#   BASE_N   基准谱节点数 64（与 result3.xlsx 的生产网格一致）
#   CRIT     烘干判据阈值 0.15                                 [kg/kg]
#   C_FLOOR  物性式求值用的含水率下限 1e-3（防 exp(-a/C) 发散） [kg/kg]
#   扰动幅度  物性与传递系数取 ±2%：题目未给经验式的标定区间，2% 足以看清各参数的相对
#             重要性，而网格与容差带来的数值误差只有 1e-5 相对量级；Tset 取 ±0.2 K。
# 扰动对象与考察指标：
#   物性侧  活化温度 E=3850 K、活化浓度指数 a、扩散系数前因子 D0；
#   传递侧  对流换热系数 h、对流传质系数 hm；
#   环境侧  烘房设定温度 T_set，以及环境口径（一阶惯性拟合式 vs 附件 1 原始序列插值）。
#   指标    问题 3 的烘干时长 t_dry，以及 t=24 h、48 h 的中心含水率；所有相对变化都以
#           本网格（BASE_N）的基准 t_dry 为分母（连续极限 205575.5 s 见论文第 11.4 节）。
# ==============================================================================
import os
from console_utf8 import fix_console
fix_console()
import numpy as np
import common as cm
from spectral import HerbSolver
import fig_core as FC

CRIT = 0.15
C_FLOOR = 1e-3
BASE_N = 64          # 灵敏度基准网格，与 result3.xlsx 的生产网格一致
OUT = os.path.join(FC.DATADIR, "sensitivity.txt")
NPZ = os.path.join(FC.DATADIR, "sensitivity.npz")
LOG = []


def log(s):
    """写一行灵敏度报告：同时进内存缓冲 LOG 与屏幕，末尾一起落盘成 txt。"""
    print(s, flush=True)
    LOG.append(str(s))


def props_pert(prob, C, T, pert):
    """按扰动因子重算物性，返回 (rho, cp, k, D)，与 common.props 的结构一致。

    prob 决定用哪套经验式（1、2/3、4 三支）；C 为干基含水率 [kg/kg]，T 为温度 [degC]；
    pert 是扰动字典，键 E / a / D0 给出对应因子的乘数（缺省 1.0，即不扰动）；
    h 与 hm 只出现在边界条件里，不参与物性式，故这里不处理。
    量纲同 common.props：rho [kg/m3]、cp [J/(kg K)]、k [W/(m K)]、D [m2/s]。
    """
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
    # 返回 (rho, cp, k, D)：h / hm 不参与物性式（它们只在边界条件里出现），
    # 由调用方按 perturb 直接改 cm.h / cm.hm，故不在此返回。
    return rho, cp, k, D


def run_case(tag, pert=None, ambient="fit", N=BASE_N, prob=3, probe1=False):
    """在给定扰动下重解一遍，返回 (t_dry [s], C0(24h), C0(48h), probe)。

    扰动只在本函数内部生效：临时代换 cm.props / cm.h / cm.hm / cm.Tset / cm.T_air /
    cm.C_air，算完在 finally 里全部还原，故各组扰动之间互不影响。扰动的口径如下：
      * E / a / D0 —— 只改物性式中的对应因子，物性关系其余部分不变；
      * h / hm     —— 只改这两个传递系数；
      * Tset       —— 只把烘房设定值平移（温差，±0.2 K）。一阶惯性式的初值锚定在实测的
                      T_air(0)=28 degC 上，故初值不受影响，被扰动的是设定值与整条升温
                      轨迹；
      * ambient    —— "fit" 用拟合曲线；"raw" 把温度与含湿量同时换回附件 1 实测序列的
                      线性插值（4 h 之后仍取各自设定值）。
    probe 是三元组 (T(0,1800), T(R,1800), C(R,1800))，仅 probe1=True 时给出，供论文
    "环境口径对结果的影响"一表直接引用；两个中心含水率取不到时返回 nan。
    tag 是调用处给这一组扰动起的名字（由调用方打印进报告），函数体内并不使用它。
    """
    import spectral
    pert = pert or {}
    base_props = cm.props
    base_T, base_C = cm.T_air, cm.C_air
    base_h, base_hm = cm.h, cm.hm
    base_Tset = cm.Tset

    def _props(p, C, T):
        rho, cp, k, D = props_pert(p, C, T, pert)
        return rho, cp, k, D

    cm.props = _props
    cm.h = base_h * pert.get("h", 1.0)
    cm.hm = base_hm * pert.get("hm", 1.0)
    if "Tset" in pert:
        # 直接平移设定值：cm.T_air 的表达式 Tset-(Tset-T_air0)exp(-t/tauT)
        # 自动把 t=0 锚在 T_air0 上，因此不会像  base_T(t)+dT  那样连初值一起改。
        cm.Tset = base_Tset + pert["Tset"]
    if ambient == "raw":
        cm.T_air = lambda t: np.where(np.asarray(t, dtype=float) <= cm.T_END,
                                      np.interp(np.asarray(t, dtype=float),
                                                cm._t_env, cm._Ta_env), cm.Tset)
        cm.C_air = lambda t: np.where(np.asarray(t, dtype=float) <= cm.T_END,
                                      np.interp(np.asarray(t, dtype=float),
                                                cm._t_env, cm._Ca_env), cm.Cset)

    try:
        # 状态按 (C_0..C_N, T_0..T_N) 排列，atol 分两段给：水分 1e-13、温度 1e-11。
        s = spectral.HerbSolver(N, prob=prob)
        at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])

        def ev(t, y):
            return np.max(y[:s.Np]) - CRIT
        ev.terminal = True
        ev.direction = -1
        e = s.solve(320000.0, rtol=1e-10, atol=at, events=ev)
        td = float(e.t_events[0][0])
        # 中心含水率取 24 h 与 48 h 两个固定时刻：事件求解不保留这两个时刻的值，
        # 故另起一次带 t_eval 的求解；只取落在 [0, min(t_dry, 48 h)] 内的时刻，
        # 取不到的就保持 nan。
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

        # 问题 1 的三个报数：供论文"环境口径对结果的影响"表直接引用，
        # 避免表中数字与脚本输出脱钩。
        probe = None
        if probe1:
            s1 = spectral.HerbSolver(64, prob=1)
            at1 = np.concatenate([np.full(s1.Np, 1e-13), np.full(s1.Np, 1e-11)])
            so1 = s1.solve(1800.0, t_eval=[1800.0], rtol=1e-11, atol=at1)
            c_c, t_c = s1.interp_u(so1.y[:, -1], np.array([0.0]))
            c_s, t_s = s1.interp_u(so1.y[:, -1], np.array([1.0]))
            probe = (float(t_c[0]), float(t_s[0]), float(c_s[0]))
    finally:
        cm.props = base_props
        cm.T_air, cm.C_air = base_T, base_C
        cm.h, cm.hm = base_h, base_hm
        cm.Tset = base_Tset
    return td, c24, c48, probe


def main():
    log("=" * 74)
    log("灵敏度分析（问题 3，烘干时长 t_dry；基准 t_dry 由 BDF 求解给出）")
    log("基准网格 N=%d，基准口径（问题 3 物性 + 一阶惯性拟合环境）" % BASE_N)
    log("注意：基准值随网格而变（N=48 得 205572.3 s，N=64 得 205574.5 s，")
    log("      连续极限 205575.5 s）。本表所有相对变化均以本基准为分母。")
    log("=" * 74)

    base = run_case("base", probe1=True)
    log("基准（无扰动）：t_dry = %.1f s = %.4f h ; C0(24h)=%.4f ; C0(48h)=%.4f"
        % (base[0], base[0] / 3600, base[1], base[2]))
    t0 = base[0]

    # 12 组扰动 = 6 个参数各两档：物性与传递系数取 ±2%，烘房设定温度取 ±0.2 K。
    # 行标签同时写清了扰动前后的数值，论文表 8 与这里的口径必须一致。
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
        td, c24, c48, _ = run_case(name, pert)
        rel = (td - t0) / t0 * 100
        rows.append((name, td, rel, td - t0, c24, c48))
        log("%-34s %12.1f %9.3f%% %10.1f" % (name, td, rel, td - t0))
        print("   C0(24h)=%.4f  C0(48h)=%.4f" % (c24, c48))

    # ---- 环境口径分叉：温度与含湿量同时换回附件 1 实测序列插值 ----
    # 这是建模选择而不是参数扰动：附件 1 的拟合残差存在系统性失配（论文图 3），
    # 故要专门看一眼"换用原始序列插值"会不会改变结论。
    log("")
    log("环境口径分叉（附件 1 原始序列插值 vs 一阶惯性拟合）：")
    log("  两种口径的 T_air 与 C_air 都同时切换；4 h 之后的平台段取各自设定值。")
    td_raw, c24r, c48r, probe_raw = run_case("raw", ambient="raw", probe1=True)
    log("  拟合口径   t_dry = %.1f s ; T(0,1800)=%.4f ; T(R,1800)=%.4f ; C(R,1800)=%.4f"
        % (t0, base[3][0], base[3][1], base[3][2]))
    log("  原始插值   t_dry = %.1f s ; T(0,1800)=%.4f ; T(R,1800)=%.4f ; C(R,1800)=%.4f"
        % (td_raw, probe_raw[0], probe_raw[1], probe_raw[2]))
    log("  差：t_dry %+.1f s (%+.4f%%) ; T(0) %+.4f K ; T(R) %+.4f K ; C(R) %+.4f"
        % (td_raw - t0, (td_raw - t0) / t0 * 100,
           probe_raw[0] - base[3][0], probe_raw[1] - base[3][1],
           probe_raw[2] - base[3][2]))
    log("  C0(24h): %.4f -> %.4f ; C0(48h): %.4f -> %.4f" % (base[1], c24r, base[2], c48r))

    # ---- 排序：按 |相对变化| 从大到小排，主导参数一眼可见 ----
    log("")
    log("按 |相对变化| 排序（主导不确定性的参数）：")
    for nm, td, rel, ab, _, _ in sorted(rows, key=lambda r: -abs(r[2])):
        log("   %-34s %+8.3f%%" % (nm, rel))

    # ---- 解析局部灵敏度对照：由 lnD = lnD0 - a/C - E/T_K 直接求偏导，与上面的全模型
    # 重算比量级。注意它是某个温度、某个含水率下的局部导数，不能代替全程重算的 t_dry
    # 变化——两者量级相符只说明扰动响应方向合理，最终结论仍以全模型重算为准。
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
             # 基准网格一并存档，便于事后核对表中数字的口径
             base_N=np.array([BASE_N]),
             probe_fit=np.array(base[3]),
             raw=np.array([td_raw, c24r, c48r]),
             probe_raw=np.array(probe_raw))
    open(OUT, "w", encoding="utf-8").write("\n".join(LOG))
    log("")
    log("已写出：" + OUT)


if __name__ == "__main__":
    main()

