# -*- coding: utf-8 -*-
# ==============================================================================
# verify_data.py —— fig22/fig23 所用验证数据的来源与一致性检查
# ==============================================================================
# 做什么：那两张"收敛性 / 积分器稳定性"图要用四类耗时复跑的结果，不可能每画一次图就重算；
#         把数字写死在绘图脚本里又会让图与代码脱钩（代码改了、图还是旧数）。故本模块负责

#         "真实复跑 → 落盘 → 校验"：recompute() 产生数据，绘图脚本只读数、且绘图前必过校验。
# 输入：common.py（物性、环境、初值）、spectral.HerbSolver、fv_theta.find_tdry 现场复跑；
#       一致性检查另读 中间数据/solutions.npz 的 p3_tdry 做交叉核对。
# 输出：<交付根>/中间数据/verification.json（带版本号 DATA_VERSION，供 fig22/fig23 使用）。
# 运行方式与依赖顺序：
#       python verify_data.py              只校验数据文件（秒级，run_all.py 与

#                                          make_figures.py 内部也是这样调用它的）；
#       python verify_data.py --recompute  真实复跑四类扫描并重写 json（约 22 min）；
#       四类扫描耗时：空间扫描约 2 min、时间扫描约 10 min、经典 FV 对照约 9 min、
#       积分器与容差约 1 min。数据文件缺失或版本不符时直接抛错并给出补救命令。
# 关键变量：
#   SPACE_N     谱节点数扫描 24/32/48/64/96/128（fig22 的自变量之一）

#   TIME_DT     定步长 theta 法的步长 40~0.25（fig22 的另一个自变量）  [s]
#   FV_N        该步长扫描用的空间网格 1600（论文表 11 的网格）
#   FV_SPACE_N  经典 FV 的空间对照 800/1600/3200/6400，固定步长 1.0   [s]
#   INTEG_CFG   积分器/容差配置 5 组（BDF×3、Radau、LSODA），加 1 s 密采样求交
#               共 6 个点，即 fig23 的 6 个刻度
#   CRIT        烘干判据阈值 0.15                                      [kg/kg]

# ==============================================================================
import json
import os
import sys

import numpy as np

from fig_core import DATADIR

DATA_PATH = os.path.join(DATADIR, "verification.json")
DATA_VERSION = "2026-09-13a"

CRIT = 0.15
# 空间扫描的网格；时间扫描的步长（s）。两者都是论文 fig22 的自变量。
SPACE_N = [24, 32, 48, 64, 96, 128]
TIME_DT = [40.0, 20.0, 10.0, 5.0, 2.0, 1.0, 0.5, 0.25]
SPACE_RTOL = 1e-11          # 与生产设置一致
FV_N = 1600                 # 定步长 theta 法的空间网格（论文表 11）
FV_SPACE_N = [800, 1600, 3200, 6400]   # fig22(a) 的经典 FV 空间对照
FV_SPACE_DT = 1.0                      # 该对照固定用的时间步（时间误差按 c*dt 扣除）
INTEG_CFG = [("BDF", 1e-9), ("BDF", 1e-10), ("BDF", 1e-12),
             ("Radau", 1e-11), ("LSODA", 1e-10)]


def _tdry(prob, N, method="BDF", rtol=1e-11, shrink=False):
    """用谱配置法求 t_dry，返回 (求解器实例, t_dry [s])。

    prob 为问题号（本模块只用到 3：附录 3 物性、半径固定）；N 是谱节点数；
    method 与 rtol 是时间积分器及其相对容差，用来考察积分器/容差带来的差异；
    shrink 只在含收缩的问题 4 下为 True。终止事件取全场最大含水率（出现在轴心 u=0）
    降到 CRIT 以下的时刻，与生产脚本同口径。
    """
    from spectral import HerbSolver
    s = HerbSolver(N, prob=prob, shrink=shrink)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])

    def ev(t, y):
        return np.max(y[:s.Np]) - CRIT
    ev.terminal = True
    ev.direction = -1
    e = s.solve(300000.0, rtol=rtol, atol=at, events=ev, method=method)
    return s, float(e.t_events[0][0])


def recompute(verbose=True):
    """真实复跑四类扫描并写进 verification.json，返回数据字典。

    verbose 控制是否边跑边打印进度。落盘的字典有四个键：space（谱网格扫描）、
    time（定步长时间扫描）、fv_space（经典 FV 空间对照）、integ（积分器与容差，
    另含一条 1 s 密采样求交的独立口径），并记录版本号与生成命令。
    """
    import common as cm
    from fv_theta import find_tdry

    def say(s):
        if verbose:
            print(s, flush=True)

    say("[verify] 1/4 空间扫描：问题 3 的 t_dry vs N")
    tdry_space = []
    for N in SPACE_N:
        _, td = _tdry(3, N, rtol=SPACE_RTOL)
        tdry_space.append(td)
        say("   N=%-4d tdry = %.4f s" % (N, td))

    say("[verify] 2/4 时间步扫描：定步长 theta 法（N=%d），最耗时" % FV_N)
    tdry_time = []
    for d in TIME_DT:
        plan = [(7200.0, min(1.0, d)), (1e18, d)]
        td = float(find_tdry(FV_N, 3, plan))
        tdry_time.append(td)
        say("   dt=%-6.2f tdry = %.4f s" % (d, td))

    say("[verify] 3/4 经典 FV 的空间对照（N=%s，dt=%.1f s）"
        % (FV_SPACE_N, FV_SPACE_DT))
    tdry_fvspace = []
    for N in FV_SPACE_N:
        plan = [(7200.0, FV_SPACE_DT), (1e18, FV_SPACE_DT)]
        td = float(find_tdry(N, 3, plan))
        tdry_fvspace.append(td)
        say("   N=%-5d tdry = %.4f s" % (N, td))

    say("[verify] 4/4 时间积分器与容差")
    rows = []
    for meth, tol in INTEG_CFG:
        _, td = _tdry(3, 64, method=meth, rtol=tol)
        rows.append({"method": meth, "rtol": tol, "tdry": td})
        say("   %-6s rtol=%.0e tdry = %.4f s" % (meth, tol, td))

    # 第 6 个点：换一种与 solve_ivp 事件机制无关的口径复核 t_dry——先在邻域内按 1 s
    # 输出求解，再对中心含水率 C(0,t) 与 CRIT 做线性求交。含水率的最大值始终出现在
    # 轴心 u=0，故事件函数取的就是 C(0,t)，这条口径与事件法等价，但"求交点"的方式
    # 完全不同，两者一致才说明 t_dry 不依赖数值求交的实现。
    s, td0 = _tdry(3, 64, rtol=1e-11)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])
    tgrid = np.arange(max(0.0, td0 - 120.0), td0 + 61.0, 1.0)
    sol = s.solve(td0 + 60.0, t_eval=list(tgrid), rtol=1e-11, atol=at)
    c0 = sol.y[0, :]
    idx = np.where(c0 <= CRIT)[0]
    k = int(idx[0])
    td_dense = float(tgrid[k - 1] + (c0[k - 1] - CRIT) / (c0[k - 1] - c0[k])
                     * (tgrid[k] - tgrid[k - 1]))
    rows.append({"method": "dense1s", "rtol": 1e-11, "tdry": td_dense})
    say("   1s 密采样求交 tdry = %.4f s" % td_dense)

    data = {
        "version": DATA_VERSION,
        "produced_by": "python verify_data.py --recompute",
        "space": {"N": SPACE_N, "tdry": tdry_space, "rtol": SPACE_RTOL},
        "time": {"dt": TIME_DT, "tdry": tdry_time, "N": FV_N},
        "fv_space": {"N": FV_SPACE_N, "tdry": tdry_fvspace, "dt": FV_SPACE_DT},
        "integ": rows,
    }
    os.makedirs(DATADIR, exist_ok=True)
    with open(DATA_PATH, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)
    say("[verify] 已写出 %s" % DATA_PATH)
    return data


def load(verbose=True):
    """读入 verification.json 并返回其字典，附带版本号校验。

    verbose 控制是否打印一行载入信息（run_all.py 用 verbose=False 静默调用）。
    文件缺失或版本号与 DATA_VERSION 不符时直接抛 RuntimeError，消息里给出重新生成
    的命令——原因是不画图也比画一张来源不明的图好。
    """
    if not os.path.isfile(DATA_PATH):
        raise RuntimeError(
            "缺少验证数据文件：%s\n"
            "  它由 fig22/fig23 使用，需要先真实复跑一次：\n"
            "      python verify_data.py --recompute      （约 22 min）\n"
            "  或运行完整复跑： python run_all.py" % DATA_PATH)
    with open(DATA_PATH, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if data.get("version") != DATA_VERSION:
        raise RuntimeError(
            "验证数据文件版本不符（文件 %r / 期望 %r）：%s\n"
            "  请重新生成： python verify_data.py --recompute"
            % (data.get("version"), DATA_VERSION, DATA_PATH))
    if verbose:
        print("[verify] 载入 %s（版本 %s，由 %s 生成）"
              % (DATA_PATH, data["version"], data["produced_by"]))
    return data


def time_slope(data):
    """由时间扫描拟合 t(dt) = L - c*dt 的斜率系数 c（s per s），供扣时间误差用。

    用最细的三个步长做最小二乘，避免粗步长段的高阶项污染一阶系数。
    """
    dt = np.asarray(data["time"]["dt"], dtype=float)
    td = np.asarray(data["time"]["tdry"], dtype=float)
    k = min(3, len(dt))
    slope, _intercept = np.polyfit(dt[-k:], td[-k:], 1)
    return float(-slope)


def check_consistency(data, tdry_cached=None, verbose=True):
    """把数据文件与本轮真实算出的量逐条对照，任何一条不成立就抛 RuntimeError。

    data 是 load() 返回的字典；tdry_cached 是解缓存里的问题 3 烘干时长 [s]，给了就
    额外核对空间扫描的 N=64 项（这是"图的数据"与"这次算的解"之间唯一的硬连接）；
    verbose 控制是否打印每条的核对结果。五条校验见下面的 (1)~(5)，都是数据自洽的
    必要条件；全部通过返回 True。
    """
    problems = []
    Ns = list(data["space"]["N"])
    ts = np.asarray(data["space"]["tdry"], dtype=float)
    dts = list(data["time"]["dt"])
    td_dt = np.asarray(data["time"]["tdry"], dtype=float)

    def note(s):
        if verbose:
            print("[verify] " + s)

    # (1) 空间扫描必须单调不减（谱解随 N 增大而收敛到同一极限）
    if not np.all(np.diff(ts) >= -1e-9):
        problems.append("空间扫描的 t_dry 随 N 非单调：%s" % np.round(ts, 3).tolist())
    # (2) 与解缓存交叉核对：N=64 项必须等于缓存里真实算出的 t_dry
    if tdry_cached is not None:
        i64 = Ns.index(64) if 64 in Ns else None
        if i64 is None:
            problems.append("空间扫描缺少 N=64 项，无法与解缓存核对")
        else:
            gap = abs(ts[i64] - float(tdry_cached))
            note("空间扫描 N=64 = %.4f s vs 解缓存 %.4f s（差 %.4f s）"
                 % (ts[i64], tdry_cached, gap))
            # 0.05 s 只是"数据是否过期"的防呆阈值：口径一致时两者逐位相等（差 0.0000 s）。
            if gap > 0.05:
                problems.append(
                    "空间扫描 N=64 的 t_dry 与解缓存相差 %.3f s，数据文件已过期"
                    % gap)
    # (3) 时间扫描必须随 dt 减小而单调增大，且增量与 dt 成正比（严格一阶）
    if not np.all(np.diff(td_dt) >= -1e-9):
        problems.append("时间扫描的 t_dry 随 dt 减小非单调：%s" % np.round(td_dt, 3).tolist())
    inc = np.diff(td_dt)                    # dt 由大到小，增量应为正且递减
    if inc.size >= 3:
        # 一阶收敛的判据：log(相邻增量) 对 log(dt 的几何中点) 的斜率应为 1。
        # 注意增量是**减半**的（比值约 0.5），不是加倍——用比值判据容易反向。
        # 容差取 0.8~1.25：既容得下离散噪声，又能把斜率约 2 的二阶数据或非单调数据挡掉。
        dts_arr = np.asarray(dts, dtype=float)
        mid = np.sqrt(dts_arr[1:] * dts_arr[:-1])
        slope = float(np.polyfit(np.log(mid), np.log(inc), 1)[0])
        if not (0.8 < slope < 1.25):
            problems.append("时间扫描不满足 O(dt) 一阶收敛：增量-步长双对数斜率 %.3f（应接近 1）"
                            % slope)
        else:
            note("时间扫描一阶收敛：增量-步长双对数斜率 %.3f（应接近 1）" % slope)
    else:
        note("时间扫描只有 %d 个点，跳过收敛阶检查" % len(dts))

    # (4) 经典 FV 的空间对照应随 N 增大单调收敛，且扣去时间误差后仍不高于连续极限
    fs = np.asarray(data["fv_space"]["tdry"], dtype=float)
    c_fv = time_slope(data)                 # 由时间扫描拟合出的 O(dt) 系数 c
    if not np.all(np.diff(fs) >= -1e-9):
        problems.append("经典 FV 空间对照的 t_dry 随 N 非单调：%s" % np.round(fs, 3).tolist())
    # 留 0.5 s 余量：扣掉时间误差后若仍比谱解极限高出 0.5 s 以上，说明系数 c 不合理，
    # 或者这批数据已经过期。
    if len(fs) and fs[-1] + c_fv * data["fv_space"]["dt"] > ts[-1] + 0.5:
        problems.append("经典 FV 扣去时间误差后超过了谱解连续极限，说明时间误差系数 c=%.3f 不合理"
                        % c_fv)
    else:
        note("经典 FV 最细网格扣时间误差 = %.3f s，谱解极限 %.3f s（差 %.3f s）"
             % (fs[-1] + c_fv * data["fv_space"]["dt"], ts[-1],
                ts[-1] - fs[-1] - c_fv * data["fv_space"]["dt"]))

    # (5) 积分器与容差：fig23 的"极差"必须真的小
    iv = np.asarray([r["tdry"] for r in data["integ"]], dtype=float)
    spread = float(iv.max() - iv.min())
    if len(iv) != len(INTEG_CFG) + 1:
        problems.append("积分器配置数 %d 与 fig23 的 6 个刻度不符" % len(iv))
    # fig23 的结论是"六种配置的极差只有 0.007 s"，故超过 0.05 s 就说明该结论不成立。
    if spread > 0.05:
        problems.append("积分器/容差极差 %.4f s 过大，fig23 的结论不成立" % spread)
    else:
        note("积分器与容差极差 %.4f s（%d 个配置）" % (spread, len(iv)))

    if problems:
        raise RuntimeError(
            "验证数据不一致，已中止绘图（避免画出与代码结果脱钩的图）：\n  - "
            + "\n  - ".join(problems)
            + "\n  请重新生成： python verify_data.py --recompute")
    note("一致性检查通过")
    return True


def _cached_tdry():
    """从解缓存里读问题 3 的 t_dry，用于与空间扫描的 N=64 项交叉核对。"""
    path = os.path.join(DATADIR, "solutions.npz")
    if not os.path.isfile(path):
        return None
    try:
        z = np.load(path)
        return float(z["p3_tdry"][0])
    except Exception as exc:                      # 缓存损坏或键缺失：跳过核对即可
        print("[verify] 读解缓存失败，跳过 N=64 交叉核对：%s" % exc)
        return None


def main(argv=None):
    """命令行入口：带 --recompute 就先重跑四类扫描，然后统一做一次一致性检查。

    argv 不给时取 sys.argv[1:]。返回 0 表示校验通过；不通过时异常抛出，
    由 sys.exit 带出非零退出码。
    """
    argv = sys.argv[1:] if argv is None else argv
    if "--recompute" in argv:
        recompute()
    check_consistency(load(), tdry_cached=_cached_tdry())
    return 0


if __name__ == "__main__":
    sys.exit(main())
