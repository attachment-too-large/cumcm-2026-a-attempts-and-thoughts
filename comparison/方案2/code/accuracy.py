# -*- coding: utf-8 -*-
# ============================================================================
# accuracy.py —— 干燥时间的精度评估：加密网格/步长 + Richardson 外推
#
# 对问题 3、问题 4 各跑 8 组配置（两种布局 · 四档网格数 · 两种时间步安排里挑出
# 的 8 个组合），在每秒采样的最大含水率 Cmax(t) 上求它首次降到 0.15 kg/kg 的
# 时刻，即干燥时间 t_dry。再用最细几档做 Richardson 外推，估计网格无限细时的
# 极限值；并把「最细与次细两档之差」当作剩余离散误差的保守估计，给出 t_dry
# 的不确定度。
#
# 输入数据：物性用附录 3 / 4 的经验关系，环境历程来自附件 1，问题 4 的收缩半径
#           来自附件 2（都由 common.py 读取）。
# 输出：data/accuracy_summary.txt（控制台报告全文），本脚本不写其它文件。
#
# 关键变量：
#   C_TARGET = 0.15 kg/kg  干燥终点判据；T_END = 2.6e5 s 模拟终止时刻
#   SCHED_FAST  前 7200 s 用 1 s 步长、之后放宽到 10 s 的两段式时间步安排
#   SCHED_1S    全程 1 s 步长，作为时间步收敛的对照
#   Cmax(t)     输出点（轴心 xi=0 与表面 xi=1）上的最大含水率。中心在整个过程中
#               都比表面湿，是最后干的点，故它降到判据以下就代表整根药材都达标
#   p_obs       由最细三档估出的观测收敛阶；t_rich2 / t_richp 按名义阶与观测阶
#               分别外推得到的极限值
# ============================================================================
"""
运行方式：在目录 A题解答 下执行 python code/accuracy.py（约需数分钟）。

每个问题两张表：8 组配置的 t_dry 与耗时、最细几档的收敛性与外推极限。
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
SCHED_FAST = [(7200.0, 1.0), (T_END, 10.0)]
SCHED_1S = [(T_END, 1.0)]


def drying_time(prop, radius, n_cell, schedule, grid="node"):
    """按指定布局、网格数与时间步安排算一次干燥时间 t_dry [s]。

    schedule 为 common.simulate 用的两段式步长安排 [(终止时刻, 步长), ...]；
    输出时刻取每 1 s 一点、采样位置为轴心 xi=0 与表面 xi=1，取两点上的较大值
    作为 Cmax(t)（中心始终比表面湿，故它就是剖面上的最大值）。
    返回 (t_dry, 求解结果, Cmax 序列)：t_dry 由跨越判据的相邻两个采样点线性插值
    得到；若 T_END 内始终没降到 C_TARGET 以下则为 nan。
    """
    # 这里第 4 个位置传的 1.0 是 simulate 的默认步长，实际推进完全由 dt_schedule
    # 决定，故该参数在本调用中没有作用
    res = cm.simulate(n_cell, prop, T_END, 1.0, cm.RoomConditions(),
                      radius=radius, grid=grid, dt_schedule=schedule,
                      t_report=np.arange(1.0, T_END + 0.5, 1.0),
                      xi_report=np.array([0.0, 1.0]))
    Cmax = res["C"].max(axis=1)
    idx = np.where(Cmax < C_TARGET)[0]
    if not len(idx):
        return np.nan, res, Cmax
    i0 = idx[0]
    t0, t1 = res["t"][i0 - 1], res["t"][i0]
    v0, v1 = Cmax[i0 - 1], Cmax[i0]
    return t0 + (C_TARGET - v0) * (t1 - t0) / (v1 - v0), res, Cmax


def main():
    """跑完两个问题的 8 组配置，做外推并把报告写入 accuracy_summary.txt。

    P(s) 既打印（带 flush，便于在后台日志里看进度）又把同一行收进 out。
    """
    rad = cm.RadiusHistory(mode="table")
    out = []
    P = lambda s: (print(s, flush=True), out.append(s))

    for tag, prop, radius in [("problem 3 (Appendix 3)", cm.props_problem23, None),
                              ("problem 4 (Appendix 4, shrinking)",
                               cm.props_problem4, rad)]:
        P("=" * 78)
        P(tag)
        P("=" * 78)
        P("  %-6s %-6s %-8s %-16s %-12s %-10s" %
          ("grid", "n", "dt", "t_dry / s", "t_dry / h", "wall / s"))
        recs = []
        # 8 组配置：(布局, 控制体数, 步长安排, 表头里的步长标签)
        # 前 6 组都是节点中心，网格数逐档翻倍；后两组换成单元中心布局作对照。
        for grid, n_cell, sched, dtlabel in [
                ("node", 200, SCHED_FAST, "1-10"),
                ("node", 400, SCHED_FAST, "1-10"),
                ("node", 800, SCHED_FAST, "1-10"),
                ("node", 1600, SCHED_FAST, "1-10"),
                ("node", 800, SCHED_1S, "1"),
                ("node", 1600, SCHED_1S, "1"),
                ("cell", 800, SCHED_FAST, "1-10"),
                ("cell", 1600, SCHED_FAST, "1-10")]:
            t0 = time.time()
            td, res, Cmax = drying_time(prop, radius, n_cell, sched, grid=grid)
            wall = time.time() - t0
            P("  %-6s %-6d %-8s %-16.2f %-12.4f %-10.1f"
              % (grid, n_cell, dtlabel, td, td / 3600.0, wall))
            recs.append(dict(grid=grid, n=n_cell, dt=dtlabel, t_dry=td))
        df = pd.DataFrame(recs)
        # 外推只针对「节点中心 + 1-10 步长」这一列：只有它是网格数逐档翻倍的
        # 加密序列，其他行（单元中心、全程 1 s 步长）另作对照
        seq = df[(df.grid == "node") & (df.dt == "1-10")].sort_values("n")
        v = seq["t_dry"].to_numpy()
        n = seq["n"].to_numpy()
        P("")
        P("  node-grid convergence (dt = 1 s -> 10 s):")
        for k in range(len(v)):
            if k == 0:
                P("    n=%5d  t_dry = %10.2f s" % (n[k], v[k]))
            else:
                # 观测收敛阶：网格数每档翻倍，故取以 2 为底的对数比；第 2 档
                # 还没有两个差值可用，记 nan（打印成 "-"）
                order = np.log2(abs(v[k - 1] - v[k - 2]) /
                                abs(v[k] - v[k - 1])) if k >= 2 else np.nan
                P("    n=%5d  t_dry = %10.2f s   diff = %7.2f s   order = %s"
                  % (n[k], v[k], abs(v[k] - v[k - 1]),
                     ("%.2f" % order) if order == order else "-"))
        # Richardson 外推。设误差 ~ C*h^p、网格加密一倍，则真值约为
        #   t_exact = t_fine + (t_fine - t_coarse) / (2^p - 1)。
        # 必须用本序列实测出的观测阶 p_obs，不能用名义阶 2，否则给出的极限值
        # 与同一张表里的收敛阶自相矛盾。
        p_obs = np.log2(abs(v[-2] - v[-3]) / abs(v[-1] - v[-2]))
        t_rich2 = v[-1] + (v[-1] - v[-2]) / 3.0            # p = 2 时的修正量
        t_richp = v[-1] + (v[-1] - v[-2]) / (2.0 ** p_obs - 1.0)
        P("  Richardson limit (nominal order 2)    : %.2f s = %.4f h"
          % (t_rich2, t_rich2 / 3600.0))
        P("  Richardson limit (observed order %.2f): %.2f s = %.4f h"
          % (p_obs, t_richp, t_richp / 3600.0))
        P("  the two differ by                     : %.2f s (%.4f %%)"
          % (abs(t_rich2 - t_richp), 100.0 * abs(t_rich2 - t_richp) / t_rich2))
        # 不确定度取 |最细档 - 次细档|/3，即 p=2 时的外推修正量：它同时也是
        # 「网格再加密一档、结果还能移动多少」的保守估计
        t_rich = t_rich2
        P("  spread over the finest configurations : %.2f s (%.4f %%)"
          % (abs(v[-1] - df[df.dt == "1"].t_dry.min()),
             100.0 * abs(v[-1] - df[df.dt == "1"].t_dry.min()) / v[-1]))
        P("  => t_dry = %.0f +/- %.0f s = %.2f +/- %.2f h"
          % (t_rich, abs(v[-1] - v[-2]) / 3.0,
             t_rich / 3600.0, abs(v[-1] - v[-2]) / 3.0 / 3600.0))
        P("")

    with open(os.path.join(cm.DATA_DIR, "accuracy_summary.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("   saved ->", os.path.join(cm.DATA_DIR, "accuracy_summary.txt"))


if __name__ == "__main__":
    main()
