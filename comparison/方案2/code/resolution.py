# -*- coding: utf-8 -*-
# ============================================================================
# resolution.py —— 分辨率研究：网格与步长加密后结果是否已经收敛
#
# 对问题 3（半径不变）与问题 4（有收缩）分别加密径向网格：直线法取 75~1200 个
# 节点；有限体积分别用节点中心与单元中心两种布局、各 200~1600 个控制体。比较
# t = 2.0e5 s 时的中心与表面含水率，按相邻两档差值的对数比算出观测收敛阶；
# 再把完整 PDE 的截面平均含水率与第一特征值降维模型逐时刻对比，量化「只保留
# 第一模态」带来的偏差。
#
# 输入数据：物性用附录 3 / 4 的经验关系，环境历程来自附件 1，收缩半径来自
#           附件 2（都由 common.py 读取）。直线法直接复用 crosscheck.mol_solve。
# 输出：data/resolution_problem3.csv、data/resolution_problem4.csv（各方法的
#       逐档结果表），data/resolution_summary.txt（控制台报告全文）。
#
# 关键变量：
#   T_TARGET = 2.0e5 s  比较时刻；C_TARGET = 0.15 kg/kg 干燥终点判据
#   order  观测收敛阶：相邻两档差值的对数比
#   XI_W   21 个归一化半径位置 xi = 0, 0.05, ..., 1.0
#   W_MEAN 对应的截面平均权重：梯形法则的 2*xi*0.05，两端点权重减半，
#          使得 sum(C*W_MEAN) = 2*int_0^1 C*xi dxi 就是截面平均含水率
# ============================================================================
"""
运行方式：在目录 A题解答 下执行 python code/resolution.py（约需数分钟）。

报告分三节：问题 3 的分辨率表与收敛阶、问题 4 的分辨率表、截面平均含水率在
完整 PDE 与降维模型之间的偏差。
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
import crosscheck as cc

T_TARGET = 2.0e5
C_TARGET = 0.15
XI_W = np.array([round(0.1 * k, 1) for k in range(21)]) / 2.0
W_MEAN = 2.0 * XI_W * 0.05
W_MEAN[0] *= 0.5
W_MEAN[-1] *= 0.5


def mean_curve(res):
    """算截面平均含水率 <C> = 2*int_0^1 C*xi dxi，两种求解器的结果都能直接用。

    有限体积在 XI_W 上采样，直接用预先算好的权重 W_MEAN 加权求和即可；
    直线法是节点中心网格（n+1 个节点），须临时组装一次节点体积权重再累加。
    返回与 res["t"] 等长的一维数组，单位 kg/kg。
    """
    C = res["C"]
    if C.shape[1] == len(XI_W):
        return C @ W_MEAN
    n = C.shape[1] - 1                      # 直线法：节点中心，n+1 个节点
    h = 1.0 / n
    vol = np.empty(n + 1)
    vol[0] = h * h / 8.0
    vol[1:n] = (np.arange(1, n) * h) * h
    vol[n] = h / 2.0 - h * h / 8.0
    return 2.0 * (C @ vol)


def main():
    """跑完三节并把报告写入 resolution_summary.txt，同时导出两张 CSV。

    P(s) 既打印到控制台又把同一行收进 out，最后一次性写文件。
    """
    room = cm.RoomConditions()
    rad = cm.RadiusHistory(mode="table")
    out = []
    P = lambda s: (print(s), out.append(s))

    P("=" * 78)
    P("RESOLUTION STUDY, PROBLEM 3 (Appendix 3, fixed radius)")
    P("=" * 78)
    # 把 T_TARGET 并进采样时刻，是为了让直线法的结果正好落在 t = 2.0e5 s 上。
    # 若只按 600 s 的均匀网格采样，最近点是 199800 s，读到的其实是另一个时刻的
    # 值，与有限体积各档（用 t_report = T_TARGET）不同步，会带来约 0.04% 的偏差，
    # 比这张表要分辨的方法间差异还大。
    t_ev = np.unique(np.concatenate([np.arange(0.0, T_TARGET + 1.0, 600.0),
                                     [T_TARGET]]))
    rows = []
    for nn in (75, 150, 300, 600, 1200):
        t0 = time.time()
        r = cc.mol_solve(cm.props_problem23, T_TARGET, n_node=nn, room=room,
                         t_eval=t_ev)
        i = int(np.argmin(np.abs(r["t"] - T_TARGET)))
        # dt=600 是输出采样间隔，不是求解步长：BDF 的步长由误差控制自适应选择，
        # 与下面有限体积各档固定 1 s 的步长含义不同
        rows.append(dict(method="MOL/BDF", n=nn, dt=600.0,
                         C_center=r["C"][i, 0], C_surf=r["C"][i, -1],
                         C_mean=float(mean_curve(r)[i]), wall=time.time() - t0))
    for grid in ("node", "cell"):
        for nc in (200, 400, 800, 1600):
            t0 = time.time()
            r = cm.simulate(nc, cm.props_problem23, T_TARGET, 1.0, room,
                            grid=grid, t_report=np.array([T_TARGET]),
                            xi_report=XI_W)
            rows.append(dict(method="FV/CN-%s" % grid, n=nc, dt=1.0,
                             C_center=r["C"][0, 0], C_surf=r["C"][0, -1],
                             C_mean=float(r["C"][0] @ W_MEAN),
                             wall=time.time() - t0))
    df = pd.DataFrame(rows)
    P(df.to_string(index=False))
    df.to_csv(os.path.join(cm.DATA_DIR, "resolution_problem3.csv"), index=False)

    P("")
    P("Convergence of C_centre at t = 2.0e5 s, per scheme")
    # 观测收敛阶：若误差 ~ C*h^p，则相邻两档的差值之比给出 p = log(前差/后差)
    # / log(网格数之比)。第 2 档（k=1）还没有两个差值可用，故记 nan。
    refs = {}
    for grid in ("node", "cell"):
        ref = df[df.method == "FV/CN-%s" % grid].sort_values("n")
        c = ref["C_center"].to_numpy()
        n = ref["n"].to_numpy()
        refs[grid] = c[-1]
        for k in range(1, len(c)):
            order = np.log(abs(c[k - 1] - c[k - 2]) / abs(c[k] - c[k - 1])) / \
                np.log(n[k] / n[k - 1]) if k >= 2 else np.nan
            P("   FV/CN-%-4s n=%5d  C=%0.7f   diff=%0.2e   observed order=%s"
              % (grid, n[k], c[k], abs(c[k] - c[k - 1]),
                 "%.2f" % order if order == order else "-"))
    m = df[df.method == "MOL/BDF"].sort_values("n")
    cm_ = m["C_center"].to_numpy()
    nm = m["n"].to_numpy()
    for k in range(1, len(cm_)):
        order = np.log(abs(cm_[k - 1] - cm_[k - 2]) / abs(cm_[k] - cm_[k - 1])) / \
            np.log(nm[k] / nm[k - 1]) if k >= 2 else np.nan
        P("   MOL    n=%5d  C=%0.7f   diff=%0.2e   observed order=%s"
          % (nm[k], cm_[k], abs(cm_[k] - cm_[k - 1]),
             "%.2f" % order if order == order else "-"))
    P("   node(N=1600) - cell(N=1600) = %.2e  (%.3f %% of the value)"
      % (refs["node"] - refs["cell"],
         100.0 * abs(refs["node"] - refs["cell"]) / refs["node"]))
    P("   node(N=1600) - MOL(n=1200)  = %.2e  (%.3f %% of the value)"
      % (refs["node"] - cm_[-1],
         100.0 * abs(refs["node"] - cm_[-1]) / refs["node"]))

    P("")
    P("=" * 78)
    P("RESOLUTION STUDY, PROBLEM 4 (Appendix 4, shrinkage)")
    P("=" * 78)
    rows = []
    for nn in (150, 300, 600, 1200):
        r = cc.mol_solve(cm.props_problem4, T_TARGET, n_node=nn, room=room,
                         radius=rad, t_eval=t_ev)
        i = int(np.argmin(np.abs(r["t"] - T_TARGET)))
        rows.append(dict(method="MOL/BDF", n=nn, C_center=r["C"][i, 0],
                         C_surf=r["C"][i, -1], C_mean=float(mean_curve(r)[i])))
    for grid in ("node", "cell"):
        for nc in (200, 400, 800, 1600):
            r = cm.simulate(nc, cm.props_problem4, T_TARGET, 1.0, room,
                            radius=rad, grid=grid,
                            t_report=np.array([T_TARGET]), xi_report=XI_W)
            rows.append(dict(method="FV/CN-%s" % grid, n=nc,
                             C_center=r["C"][0, 0], C_surf=r["C"][0, -1],
                             C_mean=float(r["C"][0] @ W_MEAN)))
    df4 = pd.DataFrame(rows)
    P(df4.to_string(index=False))
    df4.to_csv(os.path.join(cm.DATA_DIR, "resolution_problem4.csv"), index=False)

    P("")
    P("=" * 78)
    P("MEAN MOISTURE: FULL PDE vs REDUCED FIRST-EIGENVALUE MODEL")
    P("=" * 78)
    for tag, prop, radius in [("problem 3", cm.props_problem23, None),
                              ("problem 4", cm.props_problem4, rad)]:
        # 完整 PDE 用最细网格（1600 个控制体）算一遍，再与降维模型在同一时刻比较；
        # 降维模型的时间轴是 30 s 均匀的，故用 np.interp 插到 PDE 的采样时刻上
        r = cm.simulate(1600, prop, T_TARGET, 1.0, room, radius=radius,
                        grid="node", t_report=t_ev, xi_report=XI_W)
        cbar_pde = r["C"] @ W_MEAN
        ts, cs = cc.lumped_mean(prop, room=room, radius=radius)
        cbar_lum = np.interp(r["t"], ts, cs)
        # 前 6 h 还在预热与升温阶段，两种模型的差异被瞬态掩盖，只统计之后的部分
        mask = r["t"] >= 6.0 * 3600.0
        rel = np.abs(cbar_lum[mask] - cbar_pde[mask]) / cbar_pde[mask]
        P("  %s : Cbar(6h,12h,24h,48h) PDE = %s"
          % (tag, np.round(np.interp([6, 12, 24, 48], r["t"] / 3600.0,
                                     cbar_pde), 5)))
        P("  %s : Cbar(6h,12h,24h,48h) LUM = %s"
          % (tag, np.round(np.interp([6, 12, 24, 48], ts / 3600.0, cs), 5)))
        P("  %s : max relative deviation after 6 h = %.2f %%"
          % (tag, 100.0 * rel.max()))
        # 与中心判据的换算：第一模态下 C_center/C_mean = mu1/(2*J1(mu1)) ≈ 2.3
        # （见 crosscheck.shape_factor），所以截面平均总是先降到 0.15 kg/kg，
        # 中心的同一判据要更晚才满足，两个时刻不能混用。
        j = np.where(cbar_pde <= 0.15)[0]
        if len(j):
            P("      the cross-section MEAN reaches 0.15 kg/kg at %.2f h "
              "(the centre needs longer)" % (r["t"][j[0]] / 3600.0))

    with open(os.path.join(cm.DATA_DIR, "resolution_summary.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("\n   saved ->", os.path.join(cm.DATA_DIR, "resolution_summary.txt"))


if __name__ == "__main__":
    main()
