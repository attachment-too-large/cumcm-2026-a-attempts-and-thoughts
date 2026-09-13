# -*- coding: utf-8 -*-
# ============================================================================
# p23_process.py —— 问题 2、问题 3：整段烘干过程的求解与烘干时间确定
#
# 做什么：按一维轴对称热-质扩散模型求解全程（0 ~ 72 h）的温度场与水分浓度场，
#         再用「全场含水率最大值降到阈值以下」的判据求出烘干时长 t_dry。
#
# 输入：物性取题目附录 3（变物性经验式）；烘房环境为「预热平衡阶段按附件 1
#       上升，t > 14400 s 后保持恒温恒湿」。
#
# 输出：data/result2.xlsx      结果文件（1 s x 0.1 cm，覆盖 0~72 h 全过程）
#       data/result3.xlsx      结果文件（60 s x 0.1 cm，直到烘干结束时刻）
#       data/p2_table3_T.csv   论文表 3：每 0.5 h、若干半径处的温度
#       data/p2_table4_C.csv   论文表 4：同上，但为水分浓度
#       data/p3_table5_C.csv   论文表 5：每 6 h 的水分浓度与烘干结束时刻
#       data/p23_fields.npz    完整场，供绘图脚本读取
#       data/p23_summary.txt   诊断量：烘干时长、判据检验、逐时剖面
#
# 关键变量：
#   N_CELL      径向控制体数目，取 1600（节点中心格式下节点数为 N_CELL + 1）
#   T_END       模拟终了时刻 72 h                                    [s]
#   DT_SCHEDULE 时间步长安排：前 7200 s 用 1 s，其后用 10 s
#   C_TARGET    烘干判据阈值，0.15 kg/kg（干基含水率）
# ============================================================================

import os
import sys
import time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm
import compact_xlsx

N_CELL = 1600
T_END = 259200.0                      # 72 h，过程持续 2~3 天            [s]
DT_SCHEDULE = [(7200.0, 1.0), (T_END, 10.0)]   # 前 7200 s 用 1 s，其后 10 s
C_TARGET = 0.15                       # 烘干判据阈值                  [kg/kg]
HEADER = "时间\\到药材中心的距离"
FINE_RADII = [round(0.1 * k, 1) for k in range(21)]    # 结果文件半径 0~2.0 [cm]
TAB_RADII = [0.0, 0.5, 1.0, 1.5, 2.0]                  # 论文表输出半径    [cm]


def main():
    t_wall = time.time()
    print("=" * 74)
    print("PROBLEMS 2 and 3 -- whole drying process, Appendix 3 properties")
    print("=" * 74)

    room = cm.RoomConditions()           # 烘房温湿度历程，即表面的 Robin 边界
    t_fine = np.arange(1.0, T_END + 0.5, 1.0)
    xi_fine = np.array(FINE_RADII, dtype=float) / 2.0

    # track_max=True 时每个输出时刻都记录全场最大值 Cmax、轴心值 C_axis 等，
    # 问题 3 的判据必须用全场最值，而不能只在采样点上取最值。
    res = cm.simulate(N_CELL, cm.props_problem23, T_END, 1.0, room,
                      grid="node", t_report=t_fine, xi_report=xi_fine,
                      dt_schedule=DT_SCHEDULE, track_max=True, verbose=True)
    print("   simulation finished: %d steps, %.1f s wall"
          % (res["n_step"], time.time() - t_wall))

    # ------------------------------------------------------------------
    # 问题 2 的结果文件
    # ------------------------------------------------------------------
    out2 = os.path.join(cm.DATA_DIR, "result2.xlsx")
    # 距离列写成数值，与附件 3 模板一致（见 p1_preheat.py 的说明）
    hdr = [HEADER] + [float(r) for r in FINE_RADII]
    size2 = compact_xlsx.write_xlsx(out2, {
        "温度": (hdr, np.round(np.column_stack([t_fine, res["T"]]), 4)),
        "水分浓度": (hdr, np.round(np.column_stack([t_fine, res["C"]]), 4))})
    print("   result2.xlsx -> %s  (%d 行 x %d 列, %.2f MB)"
          % (out2, len(t_fine), len(hdr), size2 / 1e6))

    # 论文表 3、表 4：前 3 h 内每 0.5 h 取一个时刻
    t_tab = np.arange(1800.0, 10800.0 + 1.0, 1800.0)
    res_tab = cm.simulate(N_CELL, cm.props_problem23, 10800.0, 1.0, room,
                          grid="node", t_report=t_tab,
                          xi_report=np.array(TAB_RADII) / 2.0)
    tab_T = cm.result_frame(t_tab / 3600.0, TAB_RADII, res_tab["T"],
                            time_header="时间/h")
    tab_C = cm.result_frame(t_tab / 3600.0, TAB_RADII, res_tab["C"],
                            time_header="时间/h")
    tab_T.to_csv(os.path.join(cm.DATA_DIR, "p2_table3_T.csv"), index=False,
                 encoding="utf-8-sig")
    tab_C.to_csv(os.path.join(cm.DATA_DIR, "p2_table4_C.csv"), index=False,
                 encoding="utf-8-sig")

    # ------------------------------------------------------------------
    # 问题 3：烘干时间
    # ------------------------------------------------------------------
    C_max = res["Cmax"]                    # 全场逐时刻的最大含水率
    below = np.where(C_max < C_TARGET)[0]
    if len(below) == 0:
        raise RuntimeError("drying target not reached within %.0f s" % T_END)
    i0 = below[0]
    # 在最后一个高于阈值的采样点与第一个低于阈值的采样点之间线性插值，
    # 得到精确的穿越时刻
    t1, t0 = res["t"][i0], res["t"][i0 - 1]
    v1, v0 = C_max[i0], C_max[i0 - 1]
    t_dry = t0 + (C_TARGET - v0) * (t1 - t0) / (v1 - v0)
    # Cmax 是全场最大值、C_axis 是轴心值，故「中心最后干」这一说法等价于
    # Cmax == C_axis。而 argmax 的**下标**是另一回事：内部剖面还均匀到舍入
    # 精度时多个节点并列最大，下标只是其中任意一个，因此把下标落在轴心的
    # 次数单独报告，而不把它的比例当作该结论的强度。
    excess = (C_max - res["C_axis"]) / C_max
    on_axis = bool(np.all(excess <= 1e-12))
    n_idx_axis = int(np.sum(res["xi_at_max"] == 0))
    print("   first sample with max_r C < %.2f : t = %.0f s ; interpolated "
          "t_dry = %.1f s = %.4f h" % (C_TARGET, t1, t_dry, t_dry / 3600.0))
    print("   criterion evaluated on the FULL node field (%d nodes):"
          % (N_CELL + 1))
    print("   max_r C(r,t) - C(0,t) <= 1e-12*Cmax at every one of the %d output "
          "times: %s   (largest relative excess = %.3g)"
          % (len(C_max), on_axis, float(excess.max())))
    print("   argmax index is the axis at %d of them; at the other %d the "
          "interior is still uniform to round-off and the index is a tie"
          % (n_idx_axis, len(C_max) - n_idx_axis))

    # ---- result3：60 s x 0.1 cm，直到烘干结束时刻 ----------------------
    t_end3 = np.ceil(t_dry / 60.0) * 60.0
    t3 = np.arange(60.0, t_end3 + 0.5, 60.0)
    idx = np.searchsorted(res["t"], t3)
    C3 = res["C"][idx]
    out3 = os.path.join(cm.DATA_DIR, "result3.xlsx")
    compact_xlsx.write_xlsx(out3, {
        "Sheet1": (hdr, np.round(np.column_stack([t3, C3]), 4))})
    print("   result3.xlsx -> %s  (%d 行 x %d 列, t up to %.0f s = %.2f h)"
          % (out3, len(t3), len(hdr), t_end3, t_end3 / 3600.0))

    # ---- 论文表 5：每 6 h 一行的水分浓度 -------------------------------
    t_tab5 = list(np.arange(6.0, np.floor(t_dry / 3600.0 / 6.0) * 6.0 + 0.5, 6.0))
    t_tab5_s = np.array(t_tab5) * 3600.0
    idx5 = np.searchsorted(res["t"], t_tab5_s)
    # 输出半径共 21 个（步长 0.1 cm），取第 0/5/10/15/20 个即 0、0.5、1.0、
    # 1.5、2.0 cm 五列；末行补上烘干结束时刻所在的采样点
    C5 = res["C"][idx5][:, [0, 5, 10, 15, 20]]
    C5 = np.vstack([C5, res["C"][i0 - 1][[0, 5, 10, 15, 20]]])
    lab5 = ["%.1f" % v for v in t_tab5] + ["烘干结束(%.4f h)" % (t_dry / 3600.0)]
    tab5 = pd.DataFrame(np.round(C5, 4),
                        columns=["%.1f" % r for r in TAB_RADII])
    tab5.insert(0, "时间/h", lab5)
    tab5.to_csv(os.path.join(cm.DATA_DIR, "p3_table5_C.csv"), index=False,
                encoding="utf-8-sig")

    # ------------------------------------------------------------------
    # 保存绘图脚本需要的全部数组
    # ------------------------------------------------------------------
    np.savez_compressed(
        os.path.join(cm.DATA_DIR, "p23_fields.npz"),
        t=t_fine[::60], radii_cm=np.array(FINE_RADII),
        T=res["T"][::60].astype(np.float32),
        C=res["C"][::60].astype(np.float32),
        t_tab=t_tab, T_tab=res_tab["T"], C_tab=res_tab["C"],
        t_dry=np.array([t_dry]))

    # ---- 诊断信息：汇总写入 data/p23_summary.txt -----------------------
    lines = []
    lines.append("PROBLEMS 2-3 diagnostics")
    lines.append("cells = %d, schedule = %s, steps = %d, wall = %.1f s"
                 % (N_CELL, DT_SCHEDULE, res["n_step"], time.time() - t_wall))
    lines.append("drying time t_dry = %.4f s = %.4f h = %.4f d"
                 % (t_dry, t_dry / 3600.0, t_dry / 86400.0))
    lines.append("C at the centre and at the surface when the target is met:")
    lines.append("   C_center = %.6f, C_surface = %.6f, C_max = %.6f"
                 % (res["C"][i0, 0], res["C"][i0, -1], C_max[i0]))
    lines.append("")
    lines.append("   t/h    T_center  T_surf   C(r=0)  C(0.5)  C(1.0)  C(1.5)  "
                 "C(2.0)")
    for th in np.arange(0.0, T_END / 3600.0 + 0.5, 6.0):
        j = int(np.searchsorted(res["t"], th * 3600.0))
        j = min(j, len(res["t"]) - 1)
        lines.append("%8.1f %8.3f %8.3f  %8.4f %7.4f %7.4f %7.4f %7.4f"
                     % (res["t"][j] / 3600.0, res["T"][j, 0], res["T"][j, -1],
                        res["C"][j, 0], res["C"][j, 5], res["C"][j, 10],
                        res["C"][j, 15], res["C"][j, 20]))
    lines.append("")
    # 断面平均干基含水率 <C> = 2 * int_0^1 C(xi) * xi dxi：在 21 个输出半径上
    # 用梯形法积分，两端权重取半（梯形公式的端点处理）
    w_mean = 2.0 * xi_fine * (xi_fine[1] - xi_fine[0])
    w_mean[0] *= 0.5
    w_mean[-1] *= 0.5
    lines.append("mean moisture at t_dry = %.6f kg/kg"
                 % float(res["C"][i0] @ w_mean))
    with open(os.path.join(cm.DATA_DIR, "p23_summary.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("\n   total wall time %.1f s" % (time.time() - t_wall))


if __name__ == "__main__":
    main()

