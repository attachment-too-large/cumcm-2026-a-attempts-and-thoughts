# -*- coding: utf-8 -*-
# ============================================================================
# p4_shrinkage.py —— 问题 4：考虑收缩的烘干过程求解与烘干时间确定
#
# 做什么：药材半径由 2.000 cm 缩到 1.198 cm，求解域随之收缩。改用物质坐标
#         xi = r/R(t) 后，各向同性收缩下的对流项严格为零，方程退化为标准扩散式，
#         故与问题 1~3 共用同一数值内核。
#
# 输入：物性取题目附录 4；半径历史 R(t) 取附件 2 实测数据；烘房环境同问题 2、3。
#
# 输出：data/result4.xlsx               结果文件（60 s x 0.1 cm，物理坐标，末列为
#                                       「药材表面」，收缩后落在药材外的位置留空）
#       data/result4_material_coord.xlsx 物质坐标完整剖面（无空缺）
#       data/p4_table6_C.csv           论文表 6：每 6 h 的水分浓度
#       data/p4_fields.npz             完整场；data/p4_summary.txt 诊断量
#
# 关键变量：
#   N_CELL      径向控制体数目取 1600（节点中心格式下节点数为 N_CELL + 1）
#   T_END       模拟终了时刻 72 h                                     [s]
#   DT_SCHEDULE 步长安排：前 7200 s 用 1 s，其后用 10 s
#   C_TARGET    烘干判据阈值 0.15 kg/kg（干基含水率）
#   FINE_RADII  结果文件的物理半径列 0~2.0 [cm]；XI_MAT 是物质坐标 xi = 0~1
#   TAB_RADII   论文表 6 的物理半径列 0~1.5 [cm]
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
T_END = 259200.0                      # 72 h                            [s]
DT_SCHEDULE = [(7200.0, 1.0), (T_END, 10.0)]   # 前 7200 s 用 1 s，其后 10 s
C_TARGET = 0.15                       # 烘干判据阈值                  [kg/kg]
HEADER = "时间\\到药材中心的距离"
FINE_RADII = [round(0.1 * k, 1) for k in range(21)]   # 物理半径 0~2.0    [cm]
XI_MAT = [round(0.05 * k, 2) for k in range(21)]      # 物质坐标 xi = r/R，0~1
TAB_RADII = [0.0, 0.5, 1.0, 1.5]                      # 论文表 6 半径     [cm]


def main():
    t_wall = time.time()
    print("=" * 74)
    print("PROBLEM 4 -- shrinking herb, Appendix 4 properties")
    print("=" * 74)

    room = cm.RoomConditions()              # 烘房温湿度历程，与问题 2、3 同一台烘房
    radius = cm.RadiusHistory(mode="table")  # 附件 2 实测半径的分段线性插值
    print("   R(0) = %.4f cm, R(259200 s) = %.4f cm, ratio = %.4f"
          % (radius.R(0.0) * 100, radius.R(T_END) * 100,
             radius.R(T_END) / radius.R(0.0)))

    t_fine = np.arange(1.0, T_END + 0.5, 1.0)
    n_rad = len(FINE_RADII)

    # ------------------------------------------------------------------
    # 每个输出时刻采样一行，行内分三段（R(t) 收缩，物理半径在物质坐标上逐时刻
    # 移动）：[0:21] 固定 xi = 0.00~1.00 写物质坐标文件；[21:42] xi = r_j / R(t)
    # 写物理坐标文件；[42] xi = 1 为当前表面。落到药材之外的点把 xi 截断到 1，
    # 以免二次重构做很远的外推，这些位置随后由 r_j > R(t) 的判断置空。
    # ------------------------------------------------------------------
    xi_mat = np.asarray(XI_MAT, dtype=float)

    # 按上面的布局拼出该时刻的三段采样位置
    def xi_layout(t):
        R = float(radius.R(t))
        xi_phys = np.minimum(np.asarray(FINE_RADII, dtype=float) * 1.0e-2 / R,
                             1.0)
        return np.concatenate([xi_mat, xi_phys, [1.0]])

    res = cm.simulate(N_CELL, cm.props_problem4, T_END, 1.0, room,
                      radius=radius, grid="node", t_report=t_fine,
                      xi_report=xi_layout,
                      dt_schedule=DT_SCHEDULE, track_max=True, verbose=True)
    print("   simulation finished: %d steps, %.1f s wall"
          % (res["n_step"], time.time() - t_wall))

    C_mat = res["C"][:, :n_rad]                 # 物质坐标下的完整剖面
    C_phys = res["C"][:, n_rad:2 * n_rad]       # 在当前物理半径处采样得到的值
    C_surf = res["C"][:, -1]                    # 当前药材表面的值

    # ------------------------------------------------------------------
    # 烘干时间：判据作用在全部节点的场上
    # ------------------------------------------------------------------
    C_max = res["Cmax"]
    below = np.where(C_max < C_TARGET)[0]
    if len(below) == 0:
        raise RuntimeError("target not reached within %.0f s" % T_END)
    i0 = below[0]
    t1, t0 = res["t"][i0], res["t"][i0 - 1]
    v1, v0 = C_max[i0], C_max[i0 - 1]
    t_dry = t0 + (C_TARGET - v0) * (t1 - t0) / (v1 - v0)
    # Cmax 是全场最大值、C_axis 是轴心值，故「中心最后干」等价于 Cmax == C_axis。
    # 这里容易混两件事：excess = (Cmax - C_axis)/Cmax > 0 只在最大值偏离轴心时
    # 成立，这才是真正的论断，它在所有时刻都成立到机器精度；而 argmax 的下标
    # 是另一回事，内部剖面还均匀到舍入精度时（最初若干步）多个节点并列最大，
    # 下标只是其中任意一个，其落在轴心的比例并不能衡量该论断的强度。
    excess = (C_max - res["C_axis"]) / C_max
    on_axis = bool(np.all(excess <= 1e-12))
    n_idx_axis = int(np.sum(res["xi_at_max"] == 0))
    n_deg = len(C_max) - n_idx_axis
    print("   criterion evaluated on the FULL node field (%d nodes):"
          % (N_CELL + 1))
    print("   max_r C(r,t) - C(0,t) <= 1e-12*Cmax at every one of the %d output "
          "times: %s   (largest relative excess = %.3g)"
          % (len(C_max), on_axis, float(excess.max())))
    print("   argmax index is the axis at %d of them; at the other %d the "
          "interior is still uniform to round-off and the index is a tie"
          % (n_idx_axis, n_deg))
    print("   t_dry = %.1f s = %.4f h = %.4f d  (R = %.4f cm at that time)"
          % (t_dry, t_dry / 3600.0, t_dry / 86400.0,
             float(radius.R(t_dry)) * 100))

    # ------------------------------------------------------------------
    # result4.xlsx：各列是**当前**的物理距离 r
    # ------------------------------------------------------------------
    t_end4 = np.ceil(t_dry / 60.0) * 60.0
    t4 = np.arange(60.0, t_end4 + 0.5, 60.0)
    idx = np.searchsorted(res["t"], t4)
    Cp4, Cs4, R4 = C_phys[idx].copy(), C_surf[idx], res["R"][idx]

    block = np.column_stack([t4, Cp4, Cs4])
    # 距离列写成数值，与附件 3 模板一致（见 p1_preheat.py 的说明）；
    # 末列「药材表面」是模板给出的文字列名，保持字符串
    hdr4 = ([HEADER] + [float(r) for r in FINE_RADII] + ["药材表面"])
    for j, r_cm in enumerate(FINE_RADII):
        block[r_cm > R4 * 100.0 + 1e-12, j + 1] = np.nan   # 已超出药材范围
    out4 = os.path.join(cm.DATA_DIR, "result4.xlsx")
    compact_xlsx.write_xlsx(out4, {"Sheet1": (hdr4, block)})
    print("   result4.xlsx -> %s  (%d 行 x %d 列, 物理坐标)"
          % (out4, len(t4), len(hdr4)))

    # 附加文件：同一解在物质坐标下的剖面（21 列始终完整、无空缺）
    hdr_xi = [HEADER, "半径/cm"] + ["xi=%.2f" % v for v in XI_MAT]
    compact_xlsx.write_xlsx(
        os.path.join(cm.DATA_DIR, "result4_material_coord.xlsx"),
        {"Sheet1": (hdr_xi,
                    np.column_stack([t4, R4 * 100.0, C_mat[idx]]))})
    print("   result4_material_coord.xlsx -> 物质坐标完整剖面（21 列无空缺）")

    # ------------------------------------------------------------------
    # 论文表 6：每 6 h 一行，物理半径每 0.5 cm 一列
    # ------------------------------------------------------------------
    t6_h = np.arange(6.0,
                     np.floor(t_dry / 3600.0 / 6.0) * 6.0 + 0.5, 6.0)
    t6s = np.concatenate([t6_h * 3600.0, [t_dry]])
    idx6 = np.clip(np.searchsorted(res["t"], t6s), 0, len(res["t"]) - 1)
    R6 = res["R"][idx6] * 100.0                        # 这些时刻的半径 [cm]
    rows = []
    for k in range(len(t6s)):
        row = []
        for r_cm in TAB_RADII:
            if r_cm <= R6[k] + 1e-12:
                xi_t = r_cm / R6[k]                    # 换算成该时刻的物质坐标
                row.append(float(cm.interp_cells(C_mat[idx6[k]], xi_mat,
                                                 xi_t)[0]))
            else:
                row.append(np.nan)                     # 该点已不在药材内
        row.append(float(C_surf[idx6[k]]))             # 当前表面值
        rows.append(row)
    tab6 = pd.DataFrame(np.round(np.array(rows, dtype=float), 4),
                        columns=["%.1f" % r for r in TAB_RADII] + ["药材表面"])
    labels = ["%.1f" % v for v in t6_h] + \
             ["烘干结束(%.4f h, R=%.2f cm)" % (t_dry / 3600.0, R6[-1])]
    tab6.insert(0, "时间/h", labels)
    tab6.to_csv(os.path.join(cm.DATA_DIR, "p4_table6_C.csv"), index=False,
                encoding="utf-8-sig")

    # ------------------------------------------------------------------
    # 保存的场：所有数组共用同一时间轴（t_fine[::60]，每分钟一个采样），便于
    # 绘图脚本按下标对齐。C_phys 是在物理半径处采样得到的解（绘图时还会按已
    # 知的 R(t) 掩掉超出药材的列），而不是把 60 s 间隔的结果网格重新切片——
    # Cp4 定义在 t4 上，而 t4 到 t_dry 就结束，长度与前者不同。
    # ------------------------------------------------------------------
    np.savez_compressed(
        os.path.join(cm.DATA_DIR, "p4_fields.npz"),
        t=t_fine[::60], radii_cm=np.array(FINE_RADII),
        xi_mat=xi_mat,
        C=C_mat[::60].astype(np.float32),
        T=res["T"][::60, :n_rad].astype(np.float32),
        C_phys=res["C"][::60, n_rad:2 * n_rad].astype(np.float32),
        R=res["R"][::60].astype(np.float32),
        t_dry=np.array([t_dry]))
    assert (C_mat[::60].shape[0] == res["R"][::60].shape[0]
            == res["T"][::60, :n_rad].shape[0] == t_fine[::60].shape[0]), \
        "all saved arrays must share one time axis"

    lines = []
    lines.append("PROBLEM 4 diagnostics")
    lines.append("cells = %d, schedule = %s, steps = %d, wall = %.1f s"
                 % (N_CELL, DT_SCHEDULE, res["n_step"], time.time() - t_wall))
    lines.append("R(0) = %.4f cm, R(14400 s) = %.4f cm, R(259200 s) = %.4f cm"
                 % (radius.R(0.0) * 100, radius.R(14400.0) * 100,
                    radius.R(259200.0) * 100))
    lines.append("drying time t_dry = %.4f s = %.4f h = %.4f d"
                 % (t_dry, t_dry / 3600.0, t_dry / 86400.0))
    lines.append("criterion evaluated on the FULL node field (n+1 = %d nodes);"
                 % (N_CELL + 1))
    lines.append("max_r C(r,t) - C(0,t) <= 1e-12*Cmax at every one of the %d "
                 "output times: %s (largest relative excess %.3g)"
                 % (len(C_max), on_axis, float(excess.max())))
    lines.append("the argmax INDEX is the axis at %d of them; at the other %d "
                 "the interior is still uniform to round-off so the index is a "
                 "floating-point tie" % (n_idx_axis, n_deg))
    lines.append("surface radius at t_dry = %.4f cm" % (R6[-1]))
    lines.append("")
    lines.append("(a) MATERIAL coordinate xi = r/R(t) -- the complete profile:")
    lines.append("   t/h    R/cm   T_center T_surf  C(xi=0) C(.25)  C(.5)   "
                 "C(.75)  C(xi=1)")
    for th in np.arange(0.0, t_dry / 3600.0 + 0.5, 6.0):
        j = int(np.searchsorted(res["t"], th * 3600.0))
        j = min(j, len(res["t"]) - 1)
        lines.append("%8.1f %7.3f %8.3f %8.3f %8.4f %7.4f %7.4f %7.4f %7.4f"
                     % (res["t"][j] / 3600.0, res["R"][j] * 100,
                        res["T"][j, 0], res["T"][j, n_rad - 1],
                        C_mat[j, 0], C_mat[j, 5], C_mat[j, 10],
                        C_mat[j, 15], C_mat[j, 20]))
    lines.append("")
    lines.append("(b) PHYSICAL radius r/cm -- what result4.xlsx and Table 6 hold;")
    lines.append("    a dash means the point is already outside the material:")
    lines.append("   t/h    R/cm  " + "".join("%8.1f" % r
                                              for r in TAB_RADII) + "  surface")
    for k in range(len(t6s)):
        vals = "".join("       -" if not np.isfinite(v) else "%8.4f" % v
                       for v in rows[k][:len(TAB_RADII)])
        lines.append("%8.1f %7.3f %s %8.4f"
                     % (t6s[k] / 3600.0, R6[k], vals, rows[k][-1]))
    with open(os.path.join(cm.DATA_DIR, "p4_summary.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("\n   total wall time %.1f s" % (time.time() - t_wall))


if __name__ == "__main__":
    main()
