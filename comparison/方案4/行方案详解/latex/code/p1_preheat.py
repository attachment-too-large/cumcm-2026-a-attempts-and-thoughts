# -*- coding: utf-8 -*-
# ============================================================================
# p1_preheat.py —— 问题 1：预热平衡阶段（0 ~ 1800 s）的求解与结果输出
#
# 做什么：按一维轴对称热-质扩散模型求解药材内部的温度场 T(r,t) 与水分浓度场
#         C(r,t)。物性取常数、烘房环境取实测序列，与问题 2~4 共用同一求解内核。
#
# 输入：物性取题目附录 2（rho = 820 kg/m3、cp = 2600 J/(kg K)、k = 0.36 W/(m K)、
#       D = 7e-9 exp(-0.89/C) m2/s）；烘房温湿度取附件 1 的实测序列。
#
# 输出：data/result1.xlsx      结果文件（工作表「温度」「水分浓度」，1 s x 0.1 cm）
#       data/p1_table1_T.csv   论文表 1：若干时刻、若干半径处的温度
#       data/p1_table2_C.csv   论文表 2：同上，但为水分浓度
#       data/p1_fields.npz     完整场，供绘图脚本读取
#       data/p1_summary.txt    诊断量：Biot 数、特征时间、失水量等
#
# 关键变量：
#   N_CELL  径向控制体数目，取 1600（节点中心格式下节点数为 N_CELL + 1）
#   DT      时间步长 [s]
#   T_END   模拟终了时刻，即预热平衡阶段结束时刻 [s]
#   HEADER  结果表左上角的表头文字，与附件 3 模板一致
# ============================================================================

import os
import sys
import time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm
import compact_xlsx

N_CELL = 1600          # 径向控制体数目
DT = 1.0               # 时间步长                                          [s]
T_END = 1800.0         # 预热平衡阶段结束时刻，即模拟终了时刻              [s]
HEADER = "时间\\到药材中心的距离"
TAB_TIMES = [100, 300, 600, 900, 1200, 1500, 1800]     # 论文表输出时刻 [s]
TAB_RADII = [0.0, 0.5, 1.0, 1.5, 2.0]                  # 论文表输出半径 [cm]
FINE_RADII = [round(0.1 * k, 1) for k in range(21)]    # 结果文件半径 0~2.0 [cm]


def main():
    t_wall = time.time()
    print("=" * 74)
    print("PROBLEM 1  preheat-equilibrium stage, 0 ... 1800 s")
    print("=" * 74)

    room = cm.RoomConditions()           # 烘房温湿度历程，即表面的 Robin 边界
    diag = cm.characteristic_numbers()   # 无量纲数与特征时间
    print("   alpha_heat = %.6e m2/s   Bi_heat = %.4f   tau_heat = %.1f s"
          % (diag["alpha"], diag["bi_heat"], diag["tau_heat"]))
    print("   D(C0)      = %.6e m2/s   Bi_mass = %.4f   tau_mass = %.1f s"
          % (diag["D_ref"], diag["bi_mass"], diag["tau_mass"]))

    # ---- 论文表格：在选定的若干时刻与半径处采样 -------------------------
    t_tab = np.array(TAB_TIMES, dtype=float)
    xi_tab = np.array(TAB_RADII, dtype=float) / 2.0   # 半径列以 cm 给出，R0=2 cm
    res_tab = cm.simulate(N_CELL, cm.props_problem1, T_END, DT, room,
                          grid="node", t_report=t_tab, xi_report=xi_tab)

    # ---- 结果文件：1 s x 0.1 cm 的完整输出网格 --------------------------
    t_fine = np.arange(1.0, T_END + 0.5, 1.0)
    xi_fine = np.array(FINE_RADII, dtype=float) / 2.0
    res_fine = cm.simulate(N_CELL, cm.props_problem1, T_END, DT, room,
                           grid="node", t_report=t_fine, xi_report=xi_fine)

    # ---- 结果文件 data/result1.xlsx ------------------------------------
    out = os.path.join(cm.DATA_DIR, "result1.xlsx")
    # 距离列写成数值（0, 0.1, ..., 2）而不是字符串 "0.0"/"0.1"，与附件 3
    # 模板的第 1 行一致：模板中这些单元格是数值，字符串表头会让按数值读取
    # 结果文件的程序对不上列。
    hdr = [HEADER] + [float(r) for r in FINE_RADII]
    size1 = compact_xlsx.write_xlsx(out, {
        "温度": (hdr, np.round(np.column_stack([t_fine, res_fine["T"]]), 4)),
        "水分浓度": (hdr, np.round(np.column_stack([t_fine, res_fine["C"]]), 4))})
    print("   result1.xlsx  ->  %s   (%d 行 x %d 列, %.2f MB)"
          % (out, len(t_fine), len(hdr), size1 / 1e6))

    # ---- 论文表 1、表 2：选定时刻与半径处的温度、水分浓度 ---------------
    tab_T = cm.result_frame(t_tab, TAB_RADII, res_tab["T"])
    tab_C = cm.result_frame(t_tab, TAB_RADII, res_tab["C"])
    tab_T.to_csv(os.path.join(cm.DATA_DIR, "p1_table1_T.csv"), index=False,
                 encoding="utf-8-sig")
    tab_C.to_csv(os.path.join(cm.DATA_DIR, "p1_table2_C.csv"), index=False,
                 encoding="utf-8-sig")

    # ---- 保存完整场，供绘图脚本读取 -------------------------------------
    np.savez_compressed(
        os.path.join(cm.DATA_DIR, "p1_fields.npz"),
        t_fine=t_fine, radii_cm=np.array(FINE_RADII),
        T_fine=res_fine["T"], C_fine=res_fine["C"],
        t_tab=t_tab, radii_tab=np.array(TAB_RADII),
        T_tab=res_tab["T"], C_tab=res_tab["C"])

    # ---- 诊断信息：汇总写入 data/p1_summary.txt -------------------------
    lines = []
    lines.append("PROBLEM 1 diagnostics")
    lines.append("alpha_heat = %.6e m2/s, Bi_heat = %.6f, tau_heat = %.1f s"
                 % (diag["alpha"], diag["bi_heat"], diag["tau_heat"]))
    lines.append("D(C0) = %.6e m2/s, Bi_mass = %.6f, tau_mass = %.1f s"
                 % (diag["D_ref"], diag["bi_mass"], diag["tau_mass"]))
    lines.append("cells = %d, dt = %.3f s, wall time = %.2f s"
                 % (N_CELL, DT, time.time() - t_wall))

    C0 = res_tab["C"][:, 0]      # 各输出时刻的**中心**（xi=0）含水率
    Cs = res_tab["C"][:, -1]     # 各输出时刻的**表面**（xi=1）含水率
    T0 = res_tab["T"][:, 0]      # 中心温度
    Ts = res_tab["T"][:, -1]     # 表面温度
    lines.append("")
    lines.append("        t/s   T_center  T_surf   C_center  C_surf   "
                 "T_air    C_air")
    for j, tt in enumerate(t_tab):
        lines.append("%10.0f %9.4f %9.4f %9.4f %9.4f %8.3f %8.5f"
                     % (tt, T0[j], Ts[j], C0[j], Cs[j],
                        float(room.T_air(tt)), float(room.C_air(tt))))

    # 断面平均干基含水率 <C> = 2 * int_0^1 C(xi) * xi dxi：在 21 个输出半径上
    # 用梯形法积分，两端权重取半（梯形公式的端点处理）。
    xi_pts = np.array(FINE_RADII, dtype=float) / 2.0
    d_xi = xi_pts[1] - xi_pts[0]
    w_vol = 2.0 * xi_pts * d_xi          # 被积函数里的 2*xi*dxi
    w_vol[0] *= 0.5
    w_vol[-1] *= 0.5
    C_mean = res_fine["C"] @ w_vol
    T_mean = res_fine["T"] @ w_vol

    lines.append("")
    lines.append("surface temperature rise  : %.4f -> %.4f C" % (Ts[0], Ts[-1]))
    lines.append("centre  temperature rise  : %.4f -> %.4f C" % (T0[0], T0[-1]))
    lines.append("surface moisture drop     : %.4f -> %.4f kg/kg"
                 % (Cs[0], Cs[-1]))
    lines.append("centre  moisture drop     : %.4f -> %.4f kg/kg"
                 % (C0[0], C0[-1]))
    lines.append("cross-section mean moisture: %.4f -> %.4f kg/kg"
                 % (C_mean[0], C_mean[-1]))
    lines.append("cross-section mean temp    : %.4f -> %.4f C"
                 % (T_mean[0], T_mean[-1]))
    dry_density = 820.0 / (1.0 + cm.C_INIT)        # 单位湿体积中的干物质质量
    herb_volume = np.pi * cm.R0 ** 2 * cm.HERB_LENGTH
    water_removed = dry_density * herb_volume * (C_mean[0] - C_mean[-1]) * 1000.0
    lines.append("water removed in 1800 s    : %.3f g (per 25 cm stem)"
                 % water_removed)
    lines.append("initial water content      : %.3f g"
                 % (dry_density * herb_volume * cm.C_INIT * 1000.0))

    # 无量纲检验：温度与水分在 1800 s 内的 Fourier 数
    lines.append("Fo_heat(1800 s) = %.4f ; Fo_mass(1800 s) with D(C0) = %.4f"
                 % (diag["alpha"] * 1800.0 / cm.R0 ** 2,
                    diag["D_ref"] * 1800.0 / cm.R0 ** 2))
    lines.append("number of time steps = %d" % res_fine["n_step"])

    with open(os.path.join(cm.DATA_DIR, "p1_summary.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("\n   wall time %.2f s" % (time.time() - t_wall))


if __name__ == "__main__":
    main()

