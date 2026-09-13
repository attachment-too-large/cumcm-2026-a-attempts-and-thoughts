# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s03_solve_p23.py
# 作用  : 问题 2 与问题 3 求解(附录 3 物性, 同一求解器)。
#         问题 2: 预热平衡 + 恒温干燥全程模型, 给出 3 h 内每 0.5 h 的结果;
#         问题 3: 同一模型推进到各处含水率均低于 0.15 kg/kg, 给出干燥时间。
#         环境: t<=14400 s 用附件 1 线性插值; t>14400 s 用平台假定(基准取末段
#               3600 s 实测均值), 平台取值本身作为假设在 s06 做灵敏度分析。
#         产出:
#           results/result2.xlsx / result3.xlsx
#           data/p2_table3_T.csv / p2_table4_C.csv / p3_table5_C.csv
#           data/p23_dt_convergence.csv / p23_mesh_convergence.csv
#           data/p23_summary.txt
# 运行  : python s03_solve_p23.py
# =============================================================================
"""Problems 2 and 3: full drying model with the Appendix-3 correlations."""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drying_core as dc      # noqa: E402
import scenarios as sc        # noqa: E402

N_BASE = 1600
T_END_P2 = 10800.0            # 3 h
T_END_P3 = 400000.0           # 上限, 实际由 c_stop 提前结束
DT_P2_BASE = 0.1              # P2 基准时间步 [s]
DT_P3_BASE = 1.0              # P3 基准时间步 [s]
DT_STUDY_P2 = (5.0, 2.0, 1.0, 0.5, 0.25, 0.1)
DT_STUDY_P3 = (60.0, 30.0, 10.0, 5.0, 2.0, 1.0)
DT_MESH = 10.0
GRIDS = (200, 400, 800, 1600, 3200)

TABLE3_TIMES_H = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
TABLE5_TIMES_H = np.array([6.0, 12.0, 18.0, 24.0, 30.0, 36.0, 42.0, 48.0, 54.0])
POS_TABLE_CM = np.array([0.0, 0.5, 1.0, 1.5, 2.0])


def run_p2(n_cells, dt, env=None, xi=None, t_end=T_END_P2, record_dt=None):
    """问题 2 求解。"""
    env = env if env is not None else sc.load_env()
    return dc.simulate(n_cells, t_end, dt, dc.props_appendix3, env,
                       sc.fixed_radius(), theta=1.0, record_dt=record_dt,
                       record_xi=xi)


def run_p3(n_cells, dt, env=None, xi=None, record_dt=None):
    """问题 3 求解(以 max C < 0.15 为停止判据)。"""
    env = env if env is not None else sc.load_env()
    return dc.simulate(n_cells, T_END_P3, dt, dc.props_appendix3, env,
                       sc.fixed_radius(), theta=1.0, record_dt=record_dt,
                       record_xi=xi, c_stop=dc.C_DRY)


def main():
    lines = []
    push = lines.append
    env = sc.load_env()
    pos_cm = sc.positions_cm_fixed()
    xi = pos_cm * 1e-2 / dc.R0

    push("=" * 78)
    push("烘房环境")
    push("=" * 78)
    push(f"附件 1 覆盖 0~{env.t_last:.0f} s; 之后采用平台假定({env.T_plateau:.4f} degC, "
         f"{env.C_plateau:.6f} kg/kg = 末段 3600 s 实测均值)")
    push("该平台值是模型假设而非数据唯一推论, 其影响在 s06 灵敏度分析中量化。")

    # =====================================================================
    # 问题 2
    # =====================================================================
    push("")
    push("=" * 78)
    push("问题 2: 3 h 内结果(附录 3 物性)")
    push("=" * 78)
    dt_rows = []
    for dt in DT_STUDY_P2:
        r = run_p2(N_BASE, dt, env=env, xi=xi, record_dt=1.0)
        cg, tg = sc.read_position_series(r, pos_cm, surface_last=True)
        times_h = r.times / 3600.0
        sel = [int(np.argmin(np.abs(times_h - h))) for h in TABLE3_TIMES_H]
        dt_rows.append({"dt/s": dt,
                        "T(r=0,3h)": tg[-1, 0], "T(surf,3h)": tg[-1, -1],
                        "C(r=0,3h)": cg[-1, 0], "C(surf,3h)": cg[-1, -1],
                        "T(r=0,0.5h)": tg[sel[0], 0], "C(r=0,0.5h)": cg[sel[0], 0]})
    dtab2 = pd.DataFrame(dt_rows)
    push("时间步收敛性(3 h 结果):")
    push(dtab2.to_string(index=False, float_format=lambda v: f"{v:.6f}"))
    dtab2.to_csv(os.path.join(sc.DIR_DATA, "p2_dt_convergence.csv"),
                 index=False, encoding="utf-8-sig")
    vals = dtab2["C(r=0,3h)"].to_numpy(float)
    e_last = abs(vals[-2] - vals[-1])
    e_prev = abs(vals[-3] - vals[-2])
    dts2 = dtab2['dt/s'].to_numpy(float)
    p_obs = (np.log(e_prev / e_last) / np.log((dts2[-3] - dts2[-2]) / (dts2[-2] - dts2[-1]))
             if e_last > 0 and e_prev > 0 else np.nan)
    ext = vals[-1] - e_last / (((dts2[-3] - dts2[-2]) / (dts2[-2] - dts2[-1])) ** p_obs - 1.0)
    push(f"含水率中心值: 观测收敛阶 p = {p_obs:.2f}; dt->0 外推 = {ext:.6f}; "
         f"基准 dt={DT_P2_BASE} s 的偏差 = {abs(vals[-1]-ext):.2e}")

    res2 = run_p2(N_BASE, DT_P2_BASE, env=env, xi=xi, record_dt=1.0)
    C2, T2 = sc.read_position_series(res2, pos_cm, surface_last=True)
    t2 = res2.times
    keep = t2 >= 1.0 - 1e-9
    sc.write_result_xlsx(os.path.join(sc.DIR_RES, "result2.xlsx"), t2[keep], pos_cm,
                         [T2[keep], C2[keep]], ["温度", "水分浓度"], surface_last=True)
    push(f"已写出 results/result2.xlsx: {int(keep.sum())} 行 x {len(pos_cm)} 列 x 2 个表")

    idx_p = [int(np.argmin(np.abs(pos_cm - pp))) for pp in POS_TABLE_CM]
    idx_h = [int(np.argmin(np.abs(t2 / 3600.0 - h))) for h in TABLE3_TIMES_H]
    tab3 = pd.DataFrame(T2[np.ix_(idx_h, idx_p)], index=TABLE3_TIMES_H,
                        columns=[f"{v:g} cm" for v in POS_TABLE_CM])
    tab4 = pd.DataFrame(C2[np.ix_(idx_h, idx_p)], index=TABLE3_TIMES_H,
                        columns=[f"{v:g} cm" for v in POS_TABLE_CM])
    tab3.index.name = "时间/h"
    tab4.index.name = "时间/h"
    tab3.round(4).to_csv(os.path.join(sc.DIR_DATA, "p2_table3_T.csv"), encoding="utf-8-sig")
    tab4.round(4).to_csv(os.path.join(sc.DIR_DATA, "p2_table4_C.csv"), encoding="utf-8-sig")
    push("")
    push("表 3 3 小时内药材的温度 / degC")
    push(tab3.round(4).to_string())
    push("")
    push("表 4 3 小时内药材的水分浓度 / (kg/kg)")
    push(tab4.round(4).to_string())

    # =====================================================================
    # 问题 3
    # =====================================================================
    push("")
    push("=" * 78)
    push("问题 3: 干燥时间(max C < 0.15 kg/kg)")
    push("=" * 78)
    dt3_rows = []
    for dt in DT_STUDY_P3:
        r = run_p3(N_BASE, dt, env=env, xi=xi, record_dt=None)
        td = sc.drying_time(r)
        dt3_rows.append({"dt/s": dt, "n_steps": r.n_steps,
                         "t_dry/s": td, "t_dry/h": td / 3600.0,
                         "wall/s": r.wall_seconds})
    dtab3 = pd.DataFrame(dt3_rows)
    push("时间步收敛性:")
    push(dtab3.to_string(index=False, float_format=lambda v: f"{v:.6f}"))

    mesh3_rows = []
    for n in GRIDS:
        r = run_p3(n, DT_MESH, env=env, xi=xi, record_dt=None)
        td = sc.drying_time(r)
        mesh3_rows.append({"N": n, "t_dry/s": td, "t_dry/h": td / 3600.0,
                           "wall/s": r.wall_seconds})
    mtab3 = pd.DataFrame(mesh3_rows)
    push("网格收敛性(统一 dt = 10 s):")
    push(mtab3.to_string(index=False, float_format=lambda v: f"{v:.6f}"))
    pd.concat([dtab3.assign(kind="dt"), mtab3.rename(columns={"N": "dt/s"}).assign(kind="mesh")],
              ignore_index=True).to_csv(
        os.path.join(sc.DIR_DATA, "p23_convergence.csv"), index=False, encoding="utf-8-sig")

    res3 = run_p3(N_BASE, DT_P3_BASE, env=env, xi=xi, record_dt=60.0)
    td3 = sc.drying_time(res3)
    push(f"基准算例: N={N_BASE}, dt={DT_P3_BASE} s, 步数={res3.n_steps}, "
         f"墙钟={res3.wall_seconds:.1f} s")
    push(f"*** 问题 3 干燥时间 t3 = {td3:.1f} s = {td3/3600.0:.4f} h ***")
    push(f"末态 max C = {res3.step_maxC[-1]:.6f}, 中心 C = "
         f"{float(dc.apply_weights(res3.C_final, *dc.lagrange_weights_grid(N_BASE,[0.0]))[0]):.6f}")

    C3, T3 = sc.read_position_series(res3, pos_cm, surface_last=True)
    t3 = res3.times
    # 只保留 60 s 网格且不晚于干燥终点的行
    keep3 = (t3 >= 60.0 - 1e-9) & (t3 <= td3 + 1e-9)
    sc.write_result_xlsx(os.path.join(sc.DIR_RES, "result3.xlsx"), t3[keep3], pos_cm,
                         [C3[keep3]], ["水分浓度"], surface_last=True)
    push(f"已写出 results/result3.xlsx: {int(keep3.sum())} 行 x {len(pos_cm)} 列")

    # 表 5: 每 6 h + 烘干结束时间
    idx5 = [int(np.argmin(np.abs(t3 / 3600.0 - h))) for h in TABLE5_TIMES_H]
    idx5 = [i for i in idx5 if t3[i] <= td3]
    times5 = list(t3[idx5] / 3600.0)
    tab5 = pd.DataFrame(C3[np.ix_(idx5, idx_p)], index=[f"{h:.0f}" for h in times5],
                        columns=[f"{v:g} cm" for v in POS_TABLE_CM])
    # 烘干结束时间行: 在最后两条 60 s 记录之间线性插值
    i_last = int(np.searchsorted(t3, td3))
    i_last = min(max(i_last, 1), len(t3) - 1)
    w = (td3 - t3[i_last - 1]) / (t3[i_last] - t3[i_last - 1])
    c_end = C3[i_last - 1] * (1 - w) + C3[i_last] * w
    c_end[0] = dc.C_DRY          # 判据保证中心恰为 0.15
    tab5.loc["烘干结束时间"] = c_end[idx_p]
    tab5.index.name = "时间/h"
    tab5.round(4).to_csv(os.path.join(sc.DIR_DATA, "p3_table5_C.csv"), encoding="utf-8-sig")
    push("")
    push("表 5 药材烘干过程的水分浓度 / (kg/kg)")
    push(tab5.round(4).to_string())
    push(f"(烘干结束时间 = {td3/3600.0:.4f} h)")

    # --------------------------------------------------- 守恒与物理合理性
    push("-" * 78)
    push("守恒性与物理合理性检查(问题 3)")
    push("-" * 78)
    w0, wf, flux = res3.water[0], res3.water[-1], res3.cum_flux[-1]
    push(f"总水量 {w0:.10f} -> {wf:.10f}; 累计逸出 {flux:.10e}; "
         f"相对闭合误差 {abs((w0-wf)-flux)/(w0-wf)*100:.6f} %")
    push(f"max(C) 单调不增 = {bool(np.all(np.diff(res3.step_maxC) <= 1e-12))}")
    push(f"C 始终非负      = {bool(np.all(res3.step_minC > 0))} (最小 {res3.step_minC.min():.6f})")
    push(f"全场温度从不越过环境温度上界 {res3.step_Tair.max():.4f} degC: "
         f"{bool(res3.step_maxT.max() <= res3.step_Tair.max() + 1e-6)} "
         f"(实际最高 {res3.step_maxT.max():.4f} degC)")
    push(f"全场温度从不低于初始温度 {dc.T_INIT_C} degC: "
         f"{bool(res3.step_minT.min() >= dc.T_INIT_C - 1e-6)} "
         f"(实际最低 {res3.step_minT.min():.4f} degC)")
    push(f"末态中心温度 {res3.T_final[0]:.4f} degC, 表面温度 {res3.t_surf_final:.4f} degC")

    text = "\n".join(lines)
    with open(os.path.join(sc.DIR_DATA, "p23_summary.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
