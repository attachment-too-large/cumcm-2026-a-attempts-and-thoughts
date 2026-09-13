# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s02_solve_p1.py
# 作用  : 问题 1 求解(预热平衡阶段, 0~1800 s)。
#         模型: 附录 2 常数物性; 烘房环境取附件 1 原始数据线性插值;
#               半单元表面重构; 隐式 Euler(theta=1) + Picard 迭代。
#         产出:
#           results/result1.xlsx          温度 / 水分浓度 两表, 1 s x 0.1 cm 全网格
#           data/p1_table1_T.csv          表 1(温度, 7 个时刻 x 5 个位置)
#           data/p1_table2_C.csv          表 2(水分浓度)
#           data/p1_mesh_convergence.csv  网格收敛证据
#           data/p1_dt_convergence.csv    时间步收敛证据(含 Richardson 外推)
#           data/p1_summary.txt           结果摘要 + 守恒性检查
# 运行  : python s02_solve_p1.py
# =============================================================================
"""Problem 1: preheat stage on the fixed cylinder (0-1800 s)."""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drying_core as dc      # noqa: E402
import scenarios as sc        # noqa: E402

N_BASE = 1600                 # 基准网格(网格收敛表: N=1600 与 N=3200 偏差 ~3e-6)
DT_BASE = 0.02                # 基准时间步 [s](隐式 Euler 为一阶, 见 dt 收敛表)
T_END = 1800.0
TABLE_TIMES = np.array([100.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0])
TABLE_POS_CM = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
GRIDS = (200, 400, 800, 1600, 3200)
DTS = (10.0, 5.0, 2.0, 1.0, 0.5, 0.25, 0.1, 0.05, 0.02)
DT_MESH = 0.1                 # 网格收敛研究统一采用的时间步(差值中 dt 误差相消)


def run(n_cells, dt, record_dt=None, env=None, xi=None, t_end=T_END):
    """按给定离散配置求解问题 1(采样模式记录)。"""
    env = env if env is not None else sc.load_env()
    return dc.simulate(n_cells, t_end, dt, dc.props_appendix2, env,
                       sc.fixed_radius(), theta=1.0, record_dt=record_dt,
                       record_xi=xi)


def main():
    lines = []
    push = lines.append
    env = sc.load_env()

    pos_cm = sc.positions_cm_fixed()
    xi = pos_cm * 1e-2 / dc.R0
    idx_t = None

    # ---------------------------------------------------------------- 基准算例
    res = run(N_BASE, DT_BASE, record_dt=1.0, env=env, xi=xi)
    times = res.times
    C_grid, T_grid = sc.read_position_series(res, pos_cm, surface_last=True)
    push("=" * 78)
    push("问题 1 基准算例")
    push("=" * 78)
    push(f"网格 N = {N_BASE}, 时间步 dt = {DT_BASE} s, 步数 = {res.n_steps}, "
         f"平均 Picard 迭代 = {res.picard_iters.mean():.2f}, 墙钟 = {res.wall_seconds:.2f} s")
    push(f"输出点数 = {res.times.size}(每 1 s), 位置数 = {len(pos_cm)}(0~2 cm 每 0.1 cm)")

    keep = times >= 1.0 - 1e-9
    t_out = times[keep]
    sc.write_result_xlsx(os.path.join(sc.DIR_RES, "result1.xlsx"), t_out, pos_cm,
                         [T_grid[keep], C_grid[keep]], ["温度", "水分浓度"],
                         surface_last=True)
    push(f"已写出 results/result1.xlsx: {len(t_out)} 行 x {len(pos_cm)} 列 x 2 个表")

    # ------------------------------------------------------------ 表 1 / 表 2
    idx_t = [int(np.argmin(np.abs(times - tt))) for tt in TABLE_TIMES]
    idx_p = [int(np.argmin(np.abs(pos_cm - pp))) for pp in TABLE_POS_CM]
    tab_T = pd.DataFrame(T_grid[np.ix_(idx_t, idx_p)], index=TABLE_TIMES.astype(int),
                         columns=[f"{v:g} cm" for v in TABLE_POS_CM])
    tab_C = pd.DataFrame(C_grid[np.ix_(idx_t, idx_p)], index=TABLE_TIMES.astype(int),
                         columns=[f"{v:g} cm" for v in TABLE_POS_CM])
    tab_T.index.name = "时间/s"
    tab_C.index.name = "时间/s"
    tab_T.round(4).to_csv(os.path.join(sc.DIR_DATA, "p1_table1_T.csv"),
                          encoding="utf-8-sig")
    tab_C.round(4).to_csv(os.path.join(sc.DIR_DATA, "p1_table2_C.csv"),
                          encoding="utf-8-sig")
    push("")
    push("表 1 30 分钟内药材的温度 / degC")
    push(tab_T.round(4).to_string())
    push("")
    push("表 2 30 分钟内药材的水分浓度 / (kg/kg)")
    push(tab_C.round(4).to_string())

    # -------------------------------------------------------- 网格收敛性证据
    push("")
    push("-" * 78)
    push("网格收敛性: 表 1/表 2 的全部取值相对最细网格(N=3200)的最大偏差")
    push("-" * 78)
    store = {}
    for n in GRIDS:
        r = run(n, DT_MESH, record_dt=1.0, env=env, xi=xi)
        cg, tg = sc.read_position_series(r, pos_cm, surface_last=True)
        store[n] = (tg[np.ix_(idx_t, idx_p)], cg[np.ix_(idx_t, idx_p)])
    ref_T, ref_C = store[GRIDS[-1]]
    conv_rows = []
    for n in GRIDS:
        vT, vC = store[n]
        conv_rows.append({"N": n,
                          "max|dT|": np.max(np.abs(vT - ref_T)),
                          "max|dC|": np.max(np.abs(vC - ref_C)),
                          "C_surf(1800s)": vC[-1, -1],
                          "T_center(1800s)": vT[-1, 0]})
    conv = pd.DataFrame(conv_rows)
    push(conv.to_string(index=False, float_format=lambda v: f"{v:.3e}"))
    conv.to_csv(os.path.join(sc.DIR_DATA, "p1_mesh_convergence.csv"),
                index=False, encoding="utf-8-sig")
    push(f"N=1600 与 N=3200 的最大偏差 {conv['max|dC|'].iloc[-2]:.2e}(含水率) / "
         f"{conv['max|dT|'].iloc[-2]:.2e}(温度), 已远小于四位小数, 故基准取 N=1600。")

    # ------------------------------------------------------ 时间步收敛性证据
    push("-" * 78)
    push("时间步收敛性(隐式 Euler 一阶): 中心与表面值")
    push("-" * 78)
    dt_rows = []
    for dt in DTS:
        r = run(N_BASE, dt, record_dt=1.0, env=env, xi=xi)
        cg, tg = sc.read_position_series(r, pos_cm, surface_last=True)
        dt_rows.append({"dt/s": dt,
                        "T_center": tg[-1, 0], "T_surface": tg[-1, -1],
                        "C_center": cg[-1, 0], "C_surface": cg[-1, -1],
                        "n_steps": r.n_steps})
    dtab = pd.DataFrame(dt_rows)
    push(dtab.to_string(index=False, float_format=lambda v: f"{v:.6f}"))

    # Richardson: dt 逐次减半, 一阶格式 x(dt) = E + C*dt, 故 E = 2*x(dt/2) - x(dt)
    push("-" * 78)
    push("Richardson 外推(dt -> 0)")
    push("-" * 78)
    rich = {}
    for col in ("T_center", "T_surface", "C_center", "C_surface"):
        vals = dtab[col].to_numpy(float)
        e_last = abs(vals[-2] - vals[-1])
        e_prev = abs(vals[-3] - vals[-2])
        dts = dtab['dt/s'].to_numpy(float)
        # observed order: p = ln(e_prev/e_last)/ln(dt_prev/dt_last), no halving assumed
        p_obs = (np.log(e_prev / e_last) / np.log((dts[-3] - dts[-2]) / (dts[-2] - dts[-1]))
                 if e_last > 0 and e_prev > 0 else np.nan)
        ext = vals[-1] - e_last / (((dts[-3] - dts[-2]) / (dts[-2] - dts[-1])) ** p_obs - 1.0)
        rich[col] = {"observed_order": p_obs, "extrapolated": ext,
                     "base_dt0.02": vals[-1], "bias_at_base": abs(vals[-1] - ext)}
        push(f"{col:10s}: 观测收敛阶 p = {p_obs:.2f}, dt->0 极限 = {ext:.6f}, "
             f"基准 dt=0.02 s 的偏差 = {abs(vals[-1]-ext):.2e}")
    pd.DataFrame(rich).T.to_csv(os.path.join(sc.DIR_DATA, "p1_richardson.csv"),
                                encoding="utf-8-sig")
    dtab.to_csv(os.path.join(sc.DIR_DATA, "p1_dt_convergence.csv"),
                index=False, encoding="utf-8-sig")

    # ------------------------------------------------------------ 守恒与边界
    push("-" * 78)
    push("守恒性与物理合理性检查")
    push("-" * 78)
    w0 = res.water[0]
    wf = res.water[-1]
    flux = res.cum_flux[-1]
    push(f"初始总水量(归一化) = {w0:.10f}, 末态 = {wf:.10f}, 减少 = {w0-wf:.10e}")
    push(f"累计表面逸出量     = {flux:.10e}")
    push(f"相对闭合误差       = {abs((w0-wf)-flux)/(w0-wf)*100:.6f} %")
    push(f"max(C) 单调不增    = {bool(np.all(np.diff(res.step_maxC) <= 1e-12))}")
    push(f"min(C) >= 0        = {bool(np.all(res.step_minC > 0))}")
    push(f"末态温度区间       = [{res.T_final.min():.4f}, {res.T_final.max():.4f}] degC; "
         f"环境温度上界 = {env.T.max():.4f} degC")
    push(f"末态含水率区间     = [{res.C_final.min():.4f}, {res.C_final.max():.4f}] kg/kg "
         "(应符合 max<初始值、min>0 的单调下降规律)")

    text = "\n".join(lines)
    with open(os.path.join(sc.DIR_DATA, "p1_summary.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
