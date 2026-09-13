# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s04_solve_p4.py
# 作用  : 问题 4 求解(考虑水分流失引起的尺寸变化)。
#         模型: 附录 4 物性; 半径历史取附件 2; 采用物质坐标 xi = r/R(t), 把
#               变动域 [0,R(t)] 映射为固定区间 [0,1], 时间导数取固定材料点,
#               因此不出现"收缩速度 x 梯度"的附加对流项。
#         输出位置: 到药材中心的物理距离(0~1.1 cm, 全程位于药材内部)
#                   加一列"药材表面"(由半单元重构, r = R(t))。
#         产出:
#           results/result4.xlsx
#           data/p4_table6_C.csv
#           data/p4_convergence.csv
#           data/p4_summary.txt
# 运行  : python s04_solve_p4.py
# =============================================================================
"""Problem 4: drying of the shrinking cylinder (Appendix-4 correlations)."""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drying_core as dc      # noqa: E402
import scenarios as sc        # noqa: E402

N_BASE = 1600
DT_BASE = 1.0
T_END = 400000.0
DT_STUDY = (60.0, 30.0, 10.0, 5.0, 2.0, 1.0)
DT_MESH = 10.0
GRIDS = (200, 400, 800, 1600, 3200)

POS_OUT_CM = np.round(np.arange(0.0, 1.1 + 1e-9, 0.1), 4)   # 0~1.1 cm 每 0.1 cm
POS_TABLE_CM = np.array([0.0, 0.5, 1.0])
TABLE6_TIMES_H = np.array([6.0, 12.0, 18.0, 24.0, 30.0, 36.0, 42.0, 48.0])


def run_p4(n_cells, dt, env=None, radius=None, record_dt=None, xi=None,
           t_end=T_END):
    """问题 4 求解。"""
    env = env if env is not None else sc.load_env()
    radius = radius if radius is not None else sc.load_radius("linear")
    return dc.simulate(n_cells, t_end, dt, dc.props_appendix4, env, radius,
                       theta=1.0, record_dt=record_dt, record_xi=xi,
                       c_stop=dc.C_DRY)


def main():
    lines = []
    push = lines.append
    env = sc.load_env()
    radius = sc.load_radius("linear")
    R_end = radius(T_END)

    push("=" * 78)
    push("问题 4: 收缩药材的干燥时间(物质坐标形式)")
    push("=" * 78)
    push(f"半径轨迹: R(0) = {radius(0.0)*100:.4f} cm, R(末) = {R_end*100:.4f} cm, "
         f"总收缩 {100*(1-R_end/radius(0.0)):.2f}%")
    push(f"输出位置取物理距离 0~{POS_OUT_CM[-1]:.1f} cm(全程位于药材内部, "
         f"因最小半径 {R_end*100:.3f} cm), 另加一列药材表面 r=R(t)。")

    # ------------------------------------------------------------ 收敛性研究
    conv_rows = []
    for dt in DT_STUDY:
        r = run_p4(N_BASE, dt, env=env, radius=radius)
        td = sc.drying_time(r)
        conv_rows.append({"kind": "dt", "value": dt, "t_dry/s": td,
                          "t_dry/h": td / 3600.0, "n_steps": r.n_steps})
    for n in GRIDS:
        r = run_p4(n, DT_MESH, env=env, radius=radius)
        td = sc.drying_time(r)
        conv_rows.append({"kind": "mesh", "value": n, "t_dry/s": td,
                          "t_dry/h": td / 3600.0, "n_steps": r.n_steps})
    conv = pd.DataFrame(conv_rows)
    conv.to_csv(os.path.join(sc.DIR_DATA, "p4_convergence.csv"),
                index=False, encoding="utf-8-sig")
    push("")
    push("时间步收敛性(N=1600):")
    push(conv[conv["kind"] == "dt"].drop(columns=["kind"]).to_string(
        index=False, float_format=lambda v: f"{v:.6f}"))
    push("网格收敛性(dt=10 s):")
    push(conv[conv["kind"] == "mesh"].drop(columns=["kind"]).to_string(
        index=False, float_format=lambda v: f"{v:.6f}"))
    d = conv[conv["kind"] == "dt"]["t_dry/h"].to_numpy(float)
    e_last, e_prev = abs(d[-2] - d[-1]), abs(d[-3] - d[-2])
    p_obs = np.log2(e_prev / e_last) if e_last > 0 else np.nan
    ext = d[-1] - e_last / (2.0 ** p_obs - 1.0)
    push(f"干燥时间的时间步观测收敛阶 p = {p_obs:.2f}; dt->0 外推 = {ext:.5f} h; "
         f"基准 dt={DT_BASE} s 偏差 = {abs(d[-1]-ext)*3600:.1f} s")

    # ---------------------------------------------------------------- 基准算例
    res = run_p4(N_BASE, DT_BASE, env=env, radius=radius, record_dt=60.0)
    td = sc.drying_time(res)
    push("")
    push("=" * 78)
    push(f"*** 问题 4 干燥时间 t4 = {td:.1f} s = {td/3600.0:.4f} h ***")
    push("=" * 78)
    push(f"基准算例: N={N_BASE}, dt={DT_BASE} s, 步数={res.n_steps}, "
         f"墙钟={res.wall_seconds:.1f} s, 平均 Picard={res.picard_iters.mean():.2f}")

    # 输出序列: 物理位置 + 表面
    C_out, T_out = sc.read_position_series(res, POS_OUT_CM, surface_last=True)
    t_out = res.times
    keep = (t_out >= 60.0 - 1e-9) & (t_out <= td + 1e-9)
    sc.write_result_xlsx(os.path.join(sc.DIR_RES, "result4.xlsx"), t_out[keep],
                         POS_OUT_CM, [C_out[keep]], ["水分浓度"], surface_last=True,
                         surface_header=True)
    push(f"已写出 results/result4.xlsx: {int(keep.sum())} 行 x {len(POS_OUT_CM)} 列 "
         f"(最后一列为药材表面)")

    # 表 6
    pos_tab = np.concatenate([POS_TABLE_CM, [np.nan]])   # 末位为表面占位
    C_tab, _ = sc.read_position_series(res, pos_tab, surface_last=True)
    idx_h = [int(np.argmin(np.abs(t_out / 3600.0 - h))) for h in TABLE6_TIMES_H]
    idx_h = [i for i in idx_h if t_out[i] <= td]
    cols = [f"{v:g} cm" for v in POS_TABLE_CM] + [sc.SURFACE_LABEL]
    tab6 = pd.DataFrame(C_tab[np.ix_(idx_h, range(len(cols)))], index=[
        f"{t_out[i]/3600.0:.0f}" for i in idx_h], columns=cols)
    i_last = min(max(int(np.searchsorted(t_out, td)), 1), len(t_out) - 1)
    w = (td - t_out[i_last - 1]) / (t_out[i_last] - t_out[i_last - 1])
    c_end = C_tab[i_last - 1] * (1 - w) + C_tab[i_last] * w
    c_end[0] = dc.C_DRY
    tab6.loc["烘干结束时间"] = c_end
    tab6.index.name = "时间/h"
    tab6.round(4).to_csv(os.path.join(sc.DIR_DATA, "p4_table6_C.csv"), encoding="utf-8-sig")
    push("")
    push("表 6 药材烘干过程的水分浓度 / (kg/kg)")
    push(tab6.round(4).to_string())
    push(f"(烘干结束时间 = {td/3600.0:.4f} h; 该时刻药材半径 = {radius(td)*100:.4f} cm)")

    # --------------------------------------------------- 守恒与物理合理性
    push("-" * 78)
    push("守恒性与物理合理性检查(问题 4)")
    push("-" * 78)
    w0, wf, flux = res.water[0], res.water[-1], res.cum_flux[-1]
    push(f"归一化总水量 {w0:.10f} -> {wf:.10f} (已乘 R^2); 累计逸出 {flux:.10e}")
    push(f"相对闭合误差 {abs((w0-wf)-flux)/(w0-wf)*100:.6f} %")
    push(f"max(C) 单调不增 = {bool(np.all(np.diff(res.step_maxC) <= 1e-12))}")
    push(f"C 始终非负      = {bool(np.all(res.step_minC > 0))} (最小 {res.step_minC.min():.6f})")
    push(f"全场温度上界检查: 最高 {res.step_maxT.max():.4f} degC vs 环境最高 "
         f"{res.step_Tair.max():.4f} degC -> "
         f"{bool(res.step_maxT.max() <= res.step_Tair.max() + 1e-6)}")

    # 半径插值方式对比(线性 vs 保形三次)
    push("-" * 78)
    push("半径插值方式对比(dt=10 s, N=1600)")
    push("-" * 78)
    for mode in ("linear", "pchip"):
        r = run_p4(N_BASE, DT_MESH, env=env, radius=sc.load_radius(mode))
        push(f"  {mode:7s}: t4 = {sc.drying_time(r)/3600.0:.4f} h")

    text = "\n".join(lines)
    with open(os.path.join(sc.DIR_DATA, "p4_summary.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
