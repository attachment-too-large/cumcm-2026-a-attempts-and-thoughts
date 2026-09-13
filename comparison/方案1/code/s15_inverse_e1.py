# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s15_inverse_e1.py
# 作用  : 用"网格极小点 + 局部曲率"重估 h_m, 刷新 data/inverse_estimate.csv 与
#         data/inverse_residuals.csv(其余列保持 s14 的结果不变)。
#
# 为什么: 目标函数 SSE(h_m) 在极小点附近非常尖锐(偏离 5% 即使 SSE 增大数十倍),
#   用整条扫描曲线做全局二次拟合会把顶点拉偏(实测偏 +7%, 且 42 倍标准差之外)。
#   正确的做法是取扫描网格的极小点作为估计值, 并用该点两侧三点的二阶差分
#   估计曲率, 再由  sigma^2 * Hessian^-1  的线性化公式给出标准差。
# 运行  : python s15_inverse_e1.py
# =============================================================================
"""Robust re-estimation of h_m from the synthetic-data experiment."""

from __future__ import annotations

import io
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drying_core as dc      # noqa: E402
import scenarios as sc        # noqa: E402

N_INV, DT_INV, NOISE, SEED = 400, 2.0, 0.005, 20260911


def main():
    env = sc.load_env()
    xi = np.array([0.0, 2.0]) * 1e-2 / dc.R0

    def run_model(hm):
        def props(C, T_C):
            rho, cp, k, D = dc.props_appendix2(C, T_C)
            return rho, cp, k, D
        r = dc.simulate(N_INV, 1800.0, DT_INV, props, env, sc.fixed_radius(),
                        theta=1.0, record_dt=60.0, record_xi=xi, h_mass=hm)
        return r.times[1:], np.column_stack([r.C[1:, 0], r.c_surf[1:]])

    hm_true = dc.HM_CONV
    _, y_true = run_model(hm_true)
    rng = np.random.default_rng(SEED)
    y_obs = y_true + rng.normal(0.0, NOISE, size=y_true.shape)

    grid = np.linspace(0.4e-6, 1.6e-6, 13)
    sse = np.array([float(np.sum((run_model(h)[1] - y_obs) ** 2)) for h in grid])
    i0 = int(np.argmin(sse))
    x_fit = float(grid[i0])
    h = float(grid[1] - grid[0])
    curv = float((sse[i0 + 1] - 2.0 * sse[i0] + sse[i0 - 1]) / h ** 2)
    n_obs = y_obs.size
    sigma2 = float(sse[i0]) / (n_obs - 2)
    sd = float(np.sqrt(2.0 * sigma2 / curv)) if curv > 0 else np.nan

    _, y_fit = run_model(x_fit)
    resid = (y_obs - y_fit).ravel()
    ac1 = float(np.corrcoef(resid[:-1], resid[1:])[0, 1])

    lines = []
    push = lines.append
    push("E1 单参数反演(网格极小点 + 局部曲率):")
    push(f"  扫描网格: {grid[0]*1e6:.2f}e-6 ~ {grid[-1]*1e6:.2f}e-6, {grid.size} 点")
    push(f"  SSE 曲线: {np.round(sse, 6)}")
    push(f"  网格极小点 = {x_fit*1e7:.4f}e-7 m/s (真值 {hm_true*1e7:.4f}e-7)")
    push(f"  局部曲率 = {curv:.4e}; 残差方差 = {sigma2:.4e}")
    push(f"  标准差 = {sd*1e7:.4f}e-7 m/s; 95% 置信区间 = "
         f"[{(x_fit-1.96*sd)*1e7:.4f}, {(x_fit+1.96*sd)*1e7:.4f}]e-7 m/s")
    push(f"  真值{'落在' if abs(x_fit-hm_true) <= 1.96*sd else '不落在'}95% 置信区间内")
    push(f"  残差: 均值 {resid.mean():+.3e}, 标准差 {resid.std(ddof=1):.4e} "
         f"(噪声 {NOISE:.4e}, 比值 {resid.std(ddof=1)/NOISE:.3f}), 一阶自相关 {ac1:+.3f}")
    text = "\n".join(lines)
    print(text)

    p = os.path.join(sc.DIR_DATA, "inverse_estimate.csv")
    est = pd.read_csv(p)
    est.loc[0, "hm_fit"] = x_fit
    est.loc[0, "hm_sd"] = sd
    est.loc[0, "hm_lo95"] = x_fit - 1.96 * sd
    est.loc[0, "hm_hi95"] = x_fit + 1.96 * sd
    est.loc[0, "resid_sd"] = resid.std(ddof=1)
    est.loc[0, "resid_ac1"] = ac1
    est.to_csv(p, index=False, encoding="utf-8-sig")

    pd.DataFrame({"t_s": np.tile(run_model(x_fit)[0], 2),
                  "kind": np.repeat(["center", "surface"], y_true.shape[0]),
                  "obs": y_obs.ravel(), "fit": y_fit.ravel(),
                  "resid": resid}).to_csv(
        os.path.join(sc.DIR_DATA, "inverse_residuals.csv"),
        index=False, encoding="utf-8-sig")
    with io.open(os.path.join(sc.DIR_DATA, "inverse_summary_e1.txt"), "w",
                 encoding="utf-8", newline="\n") as f:
        f.write(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
