# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: probe_smoke.py
# 作用  : 求解器冒烟测试。用「新模型」PDF 中给出的配置(N=300, 附件1 原始插值,
#         半单元表面重构)复现问题 1(1800 s)与问题 2(3 h)的关键数字, 判断本实现
#         是否与该模型的实现等价。不做任何结果落盘以外的永久改动。
# 目标值(PDF): P1 中心 T=33.5766, 表面 T=36.7863; 中心 C=2.5500, 表面 C=1.5103
#              P2 中心 T=49.8494, 表面 T=49.9664; 中心 C=1.7662, 表面 C=1.0081
# 运行  : python probe_smoke.py
# =============================================================================
"""Smoke test against the target numbers published in the new-model document."""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drying_core as dc  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_env():
    df = pd.read_excel(os.path.join(BASE, "data", "input", "attachment1_room.xlsx"))
    df.columns = ["t", "T", "C"]
    return dc.RoomConditions(df["t"].to_numpy(float), df["T"].to_numpy(float),
                             df["C"].to_numpy(float), plateau_rule="tail_mean")


def report(name, res, props, env):
    idx0, w0 = dc.lagrange_weights_grid(res.n_cells, [0.0])
    c0 = float(dc.apply_weights(res.C_final, idx0, w0)[0])
    t0 = float(dc.apply_weights(res.T_final, idx0, w0)[0])
    cs = float(res.c_surf_final)
    ts = float(res.t_surf_final)
    print(f"{name}: 中心 T={t0:.4f} C={c0:.4f} | 表面 T={ts:.4f} C={cs:.4f} "
          f"| 步数={res.n_steps} 平均Picard={res.picard_iters.mean():.2f} "
          f"耗时={res.wall_seconds:.2f}s")


def main():
    env = load_env()
    rfix = dc.constant_radius(dc.R0)

    print("平台取值(tail_mean):", f"T={env.T_plateau:.4f} degC, C={env.C_plateau:.6f}")

    for n in (300, 800):
        res = dc.simulate(n, 1800.0, 1.0, dc.props_appendix2, env, rfix, theta=1.0,
                          record_dt=None)
        report(f"P1 N={n}", res, dc.props_appendix2, env)

    for n in (300,):
        res = dc.simulate(n, 10800.0, 1.0, dc.props_appendix3, env, rfix, theta=1.0,
                          record_dt=None)
        report(f"P2 N={n}", res, dc.props_appendix3, env)

    # 计时基准: 为问题 3/4 选择网格与时间步提供依据
    import time
    for n, dt in ((400, 60.0),):
        t0 = time.perf_counter()
        res = dc.simulate(n, 3600.0, dt, dc.props_appendix3, env, rfix, theta=1.0,
                          record_dt=None)
        dt_wall = time.perf_counter() - t0
        print(f"计时基准 N={n} dt={dt}: {res.n_steps} 步耗时 {dt_wall:.3f}s "
              f"-> 单步 {1e6*dt_wall/res.n_steps:.1f} us; "
              f"推算 57 h(205200 s) 需 {dt_wall/res.n_steps*205200/dt/60:.2f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
