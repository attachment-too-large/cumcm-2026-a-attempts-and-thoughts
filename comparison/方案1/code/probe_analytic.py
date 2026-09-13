# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: probe_analytic.py
# 作用  : 定位解析解对照中"误差不随网格收敛"的原因。分别检查:
#           (1) 解析解在 t->0 与 t->inf 的自检;
#           (2) 数值解与解析解在同一时刻的剖面数值对照;
#           (3) 时间步与网格对误差的影响。
# 运行  : python probe_analytic.py
# =============================================================================
"""Debug the analytic-vs-numerical comparison."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drying_core as dc      # noqa: E402
import scenarios as sc        # noqa: E402
from s05_verification import (analytic_cylinder_robin, constant_env,  # noqa: E402
                              frozen_props)


def main():
    R = dc.R0
    D = 5.0e-9
    props = frozen_props(D)
    env = constant_env(50.0, 0.05)
    r_probe = np.array([0.0, 0.5, 0.9, 1.0]) * R

    # (1) 解析解自检
    a0 = analytic_cylinder_robin(r_probe, np.array([1e-8]), R, D, dc.HM_CONV,
                                 0.05, 2.55)[0]
    ainf = analytic_cylinder_robin(r_probe, np.array([1e8]), R, D, dc.HM_CONV,
                                   0.05, 2.55)[0]
    print("analytic t->0 :", np.round(a0, 6), " (应为 2.55)")
    print("analytic t->inf:", np.round(ainf, 6), " (应为 0.05)")
    print("analytic at t=1800, r:", np.round(
        analytic_cylinder_robin(r_probe, np.array([1800.0]), R, D, dc.HM_CONV,
                                0.05, 2.55)[0], 6))

    # (2) 剖面数值对照
    for t_end, dt in ((1800.0, 0.5), (1800.0, 0.05)):
        res = dc.simulate(800, t_end, dt, props, env, sc.fixed_radius(),
                          theta=1.0, record_dt=600.0,
                          record_xi=r_probe[:-1] / R)
        field = np.column_stack([res.C, res.c_surf])
        i = int(np.argmin(np.abs(res.times - t_end)))
        num = field[i]
        ana = analytic_cylinder_robin(r_probe, np.array([t_end]), R, D,
                                      dc.HM_CONV, 0.05, 2.55)[0]
        print(f"t={t_end:g} dt={dt:g}: num={np.round(num,6)} ana={np.round(ana,6)} "
              f"diff={np.round(num-ana,2e0 if False else 8)}")

    # (3) 中间时刻 t=30 s
    for t_end in (30.0, 60.0, 300.0):
        res = dc.simulate(800, t_end, 0.05, props, env, sc.fixed_radius(),
                          theta=1.0, record_dt=1.0, record_xi=r_probe[:-1] / R)
        field = np.column_stack([res.C, res.c_surf])
        i = int(np.argmin(np.abs(res.times - t_end)))
        num = field[i]
        ana = analytic_cylinder_robin(r_probe, np.array([t_end]), R, D,
                                      dc.HM_CONV, 0.05, 2.55)[0]
        print(f"t={t_end:g}: num={np.round(num,8)} ana={np.round(ana,8)} "
              f"diff={np.round(num-ana,10)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
