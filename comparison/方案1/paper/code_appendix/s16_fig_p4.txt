# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s16_fig_p4.py
# 作用  : 在独立进程中生成问题 4 的论文插图(fig19~fig21)。
#
# 原因说明: s07 在同一个进程里连续跑问题 1~4 的长时间算例时, 问题 4 的隐式推进
#   会在某个时间步抛出 "array must not contain infs or NaNs"; 而完全相同的问题 4
#   算例在单独进程里稳定复现(本脚本即为此路径, 同时 s04 的基准算例、
#   s11 的守恒检验都在各自进程里正常完成)。说明该现象与进程内的累积状态有关,
#   而不是算例本身不收敛。为不阻塞交付, 这里把问题 4 的绘图单独成进程执行,
#   并把该现象如实记录在自检报告中。
# 运行  : python s16_fig_p4.py
# =============================================================================
"""Generate the problem-4 figures in a fresh process."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drying_core as dc      # noqa: E402
import plotstyle as ps        # noqa: E402
import scenarios as sc        # noqa: E402
import s07_figures as f7      # noqa: E402


def main():
    ps.setup()
    res4 = dc.simulate(400, 400000.0, 5.0, dc.props_appendix4, sc.load_env(),
                       sc.load_radius("linear"), theta=1.0, c_stop=dc.C_DRY,
                       record_dt=60.0,
                       snapshot_times=list(np.arange(1800.0, 220000.0, 1800.0)))
    print("P4 solved: t_dry = %.4f h" % (sc.drying_time(res4) / 3600.0))
    res3 = dc.simulate(400, 400000.0, 5.0, dc.props_appendix3, sc.load_env(),
                       sc.fixed_radius(), theta=1.0, c_stop=dc.C_DRY,
                       record_dt=60.0,
                       snapshot_times=list(np.arange(3600.0, 220000.0, 3600.0)))
    print("P3 solved: t_dry = %.4f h" % (sc.drying_time(res3) / 3600.0))
    out = f7.fig19_21_p4(res4, res3)
    print("OK", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
