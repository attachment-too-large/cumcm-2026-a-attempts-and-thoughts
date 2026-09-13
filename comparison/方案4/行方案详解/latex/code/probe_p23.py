# -*- coding: utf-8 -*-
# ============================================================================
# probe_p23.py —— 粗网格试算探针（不属于交付计算）
#
# 用途：用 N=400 的粗网格把问题 2/3 快速跑一遍，每隔半小时打印中心与各半径
#       处的含水率，并给出中心首次降到 0.15 kg/kg 的时刻。据此核对干燥时长的
#       量级、确定模拟终了时刻与输出采样间隔；正式结果由 p23_process.py 产出。
#
# 主要变量：
#   N        径向控制体数目（粗网格，只求量级）
#   T_END    模拟终了时刻 [s]，取 72 h（题目称烘干过程持续 2~3 天）
#   t_rep    输出时刻序列，每 1800 s 一个
#   xi_rep   采样位置（归一化半径 xi = r/R0）
# ============================================================================

import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm

N = 400
T_END = 259200.0
room = cm.RoomConditions()
t_rep = np.arange(0.0, T_END + 1.0, 1800.0)
xi_rep = np.array([0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0])

t0 = time.time()
res = cm.simulate(N, cm.props_problem23, T_END, 1.0, room,
                  t_report=t_rep, xi_report=xi_rep,
                  dt_schedule=[(7200.0, 1.0), (T_END, 10.0)])
print("wall time %.1f s, steps=%d" % (time.time() - t0, res["n_step"]))
print("   t/h   T_air  T_cen  T_sur   C_cen   C_0.5R  C_0.9R  C_sur")
for j, tt in enumerate(res["t"]):
    if tt % 3600.0 > 0.5:
        continue                        # 只打印整小时，压缩输出
    row = res["C"][j]
    print("%7.1f %6.2f %6.2f %6.2f  %7.4f %7.4f %7.4f %7.4f"
          % (tt / 3600.0, float(room.T_air(tt)), res["T"][j, 0], res["T"][j, -1],
             row[0], row[2], row[4], row[-1]))

# 中心含水率首次低于判据 0.15 kg/kg 的时刻
idx = np.where(res["C"][:, 0] < 0.15)[0]
if len(idx):
    print("centre reaches 0.15 at t = %.1f s = %.2f h"
          % (res["t"][idx[0]], res["t"][idx[0]] / 3600.0))
else:
    print("centre NOT below 0.15 within %.1f h" % (T_END / 3600.0))
