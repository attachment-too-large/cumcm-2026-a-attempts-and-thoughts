# -*- coding: utf-8 -*-
"""用 paper 自己的格式复现其表7 的问题3 序列，判断 16 s 差异是数值还是模型问题。"""
import numpy as np
import common as cm
from paper_fv import simulate, find_tdry

out = []
out.append("== 复现 paper 表7 问题3 序列（节点中心有限体积 + theta 法）==")
out.append("paper: N=200:205358.0  N=400:205488.3  N=800:205537.3  N=1600:205553.8")
out.append("       N=1600,dt=1s:205567.4 ; 单元中心N=1600:205528.5 ; 外推:205559.3")

# 时间步计划：前 7200 s 用 1 s，其后 10 s
plan_paper = [(7200.0, 1.0), (1e18, 10.0)]
plan_dt1 = [(1e18, 1.0)]

for N in [200, 400, 800, 1600]:
    td = find_tdry(N, 3, plan_paper)
    out.append("N=%-5d dt=(1s->10s)  tdry = %9.1f s = %8.4f h" % (N, td, td / 3600))
for N in [800, 1600]:
    td = find_tdry(N, 3, plan_dt1)
    out.append("N=%-5d dt=1s 全程     tdry = %9.1f s = %8.4f h" % (N, td, td / 3600))

# 时间步收敛性（N=1600）
out.append("")
out.append("== N=1600 时间步收敛性 ==")
for d in [40.0, 20.0, 10.0, 5.0, 2.0, 1.0]:
    plan = [(7200.0, min(1.0, d)), (1e18, d)]
    td = find_tdry(1600, 3, plan)
    out.append("dt=%-5.1f tdry = %9.1f s" % (d, td))

open(r"C:\Users\qing1\Desktop\数A\new_approach\step8_paper_repro.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
