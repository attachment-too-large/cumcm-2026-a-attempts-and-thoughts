# -*- coding: utf-8 -*-
"""复现 paper 表7 的问题4 序列 + Richardson，并与新解法（谱配置）对照。"""
import numpy as np
from paper_fv import find_tdry

out = []
out.append("== 问题4：复现 paper 表7 序列（节点中心FV + theta法 + 附件2 的 R(t)）==")
out.append("paper: N=200:182726.8  N=400:182753.9  N=800:182761.9  N=1600:182764.0")
out.append("       单元中心N=1600:182756.0 ; Richardson: 182764.7 ; 新解法(谱):182778.3")

plan = [(7200.0, 1.0), (1e18, 10.0)]
seq = {}
for N in [200, 400, 800, 1600]:
    td = find_tdry(N, 4, plan, shrink=True)
    seq[N] = td
    out.append("N=%-5d dt=(1s->10s)  tdry = %10.2f s = %8.4f h" % (N, td, td / 3600))

# Richardson（标称二阶）
N1, N2 = 800, 1600
t1, t2 = seq[N1], seq[N2]
rich2 = t2 + (t2 - t1) / 3.0
out.append("")
out.append("Richardson(标称二阶) = %.2f + (%.2f-%.2f)/3 = %.2f s" % (t2, t2, t1, rich2))
out.append("paper 报告值 182764.7 s -> 偏差 %+.2f s" % (rich2 - 182764.7))

out.append("")
out.append("== N=1600 时间步外推（检验一阶时间误差）==")
lad = {}
for d in [10.0, 5.0, 2.0, 1.0, 0.5]:
    td = find_tdry(1600, 4, [(7200.0, min(1.0, d)), (1e18, d)], shrink=True)
    lad[d] = td
    out.append("dt=%-5.2f tdry = %10.2f s" % (d, td))
k = (lad[2.0] - lad[1.0]) / (2.0 - 1.0)
out.append("一阶外推 Limit(N=1600, dt->0) = %.2f + %.2f = %.2f s"
           % (lad[1.0], k, lad[1.0] + k))
out.append("新解法（谱、时空均收敛）= 182778.30 s ; 差 %.2f s"
           % (182778.30 - (lad[1.0] + k)))

out.append("")
out.append("== 空间加密（dt=1s 全程）==")
for N in [800, 1600, 3200]:
    td = find_tdry(N, 4, [(1e18, 1.0)], shrink=True)
    out.append("N=%-5d dt=1s tdry = %10.2f s" % (N, td))

open(r"C:\Users\qing1\Desktop\数A\new_approach\step11_p4_repro.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
