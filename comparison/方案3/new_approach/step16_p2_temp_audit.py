# -*- coding: utf-8 -*-
"""问题2 温度场的三方对账：本文谱解 vs paper格式(充分收敛) vs 参考论文表3。

目的：判断 paper 表3 前几档比本文低 1e-4 究竟是本文的误差还是 paper 的离散误差。
"""
import numpy as np
import common as cm
from spectral import HerbSolver
from paper_fv import simulate

out = []
TR = [1800.0, 3600.0, 5400.0]
IDX_R = [0, 8, 16, 24, 32]     # N=32 时的 0,0.5,1,1.5,2 cm ; 由 N 决定

# ---------- 1) 本文谱解（加密 + 收紧容差）----------
out.append("== 1) 本文谱配置法（空间加密 + 收紧容差）==")
out.append("N      rtol     T(0,1800)     T(R,1800)     T(0,3600)     T(R,3600)")
spec = {}
for N, rt in [(64, 1e-11), (96, 1e-12), (128, 1e-12)]:
    s = HerbSolver(N, prob=2)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-12)])
    sol = s.solve(5400.0, t_eval=TR, rtol=rt, atol=at)
    v = []
    for i in range(len(TR)):
        Cn, Tn = s.interp_u(sol.y[:, i], np.array([0.0, 1.0]))
        v += [Tn[0], Tn[1]]
    spec[N] = v
    out.append("N=%-4d %.0e  %.9f  %.9f  %.9f  %.9f" % (N, rt, v[0], v[1], v[2], v[3]))
out.append("paper 表3 :      32.5656       35.7831       40.4220       42.8629")

# ---------- 2) paper 格式，但充分收敛 ----------
out.append("")
out.append("== 2) paper 格式（节点中心FV+theta法），加密到接近收敛 ==")
out.append("N      dt      T(0,1800)     T(R,1800)     T(0,3600)     T(R,3600)")
for N, dt in [(800, 1.0), (1600, 1.0), (3200, 0.5), (3200, 0.25)]:
    xi, tt, Cs, Ts = simulate(N, 2, 5400.0, [(1e18, dt)],
                              t_report=TR, report_idx=TR)
    v = []
    for i in range(len(TR)):
        Ti = Ts[i]
        v += [Ti[0], Ti[N]]                 # xi=0 与 xi=1 恰好是首末节点
    out.append("N=%-4d %-6.2f %.9f  %.9f  %.9f  %.9f" % (N, dt, v[0], v[1], v[2], v[3]))

# ---------- 3) paper 格式：时间步 / 空间步收敛趋势 ----------
out.append("")
out.append("== 3) paper 格式在 (0.5h, r=0) 的收敛趋势 ==")
out.append("N=1600 固定, 变 dt：")
for dt in [4.0, 2.0, 1.0, 0.5, 0.25]:
    xi, tt, Cs, Ts = simulate(1600, 2, 1800.0, [(1e18, dt)],
                              t_report=[1800.0], report_idx=[1800.0])
    out.append("   dt=%-5.2f  T(0,1800) = %.9f" % (dt, Ts[0][0]))
out.append("dt=0.25 固定, 变 N：")
for N in [400, 800, 1600, 3200, 6400]:
    xi, tt, Cs, Ts = simulate(N, 2, 1800.0, [(1e18, 0.25)],
                              t_report=[1800.0], report_idx=[1800.0])
    out.append("   N=%-5d  T(0,1800) = %.9f" % (N, Ts[0][0]))

open(r"C:\Users\qing1\Desktop\数A\new_approach\step16_p2_temp_audit.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
