# -*- coding: utf-8 -*-
"""验证 tdry 的数值可信度：
  (a) 谱解：换时间积分器/容差，看 tdry 是否稳定；
  (b) 谱解：用 1 s 采样的线性求交 vs 事件求交，检查事件定位是否可靠；
  (c) paper 格式 FV：N=1600 时间步向 0 外推，看是否收敛到谱解的值。
"""
import numpy as np
import common as cm
from spectral import HerbSolver
from paper_fv import find_tdry

out = []
CRIT = 0.15


def spectral_tdry(N, method, rtol, atol_C=1e-13, atol_T=1e-11):
    s = HerbSolver(N, prob=3)
    at = np.concatenate([np.full(s.Np, atol_C), np.full(s.Np, atol_T)])

    def ev(t, y):
        return np.max(y[:s.Np]) - CRIT
    ev.terminal = True
    ev.direction = -1
    sol = s.solve(300000.0, rtol=rtol, atol=at, events=ev, method=method,
                  first_step=1e-3)
    return s, sol, float(sol.t_events[0][0])


out.append("== (a) 谱解：时间积分器/容差敏感性（N=64）==")
for method, tol in [("BDF", 1e-9), ("BDF", 1e-10), ("BDF", 1e-11), ("BDF", 1e-12),
                    ("Radau", 1e-9), ("Radau", 1e-11), ("LSODA", 1e-10)]:
    try:
        s, sol, td = spectral_tdry(64, method, tol)
        out.append("%-6s rtol=%.0e  tdry = %11.3f s  nfev=%d" % (method, tol, td, sol.nfev))
    except Exception as e:
        out.append("%-6s rtol=%.0e  FAILED %s" % (method, tol, e))

out.append("")
out.append("== (b) 事件求交 vs 1 s 采样线性求交（N=64, BDF rtol=1e-11）==")
s, sol, T1 = spectral_tdry(64, "BDF", 1e-11)
at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])
grid = np.arange(T1 - 400.0, T1 + 400.0, 1.0)
sol2 = s.solve(T1 + 500.0, t_eval=grid, rtol=1e-11, atol=at)
Cmax = np.array([np.max(sol2.y[:s.Np, j]) for j in range(sol2.y.shape[1])])
k = int(np.argmax(Cmax < CRIT))
T2 = grid[k - 1] + (CRIT - Cmax[k - 1]) / (Cmax[k] - Cmax[k - 1]) * 1.0
out.append("事件求交   tdry = %.3f s" % T1)
out.append("1s采样求交 tdry = %.3f s   差 %.3f s" % (T2, T1 - T2))
out.append("该处 dC/dt = %.4e /s  (1 s 对应 %.2e kg/kg)" % (Cmax[k] - Cmax[k - 1], abs(Cmax[k] - Cmax[k - 1])))

out.append("")
out.append("== (c) paper 格式 FV N=1600：时间步外推 ==")
for d in [5.0, 2.0, 1.0, 0.5, 0.25]:
    plan = [(7200.0, min(1.0, d)), (1e18, d)]
    td = find_tdry(1600, 3, plan)
    out.append("dt=%-5.2f tdry = %11.2f s" % (d, td))

out.append("")
out.append("== (d) paper 格式 FV：空间加密（dt=1s 全程）==")
for N in [800, 1600, 3200, 6400]:
    td = find_tdry(N, 3, [(1e18, 1.0)])
    out.append("N=%-5d dt=1s tdry = %11.2f s" % (N, td))

open(r"C:\Users\qing1\Desktop\数A\new_approach\step10_verify_tdry.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
