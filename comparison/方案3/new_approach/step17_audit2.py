# -*- coding: utf-8 -*-
"""两个诊断：
(a) 问题2 温度：paper 格式外推 dt->0 是否收敛到本文谱解的值（判定谁的 1e-4 在误差内）；
(b) 环境口径分叉：用附件1 原始插值代替一阶惯性拟合，看答案会怎样变
    —— 这是模型层真正的自由度，paper 自己报告过该分叉。
"""
import numpy as np
import common as cm
from spectral import HerbSolver
from paper_fv import simulate

out = []

# ---------- (a) ----------
out.append("== (a) 问题2 (0.5h, r=0) 的收敛值 ==")
out.append("本文谱解 N=64/96/128, rtol=1e-11~1e-12 :")
for N, rt in [(64, 1e-11), (96, 1e-12), (128, 1e-12)]:
    s = HerbSolver(N, prob=2)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-12)])
    sol = s.solve(1800.0, t_eval=[1800.0], rtol=rt, atol=at)
    Cn, Tn = s.interp_u(sol.y[:, -1], np.array([0.0]))
    out.append("   N=%-4d %.9f" % (N, Tn[0]))
out.append("paper 格式 N=1600 的 dt 序列（一阶）：")
lad = []
for dt in [4.0, 2.0, 1.0, 0.5, 0.25, 0.125]:
    xi, tt, Cs, Ts = simulate(1600, 2, 1800.0, [(1e18, dt)],
                              t_report=[1800.0], report_idx=[1800.0])
    lad.append((dt, Ts[0][0]))
    out.append("   dt=%-6.3f %.9f" % (dt, Ts[0][0]))
c = (lad[-2][1] - lad[-1][1]) / (lad[-2][0] - lad[-1][0])
out.append("   一阶外推 dt->0 : %.9f" % (lad[-1][1] - c * lad[-1][0]))
out.append("paper 表3 印制值      : 32.5656")

# ---------- (b) 环境口径分叉 ----------
out.append("")
out.append("== (b) 环境口径分叉（问题1，其余一切不变）==")
_T_air_fit, _C_air_fit = cm.T_air, cm.C_air


def T_air_raw(t):
    t = np.asarray(t, dtype=float)
    v = np.interp(t, cm._t_env, cm._Ta_env)
    return np.where(t <= cm.T_END, v, cm.Tset)


def C_air_raw(t):
    t = np.asarray(t, dtype=float)
    v = np.interp(t, cm._t_env, cm._Ca_env)
    return np.where(t <= cm.T_END, v, cm.Cset)


for tag, fT, fC in [("一阶惯性拟合（本文=paper 口径）", _T_air_fit, _C_air_fit),
                    ("附件1 原始序列插值", T_air_raw, C_air_raw)]:
    cm.T_air, cm.C_air = fT, fC
    s = HerbSolver(64, prob=1)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])
    sol = s.solve(1800.0, t_eval=[1800.0], rtol=1e-11, atol=at)
    Cn, Tn = s.interp_u(sol.y[:, -1], np.array([0.0, 1.0]))
    out.append("%-28s T(0,1800)=%.4f  T(R,1800)=%.4f  C(R,1800)=%.4f"
               % (tag, Tn[0], Tn[1], Cn[1]))
cm.T_air, cm.C_air = _T_air_fit, _C_air_fit

out.append("")
out.append("目标答案: T(0)=34.0425  T(R)=37.1989  C(R)=1.5109")

open(r"C:\Users\qing1\Desktop\数A\new_approach\step17_audit2.txt", "w",
     encoding="utf-8").write("\n".join(out))
print("done")
