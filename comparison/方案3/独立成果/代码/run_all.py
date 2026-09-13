# -*- coding: utf-8 -*-
"""一键复跑：问题 1-4 全部结果 + 全部验证。输出 复跑报告.txt。

用法：  python run_all.py            （约 10-20 分钟）
      python run_all.py --fast     （跳过最重的 t_dry 稳定性扫描）
"""
import sys
import time
import numpy as np
import common as cm
import analytic_p1 as an
from spectral import HerbSolver
from fv_theta import simulate as fv_theta_run, find_tdry as fv_theta_tdry
from fv_uniform import solve_fv, sample_xi

FAST = "--fast" in sys.argv
OUT = []
T0 = time.time()


def w(s=""):
    OUT.append(str(s))
    print(str(s), flush=True)


def hdr(t):
    w("")
    w("=" * 78)
    w(t)
    w("=" * 78)


RCOLS = [0.0, 0.005, 0.010, 0.015, 0.020]

# ======================================================================
hdr("0. 烘房环境标定（附件1，一阶惯性拟合，初值固定为实测值）")
_te, _Ta, _Ca = cm._t_env, cm._Ta_env, cm._Ca_env


def _mT(t, Ts, tau):
    return Ts - (Ts - cm.T_air0) * np.exp(-t / tau)


def _mC(t, Cs, tau):
    return Cs - (Cs - cm.C_air0) * np.exp(-t / tau)


def _stats(f, p, t, y):
    r = y - f(t, *p)
    rmse = float(np.sqrt(np.mean(r ** 2)))
    r2 = 1.0 - float(np.sum(r ** 2) / np.sum((y - y.mean()) ** 2))
    J = np.column_stack([(f(t, *(p + np.eye(len(p))[i] * (abs(p[i]) * 1e-7 + 1e-12)))
                          - f(t, *p)) / (abs(p[i]) * 1e-7 + 1e-12)
                         for i in range(len(p))])
    s2 = float(np.sum(r ** 2) / (len(t) - len(p)))
    se = np.sqrt(np.diag(s2 * np.linalg.inv(J.T @ J)))
    return rmse, r2, se


_pT = np.array([cm.Tset_fit, cm.tauT_fit])
_pC = np.array([cm.Cset_fit, cm.tauC_fit])
rmseT, r2T, seT = _stats(_mT, _pT, _te, _Ta)
rmseC, r2C, seC = _stats(_mC, _pC, _te, _Ca)
w("T_set = %.6f degC  (标准误 %.4f)   tau_T = %.4f s  (标准误 %.2f)"
  % (cm.Tset_fit, seT[0], cm.tauT_fit, seT[1]))
w("C_set = %.6f kg/kg (标准误 %.2e)   tau_C = %.4f s  (标准误 %.2f)"
  % (cm.Cset_fit, seC[0], cm.tauC_fit, seC[1]))
w("温度拟合: RMSE = %.4f degC,  R^2 = %.6f" % (rmseT, r2T))
w("湿度拟合: RMSE = %.3e kg/kg, R^2 = %.6f" % (rmseC, r2C))
w("初值固定为附件1 实测首点: T=%.2f degC, C=%.6f kg/kg" % (cm.T_air0, cm.C_air0))
w("生产取值（按有效位数报告）: Tset=%.3f  tauT=%.1f  Cset=%.6f  tauC=%.1f"
  % (cm.Tset, cm.tauT, cm.Cset, cm.tauC))
w("两阶段工艺：0<=t<=14400 s 按上式，t>14400 s 恒定为 (Tset, Cset)。")

# ======================================================================
hdr("1. 问题1 —— 闭式解自检（Bessel-Duhamel）")
w(an._selftest())

# ======================================================================
hdr("2. 问题1 —— 谱配置法 vs 闭式解（只加密空间，rtol=1e-12）")
w("N      全场max误差/K    T(0,1800)误差/K   T(R,1800)误差/K")
for N in [8, 12, 16, 24, 32, 64]:
    s = HerbSolver(N, prob=1)
    at = np.concatenate([np.full(s.Np, 1e-14), np.full(s.Np, 1e-12)])
    sol = s.solve(1800.0, t_eval=[1800.0], rtol=1e-12, atol=at)
    rc = np.linspace(0, cm.R0, 201)
    Ti = s._bary(sol.y[s.Np:, -1], (rc / cm.R0) ** 2)
    Te = an.T_exact(rc, 1800.0)
    w("N=%-5d %13.3e %16.3e %17.3e"
      % (N, np.max(np.abs(Ti - Te)), Ti[0] - Te[0], Ti[-1] - Te[-1]))

# ======================================================================
hdr("3. 问题1 —— 表1（温度 / degC）")
s1 = HerbSolver(64, prob=1)
at1 = np.concatenate([np.full(s1.Np, 1e-13), np.full(s1.Np, 1e-11)])
tt1 = np.array([100., 300., 600., 900., 1200., 1500., 1800.])
sol1 = s1.solve(1800.0, t_eval=list(tt1), rtol=1e-11, atol=at1)
U1 = (np.array(RCOLS) / cm.R0) ** 2
w("t/s      " + "".join("%10s" % ("r=%.1fcm" % (x * 100)) for x in RCOLS))
tabT = np.zeros((7, 5))
tabC = np.zeros((7, 5))
for i in range(7):
    Cn, Tn = s1.interp_u(sol1.y[:, i], U1)
    tabT[i] = Tn
    tabC[i] = Cn
    w("%-8.0f " % tt1[i] + "".join("%10.4f" % v for v in Tn))

hdr("4. 问题1 —— 表2（水分浓度 / (kg/kg)）")
w("t/s      " + "".join("%10s" % ("r=%.1fcm" % (x * 100)) for x in RCOLS))
for i in range(7):
    w("%-8.0f " % tt1[i] + "".join("%10.4f" % v for v in tabC[i]))

# ======================================================================
hdr("5. 问题1 —— 谱解网格收敛性（水分场）")
w("N      C(0,1800)      C(R,1800)      C(R,1200)")
for N in [24, 32, 48, 64, 96, 128]:
    s = HerbSolver(N, prob=1)
    at = np.concatenate([np.full(s.Np, 1e-14), np.full(s.Np, 1e-12)])
    sol = s.solve(1800.0, t_eval=[1200.0, 1800.0], rtol=1e-12, atol=at)
    a = s.interp_u(sol.y[:, 0], (0.02 / cm.R0) ** 2)[0][0]
    b = s.interp_u(sol.y[:, 1], U1)[0]
    w("N=%-5d %.9f   %.9f   %.9f" % (N, b[0], b[-1], a))

# ======================================================================
hdr("6. 问题1 —— 第三种离散族互证（均匀 r 网格守恒型有限体积 + 后向Euler）")
w("N     dt      C(R,1800)      C(R,1200)      T(R,1800)")
for N, dt in [(200, 1.0), (200, 0.5), (400, 0.5), (400, 0.25)]:
    xi, Ctab, Ttab = solve_fv(N, dt, 1800.0, prob=1, t_report=[1200.0, 1800.0])
    c12 = float(np.atleast_1d(sample_xi(xi, Ctab[0], 0.02 / cm.R0))[0])
    c18 = float(np.atleast_1d(sample_xi(xi, Ctab[1], 0.02 / cm.R0))[0])
    t18 = float(np.atleast_1d(sample_xi(xi, Ttab[1], 0.02 / cm.R0))[0])
    w("N=%-4d %-6.2f %12.8f  %12.8f  %12.8f" % (N, dt, c18, c12, t18))
w("谱配置法(N=128):           1.51091810      1.65914239      37.198922")

# ======================================================================
hdr("7. 问题2 —— 表3（温度 / degC）与表4（水分浓度 / (kg/kg)）")
tt2 = np.array([1800., 3600., 5400., 7200., 9000., 10800.])
s2 = HerbSolver(64, prob=2)
at2 = np.concatenate([np.full(s2.Np, 1e-12), np.full(s2.Np, 1e-10)])
sol2 = s2.solve(10800.0, t_eval=list(tt2), rtol=1e-11, atol=at2)
w("表3  温度")
w("t/h      " + "".join("%10s" % ("r=%.1fcm" % (x * 100)) for x in RCOLS))
_p2T = np.zeros((6, 5))
_p2C = np.zeros((6, 5))
for i in range(6):
    Cn, Tn = s2.interp_u(sol2.y[:, i], U1)
    _p2T[i] = Tn
    _p2C[i] = Cn
    w("%-8.1f " % (tt2[i] / 3600) + "".join("%10.4f" % v for v in Tn))
w("")
w("表4  水分浓度")
w("t/h      " + "".join("%10s" % ("r=%.1fcm" % (x * 100)) for x in RCOLS))
for i in range(6):
    w("%-8.1f " % (tt2[i] / 3600) + "".join("%10.4f" % v for v in _p2C[i]))

hdr("8. 问题2 —— 两种离散族在 (0.5h, r=0) 的互证")
w("谱配置法 N=64/96/128 :")
for N in [64, 96, 128]:
    s = HerbSolver(N, prob=2)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-12)])
    sol = s.solve(1800.0, t_eval=[1800.0], rtol=1e-12, atol=at)
    w("   N=%-4d %.9f" % (N, s.interp_u(sol.y[:, -1], np.array([0.0]))[1][0]))
w("经典 FV+theta 格式 N=1600，时间步一阶外推：")
lad = []
for dt in [2.0, 1.0, 0.5, 0.25]:
    xi, tt, Cs, Ts = fv_theta_run(1600, 2, 1800.0, [(1e18, dt)],
                                  t_report=[1800.0], report_idx=[1800.0])
    lad.append((dt, Ts[0][0]))
    w("   dt=%-6.3f %.9f" % (dt, Ts[0][0]))
c = (lad[-2][1] - lad[-1][1]) / (lad[-2][0] - lad[-1][0])
w("   一阶外推 dt->0 : %.9f" % (lad[-1][1] - c * lad[-1][0]))

# ======================================================================
CRIT = 0.15


def dry_time(N, prob, shrink=False, method="BDF", rtol=1e-11, ac=1e-13, at_=1e-11):
    s = HerbSolver(N, prob=prob, shrink=shrink)
    at = np.concatenate([np.full(s.Np, ac), np.full(s.Np, at_)])

    def ev(t, y):
        return np.max(y[:s.Np]) - CRIT
    ev.terminal = True
    ev.direction = -1
    sol = s.solve(300000.0, rtol=rtol, atol=at, events=ev, method=method)
    return s, sol, float(sol.t_events[0][0])


hdr("9. 问题3 —— 烘干时间 t_dry")
w("谱配置法网格收敛性：")
for N in [24, 32, 48, 64, 96, 128]:
    s, sol, td = dry_time(N, 3, rtol=1e-10)
    w("N=%-5d tdry = %11.3f s = %8.4f h" % (N, td, td / 3600))
s3, solE3, td3 = dry_time(64, 3)
w("生产值 N=64: tdry = %.3f s = %.4f h = %.4f d" % (td3, td3 / 3600, td3 / 86400))
w("终点 C(0)=%.8f, C(表面)=%.8f, 最大值位于 u=%.4g"
  % (solE3.y[0, -1], solE3.y[s3.N, -1],
     s3.u[int(np.argmax(solE3.y[:s3.Np, -1]))]))

if not FAST:
    w("")
    w("时间积分器/容差敏感性（验证 tdry 不是数值噪声）：")
    for meth, tol in [("BDF", 1e-9), ("BDF", 1e-12), ("Radau", 1e-11), ("LSODA", 1e-10)]:
        _, _, td = dry_time(64, 3, method=meth, rtol=tol)
        w("   %-6s rtol=%.0e  tdry = %.3f s" % (meth, tol, td))

w("")
w("表5  药材烘干过程的水分浓度（每 6 h）")
hours = [6, 12, 18, 24, 30, 36, 42, 48, 54]
th = np.array([h * 3600.0 for h in hours] + [td3])
at3 = np.concatenate([np.full(s3.Np, 1e-13), np.full(s3.Np, 1e-11)])
sol3 = s3.solve(td3, t_eval=list(th), rtol=1e-11, atol=at3)
w("t/h      " + "".join("%10s" % ("r=%.1fcm" % (x * 100)) for x in RCOLS))
for i, h in enumerate(hours):
    Cn, _ = s3.interp_u(sol3.y[:, i], U1)
    w("%-8d " % h + "".join("%10.4f" % v for v in Cn))
Cn, _ = s3.interp_u(sol3.y[:, -1], U1)
w("结束     " + "".join("%10.4f" % v for v in Cn))

hdr("10. 问题4 —— 含收缩的烘干时间 t_dry")
w("R(0)=%.4f cm, R(tdry)=%.4f cm, R(259200)=%.4f cm"
  % (cm.R_of(0) * 100, cm.R_of(td3) * 100, cm.R_of(259200) * 100))
w("谱配置法网格收敛性：")
for N in [24, 32, 48, 64, 96]:
    s, sol, td = dry_time(N, 4, shrink=True, rtol=1e-10)
    w("N=%-5d tdry = %11.3f s = %8.4f h = %.4f d" % (N, td, td / 3600, td / 86400))
s4, solE4, td4 = dry_time(64, 4, shrink=True)
w("生产值 N=64: tdry = %.3f s = %.4f h = %.4f d" % (td4, td4 / 3600, td4 / 86400))

w("")
w("表6  药材烘干过程的水分浓度（列为当前物理距离, r_j>R(t) 留空）")
hours4 = [h for h in hours if h * 3600.0 < td4]      # 问题4 的 tdry < 54 h
th4 = np.array([h * 3600.0 for h in hours4] + [td4])
at4 = np.concatenate([np.full(s4.Np, 1e-13), np.full(s4.Np, 1e-11)])
sol4 = s4.solve(td4, t_eval=list(th4), rtol=1e-11, atol=at4)
w("t/h      R/cm   r=0.0     r=0.5     r=1.0     r=1.5     表面")
for i, h in enumerate(hours4):
    Rt = s4.R_of(th4[i])
    row = "%-8d %.3f " % (h, Rt * 100)
    for r in RCOLS[:4]:
        if r <= Rt + 1e-12:
            Cn, _ = s4.interp_u(sol4.y[:, i], (r / Rt) ** 2)
            row += "%10.4f" % Cn[0]
        else:
            row += "%10s" % "—"
    row += "%10.4f" % sol4.y[s4.N, i]
    w(row)

# ======================================================================
hdr("11. 烘干时长对时间离散的敏感性（对照格式为何不可直接引用）")
w("经典 FV+theta 格式，N=1600，长时间段步长放宽到 Δt：")
for d in [10.0, 5.0, 2.0, 1.0, 0.5]:
    plan = [(7200.0, min(1.0, d)), (1e18, d)]
    td = fv_theta_tdry(1600, 3, plan)
    w("   Δt=%-5.2f  tdry = %10.2f s   (相对谱解 %+.2f s)" % (d, td, td - td3))
w("")
w("结论：该格式的 tdry 随 Δt 近似线性变化（每减半步长增量减半），属一阶时间误差；")
w("因此任何在固定粗步长序列上做网格外推得到的值都不是连续极限。")

hdr("12. 关键答案汇总")
w("问题1 中心温度 T(0,1800)   = %.6f  -> %.4f" % (tabT[6][0], tabT[6][0]))
w("问题1 表面温度 T(R,1800)   = %.6f  -> %.4f" % (tabT[6][4], tabT[6][4]))
w("问题1 表面含水 C(R,1800)   = %.6f  -> %.4f" % (tabC[6][4], tabC[6][4]))
w("问题2 中心温度 T(0,3h)     = %.6f  -> %.4f" % (_p2T[5][0], _p2T[5][0]))
w("问题2 中心含水 C(0,3h)     = %.6f  -> %.4f" % (_p2C[5][0], _p2C[5][0]))
w("问题2 表面含水 C(R,3h)     = %.6f  -> %.4f" % (_p2C[5][4], _p2C[5][4]))
w("问题3 t_dry = %.1f s = %.4f h" % (td3, td3 / 3600))
w("问题4 t_dry = %.1f s = %.4f h" % (td4, td4 / 3600))
w("")
w("总耗时 %.1f s" % (time.time() - T0))

open("复跑报告.txt", "w", encoding="utf-8").write("\n".join(OUT))
print("written 复跑报告.txt")
