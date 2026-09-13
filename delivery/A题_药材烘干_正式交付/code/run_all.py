# -*- coding: utf-8 -*-
# ==============================================================================
# run_all.py —— 一键复跑：问题 1~4 的全部结果与全部验证
# ==============================================================================
# 做什么：把四个问题的求解、三族离散格式互证、网格与时间收敛性、烘干时间 t_dry
#         与关键答案汇总串成一趟跑完，过程与结果写成一份文本报告；只算与打印，

#         不生成 xlsx、不画图。
# 输入：common.py（物性、烘房环境、初值；附件 1/2 由它定位）、analytic_p1.py
#       （问题 1 闭式解）、spectral.py（主求解器）、fv_theta.py 与 fv_uniform.py
#       （两个对照离散族）、verify_data.py（只读验证数据，做一致性校验）。
# 输出：复跑报告.txt，写到当前工作目录；xlsx 与插图由别的脚本生成。
# 运行方式与依赖顺序：

#       1) python verify_data.py --recompute   重建验证数据（约 22 min，可跳过）
#       2) python run_all.py                   全部结果与验证（实测约 6.5 min）
#          python run_all.py --fast            再跳过时间步/积分器敏感性扫描
#       3) python make_results.py 生成 result/*.xlsx，make_figures.py 出图
#       出错时异常直接抛出、停在出错那一节；把上面的模块逐个单独跑即可定位。
# 关键变量：

#   CRIT    烘干判据的干基含水率阈值 0.15                      [kg/kg]
#   RCOLS   报告表用的 5 个物理半径 0.0/0.5/1.0/1.5/2.0        [cm]
#   N       谱节点数，扫描 24/32/48/64/96/128，生产值取 64
#   atol    时间积分绝对容差按状态排列 (C_0..C_N, T_0..T_N) 分两段给：水分
#           1e-13、温度 1e-11；rtol 取 1e-9~1e-12
#   FAST    命令行带 --fast 时为 True，跳过最耗时的敏感性扫描

# ==============================================================================

# ---------------- 控制台编码（必须在任何 print 之前） ----------------
# 中文进度信息只有先把 stdout 切到 UTF-8 才能正确输出，故 fix_console() 必须排在
# 本文件所有 print 之前；它只改编码，不参与任何计算。
from console_utf8 import fix_console
fix_console()

# ---------------- 误 import 的防护 ----------------
# 本脚本没有 main()，复跑逻辑直接写在模块顶层，import 即整段执行并改写 复跑报告.txt，
# 故拦一道：误 import 立刻报错，而不是悄悄跑十几分钟、顺手覆盖掉已有报告。
if __name__ != "__main__":
    raise RuntimeError(
        "run_all 的复跑逻辑在模块顶层，import 会立刻重跑全流程并覆盖 复跑报告.txt。"
        "请用 `python run_all.py` 运行。")

import sys
import time
import numpy as np
import common as cm
import analytic_p1 as an
from spectral import HerbSolver
from fv_theta import simulate as fv_theta_run, find_tdry as fv_theta_tdry
from fv_uniform import solve_fv, sample_xi

FAST = "--fast" in sys.argv
OUT = []                  # 报告全文先攒在内存，末尾一次性落盘，中途报错不会留半份文件
T0 = time.time()          # 计时起点，末尾打印总耗时


def w(s=""):
    """写一行报告：同时进内存缓冲与屏幕（flush=True 便于重定向时实时看到进度）。"""
    OUT.append(str(s))
    print(str(s), flush=True)


def hdr(t):
    """打一节标题：空行 + 上下两条 78 字符分隔线，其余与 w 相同。"""
    w("")
    w("=" * 78)
    w(t)
    w("=" * 78)


# 报告表统一用这 5 个物理半径 [m]：0、0.5、1.0、1.5、2.0 cm，即论文表 1~表 6 的 5 列。
RCOLS = [0.0, 0.005, 0.010, 0.015, 0.020]

# ---------------- 第 0 节：烘房环境标定（附件 1） ----------------
# 这里用 common.py 里那套拟合参数把 RMSE、R²、标准误本地重算一遍做独立复核：
# 参数取自 curve_fit 的拟合结果，误差统计自己写（见 _stats），不引用 curve_fit 的 pcov。
hdr("0. 烘房环境标定（附件1，一阶惯性拟合，初值固定为实测值）")
_te, _Ta, _Ca = cm._t_env, cm._Ta_env, cm._Ca_env


def _mT(t, Ts, tau):
    """温度惯性式（与 common._fT 同形，供 _stats 复算残差）：t [s]，返回 [degC]。"""
    return Ts - (Ts - cm.T_air0) * np.exp(-t / tau)


def _mC(t, Cs, tau):
    """含湿量惯性式（与 common._fC 同形）：t [s]，返回空气含湿量 [kg/kg]。"""
    return Cs - (Cs - cm.C_air0) * np.exp(-t / tau)


def _stats(f, p, t, y):
    """算一次拟合的 RMSE、决定系数 R² 与两个参数的近似标准误，返回 (rmse, r2, se)。

    f 是模型函数，p 是参数向量，t 是自变量 [s]，y 是实测值。
    RMSE 用 ÷n 口径：把标定值当作对整段序列的整体描述，与 common.py 报告的口径一致。
    标准误另用 ÷(n-len(p)) 的自由度口径算，因为要的是参数自身的统计不确定度。
    Jacobian 用前向差分数值求，步长取 |p_i|*1e-7 + 1e-12（1e-12 防止参数为 0 时步长
    归零）；标准误取 (J^T J)^{-1} 对角元乘 s2 后开方，即线性化后的参数标准差。
    """
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

# ---------------- 第 1 节：闭式解自检 ----------------
# an._selftest() 返回一段文本：特征方程残差 max|beta*J1-Bi*J0|、t=0 的全场极差
# （即瞬态级数在初值处的截断残差）与 nterm 敏感性表。它不做自动断言，需人眼判读：
# 残差应在机器精度量级（1e-12 上下），t=0 极差应远小于 1e-6，nterm 加到 2 以上后
# T(0,1800) 不应再变。
hdr("1. 问题1 —— 闭式解自检（Bessel-Duhamel）")
w(an._selftest())

# ---------------- 第 2 节：问题 1 谱解的空间收敛性 ----------------
# 只加密空间（N=8~64）：问题 1 是常物性线性问题，谱方法误差随 N 指数下降，故从很粗的
# N=8 起扫也能看出谱收敛；时间容差固定在最紧的 rtol=1e-12（atol 水分 1e-14、温度
# 1e-12），使残差完全由空间离散主导。对照量是 analytic_p1.T_exact 的闭式解。
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

# ---------------- 第 3、4 节：问题 1 的表 1 与表 2 ----------------
# 报告 7 个时刻（100~1800 s，每 300 s 一个点）× 5 个半径，即论文表 1、表 2 的排版；
# 生产网格 N=64，atol 按状态排列 (C_0..C_N, T_0..T_N) 分两段给：水分 1e-13、温度 1e-11。
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

# ---------------- 第 5 节：水分场的网格收敛性 ----------------
# 与第 2 节同一套 N（24~128），但看的是水分场：扩散系数含 exp(-a/C)，含水率一降低该
# 指数项就变得很陡，水分场比温度场难解，故 N 从 24 起而不是从 8 起。同时报 t=1200 s
# 的表面值，用来核对"不同时刻的收敛行为是否一致"。
hdr("5. 问题1 —— 谱解网格收敛性（水分场）")
w("N      C(0,1800)      C(R,1800)      C(R,1200)")
for N in [24, 32, 48, 64, 96, 128]:
    s = HerbSolver(N, prob=1)
    at = np.concatenate([np.full(s.Np, 1e-14), np.full(s.Np, 1e-12)])
    sol = s.solve(1800.0, t_eval=[1200.0, 1800.0], rtol=1e-12, atol=at)
    a = s.interp_u(sol.y[:, 0], (0.02 / cm.R0) ** 2)[0][0]
    b = s.interp_u(sol.y[:, 1], U1)[0]
    w("N=%-5d %.9f   %.9f   %.9f" % (N, b[0], b[-1], a))

# ---------------- 第 6 节：第三族离散的互证（有限体积 + 后向 Euler） ----------------
# 4 组 (N, dt) = (200,1.0)、(200,0.5)、(400,0.5)、(400,0.25)：同一空间网格配两档时间步，
# 后向 Euler 是一阶格式，步长减半时误差也应减半，这是判"离散误差已压到目标精度以下"
# 的依据；三种离散族给出同一个数，才说明它们解的是同一个 PDE。
hdr("6. 问题1 —— 第三种离散族互证（均匀 r 网格守恒型有限体积 + 后向Euler）")
w("N     dt      C(R,1800)      C(R,1200)      T(R,1800)")
for N, dt in [(200, 1.0), (200, 0.5), (400, 0.5), (400, 0.25)]:
    xi, Ctab, Ttab = solve_fv(N, dt, 1800.0, prob=1, t_report=[1200.0, 1800.0])
    c12 = float(np.atleast_1d(sample_xi(xi, Ctab[0], 0.02 / cm.R0))[0])
    c18 = float(np.atleast_1d(sample_xi(xi, Ctab[1], 0.02 / cm.R0))[0])
    t18 = float(np.atleast_1d(sample_xi(xi, Ttab[1], 0.02 / cm.R0))[0])
    w("N=%-4d %-6.2f %12.8f  %12.8f  %12.8f" % (N, dt, c18, c12, t18))
w("谱配置法(N=128):           1.51091810      1.65914239      37.198922")

# ---------------- 第 7 节：问题 2 的表 3 与表 4 ----------------
# 报告 6 个时刻（0.5~3 h，每 0.5 h 一个点）× 5 个半径，即论文表 3、表 4 的排版；
# 问题 2 的物性随 C、T 变化，atol 取水分 1e-12 / 温度 1e-10（比第 3 节松一档）。
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
# FV+theta 的时间方向只有一阶精度，单看某个 dt 的值无法与谱解相比；故取 4 档步长
# （2/1/0.5/0.25 s），用最细两档估一阶误差斜率 c [s/s]，再外推到 dt -> 0 的极限
# （t(dt=0) = t(dt_min) - c*dt_min），之后才与谱解对照。
lad = []
for dt in [2.0, 1.0, 0.5, 0.25]:
    xi, tt, Cs, Ts = fv_theta_run(1600, 2, 1800.0, [(1e18, dt)],
                                  report_idx=[1800.0])
    lad.append((dt, Ts[0][0]))
    w("   dt=%-6.3f %.9f" % (dt, Ts[0][0]))
c = (lad[-2][1] - lad[-1][1]) / (lad[-2][0] - lad[-1][0])
w("   一阶外推 dt->0 : %.9f" % (lad[-1][1] - c * lad[-1][0]))

# ---------------- 问题 3、4 共用的烘干判据与终止事件 ----------------
# CRIT 是题目给的烘干判据：干基含水率处处不高于 0.15 kg/kg。含水率沿半径由中心向外
# 递减，最大值在轴心 u=0（表面先干、中心最后达标），故事件函数取水分段的最大值减 CRIT：
# 它穿过 0 的时刻就是"处处达标"的首达时刻；配 direction=-1 表示只在下降方向触发。
CRIT = 0.15


def dry_time(N, prob, shrink=False, method="BDF", rtol=1e-11, ac=1e-13, at_=1e-11):
    """求烘干时长 t_dry，返回 (求解器实例, 解对象, t_dry [s])。

    N      谱节点数（节点总数为 N+1）；
    prob   问题号，3 或 4，决定用哪套物性经验式；
    shrink 是否用附件 2 的实测半径 R(t)（问题 4 为 True）；
    method 时间积分器名（BDF/Radau/LSODA），rtol 是它的相对容差；
    ac、at_ 分别是水分段与温度段的绝对容差（状态按 (C_0..C_N, T_0..T_N) 排列）。

    积分上限 300000 s 只是"肯定够长"的松上界：事件是终止事件，达标即停，不会白算。
    """
    s = HerbSolver(N, prob=prob, shrink=shrink)
    at = np.concatenate([np.full(s.Np, ac), np.full(s.Np, at_)])

    def ev(t, y):
        return np.max(y[:s.Np]) - CRIT
    ev.terminal = True
    ev.direction = -1
    sol = s.solve(300000.0, rtol=rtol, atol=at, events=ev, method=method)
    return s, sol, float(sol.t_events[0][0])


# ---------------- 第 9 节：问题 3 的烘干时间 ----------------
# 扫描 N=24~128 看 t_dry 是否随网格收敛；生产值取 N=64（与 result3.xlsx 的网格一致），
# 下面各节以及其它脚本报的 t_dry 都以它为准（N=128 的连续极限见 verify_data.py）。
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

# fig22/fig23 的验证数据由 verify_data.py 单独复跑（约 22 min，不适合塞进一键复跑），
# 这里只做一致性校验：数据文件里的 N=64 项必须与上面刚算出的生产值相符。
# 校验不通过也不中断本节，而是把原因打印进报告（见下面的 except 分支），
# 让一次复跑就能看出"图的数据"与"这次算的解"到底哪边过期了。
w("")
w("fig22/fig23 验证数据一致性校验（中间数据/verification.json）：")
try:
    import verify_data as VD
    VD.check_consistency(VD.load(verbose=False), tdry_cached=td3, verbose=False)
    w("   通过：数据文件与本次复跑一致（N=64 项 %.4f s）" % td3)
except Exception as exc:                     # 数据缺失或不一致都要显式报出来
    w("   ！！未通过：%s" % exc)

if not FAST:
    w("")
    w("时间积分器/容差敏感性（验证 tdry 不是数值噪声）：")
    # 4 个配置跨 3 个积分器族（BDF 两档容差、Radau、LSODA）：若它们的 t_dry 只差
    # 百分之一秒，说明这个时刻由方程本身决定，不是积分器或容差挑出来的。
    for meth, tol in [("BDF", 1e-9), ("BDF", 1e-12), ("Radau", 1e-11), ("LSODA", 1e-10)]:
        _, _, td = dry_time(64, 3, method=meth, rtol=tol)
        w("   %-6s rtol=%.0e  tdry = %.3f s" % (meth, tol, td))

w("")
w("表5  药材烘干过程的水分浓度（每 6 h）")
# 论文表 5 的行：6~54 h 每 6 h 一行，末行再补 t_dry 本身（刚好达标的那一行）。
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

# ---------------- 第 10 节：问题 4 的烘干时间（含实测收缩） ----------------
# 扫描 N=24~96（比问题 3 少最细的 N=128 一档）；半径历史来自附件 2，由 R_of 分段线性
# 插值取值，故这里的 t_dry 同时含物性差异与几何收缩，不能与问题 3 直接相减当收缩效应。
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

# ---------------- 第 11 节：对照格式的时间步敏感性 ----------------
# 步长计划 [(7200 s, min(1, Δt)), (1e18 s, Δt)]：前 7200 s 恒定用不超过 1 s 的小步长，
# 之后才用 Δt，这样扫到的是恒温干燥段的步长误差，不与升温段混在一起。
# Δt 从 10 s 减到 0.5 s；若 t_dry 随步长近似线性变化，该格式就是 O(Δt) 一阶，
# 任何在固定粗步长序列上做的外推都不代表连续极限。
hdr("11. 烘干时长对时间离散的敏感性（对照格式为何不可直接引用）")
w("经典 FV+theta 格式，N=1600，长时间段步长放宽到 Δt：")
for d in [10.0, 5.0, 2.0, 1.0, 0.5]:
    plan = [(7200.0, min(1.0, d)), (1e18, d)]
    td = fv_theta_tdry(1600, 3, plan)
    w("   Δt=%-5.2f  tdry = %10.2f s   (相对谱解 %+.2f s)" % (d, td, td - td3))
w("")
w("结论：该格式的 tdry 随 Δt 近似线性变化（每减半步长增量减半），属一阶时间误差；")
w("因此任何在固定粗步长序列上做网格外推得到的值都不是连续极限。")

# ---------------- 第 12 节：关键答案汇总 ----------------
# 把前面各节算出的、论文要引用的数字集中重打一遍，便于与正文表格逐个对位。
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

