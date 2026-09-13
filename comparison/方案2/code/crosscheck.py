# -*- coding: utf-8 -*-
# ============================================================================
# crosscheck.py —— 用不同实现互证干燥模型，排查「同一份代码里的同一个错误」
#
# 三条路线的独立性并不相同，看结论时要分清：
#   (A) 节点中心有限体积 + scipy.integrate.solve_ivp(BDF)：控制体划分、算术平均
#       面扩散系数与表面通量表达式都与生产代码相同，只有时间积分器不同（自适应
#       BDF vs theta 法），故它排除不了空间离散的系统性错误；生产代码另有的
#       单元中心布局（grid='cell'）才是真正不同的空间离散，由 resolution.py 对照。
#   (B) 第一特征值降维模型：dCbar/dt = -lambda1(Cbar)(Cbar - Cair)，
#       lambda1 = mu1^2 D(Cbar,T)/R^2，mu1 为 mu*J1(mu)=Bi_m*J0(mu) 的最小正根。
#       它完全不含空间离散，用截面平均含水与完整 PDE 对照，可量化「拿 D(Cbar)
#       代替剖面上平均的 <D>」带来的偏差。
#   (C) 生产有限体积解自身的网格/步长加密序列，看观测收敛阶。
#
# 输入数据：环境与烘房历程来自附件 1（common.RoomConditions 已拟合成一阶惯性式），
#           收缩半径来自附件 2（RadiusHistory），物性用附录 2 / 3 / 4 的经验关系。
# 输出：data/crosscheck_summary.txt（控制台报告全文）、
#       data/crosscheck_convergence.csv（各布局/网格/步长的收敛表）。
#
# 关键变量：
#   C_TARGET = 0.15 kg/kg  干燥终点判据：全场含水率首次全部低于该值
#   n_node / n_cell  直线法节点数与有限体积控制体数；t_eval / t_report 采样时刻 [s]
# ============================================================================
"""
运行方式：在目录 A题解答 下执行 python code/crosscheck.py。

报告分三节：A 直线法（BDF）结果、B 降维模型结果、C 生产求解器的网格/步长收敛表。
前两节的抽样数值只是给人看的中间量，真正的判据是同一时刻两套实现的差值。
"""

import os
import sys
import time
import io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import j0, j1
from scipy.optimize import brentq
import pandas as pd

if sys.stdout.encoding is None or "utf" not in sys.stdout.encoding.lower():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm

C_TARGET = 0.15


# ==========================================================================
# A) 路线 A：节点中心有限体积 + scipy BDF 自适应时间积分
# ==========================================================================
def _face_mean(a, b):
    """相邻两节点的面扩散系数，取算术平均 (a+b)/2。

    与生产求解器 common._face_mean 用同一种面系数，这条路线才构成「同一空间
    离散、不同时间积分器」的对照；取算术平均的理由见该函数的说明。
    """
    return 0.5 * (a + b)


def mol_solve(prop, t_end, n_node=300, room=None, radius=None,
              rtol=1e-9, atol=1e-11, t_eval=None):
    """直线法（MOL）求解：径向离散成常微分方程组，再交给 BDF 自适应积分。

    n_node 为节点数（节点中心布局，含 xi=0 与 xi=1 两端），t_end 为终止时刻 [s]，
    t_eval 为输出时刻数组 [s]（默认每 60 s 一点），rtol/atol 控制 BDF 的误差。
    状态向量前半是各节点含水率、后半是各节点温度。返回 dict：t 时刻、xi 节点
    位置、C 与 T（形状均为 (len(t), n+1)）、nfev 右端函数调用次数。
    """
    room = room if room is not None else cm.RoomConditions()
    n = int(n_node)
    h = 1.0 / n
    nv = n + 1
    xi = np.linspace(0.0, 1.0, nv)
    xi_f = (np.arange(n) + 0.5) * h                # 内部面位置 xi=(j+0.5)h

    vol = np.empty(nv)                             # 控制体体积 int xi dxi
    vol[0] = h * h / 8.0
    vol[1:n] = xi[1:n] * h
    vol[n] = h / 2.0 - h * h / 8.0

    def R_of(t):
        return cm.R0 if radius is None else float(radius.R(t))

    def rhs(t, y):
        C = y[:nv]
        T = y[nv:]
        rho, cp, k, D = prop(np.maximum(C, cm.C_FLOOR), T + cm.T_KELVIN)
        R = R_of(t)
        t_air = float(room.T_air(t))
        c_air = float(room.C_air(t))

        def field_rhs(phi, a, cap, coef, env):
            af = _face_mean(a[:-1], a[1:])
            gf = xi_f * af * (phi[1:] - phi[:-1]) / h
            g_s = -R * coef * (phi[n] - env)
            net = np.empty(nv)
            net[0] = gf[0]
            net[1:n] = gf[1:] - gf[:-1]
            net[n] = g_s - gf[n - 1]
            return net / (R * R * vol * cap)

        dC = field_rhs(C, D, np.ones(nv), cm.HM_CONV, c_air)
        dT = field_rhs(T, k, rho * cp, cm.H_CONV, t_air)
        return np.concatenate([dC, dT])

    y0 = np.concatenate([np.full(nv, cm.C_INIT), np.full(nv, cm.T_INIT_C)])
    if t_eval is None:
        t_eval = np.arange(0.0, t_end + 1.0, 60.0)
    sol = solve_ivp(rhs, (0.0, t_end), y0, method="BDF", t_eval=t_eval,
                    rtol=rtol, atol=atol, jac_sparsity=_sparsity(nv))
    if not sol.success:
        raise RuntimeError(sol.message)
    return dict(t=sol.t, xi=xi, C=sol.y[:nv].T, T=sol.y[nv:].T,
                nfev=sol.nfev)


def _sparsity(nv):
    """列出状态向量 (C_0..C_n, T_0..T_n) 的块三对角稀疏结构。

    BDF 是隐式方法，每步都要解以 Jacobian 为系数阵的线性系统。只有相邻节点
    之间有耦合（每个未知量连着同场的左右邻居与另一场的同号节点），把这一结构
    告诉 solve_ivp 就能用稀疏分解代替稠密分解，长时间积分才跑得动。
    """
    import scipy.sparse as sp
    rows, cols = [], []
    for i in range(nv):
        for j in (i - 1, i, i + 1):
            if 0 <= j < nv:
                for a in (0, nv):
                    for b in (0, nv):
                        rows.append(a + i)
                        cols.append(b + j)
    S = sp.lil_matrix((2 * nv, 2 * nv))
    S[np.array(rows), np.array(cols)] = 1
    return sp.csr_matrix(S)


def mol_drying_time(prop, room=None, radius=None, n_node=300, t_end=3.6e5):
    """用直线法结果求全场含水率首次降到 C_TARGET 以下的时刻 [s]。

    采样间隔 300 s，在跨越目标值的相邻两点间线性插值以提高时间分辨率。
    返回 (t_dry, 直线法结果)；若 t_end 内始终没有降到目标值，则 t_dry 为 nan。
    """
    res = mol_solve(prop, t_end, n_node=n_node, room=room, radius=radius,
                    t_eval=np.arange(0.0, t_end + 1.0, 300.0))
    Cmax = res["C"].max(axis=1)
    idx = np.where(Cmax < C_TARGET)[0]
    if not len(idx):
        return np.nan, res
    i0 = idx[0]
    t0, t1 = res["t"][i0 - 1], res["t"][i0]
    v0, v1 = Cmax[i0 - 1], Cmax[i0]
    return t0 + (C_TARGET - v0) * (t1 - t0) / (v1 - v0), res


# ==========================================================================
# B) 路线 B：第一特征值降维模型（只看截面平均含水率）
# ==========================================================================
def first_root(bi):
    """求 mu*J1(mu) = bi*J0(mu) 的最小正根 mu1（第一特征值）。

    bi 为传质 Biot 数。该根必落在 (0, 2.4048) 内（2.4048 是 J0 的第一个零点），
    故以此为搜索区间右端用 brentq 求根；bi 越小根越靠近 0。
    """
    f = lambda mu: mu * j1(mu) - bi * j0(mu)
    hi = 2.404825557695773 - 1e-9
    return brentq(f, 1e-10, hi, xtol=1e-15)


def lumped_mean(prop, room=None, radius=None, t_end=3.0e5, dt=30.0):
    """积分第一特征值模型，得到截面平均含水率随时间变化的曲线。

    prop 为物性函数，radius 为收缩半径历史（None 表示半径恒为 R0），
    t_end 为终止时刻 [s]，dt 为显式 Euler 步长 [s]（30 s 远小于 1/lambda1 的
    1e4 s 量级，稳定性没有问题）。返回 (ts, cs)：时刻 [s] 与平均含水率 [kg/kg]。
    """
    room = room if room is not None else cm.RoomConditions()
    c = cm.C_INIT
    ts, cs = [0.0], [c]
    for _ in range(int(t_end / dt)):
        t = ts[-1] + dt
        R = cm.R0 if radius is None else float(radius.R(t))
        # 准瞬时热平衡：药材很细，热扩散比水分扩散快一个多数量级，故直接把环境
        # 温度当作药材温度来算 D；时间下限 1200 s 是为了跳过预热最初空气温度
        # 尚未建立的一小段，避免用明显偏低的环境温度算 D。
        T_herb = float(room.T_air(max(t, 1200.0)))
        _, _, _, D = prop(np.array([c]), np.array([T_herb + cm.T_KELVIN]))
        D = float(D[0])
        # Bi_m = h_m R / D，下限 1e-6 只是防止 D 极小时除零
        bi = max(cm.HM_CONV * R / max(D, 1e-30), 1e-6)
        mu = first_root(bi)
        lam = mu ** 2 * D / (R * R)      # 第一特征值 lambda1 = mu1^2 D / R^2 [1/s]
        c_eq = float(room.C_air(t))      # 平衡含水率取环境含湿量
        c = c + dt * (-lam * (c - c_eq))
        if c < 0.0:
            c = 0.0                      # 显式 Euler 可能过冲到负值，截断到 0
        ts.append(t)
        cs.append(c)
    return np.array(ts), np.array(cs)


def shape_factor(bi):
    """圆柱第一径向模态的形函数 C_center / C_mean = mu/(2*J1(mu))。

    只保留第一模态时可写 C(r) = C_center*J0(mu*r/R)，对截面做加权平均后即得
    该比值，用于把「中心含水率判据」与「截面平均含水率判据」互相换算。
    """
    mu = first_root(bi)
    return mu / (2.0 * j1(mu))


# ==========================================================================
def main():
    """跑完三节，把报告写入 crosscheck_summary.txt 并导出收敛表。

    P(s) 是打印的包装：既在控制台输出，又把同一行收进 out，最后一次性写文件。
    """
    out = []
    P = lambda s: (print(s), out.append(s))
    room = cm.RoomConditions()
    rad = cm.RadiusHistory(mode="table")

    P("=" * 78)
    P("A) INDEPENDENT METHOD: node-centred finite volumes + scipy BDF")
    P("=" * 78)
    # 三组算例：(标签, 物性函数, 收缩半径历史, 终止时刻 [s], 采样间隔 [s])
    cases = [("problem 1   (Appendix 2, 0-1800 s)", cm.props_problem1, None,
              1800.0, 300.0),
             ("problem 3   (Appendix 3)", cm.props_problem23, None,
              2.6e5, 1800.0),
             ("problem 4   (Appendix 4, shrinkage)", cm.props_problem4, rad,
              2.6e5, 1800.0)]
    mol_res = {}
    for tag, prop, radius, t_end, dt_ev in cases:
        t0 = time.time()
        res = mol_solve(prop, t_end, n_node=300, room=room, radius=radius,
                        t_eval=np.arange(0.0, t_end + 1.0, dt_ev))
        mol_res[tag] = res
        P("  %-38s nodes=301  wall %.1f s  nfev=%d"
          % (tag, time.time() - t0, res["nfev"]))
        # 只打印 0、1/4、1/2、终点四个时刻，用来看整体趋势
        pick = np.array([0.0, t_end / 4.0, t_end / 2.0, t_end])
        idxs = [int(np.argmin(np.abs(res["t"] - p))) for p in pick]
        P("      t/s        : %s" % np.round(res["t"][idxs], 0))
        P("      C(centre)  : %s" % np.round(res["C"][idxs, 0], 5))
        P("      C(surface) : %s" % np.round(res["C"][idxs, -1], 5))
        if t_end > 1e5:
            # 取全场最大含水率（此处为轴上节点）首次降到 C_TARGET 的时刻
            Cmax = res["C"].max(axis=1)
            idx = np.where(Cmax < C_TARGET)[0]
            if len(idx):
                i0 = idx[0]
                tt = res["t"][i0 - 1] + (C_TARGET - Cmax[i0 - 1]) * \
                    (res["t"][i0] - res["t"][i0 - 1]) / (Cmax[i0] - Cmax[i0 - 1])
                P("      MOL drying time : %.1f s = %.4f h = %.4f d"
                  % (tt, tt / 3600.0, tt / 86400.0))

    P("")
    P("=" * 78)
    P("B) REDUCED FIRST-EIGENVALUE MODEL (mean moisture)")
    P("=" * 78)
    for tag, prop, radius in [("problem 3 (Appendix 3)", cm.props_problem23, None),
                              ("problem 4 (Appendix 4)", cm.props_problem4, rad)]:
        ts, cs = lumped_mean(prop, room=room, radius=radius)
        # 降维模型给的是截面平均含水率，阈值同样取 C_TARGET=0.15（偏保守，
        # 因为平均值总是先于中心点降到该值）
        idx = np.where(cs <= C_TARGET)[0]
        t_lumped = ts[idx[0]] if len(idx) else np.nan
        P("  %-30s t(Cbar=0.15) = %.1f s = %.4f h" % (tag, t_lumped,
                                                      t_lumped / 3600.0))
        for tn in (6.0, 12.0, 24.0, 48.0):
            j = int(np.argmin(np.abs(ts - tn * 3600.0)))
            P("      Cbar(%.0f h) = %.5f" % (tn, cs[j]))

    P("")
    P("=" * 78)
    P("C) MESH / TIME-STEP CONVERGENCE OF THE PRODUCTION SOLVER")
    P("=" * 78)
    # 每个问题都跑同一组 (布局, 控制体数, 步长)：加密网格与减小步长时，中心/表面
    # 含水率与干燥时间若稳定下来，说明离散误差已不影响结论。
    rows = []
    for tag, prop, radius in [("P1", cm.props_problem1, None),
                              ("P3", cm.props_problem23, None),
                              ("P4", cm.props_problem4, rad)]:
        for grid, n_cell, dt in [("node", 200, 1.0), ("node", 400, 1.0),
                                 ("node", 800, 1.0), ("node", 1600, 1.0),
                                 ("node", 800, 2.0), ("node", 800, 5.0),
                                 ("cell", 800, 1.0), ("cell", 1600, 1.0)]:
            # 问题 1 只算到 1800 s，远没到干燥终点，故只记中心与表面值；
            # 问题 3、4 还要在 300 s 采样序列上插值出 t_dry
            if tag == "P1":
                res = cm.simulate(n_cell, prop, 1800.0, dt, room, grid=grid,
                                  t_report=np.array([1800.0]),
                                  xi_report=np.array([0.0, 0.5, 1.0]))
                rows.append(dict(problem=tag, grid=grid, n_cell=n_cell, dt=dt,
                                 T_center=res["T"][0, 0],
                                 C_center=res["C"][0, 0],
                                 C_surf=res["C"][0, 2], t_dry=np.nan))
            else:
                t_end = 2.6e5
                res = cm.simulate(n_cell, prop, t_end, dt, room, radius=radius,
                                  grid=grid,
                                  t_report=np.arange(0.0, t_end + 1.0, 300.0),
                                  xi_report=np.array([0.0, 0.5, 1.0]))
                Cmax = res["C"].max(axis=1)
                idx = np.where(Cmax < C_TARGET)[0]
                t_dry = np.nan
                if len(idx):
                    i0 = idx[0]
                    t_dry = res["t"][i0 - 1] + (C_TARGET - Cmax[i0 - 1]) * \
                        (res["t"][i0] - res["t"][i0 - 1]) / \
                        (Cmax[i0] - Cmax[i0 - 1])
                rows.append(dict(problem=tag, grid=grid, n_cell=n_cell, dt=dt,
                                 T_center=np.nan, C_center=Cmax[-1],
                                 C_surf=res["C"][-1, 2], t_dry=t_dry))
    df = pd.DataFrame(rows)
    P(df.to_string(index=False))
    df.to_csv(os.path.join(cm.DATA_DIR, "crosscheck_convergence.csv"),
              index=False)

    with open(os.path.join(cm.DATA_DIR, "crosscheck_summary.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("\n   saved ->", os.path.join(cm.DATA_DIR, "crosscheck_summary.txt"))


if __name__ == "__main__":
    main()
