# -*- coding: utf-8 -*-
# ============================================================================
# verify_solver.py —— 数值格式的正确性验证：解析解对照 + 离散守恒性核对
#
# (A) 与解析解对照。常物性、恒环境条件下，无限长圆柱带 Robin 边界的非稳态
#     扩散问题存在 Bessel 级数解：
#        (phi-phi_inf)/(phi0-phi_inf) = sum_n A_n J0(mu_n r/R) exp(-mu_n^2 a t/R^2)
#     其中 mu_n 为 mu*J1(mu) = Bi*J0(mu) 的正根（传热 Bi = hR/k，传质用 h_m R/D），
#     A_n = 2 J1(mu_n) / (mu_n (J0(mu_n)^2 + J1(mu_n)^2))。
#     把有限体积解与级数解逐点相减，看误差随位置与网格的变化。
# (B) 离散守恒性核对。把格式对所有控制体求和，理论上应满足「内部储量变化率
#     = 表面通量」；脚本把冻结 D 算例的表面通量沿时间积分，再与截面储量的实际
#     减少量比较，用于检验离散格式是否漏掉通量。
#
# 输入数据：本脚本自己不读附件（import common 时会先确认题目附件目录存在）。
#           参数取自 common.py 中题目给定的常量：温度场用问题 1 的常数物性，
#           水分场把 D 冻结成常数，环境取恒定的 60 degC / 0.02 kg/kg —— 三者都
#           取常数，才能构造出与解析解同样的常系数问题。
# 输出：只打印控制台报告，本脚本不写文件。
#
# 关键变量：
#   alpha = k/(rho*cp) 热扩散系数 [m2/s]；bi 传热 Biot 数；roots 前几个 Bessel 根
#   times / xi_out  采样时刻 [s] 与归一化半径位置；n_root 级数保留的项数
#   D_frozen  水分扩散系数冻结值 [m2/s]；flux_int 表面水分累计流出量
# ============================================================================
"""
运行方式：在目录 A题解答 下执行 python code/verify_solver.py。

控制台报告分三段：解析解逐点误差、冻结 D 的水分守恒核对、问题 1 模型的
网格与步长收敛。
"""

import os
import sys
import numpy as np
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import brentq
import io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm


# --------------------------------------------------------------------------
def robin_roots(bi, n_root=40):
    """求 mu*J1(mu) - bi*J0(mu) = 0 的正根（bi 为 Biot 数）。

    每个根恰好夹在 J0 的相邻两个零点之间，且 J0 的第一个零点之前还有一根，
    故先取 n_root+2 个 J0 零点作为搜索区间端点，再逐区间用 brentq 求根。
    实际只取到 n_root 个根，返回按从小到大排好的一维数组。
    """
    f = lambda mu: mu * j1(mu) - bi * j0(mu)
    z = jn_zeros(0, n_root + 2)
    edges = np.concatenate(([0.0], z))
    roots = []
    for a, b in zip(edges[:-1], edges[1:]):
        a_ = a + 1.0e-10
        b_ = b - 1.0e-10
        fa, fb = f(a_), f(b_)
        if fa * fb < 0:
            roots.append(brentq(f, a_, b_, xtol=1e-15, rtol=1e-15))
        if len(roots) >= n_root:
            break
    return np.array(roots)


def analytical_cylinder(r, t, phi0, phi_inf, alpha, bi, n_root=80):
    """抛 Bessel 级数解析解，返回形状 (len(t), len(r)) 的数组。

    r 为物理半径 [m]（内部除以 R0 换成归一化半径），t 为时刻 [s]，phi0 为初始
    均匀值，phi_inf 为环境平衡值，alpha 为扩散系数 [m2/s]，bi 为 Biot 数。
    三项下标 r/t/n 靠广播一次算完所有分量的和。
    """
    r = np.atleast_1d(np.asarray(r, dtype=float))
    t = np.atleast_1d(np.asarray(t, dtype=float))
    mu = robin_roots(bi, n_root)
    A = 2.0 * j1(mu) / (mu * (j0(mu) ** 2 + j1(mu) ** 2))
    rr = (r / cm.R0)[None, :, None]                    # (1, nr, 1)
    tt = t[:, None, None]                              # (nt, 1, 1)
    terms = A[None, None, :] * j0(mu[None, None, :] * rr) \
        * np.exp(-(mu ** 2)[None, None, :] * alpha * tt / cm.R0 ** 2)
    return phi_inf + (phi0 - phi_inf) * terms.sum(axis=-1)


# --------------------------------------------------------------------------
def constant_heat_props(D_const_value):
    """返回一个物性函数 prop(C, T_K) -> (rho, cp, k, D)。

    rho、cp、k 固定取问题 1 的常数值，D 冻结为参数给定的 D_const_value [m2/s]，
    这样水分方程也变成常系数，才可能和解析解对照。仅供本脚本的验证算例使用。
    """

    def prop(C, T_K):
        C = np.maximum(np.asarray(C, dtype=float), 1e-9)
        rho = np.full_like(C, 820.0)
        cp = np.full_like(C, 2600.0)
        k = np.full_like(C, 0.36)
        D = np.full_like(C, D_const_value)
        return rho, cp, k, D

    return prop


class ConstantRoom:
    """把烘房环境固定成常数的替身类，接口与 RoomConditions 一致。

    只提供 T_air(t)、C_air(t) 两个查询，且都返回与时间无关的常数：
    t_air 单位 degC，c_air 单位 kg/kg。生产用的 RoomConditions 环境随时间变化，
    而解析解只在恒环境假设下成立，故验证时换成本类。
    """

    def __init__(self, t_air, c_air):
        self.t_air = t_air
        self.c_air = c_air

    def T_air(self, t):
        return np.full_like(np.asarray(t, dtype=float), self.t_air)

    def C_air(self, t):
        return np.full_like(np.asarray(t, dtype=float), self.c_air)


def main():
    """跑完三段验证并打印报告：解析解对照、离散守恒核对、问题 1 的收敛性。"""
    print("=" * 74)
    print("A) ANALYTICAL VALIDATION : constant properties, constant ambient")
    print("=" * 74)

    alpha = 0.36 / (820.0 * 2600.0)
    bi = cm.H_CONV * cm.R0 / 0.36
    roots = robin_roots(bi, 8)
    print("   alpha_heat = %.6e m2/s      Bi_heat = %.6f" % (alpha, bi))
    print("   first 5 roots mu_n :", np.round(roots[:5], 6))

    # ---- 温度场：常物性 + 恒环境，与级数解逐点比较 ----
    room = ConstantRoom(60.0, 0.02)
    D_frozen = 5.0e-9
    times = np.array([60.0, 300.0, 900.0, 1800.0, 3600.0])
    xi_out = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0])
    print("\n  --- heat field (T_air = 60 C, T0 = 28 C) ---")
    for n_cell, dt in [(100, 1.0), (200, 1.0), (400, 1.0), (800, 1.0),
                       (400, 2.0), (400, 0.5), (400, 0.25)]:
        res = cm.simulate(n_cell, cm.props_problem1, 3600.0, dt, room,
                          t_report=times, xi_report=xi_out)
        # 级数解的自变量是物理半径，故先把归一化位置 xi 换算回 r = xi*R0
        exact = analytical_cylinder(xi_out * cm.R0, times, cm.T_INIT_C, 60.0,
                                    alpha, bi)
        err = np.abs(res["T"] - exact)
        i_t, i_r = np.unravel_index(np.argmax(err), err.shape)
        print("   N=%4d dt=%5.2f s  max|err| = %.3e K at (t=%.0f s, xi=%.2f)"
              "   centre@1800s err = %.2e"
              % (n_cell, dt, err.max(), times[i_t], xi_out[i_r], err[3, 0]))

    # ---- 水分场：把 D 冻结成常数，否则 D 随 C 变化就没有解析解可用 ----
    D = D_frozen
    bi_m = cm.HM_CONV * cm.R0 / D
    print("\n  --- moisture field with D frozen at %.1e m2/s "
          "(Bi_m = %.4f) ---" % (D, bi_m))
    prop_const = constant_heat_props(D)
    for n_cell, dt in [(100, 1.0), (200, 1.0), (400, 1.0), (800, 1.0),
                       (400, 0.25)]:
        res = cm.simulate(n_cell, prop_const, 3600.0, dt, room,
                          t_report=times, xi_report=xi_out)
        # 水分侧同样要与解析解对齐：初值 C_INIT，环境含湿量就是 room 的 0.02
        exact = analytical_cylinder(xi_out * cm.R0, times, cm.C_INIT, 0.02,
                                    D, bi_m)
        err = np.abs(res["C"] - exact)
        print("   N=%4d dt=%5.2f s  max|err| = %.3e kg/kg" % (n_cell, dt, err.max()))

    # ---- 离散守恒性：表面通量的时间积分 vs 截面储量的实际减少量 ----
    print("\n  --- discrete water balance, D frozen (flux integral vs storage) ---")
    # 节点中心与单元中心两种布局都查：两者的表面值取法不同，后者的表面通量
    # 只有二阶外推，守恒误差理应更大，正好可以互相印证。
    for n_cell, grid in ((100, "node"), (200, "node"), (400, "node"),
                         (800, "node"), (200, "cell"), (800, "cell")):
        slv = cm.DryingSolver(n_cell, prop_const, theta=0.5, room=room,
                              grid=grid)
        slv.reset()
        t, dt = 0.0, 1.0
        flux_int = 0.0
        while t < 3600.0 - 1e-9:
            # 节点中心格式下最后一个节点就落在表面 xi=1 上，值是精确的；
            # 单元中心格式没有表面节点，只能用最后两点的线性外推
            C_s = slv.C[-1] if grid == "node" else \
                1.5 * slv.C[-1] - 0.5 * slv.C[-2]
            flux_int += 2.0 * np.pi * cm.R0 * cm.HM_CONV * (C_s - 0.02) * dt
            slv.step(dt)
            t += dt
        # W 为单位长度圆柱内的水量：2*pi*R0^2 乘上积分 int C*xi dxi（离散后即
        # sum(C*vol)）；W0 是同一积分在初始均匀场 C_INIT 下的值。
        W = 2.0 * np.pi * cm.R0 ** 2 * float((slv.C * slv.vol).sum())
        W0 = np.pi * cm.R0 ** 2 * cm.C_INIT
        print("   N=%4d %-4s  flux integral = %.10e   storage change = %.10e"
              "   rel.diff = %.2e" % (n_cell, grid, flux_int, W0 - W,
                                     abs(flux_int - (W0 - W)) / W0))

    print("\n" + "=" * 74)
    print("B) CONVERGENCE OF THE REAL PROBLEM 1 MODEL (node-centred grid)")
    print("=" * 74)
    # 这一段换回真实的烘房历程（预热平衡 + 恒温干燥两阶段）与问题 1 的实际
    # 物性，不再有解析解可比，只能看网格加密与步长减半时结果是否稳定。
    room1 = cm.RoomConditions()
    xi_out = np.array([0.0, 0.25, 0.5, 0.75, 0.875, 1.0])
    print("   %-8s %-8s %-14s %-14s %-14s" %
          ("N", "dt/s", "T(r=0,1800)", "C(r=0,1800)", "C(r=R,1800)"))
    ref = None   # 最细网格（N=1600）的解，留作下面衡量 N=800 偏差的参考
    for n_cell, dt in [(100, 1.0), (200, 1.0), (400, 1.0), (800, 1.0),
                       (1600, 1.0), (800, 0.5)]:
        res = cm.simulate(n_cell, cm.props_problem1, 1800.0, dt, room1,
                          t_report=np.array([1800.0]), xi_report=xi_out)
        line = "   %-8d %-8.2f %-14.6f %-14.6f %-14.6f" % (
            n_cell, dt, res["T"][0, 0], res["C"][0, 0], res["C"][0, -1])
        if n_cell == 1600:
            ref = res
        print(line)
    if ref is not None:
        res = cm.simulate(800, cm.props_problem1, 1800.0, 1.0, room1,
                          t_report=np.array([1800.0]), xi_report=xi_out)
        print("   difference N=800 vs N=1600 : max|dT| = %.2e K, max|dC| = %.2e"
              % (np.abs(res["T"] - ref["T"]).max(),
                 np.abs(res["C"] - ref["C"]).max()))


if __name__ == "__main__":
    main()
