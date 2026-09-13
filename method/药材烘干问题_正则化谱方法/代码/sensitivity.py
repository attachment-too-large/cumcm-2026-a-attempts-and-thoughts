# -*- coding: utf-8 -*-
"""参数扰动灵敏度分析（问题3）。

做法：把物性经验式里的某个常数乘上一个扰动因子，重新求解到烘干判据触发，
记录烘干时长相对于基准值的变化。

实现：主求解器内部调用的是 common.props 与模块级常数 h、hm、Tset。
本脚本先把它们保存下来，换成带扰动的版本，求解完再恢复。
主求解器一行都不用改。这种做法常被称为猴子补丁。

用法：python sensitivity.py
输出：灵敏度结果.txt 与本文件的打印
"""
import io
import os
import time
import numpy as np
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm
from spectral import HerbSolver

N = 64            # 空间节点参数
PROB = 3          # 问题3
CRIT = 0.15       # 烘干判据
OUT = []


def w(s=""):
    OUT.append(str(s))
    print(s, flush=True)


def tdry(N=N, prob=PROB, rtol=1e-11):
    """返回烘干时长（秒）。"""
    s = HerbSolver(N, prob=prob)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])

    def ev(t, y):
        return np.max(y[:s.Np]) - CRIT
    ev.terminal = True
    ev.direction = -1
    sol = s.solve(300000.0, rtol=rtol, atol=at, events=ev)
    return float(sol.t_events[0][0])


_orig_props = cm.props


def make_props(D0_f=1.0, a_f=1.0, E_f=1.0):
    """返回一个带扰动的物性函数。

    D0_f 前因子倍数；a_f 浓度指数倍数；E_f 活化温度倍数。
    只有含 Arrhenius 项的问题（2 和 3）才受影响。
    """
    def props(prob, C, T):
        rho, cp, k, D = _orig_props(prob, C, T)
        if prob == 1:
            return rho, cp, k, D
        C = np.asarray(C, dtype=float)
        T = np.asarray(T, dtype=float)
        Cf = np.maximum(C, cm.C_FLOOR)
        TK = T + 273.15
        a = 0.45 if prob in (2, 3) else 0.30
        E = 3850.0
        fac = np.exp(-a * (a_f - 1.0) / Cf) * np.exp(-E * (E_f - 1.0) / TK) * D0_f
        return rho, cp, k, D * fac
    return props


def run_case(label, D0_f=1.0, a_f=1.0, E_f=1.0, hm_f=1.0, h_f=1.0, Tset_d=0.0):
    """扰动若干参数，求解一次，返回烘干时长。"""
    cm.props = make_props(D0_f, a_f, E_f)
    old_hm, old_h, old_Tset = cm.hm, cm.h, cm.Tset
    cm.hm = old_hm * hm_f
    cm.h = old_h * h_f
    cm.Tset = old_Tset + Tset_d
    t0 = time.time()
    td = tdry()
    cm.props = _orig_props
    cm.hm, cm.h, cm.Tset = old_hm, old_h, old_Tset
    w("  %-30s tdry = %11.3f s   用时 %.1f s" % (label, td, time.time() - t0))
    return td


def main():
    w("参数扰动灵敏度（问题3，N=%d，BDF，rtol=1e-11）" % N)
    w("=" * 74)
    base = run_case("基准")
    w("")
    rows = [
        ("扩散系数前因子 D0 加百分之二", run_case("D0 +2%", D0_f=1.02)),
        ("扩散系数前因子 D0 减百分之二", run_case("D0 -2%", D0_f=0.98)),
        ("浓度指数 a 加百分之二", run_case("a +2%", a_f=1.02)),
        ("浓度指数 a 减百分之二", run_case("a -2%", a_f=0.98)),
        ("活化温度 E 加百分之二", run_case("E +2%", E_f=1.02)),
        ("活化温度 E 减百分之二", run_case("E -2%", E_f=0.98)),
        ("传质系数 hm 加百分之二", run_case("hm +2%", hm_f=1.02)),
        ("传质系数 hm 减百分之二", run_case("hm -2%", hm_f=0.98)),
        ("换热系数 h 加百分之二", run_case("h +2%", h_f=1.02)),
        ("换热系数 h 减百分之二", run_case("h -2%", h_f=0.98)),
        ("烘房设定温度 Tset 加 0.2 度", run_case("Tset +0.2", Tset_d=0.2)),
        ("烘房设定温度 Tset 减 0.2 度", run_case("Tset -0.2", Tset_d=-0.2)),
    ]
    w("")
    w("=" * 74)
    w("%-30s %13s %12s" % ("扰动", "tdry / s", "相对变化"))
    w("-" * 74)
    for label, td in sorted(rows, key=lambda r: -abs(r[1] - base)):
        w("%-30s %13.3f %11.4f%%" % (label, td, (td - base) / base * 100.0))
    w("-" * 74)
    w("基准 tdry = %.3f s" % base)

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "灵敏度结果.txt")
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT))
    print("written", path)


if __name__ == "__main__":
    main()
