# -*- coding: utf-8 -*-
"""烘房环境两种口径的对比实验（附件1 数据的两种用法）。

口径甲：把附件1 标定成一阶惯性式，用光滑曲线作为烘房环境。
口径乙：0 到 14400 秒直接用附件1 的表格数据线性插值，
        14400 秒以后两种口径都取拟合设定值。

两种口径只在预热平衡段不同，因此这个实验隔离出的正是
用光滑拟合曲线代替原始表格数据所带来的影响。

用法：python env_compare.py
输出：环境口径对比.txt 与本文件的打印
"""
import io
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm
from spectral import HerbSolver

OUT = []


def w(s=""):
    OUT.append(str(s))
    print(s, flush=True)


_orig_T, _orig_C = cm.T_air, cm.C_air
_t, _Ta, _Ca = cm._t_env, cm._Ta_env, cm._Ca_env


def T_tab(t):
    """附件1 表格线性插值口径。"""
    t = np.asarray(t, dtype=float)
    v = np.interp(t, _t, _Ta)
    return np.where(t <= cm.T_END, v, cm.Tset)


def C_tab(t):
    """附件1 表格线性插值口径。"""
    t = np.asarray(t, dtype=float)
    v = np.interp(t, _t, _Ca)
    return np.where(t <= cm.T_END, v, cm.Cset)


def run(prob, need_tdry=False):
    """求解一次，返回 (T(0,1800), T(R,1800), C(R,1800), tdry)。

    问题一的物性下水分极难扩散，烘干时间远超合理积分区间，
    因此只在需要时才对问题三求烘干时长。
    """
    s = HerbSolver(64, prob=prob)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])
    sol = s.solve(1800.0, t_eval=[1800.0], rtol=1e-11, atol=at)
    C, T = s.interp_u(sol.y[:, -1], np.array([0.0, 1.0]))
    t0, tR, cR = float(T[0]), float(T[1]), float(C[1])
    if not need_tdry:
        return t0, tR, cR, float("nan")

    def ev(t, y):
        return np.max(y[:s.Np]) - 0.15
    ev.terminal = True
    ev.direction = -1
    e = s.solve(300000.0, rtol=1e-11, atol=at, events=ev)
    return t0, tR, cR, float(e.t_events[0][0])


def main():
    w("烘房环境两种口径的对比")
    w("=" * 74)
    w("口径甲：T_air 与 C_air 取一阶惯性拟合式")
    w("口径乙：0 到 14400 秒取附件1 表格线性插值，其后取拟合设定值")
    w("")
    i = int(np.argmin(np.abs(_t - 1800)))
    w("空气在 t = 1800 秒处的取值：")
    w("  附件1 表格       T_air = %.3f degC    C_air = %.5f kg/kg" % (_Ta[i], _Ca[i]))
    w("  拟合曲线         T_air = %.4f degC    C_air = %.6f kg/kg"
      % (float(_orig_T(1800.0)), float(_orig_C(1800.0))))
    w("  差值             %+.4f degC          %+.6f kg/kg"
      % (float(_orig_T(1800.0)) - _Ta[i], float(_orig_C(1800.0)) - _Ca[i]))
    w("")

    cm.T_air, cm.C_air = _orig_T, _orig_C
    a = run(1)
    b = run(3, need_tdry=True)
    cm.T_air, cm.C_air = T_tab, C_tab
    c = run(1)
    d = run(3, need_tdry=True)
    cm.T_air, cm.C_air = _orig_T, _orig_C

    w("=" * 74)
    w("%-26s %14s %14s %12s" % ("", "口径甲", "口径乙", "差值"))
    w("-" * 74)
    w("%-26s %14.4f %14.4f %+12.4f" % ("问题1 T(0, 1800 s) / degC", a[0], c[0], c[0] - a[0]))
    w("%-26s %14.4f %14.4f %+12.4f" % ("问题1 T(R, 1800 s) / degC", a[1], c[1], c[1] - a[1]))
    w("%-26s %14.4f %14.4f %+12.4f" % ("问题1 C(R, 1800 s) / (kg/kg)", a[2], c[2], c[2] - a[2]))
    w("%-26s %14.1f %14.1f %+12.1f" % ("问题3 tdry / s", b[3], d[3], d[3] - b[3]))
    w("-" * 74)
    w("问题1 中心温度的相对差别：%.4f%%" % ((c[0] - a[0]) / a[0] * 100.0))
    w("问题3 烘干时长的相对差别：%.4f%%" % ((d[3] - b[3]) / b[3] * 100.0))

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "环境口径对比.txt")
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT))
    print("written", path)


if __name__ == "__main__":
    main()
