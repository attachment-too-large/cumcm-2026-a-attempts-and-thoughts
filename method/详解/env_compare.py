# -*- coding: utf-8 -*-
"""烘房环境的两种口径对比：一阶惯性拟合曲线 与 附件1表格直接插值。

两种口径只在预热平衡段（0 到 14400 秒）不同，
14400 秒以后两者都取设定值，因此这一实验隔离出的正是
用光滑拟合曲线代替原始表格数据所带来的影响。

输出：
  问题一 1800 秒的温度与表面水分浓度
  问题三 烘干时长
"""
import io
import sys
import numpy as np

CODEDIR = r"C:\Users\qing1\Desktop\bestway\药材烘干问题_正则化谱方法\代码"
sys.path.insert(0, CODEDIR)
import common as cm
from spectral import HerbSolver

OUT = []


def w(s=""):
    OUT.append(str(s))
    print(s, flush=True)


_orig_T, _orig_C = cm.T_air, cm.C_air
_t, _Ta, _Ca = cm._t_env, cm._Ta_env, cm._Ca_env


def T_tab(t):
    t = np.asarray(t, dtype=float)
    v = np.interp(t, _t, _Ta)
    return np.where(t <= cm.T_END, v, cm.Tset)


def C_tab(t):
    t = np.asarray(t, dtype=float)
    v = np.interp(t, _t, _Ca)
    return np.where(t <= cm.T_END, v, cm.Cset)


def run(prob, need_tdry=False):
    """返回 (T(0,1800), T(R,1800), C(R,1800), tdry)。

    只有 need_tdry 为真时才求烘干时长。问题一的物性下水分极难扩散，
    烘干时间远超任何合理的积分区间，因此不对问题一求烘干时长。
    """
    s = HerbSolver(64, prob=prob)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])
    sol = s.solve(1800.0, t_eval=[1800.0], rtol=1e-11, atol=at)
    C, T = s.interp_u(sol.y[:, -1], np.array([0.0, 1.0]))
    t0, tR = float(T[0]), float(T[1])
    cR = float(C[1])
    if not need_tdry:
        return t0, tR, cR, float("nan")

    def ev(t, y):
        return np.max(y[:s.Np]) - 0.15
    ev.terminal = True
    ev.direction = -1
    e = s.solve(300000.0, rtol=1e-11, atol=at, events=ev)
    return t0, tR, cR, float(e.t_events[0][0])


w("烘房环境两种口径的对比")
w("=" * 74)
w("拟合曲线口径：T_air(t) 用一阶惯性式，C_air(t) 同理")
w("表格插值口径：0 到 14400 秒直接用附件1 数据线性插值，")
w("              14400 秒以后仍取拟合设定值，两种口径在恒温段完全一致")
w("")
w("空气在 t = 1800 秒处的取值：")
i = int(np.argmin(np.abs(_t - 1800)))
w("  附件1 表格实测   T_air = %.3f degC    C_air = %.5f kg/kg" % (_Ta[i], _Ca[i]))
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
w("%-26s %14s %14s %12s" % ("", "拟合曲线", "表格插值", "差值"))
w("-" * 74)
w("%-26s %14.4f %14.4f %+12.4f" % ("问题1 T(0, 1800 s) / degC", a[0], c[0], c[0] - a[0]))
w("%-26s %14.4f %14.4f %+12.4f" % ("问题1 T(R, 1800 s) / degC", a[1], c[1], c[1] - a[1]))
w("%-26s %14.4f %14.4f %+12.4f" % ("问题1 C(R, 1800 s) / (kg/kg)", a[2], c[2], c[2] - a[2]))
w("%-26s %14.1f %14.1f %+12.1f" % ("问题3 tdry / s", b[3], d[3], d[3] - b[3]))
w("-" * 74)
w("问题3 烘干时长的相对差别：%.4f%%" % ((d[3] - b[3]) / b[3] * 100.0))
w("问题1 中心温度的相对差别：%.4f%%" % ((c[0] - a[0]) / a[0] * 100.0))

with io.open(r"C:\Users\qing1\Desktop\bestway\详解\环境口径对比.txt", "w",
             encoding="utf-8") as fh:
    fh.write("\n".join(OUT))
print("written 环境口径对比.txt")
