# -*- coding: utf-8 -*-
"""论文用图统一绘图设置与数据缓存。

约定：所有图保存到 <交付根>/图表/，分辨率 200 dpi 以上，中文标签用 Microsoft YaHei。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
FIGDIR = os.path.join(ROOT, "图表")
DATADIR = os.path.join(ROOT, "中间数据")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# ---------- 中文字体 ----------
for _fp in [r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf"]:
    if os.path.exists(_fp):
        fm.fontManager.addfont(_fp)
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10.5
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["axes.labelsize"] = 10.5
plt.rcParams["legend.fontsize"] = 9
plt.rcParams["xtick.labelsize"] = 9.5
plt.rcParams["ytick.labelsize"] = 9.5
plt.rcParams["figure.dpi"] = 200
plt.rcParams["savefig.dpi"] = 200
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.3
plt.rcParams["grid.linewidth"] = 0.6

DPI = 200


def save(fig, name):
    p = os.path.join(FIGDIR, name)
    fig.savefig(p, dpi=DPI)
    plt.close(fig)
    print("  [fig] %s" % name)
    return p


# ======================================================================
#  数据缓存
# ======================================================================
def build_cache():
    """计算并缓存四问的解，避免重复求解。"""
    import common as cm
    from spectral import HerbSolver

    cache = os.path.join(DATADIR, "solutions.npz")
    if os.path.exists(cache):
        print("[cache] 已存在，跳过重算：", cache)
        return
    print("[cache] 开始计算四问解 ...")
    N = 64
    d = {}

    # ---- 问题1 ----
    s1 = HerbSolver(N, prob=1)
    t1 = np.arange(0.0, 1800.1, 10.0)
    at1 = np.concatenate([np.full(s1.Np, 1e-13), np.full(s1.Np, 1e-11)])
    sol1 = s1.solve(1800.0, t_eval=t1, rtol=1e-11, atol=at1)
    d["p1_t"], d["p1_y"] = sol1.t, sol1.y

    # ---- 问题2 ----
    s2 = HerbSolver(N, prob=2)
    t2 = np.arange(0.0, 10800.1, 60.0)
    at2 = np.concatenate([np.full(s2.Np, 1e-12), np.full(s2.Np, 1e-10)])
    sol2 = s2.solve(10800.0, t_eval=t2, rtol=1e-11, atol=at2)
    d["p2_t"], d["p2_y"] = sol2.t, sol2.y

    # ---- 问题3 ----
    CRIT = 0.15
    s3 = HerbSolver(N, prob=3)
    at3 = np.concatenate([np.full(s3.Np, 1e-13), np.full(s3.Np, 1e-11)])

    def ev(t, y):
        return np.max(y[:s3.Np]) - CRIT
    ev.terminal = True
    ev.direction = -1
    e3 = s3.solve(300000.0, rtol=1e-11, atol=at3, events=ev)
    td3 = float(e3.t_events[0][0])
    t3 = np.linspace(0.0, td3, 1400)
    sol3 = s3.solve(td3, t_eval=list(t3), rtol=1e-11, atol=at3)
    d["p3_t"], d["p3_y"], d["p3_tdry"] = sol3.t, sol3.y, np.array([td3])

    # ---- 问题4 ----
    s4 = HerbSolver(N, prob=4, shrink=True)
    at4 = np.concatenate([np.full(s4.Np, 1e-13), np.full(s4.Np, 1e-11)])

    def ev4(t, y):
        return np.max(y[:s4.Np]) - CRIT
    ev4.terminal = True
    ev4.direction = -1
    e4 = s4.solve(300000.0, rtol=1e-11, atol=at4, events=ev4)
    td4 = float(e4.t_events[0][0])
    t4 = np.linspace(0.0, td4, 1400)
    sol4 = s4.solve(td4, t_eval=list(t4), rtol=1e-11, atol=at4)
    d["p4_t"], d["p4_y"], d["p4_tdry"] = sol4.t, sol4.y, np.array([td4])

    np.savez_compressed(cache, **d)
    print("[cache] 完成：", cache)


def load_cache():
    import common as cm
    from spectral import HerbSolver
    z = np.load(os.path.join(DATADIR, "solutions.npz"))
    return z, HerbSolver, cm
