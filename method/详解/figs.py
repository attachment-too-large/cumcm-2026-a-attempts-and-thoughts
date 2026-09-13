# -*- coding: utf-8 -*-
"""为详解文档生成教学插图。全部图输出到 详解/fig/。

图清单：
  fig_cv     控制体与守恒定律的几何（圆柱环壳 + 节点中心有限体积网格）
  fig_map    坐标变换 r -> u=(r/R)^2 的几何意义与节点分布对照
  fig_cheb   Chebyshev-Gauss-Lobatto 节点的来源（单位圆投影）与分布
  fig_flux   通量形式的谱配置（节点通量 + 表面 Robin 通量）
  fig_p1     问题1 温度与水分剖面演化
  fig_p3     问题3 干燥曲线与终点放大
  fig_conv   空间收敛性对比（谱方法 vs 有限体积）
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, Wedge

HERE = os.path.dirname(os.path.abspath(__file__))
CODEDIR = os.path.abspath(os.path.join(HERE, "..", "药材烘干问题_正则化谱方法", "代码"))
sys.path.insert(0, CODEDIR)

FP = r"C:\Windows\Fonts\msyh.ttc"
if os.path.isfile(FP):
    fm.fontManager.addfont(FP)
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10.5
plt.rcParams["savefig.dpi"] = 220
plt.rcParams["figure.dpi"] = 220

OUT = os.path.join(HERE, "fig")
os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)
    print("saved", name)


# ======================================================================
# 图 1：控制体
# ======================================================================
def fig_cv():
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.7))

    a = ax[0]
    a.set_aspect("equal")
    a.add_patch(Rectangle((-1.25, -0.62), 2.5, 1.24, fc="#eef3fa", ec="#4a6fa5", lw=1.2))
    a.add_patch(Rectangle((-1.25, -0.62), 0.30, 1.24, fc="#c9d8ee", ec="#4a6fa5", lw=1.0))
    a.add_patch(Rectangle((0.95, -0.62), 0.30, 1.24, fc="#c9d8ee", ec="#4a6fa5", lw=1.0))
    a.plot([-1.25, 1.25], [0, 0], color="#888", ls="-.", lw=1.0)
    a.annotate("", xy=(-1.25, 0.86), xytext=(1.25, 0.86),
               arrowprops=dict(arrowstyle="<->", color="#333", lw=1.1))
    a.text(0.0, 0.92, "长度 L = 0.25 m", ha="center", fontsize=9.5)
    a.annotate("", xy=(-1.25, -0.80), xytext=(-0.95, -0.80),
               arrowprops=dict(arrowstyle="<->", color="#333", lw=1.1))
    a.text(-1.10, -0.98, "dr", ha="center", fontsize=10)
    a.annotate("半径 R = 2 cm", xy=(1.10, 0.30), xytext=(1.35, 0.55),
               fontsize=9.5, arrowprops=dict(arrowstyle="->", color="#333", lw=1.0))
    a.text(0.0, 0.18, "环壳控制体：\n半径 r，厚 dr，侧面积 2*pi*r*L", ha="center",
           fontsize=9.5, color="#204060")
    a.text(-1.10, 0.60, "剖面", ha="center", fontsize=9, color="#204060")
    a.set_xlim(-1.6, 2.3)
    a.set_ylim(-1.25, 1.15)
    a.axis("off")
    a.set_title("(a) 圆柱与环壳控制体", fontsize=10)

    b = ax[1]
    b.set_aspect("equal")
    N = 8
    h = 1.0 / N
    for i in range(N + 1):
        x = i * h
        if i == 0:
            w = h * h / 8.0
            b.add_patch(Rectangle((0.0, -w / 2), h / 2, w, fc="#f6d6d6", ec="#a03a3a", lw=1.0))
        elif i == N:
            w = h / 2 - h * h / 8
            b.add_patch(Rectangle((1 - h / 2, -w / 2), h / 2, w, fc="#d6e8f6", ec="#3a6ea0", lw=1.0))
        else:
            w = x * h
            b.add_patch(Rectangle((x - h / 2, -w / 2), h, w, fc="#e8f0e8", ec="#4a7a4a", lw=1.0))
        b.plot([x], [0], "ko", ms=2.6)
    b.annotate("第一个控制体 [0, h/2]\n测度 = h^2/8\n比中间格子薄十六倍", xy=(h / 4, 0.012),
               xytext=(0.20, 0.30),
               fontsize=9, arrowprops=dict(arrowstyle="->", color="#a03a3a", lw=1.0),
               color="#a03a3a")
    b.annotate("表面控制体 [R-h/2, R]\n测度 = h/2 - h^2/8", xy=(1 - h / 4, 0.03),
               xytext=(0.40, -0.36),
               fontsize=9, arrowprops=dict(arrowstyle="->", color="#3a6ea0", lw=1.0),
               color="#3a6ea0")
    b.plot([0, 1], [0, 0], color="#555", lw=1.0)
    b.text(0.5, -0.50, "盒子高度就是该控制体的体积测度，按真实比例绘制", ha="center",
           fontsize=9.5)
    b.set_xlim(-0.05, 1.12)
    b.set_ylim(-0.62, 0.45)
    b.axis("off")
    b.set_title("(b) 常规 r 坐标的节点中心有限体积", fontsize=10)
    save(fig, "fig_cv.png")


# ======================================================================
# 图 2：坐标变换 r -> u
# ======================================================================
def fig_map():
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.6))
    R = 1.0
    r = np.linspace(0, R, 400)
    a = ax[0]
    a.plot(r / R, (r / R) ** 2, lw=2.0, color="#1f4e79", label="u = (r/R)^2")
    a.plot(r / R, r / R, lw=1.2, ls="--", color="#999", label="u = r/R（对照：恒等映射）")
    for rr in [0.0, 0.5, 0.7071, 1.0]:
        a.plot([rr], [rr ** 2], "o", color="#c00000", ms=5)
        a.annotate("", xy=(rr, rr ** 2), xytext=(rr, 0),
                   arrowprops=dict(arrowstyle="->", color="#c00000", lw=0.9, ls=":"))
    a.text(0.52, 0.16, "等距的 r\n映射成不等距的 u", fontsize=9, color="#c00000")
    a.set_xlabel("归一化半径 xi = r/R")
    a.set_ylabel("正则化坐标 u")
    a.set_title("(a) 坐标变换把拉伸与压缩分开", fontsize=10)
    a.grid(alpha=.3)
    a.legend(fontsize=8.5, loc="upper left")

    b = ax[1]
    N = 12
    j = np.arange(N + 1)
    x = np.cos(np.pi * j / N)
    u = (1 - x) / 2
    rnodes = np.sqrt(u)
    b.plot(u, np.zeros_like(u), "|", color="#1f4e79", ms=14, mew=1.6)
    b.plot(rnodes, np.ones_like(u) * 0.55, "|", color="#c00000", ms=14, mew=1.6)
    b.text(0.5, 0.16, "u 节点（u 上向两端聚集）", ha="center", fontsize=9, color="#1f4e79")
    b.text(0.62, 0.71, "对应的物理半径 r = R*sqrt(u)", ha="center", fontsize=9, color="#c00000")
    b.plot([0, 1], [1.12, 1.12], color="#555", lw=1.0)
    b.text(0.5, 1.20, "表面 u = 1（r = R）", ha="center", fontsize=9)
    b.plot([0, 0], [1.12, 0.0], color="#aaa", ls=":", lw=0.9)
    b.text(-0.04, -0.10, "轴心 u = 0", ha="left", fontsize=9)
    b.set_xlim(-0.12, 1.12)
    b.set_ylim(-0.25, 1.45)
    b.axis("off")
    b.set_title("(b) 同一批节点在两种坐标下的位置", fontsize=10)
    save(fig, "fig_map.png")


# ======================================================================
# 图 3：Chebyshev 节点的来源
# ======================================================================
def fig_cheb():
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.7))
    N = 10
    j = np.arange(N + 1)
    th = j * np.pi / N
    x = np.cos(th)
    a = ax[0]
    a.set_aspect("equal")
    t = np.linspace(0, 2 * np.pi, 400)
    a.plot(np.cos(t), np.sin(t), color="#8aa8c8", lw=1.2)
    for jj in range(N + 1):
        a.plot([0, np.cos(th[jj])], [0, np.sin(th[jj])], color="#cccccc", lw=0.7)
        a.plot(np.cos(th[jj]), np.sin(th[jj]), "o", color="#1f4e79", ms=4)
        a.plot([np.cos(th[jj]), np.cos(th[jj])], [0, np.sin(th[jj])], ls=":", color="#c00000", lw=0.8)
    a.axhline(0, color="#555", lw=1.0)
    a.text(0.0, -0.22, "把圆上等角距的点\n垂直投影到直径上", ha="center", fontsize=9)
    a.set_xlim(-1.35, 1.35)
    a.set_ylim(-1.30, 1.30)
    a.axis("off")
    a.set_title("(a) 节点来源：单位圆等角投影", fontsize=10)

    b = ax[1]
    b.plot(x, np.zeros_like(x), "|", color="#1f4e79", ms=16, mew=1.8)
    for jj in range(N + 1):
        b.plot([x[jj], x[jj]], [0, 0.4], ls=":", color="#c00000", lw=0.8)
    b.text(0.0, 0.52, "两端密、中间疏（每端约 1/N^2 量级）", ha="center", fontsize=9.5)
    b.text(-1.15, -0.16, "x = -1\n(轴心 u=0)", ha="center", fontsize=9)
    b.text(1.15, -0.16, "x = +1\n(表面 u=1)", ha="center", fontsize=9)
    b.set_xlim(-1.45, 1.45)
    b.set_ylim(-0.35, 0.75)
    b.axis("off")
    b.set_title("(b) 节点在 [-1,1] 上的分布", fontsize=10)
    save(fig, "fig_cheb.png")


# ======================================================================
# 图 4：通量形式的谱配置
# ======================================================================
def fig_flux():
    fig, ax = plt.subplots(figsize=(9.6, 2.9))
    N = 8
    j = np.arange(N + 1)
    x = np.cos(np.pi * j / N)
    u = (1 - x) / 2
    ax.plot([0, 1], [0, 0], color="#555", lw=1.4)
    for k, uu in enumerate(u):
        if k == 0:
            ax.plot([uu], [0], "o", color="#c00000", ms=8)
            ax.text(uu, -0.22, "节点 0\nu=0", ha="center", fontsize=8.5, color="#c00000")
        elif k == N:
            ax.plot([uu], [0], "o", color="#1f4e79", ms=8)
            ax.text(uu, -0.22, "节点 N\nu=1", ha="center", fontsize=8.5, color="#1f4e79")
        else:
            ax.plot([uu], [0], "o", color="#4a7a4a", ms=6)
            ax.text(uu, -0.22, "u_%d" % k, ha="center", fontsize=8.5, color="#4a7a4a")
        if k < N:
            ax.annotate("", xy=(u[k] + (u[k + 1] - u[k]) * 0.25, 0.30),
                        xytext=(u[k], 0.30),
                        arrowprops=dict(arrowstyle="->", color="#4a7a4a", lw=1.0))
    ax.text(0.5, 0.44, "内部节点通量 g_j = u_j * D_j * (dC/du)_j 由谱微分矩阵算出",
            ha="center", fontsize=9.5, color="#4a7a4a")
    ax.annotate("", xy=(0.985, 0.05), xytext=(0.66, 0.85),
                arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.2))
    ax.text(0.60, 0.92, "表面通量不由微分矩阵给出，\n直接取 Robin 条件的物理值",
            ha="left", fontsize=9.5, color="#1f4e79")
    ax.text(0.10, -0.55, "轴心通量 g_0 = u_0 * D_0 * (dC/du)_0 = 0，因为 u_0 = 0：对称条件自动满足",
            ha="left", fontsize=9.5, color="#c00000")
    ax.set_xlim(-0.08, 1.20)
    ax.set_ylim(-0.75, 1.15)
    ax.axis("off")
    save(fig, "fig_flux.png")


# ======================================================================
# 图 5/6：物理剖面（用主求解器现算）
# ======================================================================
def fig_physics():
    import common as cm
    from spectral import HerbSolver
    RCOLS = np.array([0.0, 0.005, 0.010, 0.015, 0.020])
    CRIT = 0.15

    s1 = HerbSolver(64, prob=1)
    t1 = np.arange(0.0, 1800.1, 60.0)
    at1 = np.concatenate([np.full(s1.Np, 1e-13), np.full(s1.Np, 1e-11)])
    sol1 = s1.solve(1800.0, t_eval=list(t1), rtol=1e-11, atol=at1)
    U1 = (RCOLS / cm.R0) ** 2
    T1 = np.array([s1.interp_u(sol1.y[:, i], U1)[1] for i in range(len(t1))])
    C1 = np.array([s1.interp_u(sol1.y[:, i], U1)[0] for i in range(len(t1))])

    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.5))
    rr = RCOLS * 100
    for tt in [100, 600, 1200, 1800]:
        i = int(np.argmin(np.abs(t1 - tt)))
        ax[0].plot(rr, T1[i], "o-", ms=3.5, lw=1.4, label="t=%d s" % tt)
        ax[1].plot(rr, C1[i], "s-", ms=3.5, lw=1.4, label="t=%d s" % tt)
    ax[0].set_xlabel("到中心的距离 r / cm")
    ax[0].set_ylabel("温度 T / °C")
    ax[0].set_title("(a) 问题1 温度剖面：热扩散快", fontsize=10)
    ax[1].set_xlabel("到中心的距离 r / cm")
    ax[1].set_ylabel("水分浓度 C / (kg/kg)")
    ax[1].set_title("(b) 问题1 水分剖面：只有表层在失水", fontsize=10)
    for a in ax:
        a.grid(alpha=.3)
        a.legend(fontsize=8.5)
    save(fig, "fig_p1.png")

    s3 = HerbSolver(64, prob=3)
    at3 = np.concatenate([np.full(s3.Np, 1e-13), np.full(s3.Np, 1e-11)])

    def ev3(t, y):
        return np.max(y[:s3.Np]) - CRIT

    ev3.terminal = True
    ev3.direction = -1
    solE = s3.solve(300000.0, rtol=1e-10, atol=at3, events=ev3)
    tdry3 = float(solE.t_events[0][0])

    tv = np.linspace(0.0, tdry3 * 1.06, 1500)
    sol3 = s3.solve(tdry3 * 1.06, t_eval=list(tv), rtol=1e-10, atol=at3)
    C3 = np.array([s3.interp_u(sol3.y[:, i], U1)[0] for i in range(len(tv))])
    hh = tv / 3600.0

    # 截面平均要按面积加权：每一圈的面积与半径成正比，权重是半径本身。
    # 定义是 Cbar = (2/R^2) * 积分 C(r) r dr。
    # 注意不能把切比雪夫节点上的值直接加权平均，节点在表面附近极密，
    # 那样会把平均值拉低约 0.2，达标时刻提前十个多小时。
    # 这里先把解插值到均匀半径网格，再按面积加权做梯形积分。
    rr = np.linspace(0.0, cm.R0, 401)
    wb = (-1.0) ** np.arange(s3.Np)
    wb[0] *= 0.5
    wb[-1] *= 0.5
    dif = (rr / cm.R0)[:, None] ** 2 - s3.u[None, :]
    den = (wb[None, :] / np.where(np.abs(dif) < 1e-14, 1.0, dif)).sum(axis=1, keepdims=True)
    Bm = (wb[None, :] / np.where(np.abs(dif) < 1e-14, 1.0, dif)) / den
    Cf = Bm @ sol3.y[:s3.Np, :]
    meanC = np.trapezoid(Cf * rr[:, None], rr, axis=0) * 2.0 / cm.R0 ** 2

    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.5))
    for jj, r in enumerate(RCOLS):
        ax[0].semilogy(hh, C3[:, jj], lw=1.4,
                       label=("r=%.1f cm" % (r * 100)) if jj in (0, 2, 4) else None)
    ax[0].axhline(CRIT, color="r", ls="--", lw=1.1)
    ax[0].axvline(tdry3 / 3600.0, color="k", ls=":", lw=1.1)
    ax[0].text(2, 0.165, "判据 0.15 kg/kg", color="r", fontsize=9)
    ax[0].text(tdry3 / 3600.0 - 26, 1.4, "tdry = %.2f h" % (tdry3 / 3600.0), fontsize=9)
    ax[0].set_xlabel("时间 t / h")
    ax[0].set_ylabel("水分浓度 C / (kg/kg)")
    ax[0].set_title("(a) 问题3 各半径处干燥曲线（半对数）", fontsize=10)
    ax[0].set_xlim(0, 60)
    ax[0].set_ylim(0.04, 3.2)
    ax[0].grid(alpha=.3, which="both")
    ax[0].legend(fontsize=8.5, loc="lower left")

    zoom = (hh > 48) & (hh < 60)
    ax[1].plot(hh[zoom], C3[zoom, 0], lw=1.8, label="中心 C(0,t)")
    ax[1].plot(hh[zoom], C3[zoom, -1], lw=1.8, label="表面 C(R,t)")
    ax[1].plot(hh[zoom], meanC[zoom], lw=1.5, ls="--", label="截面面积加权平均")
    ax[1].axhline(CRIT, color="r", ls="--", lw=1.2)
    ax[1].axvline(tdry3 / 3600.0, color="k", ls=":", lw=1.2)
    ax[1].plot([tdry3 / 3600.0], [CRIT], "r*", ms=12)
    ax[1].set_xlabel("时间 t / h")
    ax[1].set_ylabel("水分浓度 C / (kg/kg)")
    ax[1].set_title("(b) 终点附近放大：平均值先达标，中心最后达标", fontsize=10)
    ax[1].grid(alpha=.3)
    ax[1].legend(fontsize=8.5)
    save(fig, "fig_p3.png")
    print("   tdry = %.3f s = %.4f h" % (tdry3, tdry3 / 3600.0))
    below = np.where(meanC < CRIT)[0]
    if below.size:
        print("   截面平均首次低于 0.15 的时刻 = %.3f h" % hh[below[0]])


# ======================================================================
# 图 7：空间收敛性
# ======================================================================
def fig_conv():
    import common as cm
    from spectral import HerbSolver
    import analytic_p1 as an

    Ns = [4, 6, 8, 12, 16, 24, 32, 64]
    err = []
    for N in Ns:
        s = HerbSolver(N, prob=1)
        at = np.concatenate([np.full(s.Np, 1e-14), np.full(s.Np, 1e-12)])
        sol = s.solve(1800.0, t_eval=[1800.0], rtol=1e-12, atol=at)
        rc = np.linspace(0, cm.R0, 201)
        Ti = s._bary(sol.y[s.Np:, -1], (rc / cm.R0) ** 2)
        Te = an.T_exact(rc, 1800.0)
        err.append(float(np.max(np.abs(Ti - Te))))
        print("   N=%3d err=%.3e" % (N, err[-1]))

    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.5))
    ax[0].semilogy(Ns, err, "^-", color="#1f4e79", lw=1.6, ms=6, label="Chebyshev 谱配置法")
    ax[0].axhline(1e-12, color="#c00000", ls="--", lw=1.1, label="时间积分容差平台")
    ax[0].set_xlabel("空间节点参数 N")
    ax[0].set_ylabel("全场最大温度误差 / K")
    ax[0].set_title("(a) 谱收敛：N=8 就落到 1e-10", fontsize=10)
    ax[0].grid(alpha=.3, which="both")
    ax[0].legend(fontsize=8.5)

    Nfv = np.array([800, 1600, 3200, 6400])
    t3 = np.array([205550.98, 205567.49, 205572.43, 205573.77])
    ax[1].semilogx(Nfv, np.abs(t3 - 205574.2), "s-", color="#a03a3a", lw=1.6, ms=6,
                   label="有限体积法（二阶）")
    ax[1].semilogx(Nfv, Nfv.astype(float) ** -2 * 1.2e7, "k--", lw=1.0, label="斜率 -2 参考线")
    ax[1].set_xlabel("空间网格数 N")
    ax[1].set_ylabel("烘干时长误差 / s")
    ax[1].set_title("(b) 对照格式只达到二阶", fontsize=10)
    ax[1].grid(alpha=.3, which="both")
    ax[1].legend(fontsize=8.5)
    save(fig, "fig_conv.png")


if __name__ == "__main__":
    fig_cv()
    fig_map()
    fig_cheb()
    fig_flux()
    fig_physics()
    fig_conv()
    print("ALL FIGURES DONE ->", OUT)
