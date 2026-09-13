# -*- coding: utf-8 -*-
"""
make_two_figures.py —— 只画论文里这两张图，单文件、不依赖本包其他脚本、路径自适应。

    (1) fig09_p1_surface.png   图 10：近表面水分薄层与 Robin 条件核对（a 剖面 + b 柱状图）
    (2) fig13_p2_spacetime.png 图 12：问题 2 温度/水分时空云图 + 干壳薄层放大（三面板）

为什么单独做一个版本
--------------------
原绘图链是 `make_figures.py` → `fig_core.py` / `common.py` / `spectral.py` / `verify_data.py`，
其中含有对目录结构和机器环境的硬依赖，换一台电脑容易报路径错：

  * fig_core.py 写死 `C:\\Windows\\Fonts\\msyh.ttc` 这类绝对字体路径；
  * 数据目录名写死为中文的 `中间数据/`，输出目录写死为 `<脚本上级>/figure/`；
  * common.py 按 `../题目/`、`../`、`../../` 三层去找中文目录 `附件/`；
  * 控制台若为 GBK，打印中文会 UnicodeEncodeError 直接中断。

本脚本把这些依赖全部去掉：
  * 只 import numpy 与 matplotlib，不 import 本包任何模块；
  * 输入 `solutions.npz` 由命令行指定，或在脚本所在目录及上下若干层内自动搜索（含中英文候选目录名）；
  * 字体按"系统字体目录 + matplotlib 自带字体库"依次尝试，找不到 CJK 字体也不报错，只提示；
  * 输出目录可用 `--outdir` 指定，默认写到脚本同级的 `figure/`（自动创建）；
  * 控制台只打印 ASCII，任何编码的控制台都不会崩。

用法
----
    python make_two_figures.py                          # 自动找 solutions.npz，输出到 ./figure/
    python make_two_figures.py --solutions 路径/solutions.npz
    python make_two_figures.py --outdir 路径/out --dpi 240

依赖：Python 3.7+、numpy、matplotlib（无需 scipy）。
"""

import os
import sys

# 控制台可能是 GBK，打印到中文路径时会把整个脚本崩掉；这里把标准输出统一成 UTF-8。
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib import font_manager as fm

# ======================================================================
# 0. 题目给定常数（与 common.py 完全一致，写在这里以免依赖该模块）
# ======================================================================
R0 = 0.02          # m       药材半径
HM = 8.0e-7        # m/s     对流传质系数
C_AIR0 = 0.01963   # kg/kg   烘房含湿量初值（附件 1 实测初值）
C_SET = 0.050908   # kg/kg   烘房含湿量设定值（一阶惯性标定）
TAU_C = 2776.3     # s       含湿量一阶惯性时间常数
T_END = 14400.0    # s       环境分段点
D_DPI_DEFAULT = 240


def c_air(t):
    """烘房空气含湿量 [kg/kg]：t<=14400 s 走一阶惯性式，其后取设定值。"""
    t = np.asarray(t, dtype=float)
    ramp = C_SET - (C_SET - C_AIR0) * np.exp(-t / TAU_C)
    return np.where(t <= T_END, ramp, C_SET)


def d_problem1(c):
    """问题 1（附录 2）的有效扩散系数 D(C) = 7e-9 * exp(-0.89/C) [m^2/s]。"""
    return 7.0e-9 * np.exp(-0.89 / np.maximum(c, 1.0e-3))


# ======================================================================
# 1. 路径自适应
# ======================================================================
DATA_NAMES = ["solutions.npz"]
DATA_DIRS = ["中间数据", "data", "Data", "interim", "cache", "outputs", "."]


def find_solutions(explicit=None, script_dir=None):
    """定位 solutions.npz：先看命令行，再在脚本上下若干层里按候选目录名找。

    返回绝对路径；找不到抛 FileNotFoundError，并列出搜索过的位置。
    """
    script_dir = os.path.abspath(script_dir or os.path.dirname(os.path.abspath(__file__)))
    if explicit:
        p = os.path.abspath(explicit)
        if os.path.isdir(p):
            for nm in DATA_NAMES:
                q = os.path.join(p, nm)
                if os.path.exists(q):
                    return q
        if os.path.exists(p):
            return p
        raise FileNotFoundError("--solutions 指向的路径不存在: %s" % p)

    tried = []
    bases = [script_dir]
    cur = script_dir
    for _ in range(3):                      # 向上最多 3 层
        cur = os.path.dirname(cur)
        if cur and cur not in bases:
            bases.append(cur)
    for base in bases:
        for sub in DATA_DIRS:
            cand = os.path.join(base, sub)
            for nm in DATA_NAMES:
                q = os.path.abspath(os.path.join(cand, nm))
                tried.append(q)
                if os.path.exists(q):
                    return q
    # 兜底：在脚本目录及其上下两层里做一次限深遍历
    for base in bases[:2]:
        for root, dirs, files in os.walk(base):
            if root.count(os.sep) - base.count(os.sep) > 2:
                dirs[:] = []
                continue
            for nm in DATA_NAMES:
                if nm in files:
                    return os.path.abspath(os.path.join(root, nm))
    raise FileNotFoundError(
        "找不到 solutions.npz。可用 --solutions 直接指定。已搜索:\n  " + "\n  ".join(tried[:20]))


def resolve_outdir(explicit=None, script_dir=None):
    script_dir = os.path.abspath(script_dir or os.path.dirname(os.path.abspath(__file__)))
    out = os.path.abspath(explicit) if explicit else os.path.join(script_dir, "figure")
    os.makedirs(out, exist_ok=True)
    return out


# ======================================================================
# 2. 字体：不写死任何绝对路径
# ======================================================================
def setup_cjk_font(verbose=True):
    """依次尝试系统字体目录与 matplotlib 自带字体库里的 CJK 字体；失败也不中断。"""
    windir = os.environ.get("WINDIR") or os.environ.get("SystemRoot") or ""
    file_candidates = []
    if windir:
        fdir = os.path.join(windir, "Fonts")
        file_candidates += [os.path.join(fdir, n) for n in
                            ("msyh.ttc", "msyh.ttf", "simhei.ttf", "simsun.ttc", "Deng.ttf")]
    file_candidates += [
        "/System/Library/Fonts/PingFang.ttc",
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
    ]
    registered = []
    for fp in file_candidates:
        try:
            if os.path.exists(fp):
                fm.fontManager.addfont(fp)
                registered.append(os.path.basename(fp))
        except Exception:
            pass

    families = ["Microsoft YaHei", "SimHei", "SimSun", "DengXian",
                "Noto Sans CJK SC", "Source Han Sans SC", "PingFang SC",
                "WenQuanYi Zen Hei", "Arial Unicode MS", "DejaVu Sans"]
    have = {f.name for f in fm.fontManager.ttflist}
    picked = [f for f in families if f in have] or ["DejaVu Sans"]
    plt.rcParams["font.sans-serif"] = picked
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.unicode_minus"] = False     # CJK 字体多缺 U+2212，负号退回 ASCII
    if verbose:
        msg = "font families in use: %s" % ", ".join(picked)
        if picked == ["DejaVu Sans"]:
            msg += "  (WARNING: no CJK font found; Chinese labels may render as boxes)"
        print(msg)
    return picked


def apply_style():
    """与 fig_core.py 的全局样式一致，保证成图外观不变。"""
    plt.rcParams.update({
        "font.size": 10.5,
        "axes.titlesize": 11,
        "axes.labelsize": 10.5,
        "legend.fontsize": 9,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "figure.dpi": 200,
        "savefig.dpi": D_DPI_DEFAULT,
        "savefig.bbox": "tight",
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linewidth": 0.6,
    })


def save(fig, outdir, name, dpi):
    p = os.path.join(outdir, name)
    fig.savefig(p, dpi=dpi)
    plt.close(fig)
    print("  [fig] %s" % p)
    return p


# ======================================================================
# 3. CGL 谱解的插值（等价于 spectral.py 的 interp_u，但不依赖它）
# ======================================================================
def cgl_nodes(N):
    """Chebyshev--Gauss--Lobatto 节点映射到 u=(1-x)/2 ∈ [0,1]，u_0=0、u_N=1。"""
    j = np.arange(N + 1)
    return (1.0 - np.cos(j * np.pi / N)) / 2.0


def bary_weights(N):
    """重心插值权 w_j = (-1)^j δ_j，δ_0 = δ_N = 1/2。"""
    w = (-1.0) ** np.arange(N + 1)
    w[0] *= 0.5
    w[-1] *= 0.5
    return w


def interp_to(u_nodes, w, F, uu):
    """把 CGL 节点上的场 F（形状 (N+1, nt)）重心插值到任意 uu（形状 (nu,)）。

    返回 (nu, nt)。落在节点上的目标点直接取节点值。
    注意：重心公式的分母在目标点靠近端点时可以取负值，判"有效"只能看它是否为零。
    """
    F = np.atleast_2d(F)
    if F.shape[0] != u_nodes.size:
        F = F.T
    uu = np.asarray(uu, dtype=float)
    D = uu[:, None] - u_nodes[None, :]
    hit = np.abs(D) < 1e-13
    with np.errstate(divide="ignore", invalid="ignore"):
        T = w[None, :] / D
    T = np.where(hit, 0.0, T)
    den = T.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = (T @ F) / np.where(den == 0.0, np.nan, den)
    for r in np.where(hit.any(axis=1))[0]:        # 落在节点上的目标点：直接取节点值
        out[r, :] = F[np.argmax(hit[r]), :]
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


# ======================================================================
# 4. 两张图
# ======================================================================
def draw_fig09(p1_t, p1_y, outdir, dpi):
    """图 10：近表面水分薄层与 Robin 条件核对。"""
    N = p1_y.shape[0] // 2 - 1
    un = cgl_nodes(N)
    w = bary_weights(N)
    uu = np.linspace(0, 1, 220)
    rr = np.sqrt(uu) * R0 * 100.0            # cm

    Cn_all = interp_to(un, w, p1_y[:N + 1, :], uu)          # (nu, nt)
    Cn = Cn_all[:, -1]                                       # 最后一条剖面（t = 1800 s）
    t_last = float(np.ravel(p1_t)[-1])

    Ds = float(d_problem1(np.array([Cn[-1]]))[0])
    grad_num = (Cn[-1] - Cn[-2]) / ((rr[-1] - rr[-2]) / 100.0)
    grad_robin = -HM * (Cn[-1] - float(c_air(t_last))) / Ds
    rel = abs(abs(grad_num) - abs(grad_robin)) / abs(grad_robin) * 100.0

    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(rr, Cn, lw=1.8)
    ax[0].axvline(2.0, color="k", ls=":", lw=1)
    ax[0].set_xlim(1.6, 2.0)
    ax[0].set_xlabel("$r$ / cm")
    ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_title("(a) %d s 近表面水分剖面（薄层）" % int(t_last))

    # 两根柱值几乎相等（相对差 <1%），柱宽取 0.32：默认 0.8 会把绘图区填成两块实心
    # 色块，过细则两根之间空出一大片。同时把组中心间距从默认的 1.0 收到 0.5 并收紧
    # xlim，让数据占满画幅；柱顶标出数值，便于直接核对。
    v9 = [abs(grad_num), abs(grad_robin)]
    x9 = np.array([0.0, 0.5])
    ax[1].bar(x9, v9, width=0.32, color=["C0", "C1"])
    ax[1].set_xticks(x9)
    ax[1].set_xticklabels(["数值差分\n(Robin 左端)", "Robin 条件\n右端"])
    ax[1].set_xlim(-0.30, 0.80)
    for x, v in zip(x9, v9):
        ax[1].text(x, v, "%.1f" % v, ha="center", va="bottom", fontsize=10.5)
    ax[1].set_ylim(0, max(v9) * 1.14)
    ax[1].set_ylabel("$|\\partial C/\\partial r|$ / (kg/kg/m)")
    ax[1].set_title("(b) 表面梯度：差分值 vs Robin 预测\n(相对差 %.1f%%)" % rel)

    print("    grad_num=%.4f  grad_robin=%.4f  rel=%.3f%%" % (abs(grad_num), abs(grad_robin), rel))
    return save(fig, outdir, "fig09_p1_surface.png", dpi)


def draw_fig13(p2_t, p2_y, outdir, dpi):
    """图 12：问题 2 温度/水分时空云图，并把干壳薄层放大独立成 (c) 面板。"""
    N = p2_y.shape[0] // 2 - 1
    un = cgl_nodes(N)
    w = bary_weights(N)
    uu = np.linspace(0, 1, 220)
    rr = np.sqrt(uu) * R0 * 100.0
    t2h = np.ravel(p2_t) / 3600.0

    # pcolormesh 要求场的形状为 (len(Y), len(X)) = (nt, nu)，故按时间在前转置
    Cg2 = interp_to(un, w, p2_y[:N + 1, :], uu).T
    Tg2 = interp_to(un, w, p2_y[N + 1:, :], uu).T

    # 三幅并排：(c) 把原先压在 (b) 左上角的内嵌放大图拉出来独立成面板——小图不再遮挡
    # 主图数据，放大后的干壳薄层也看得清。画幅加宽到 13 in，按 \textwidth 排入时高度反而
    # 比原来的 9.6x3.4 更矮，不会挤占正文版面；字号相应放大以补偿缩印比例。
    fig, ax = plt.subplots(1, 3, figsize=(13.0, 3.5))
    fig.subplots_adjust(wspace=0.42)

    m0 = ax[0].pcolormesh(rr, t2h, Tg2, shading="auto", cmap="inferno")
    plt.colorbar(m0, ax=ax[0], label="温度 / $^\\circ$C")
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_ylabel("时间 / h")
    ax[0].set_title("(a) 问题2 温度时空云图")

    m1 = ax[1].pcolormesh(rr, t2h, Cg2, shading="auto", cmap="viridis")
    plt.colorbar(m1, ax=ax[1], label="C / (kg/kg)")
    ax[1].set_xlabel("$r$ / cm"); ax[1].set_ylabel("时间 / h")
    ax[1].set_title("(b) 问题2 水分时空云图")

    m2 = ax[2].pcolormesh(rr, t2h, Cg2, shading="auto", cmap="viridis")
    plt.colorbar(m2, ax=ax[2], label="C / (kg/kg)")
    ax[2].set_xlim(1.7, 2.0)
    ax[2].set_xlabel("$r$ / cm"); ax[2].set_ylabel("时间 / h")
    ax[2].set_title("(c) 干壳薄层放大（$r>1.7$ cm）")

    for a in ax:
        a.grid(False)

    lv_T2 = [32, 38, 44, 48]
    lv_C2 = [1.2, 1.6, 2.0, 2.3]
    for a, Z, lv in ((ax[0], Tg2, lv_T2), (ax[1], Cg2, lv_C2), (ax[2], Cg2, lv_C2)):
        cs = a.contour(rr, t2h, Z, levels=lv, colors="white", linewidths=0.9)
        for t in a.clabel(cs, fmt="%g", fontsize=11.5, inline=True, inline_spacing=2):
            t.set_color("#1a1a1a")
            t.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white")])
    for a in ax:
        a.tick_params(labelsize=11.5)
        a.xaxis.label.set_size(12); a.yaxis.label.set_size(12)
        a.title.set_size(12.5)

    return save(fig, outdir, "fig13_p2_spacetime.png", dpi)


# ======================================================================
# 5. 入口
# ======================================================================
def parse_args(argv):
    opts = {"solutions": None, "outdir": None, "dpi": D_DPI_DEFAULT}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("--solutions", "-s") and i + 1 < len(argv):
            opts["solutions"] = argv[i + 1]; i += 2
        elif a in ("--outdir", "-o") and i + 1 < len(argv):
            opts["outdir"] = argv[i + 1]; i += 2
        elif a == "--dpi" and i + 1 < len(argv):
            opts["dpi"] = int(argv[i + 1]); i += 2
        elif a in ("-h", "--help"):
            print(__doc__); sys.exit(0)
        else:
            print("ignore unknown arg: %s" % a); i += 1
    return opts


def main(argv=None):
    opts = parse_args(sys.argv[1:] if argv is None else argv)
    apply_style()
    setup_cjk_font()
    outdir = resolve_outdir(opts["outdir"])
    npz_path = find_solutions(opts["solutions"])
    print("solutions.npz : %s" % npz_path)
    print("outdir        : %s" % outdir)

    z = np.load(npz_path, allow_pickle=True)
    missing = [k for k in ("p1_t", "p1_y", "p2_t", "p2_y") if k not in z.files]
    if missing:
        raise KeyError("solutions.npz 缺少键: %s（现有: %s）" % (missing, list(z.files)))

    print("[1/2] fig09_p1_surface.png")
    draw_fig09(z["p1_t"], z["p1_y"], outdir, opts["dpi"])
    print("[2/2] fig13_p2_spacetime.png")
    draw_fig13(z["p2_t"], z["p2_y"], outdir, opts["dpi"])
    print("done. files written to: %s" % outdir)


if __name__ == "__main__":
    main()
