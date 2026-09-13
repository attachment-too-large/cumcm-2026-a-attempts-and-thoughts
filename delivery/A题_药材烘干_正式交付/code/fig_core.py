# -*- coding: utf-8 -*-
# ==============================================================================
# fig_core.py —— 论文插图的统一绘图设置、出图入口与四问解缓存
# ==============================================================================
# 做什么：集中设定 matplotlib 的中文字体、字号、网格与出图分辨率，并对外只留一个
#         出图入口 save()，使全部插图的大小字号风格一致；同时缓存四问的解，让
#         make_figures.py 等脚本不必每次重新求解。
# 输入：common.py、spectral.py 与本文件的源码（用于算缓存指纹）；
#       中间数据/solutions.npz（存在且版本/指纹相符时直接复用，否则重算）。
# 输出：figure/*.png                由 save() 逐张写出，240 dpi
#       中间数据/solutions.npz      由 build_cache() 写出四问解 + 版本号 + 源码指纹
# 关键变量：
#   FIGDIR / DATADIR  出图目录 figure/ 与中间数据/（相对交付根，不存在则自动创建）
#   DPI    出图分辨率 240 dpi：插图排版时会缩到 6.30 in 版心内，缩放系数约 0.76，
#          按 200 dpi 渲染印后只剩约 260 ppi，240 dpi 才够 310 ppi 以上
#   CACHE_VERSION  缓存版本号字符串；改动物理模型、参数、离散或容差必须同步递增
#   N      四问解统一使用的谱节点数 64（Np = 65 个 CGL 节点，两场共 130 个状态量）
#   CRIT   烘干判据阈值 0.15 kg/kg（干基含水率），仅用于问题 3、4 的事件求交
# ==============================================================================
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
FIGDIR = os.path.join(ROOT, "figure")
DATADIR = os.path.join(ROOT, "中间数据")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# ---------------- 中文字体与字号 ----------------
# 不把字体文件注册进 fontManager，中文标签就渲染成一片方框：matplotlib 自带的
# 字体表里没有 CJK 字体。axes.unicode_minus=False 同理——CJK 字体缺 U+2212 字形，
# 负号不退回 ASCII 连字符也会显示成方框。
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
plt.rcParams["savefig.dpi"] = 240
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.3
plt.rcParams["grid.linewidth"] = 0.6

# 出图分辨率：插图按 \textwidth 排版后会被缩到 6.30 in 版心内，缩放系数约 0.76，
# 按 200 dpi 渲染到纸面上的有效像素密度只剩约 260 ppi，240 dpi 才有 310 ppi 以上。
# 注意 dpi 只决定像素密度，不改变图内字号的磅值（磅值由 figsize 与版心宽度之比
# 决定），故提高 dpi 不会动到任何一张图的版面与版式。
DPI = 240


def save(fig, name):
    """把画布 fig 存成 figure/<name>，关闭画布并打印一行日志；返回写出的绝对路径。

    所有插图都收敛到这一个入口，是为了把分辨率和紧裁剪钉死在一处——逐图手写
    savefig 很容易漏掉 dpi。注意本函数会 close(fig)，调用方不得再复用该画布。
    """
    p = os.path.join(FIGDIR, name)
    fig.savefig(p, dpi=DPI)
    plt.close(fig)
    print("  [fig] %s" % name)
    return p


# ---------------- 四问解缓存 ----------------
# 缓存版本号：任何改变解的内容的改动（物理模型、参数、离散、容差）都必须同时
# 递增此号——否则旧的 solutions.npz 会被静默沿用，图上画的是旧解。
CACHE_VERSION = "2026-09-13.N64.rtol1e-11"

import hashlib


def _code_fingerprint():
    """对求解相关源文件（common.py、spectral.py、fig_core.py）整体取 SHA-256，
    返回前 16 位十六进制字符串。

    把它随解一起落盘，下次启动时对比即可发现"改了代码但缓存没重算"；只靠版本号
    字符串要靠人记得手工递增，漏改一次就会画出旧解。
    """
    h = hashlib.sha256()
    for name in ("common.py", "spectral.py", "fig_core.py"):
        p = os.path.join(HERE, name)
        if os.path.exists(p):
            with open(p, "rb") as fh:
                h.update(fh.read())
    return h.hexdigest()[:16]


CACHE_FINGERPRINT = _code_fingerprint()


def build_cache(force=False):
    """求解问题 1~4 并写入 中间数据/solutions.npz，供出图脚本复用，避免重复求解。

    缓存中的键（时间一律 [s]；y 的前 Np 行是水分浓度 C [kg/kg]、后 Np 行是温度
    T [degC]，按 CGL 节点 u = 0..1 排列）：
      p1_t, p1_y              问题 1：0~1800 s，报告点间隔 10 s
      p2_t, p2_y              问题 2：0~10800 s，报告点间隔 60 s
      p3_t, p3_y, p3_tdry     问题 3：0~t_dry，均匀 1400 点；t_dry [s]
      p4_t, p4_y, p4_tdry     问题 4：含收缩，同上
      version, fingerprint    版本号与源码指纹，用于判定缓存是否对应当前代码

    force=True（或环境变量 A_FORCE_CACHE=1）时无条件重算；版本号或源码指纹与当前
    代码不符时同样重算，故缓存不会悄悄沿用旧解。
    """
    import common as cm
    from spectral import HerbSolver

    cache = os.path.join(DATADIR, "solutions.npz")
    if os.environ.get("A_FORCE_CACHE") == "1":
        force = True
    if os.path.exists(cache) and not force:
        try:
            z = np.load(cache)
            ver = str(z["version"]) if "version" in z else "<无版本号>"
            fp = str(z["fingerprint"]) if "fingerprint" in z else "<无指纹>"
        except Exception as exc:                      # 缓存损坏：重算而不是崩溃
            print("[cache] 缓存无法读取（%s），改为重算" % exc)
        else:
            if ver == CACHE_VERSION and fp == CACHE_FINGERPRINT:
                print("[cache] 命中：%s（版本 %s）" % (cache, ver))
                return
            print("[cache] 已存在的缓存与当前代码不符，重算：")
            print("        缓存版本 %s / 当前版本 %s" % (ver, CACHE_VERSION))
            print("        缓存指纹 %s / 当前指纹 %s" % (fp, CACHE_FINGERPRINT))
            print("        若确认要沿用旧解，只能先把 solutions.npz 备份回来，"
                  "再跳过 build_cache() 直接调用 load_cache()；")
            print("        A_FORCE_CACHE=0 并不代表沿用旧解（版本/指纹不符时仍会重算）。")
    print("[cache] 开始计算四问解 ...")
    N = 64                                # 谱节点数：Np = 65 个 CGL 节点，两场共 130 个状态量
    d = {}

    # ---- 问题1：0~1800 s，报告点间隔 10 s ----
    s1 = HerbSolver(N, prob=1)
    t1 = np.arange(0.0, 1800.1, 10.0)
    # atol 按状态块分开给：前 Np 个是水分浓度 C [kg/kg]（量级 1），后 Np 个是温度
    # T [degC]（量级 3e2），两块量级差两个数量级，同一个绝对容差对两者并不等价。
    at1 = np.concatenate([np.full(s1.Np, 1e-13), np.full(s1.Np, 1e-11)])
    sol1 = s1.solve(1800.0, t_eval=t1, rtol=1e-11, atol=at1)
    d["p1_t"], d["p1_y"] = sol1.t, sol1.y

    # ---- 问题2：0~10800 s，报告点间隔 60 s ----
    s2 = HerbSolver(N, prob=2)
    t2 = np.arange(0.0, 10800.1, 60.0)
    at2 = np.concatenate([np.full(s2.Np, 1e-12), np.full(s2.Np, 1e-10)])
    sol2 = s2.solve(10800.0, t_eval=t2, rtol=1e-11, atol=at2)
    d["p2_t"], d["p2_y"] = sol2.t, sol2.y

    # ---- 问题3：先只求解到烘干时刻，再按均匀时刻表回放，供画光滑曲线 ----
    CRIT = 0.15                            # 烘干判据：中心含水率降到 0.15 kg/kg（干基）
    s3 = HerbSolver(N, prob=3)
    at3 = np.concatenate([np.full(s3.Np, 1e-13), np.full(s3.Np, 1e-11)])

    def ev(t, y):
        # y[:Np] 是水分浓度场，其最大值出现在干燥最慢的中心点。terminal 让积分在事件处
        # 停下，不必白算到 t=300000 s；direction=-1 只认由大变小穿越（初值 2.55 kg/kg
        # 远高于阈值，不会误触发）。
        return np.max(y[:s3.Np]) - CRIT
    ev.terminal = True
    ev.direction = -1
    e3 = s3.solve(300000.0, rtol=1e-11, atol=at3, events=ev)
    td3 = float(e3.t_events[0][0])
    t3 = np.linspace(0.0, td3, 1400)        # 1400 个点足够把 60 h 量级的曲线画光滑
    sol3 = s3.solve(td3, t_eval=list(t3), rtol=1e-11, atol=at3)
    d["p3_t"], d["p3_y"], d["p3_tdry"] = sol3.t, sol3.y, np.array([td3])

    # ---- 问题4：物性换成附录 4，并开启附件 2 的实测收缩 R(t) ----
    s4 = HerbSolver(N, prob=4, shrink=True)
    at4 = np.concatenate([np.full(s4.Np, 1e-13), np.full(s4.Np, 1e-11)])

    def ev4(t, y):
        # 事件定义同问题 3 的 ev；问题 4 的判据仍取当前半径下的中心含水率
        return np.max(y[:s4.Np]) - CRIT
    ev4.terminal = True
    ev4.direction = -1
    e4 = s4.solve(300000.0, rtol=1e-11, atol=at4, events=ev4)
    td4 = float(e4.t_events[0][0])
    t4 = np.linspace(0.0, td4, 1400)
    sol4 = s4.solve(td4, t_eval=list(t4), rtol=1e-11, atol=at4)
    d["p4_t"], d["p4_y"], d["p4_tdry"] = sol4.t, sol4.y, np.array([td4])

    # 版本与指纹随解一起落盘，供下次启动时判定缓存是否仍然对应当前代码
    d["version"] = np.array(CACHE_VERSION)
    d["fingerprint"] = np.array(CACHE_FINGERPRINT)
    np.savez_compressed(cache, **d)
    print("[cache] 完成：", cache)
    print("[cache] 版本 %s，指纹 %s" % (CACHE_VERSION, CACHE_FINGERPRINT))


def load_cache():
    """读出 中间数据/solutions.npz，返回 (数据字典 z, HerbSolver 类, common 模块)。

    连带返回后两者，是为了让出图脚本一次拿到画图所需的全部入口；本函数只做读取，
    版本号与源码指纹的校验在 build_cache() 里，故调用前必须先跑一次 build_cache()。
    """
    import common as cm
    from spectral import HerbSolver
    z = np.load(os.path.join(DATADIR, "solutions.npz"))
    return z, HerbSolver, cm
