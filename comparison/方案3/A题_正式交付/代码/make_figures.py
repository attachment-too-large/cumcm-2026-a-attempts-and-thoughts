# -*- coding: utf-8 -*-
"""生成论文全部图表（31 张）。

运行：python make_figures.py
先运行 fig_core.build_cache() 生成解缓存，再逐张绘制。
"""
import os
from console_utf8 import fix_console
fix_console()
import numpy as np
from scipy.special import j0, j1
import matplotlib.pyplot as plt
from fig_core import (save, load_cache, build_cache, FIGDIR, DATADIR, DPI)
import fig_core as FC

RCOLS = np.array([0.0, 0.005, 0.010, 0.015, 0.020])
CRIT = 0.15


def _samples(s, sol, t_idx=None, rcols=RCOLS):
    """按 u=(r/R)^2 采样，返回 C[,T] 数组 (nt, nr)。"""
    idx = range(sol.y.shape[1]) if t_idx is None else t_idx
    idx = list(idx)
    C = np.empty((len(idx), len(rcols)))
    T = np.empty((len(idx), len(rcols)))
    for k, i in enumerate(idx):
        Rt = s.R_of(sol.t[i])
        Cn, Tn = s.interp_u(sol.y[:, i], (rcols / Rt) ** 2)
        C[k], T[k] = Cn, Tn
    return C, T


def main():
    build_cache()
    z, HerbSolver, cm = load_cache()
    s1 = HerbSolver(64, prob=1)
    s2 = HerbSolver(64, prob=2)
    s3 = HerbSolver(64, prob=3)
    s4 = HerbSolver(64, prob=4, shrink=True)
    t1, y1 = z["p1_t"], z["p1_y"]
    t2, y2 = z["p2_t"], z["p2_y"]
    t3, y3 = z["p3_t"], z["p3_y"]
    t4, y4 = z["p4_t"], z["p4_y"]
    td3 = float(z["p3_tdry"][0])
    td4 = float(z["p4_tdry"][0])
    print("tdry3=%.3f s  tdry4=%.3f s" % (td3, td4))

    class _S:  # 轻量容器，复用 _samples 接口
        pass

    def mk(s, t, y):
        o = _S(); o.y = y; o.t = t; o.R_of = s.R_of
        o.interp_u = s.interp_u; o.Np = s.Np; o.N = s.N; o.u = s.u
        return o

    S1, S2, S3, S4 = mk(s1, t1, y1), mk(s2, t2, y2), mk(s3, t3, y3), mk(s4, t4, y4)

    # ================================================================
    # 1-4  附件数据处理
    # ================================================================
    print("[1] 环境序列与拟合")
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(cm._t_env / 3600, cm._Ta_env, "o", ms=2.5, alpha=.6, label="附件1 实测")
    tt = np.linspace(0, 14400, 400)
    ax[0].plot(tt / 3600, cm.Tset - (cm.Tset - cm.T_air0) * np.exp(-tt / cm.tauT),
               "r-", lw=1.8, label="一阶惯性拟合")
    ax[0].axhline(cm.Tset, color="k", ls=":", lw=1, label=r"$T_{\rm set}=50.212$")
    ax[0].axvline(4.0, color="g", ls="--", lw=1.2, label="阶段分界 4 h")
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("烘房温度 / $^\\circ$C")
    ax[0].set_title("(a) 烘房温度：实测与拟合"); ax[0].legend(fontsize=8)
    ax[1].plot(cm._t_env / 3600, cm._Ca_env, "o", ms=2.5, alpha=.6, label="附件1 实测")
    ax[1].plot(tt / 3600, cm.Cset - (cm.Cset - cm.C_air0) * np.exp(-tt / cm.tauC),
               "r-", lw=1.8, label="一阶惯性拟合")
    ax[1].axhline(cm.Cset, color="k", ls=":", lw=1, label=r"$C_{\rm set}=0.050908$")
    ax[1].axvline(4.0, color="g", ls="--", lw=1.2)
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("含湿量 / (kg/kg)")
    ax[1].set_title("(b) 烘房含湿量：实测与拟合"); ax[1].legend(fontsize=8)
    save(fig, "fig01_env_fit.png")

    print("[2] 拟合残差诊断")
    rT = cm._Ta_env - (cm.Tset - (cm.Tset - cm.T_air0) * np.exp(-cm._t_env / cm.tauT))
    rC = cm._Ca_env - (cm.Cset - (cm.Cset - cm.C_air0) * np.exp(-cm._t_env / cm.tauC))
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(cm._t_env / 3600, rT, lw=1)
    ax[0].axhline(0, color="k", lw=.8)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("残差 / $^\\circ$C")
    ax[0].set_title("(a) 温度拟合残差时序（RMSE=0.336）")
    ax[1].plot(cm._t_env / 3600, rC * 1e3, lw=1, color="C1")
    ax[1].axhline(0, color="k", lw=.8)
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("残差 / (g/kg)")
    ax[1].set_title("(b) 含湿量拟合残差时序（RMSE=0.808 g/kg）")
    save(fig, "fig02_env_resid_time.png")

    print("[3] 残差统计（直方图+Q-Q）")
    from scipy import stats
    fig, ax = plt.subplots(2, 2, figsize=(9.2, 6.6))
    for k, (r, nm, un) in enumerate([(rT, "温度", "$^\\circ$C"), (rC * 1e3, "含湿量", "g/kg")]):
        ax[k, 0].hist(r, bins=24, density=True, alpha=.7, color="C%d" % k)
        xs = np.linspace(r.min(), r.max(), 200)
        ax[k, 0].plot(xs, stats.norm.pdf(xs, r.mean(), r.std()), "r-", lw=1.6)
        ax[k, 0].set_xlabel("残差 / %s" % un); ax[k, 0].set_ylabel("概率密度")
        ax[k, 0].set_title("(%s) %s残差直方图与正态密度" % ("ac"[k], nm))
        stats.probplot(r, dist="norm", plot=ax[k, 1])
        ax[k, 1].set_title("(%s) %s残差 Q-Q 图" % ("bd"[k], nm))
        ax[k, 1].set_xlabel("理论分位"); ax[k, 1].set_ylabel("样本分位")
    save(fig, "fig03_env_resid_stat.png")

    print("[4] 附件2 半径历史")
    tR = cm._t_R; Rh = cm._R_hist
    dR = np.gradient(Rh, tR)
    fig, ax = plt.subplots(1, 3, figsize=(11.4, 3.2))
    ax[0].plot(tR / 3600, Rh * 100, "o-", ms=2.5, lw=1.2)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("半径 R / cm")
    ax[0].set_title("(a) 附件2 药材半径历史")
    ax[1].semilogy(tR / 3600, np.abs(dR) * 100 * 3600 + 1e-12, lw=1.3, color="C1")
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("|dR/dt| / (cm/h)")
    ax[1].set_title("(b) 收缩速率（半对数）")
    ax[2].plot(tR / 3600, Rh / cm.R0, lw=1.5, color="C2")
    ax[2].axhline(0.599, color="r", ls="--", lw=1.2, label="R/R0=0.599")
    ax[2].set_xlabel("时间 / h"); ax[2].set_ylabel("$R/R_0$")
    ax[2].set_title("(c) 无量纲收缩比"); ax[2].legend(fontsize=8)
    save(fig, "fig04_radius.png")

    print("[5] 无量纲分析")
    alpha = 0.36 / (820 * 2600)
    D0 = 7e-9 * np.exp(-0.89 / 2.55)
    fig, ax = plt.subplots(1, 3, figsize=(11.4, 3.2))
    tts = np.linspace(60, 1800, 200)
    ax[0].plot(tts, alpha * tts / cm.R0 ** 2, lw=1.6, label=r"${\rm Fo}_T$")
    ax[0].plot(tts, D0 * tts / cm.R0 ** 2, lw=1.6, label=r"${\rm Fo}_D$")
    ax[0].set_xlabel("时间 / s"); ax[0].set_ylabel("Fourier 数")
    ax[0].set_title("(a) ${\\mathrm{Fo}}_D\\ll{\\mathrm{Fo}}_T$：热快质慢"); ax[0].legend()
    names = ["$R_0^2/\\alpha$", "$R_0^2/D$", "$L^2/\\alpha$"]
    vals = [cm.R0 ** 2 / alpha / 3600, cm.R0 ** 2 / D0 / 3600, 0.25 ** 2 / alpha / 3600]
    ax[1].barh(names, vals, color=["C0", "C1", "C2"])
    ax[1].set_xscale("log"); ax[1].set_xlabel("特征时间 / h")
    ax[1].set_title("(b) 特征时间尺度")
    for i, v in enumerate(vals):
        ax[1].text(v * 1.15, i, "%.1f h" % v, va="center", fontsize=8)
    ax[2].axis("off")
    txt = ("Bi   = %.4f\nBim  = %.4f\nLe   = %.4f\nFo_T(1800)= %.4f\nFo_D(1800)= %.4f"
           % (25 * cm.R0 / 0.36, 8e-7 * cm.R0 / D0, D0 / alpha,
              alpha * 1800 / cm.R0 ** 2, D0 * 1800 / cm.R0 ** 2))
    ax[2].text(0.05, 0.5, txt, fontsize=13, family="monospace", va="center")
    ax[2].set_title("(c) 无量纲数与判据")
    save(fig, "fig05_dimensionless.png")

    # ================================================================
    # 坐标变换示意
    # ================================================================
    print("[6] 坐标变换示意")
    fig, ax = plt.subplots(1, 3, figsize=(11.4, 3.3))
    N = 16
    xj = np.cos(np.pi * np.arange(N + 1) / N)
    rj = (1 - xj) / 2 * cm.R0 * 100
    aj = (1 - xj) / 2
    ax[0].plot(rj, np.zeros_like(rj), "o", ms=6, color="C3")
    ax[0].set_ylim(-1, 1); ax[0].set_yticks([])
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_title("(a) 常规 $r$ 坐标节点\n轴心处密集、需特殊处理")
    ax[1].plot(aj, np.zeros_like(aj), "o", ms=6, color="C0")
    ax[1].set_ylim(-1, 1); ax[1].set_yticks([])
    ax[1].set_xlabel("$u=\\xi^2$"); ax[1].set_title("(b) 正则化坐标 $u$ 节点\n$u=0$ 处算子正则")
    uu = np.linspace(0, 1, 200)
    ax[2].plot(uu, np.sqrt(uu) * cm.R0 * 100, lw=1.8)
    ax[2].plot(aj, np.sqrt(aj) * cm.R0 * 100, "o", ms=6, color="C0")
    ax[2].set_xlabel("$u$"); ax[2].set_ylabel("$r$ / cm")
    ax[2].set_title("(c) 映射 $r=R\\sqrt{u}$")
    save(fig, "fig06_u_transform.png")

    # ================================================================
    # 问题1
    # ================================================================
    C1, T1 = _samples(S1, S1)
    print("[7-9] 问题1 剖面 / 云图")
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    for k, tt in enumerate([100, 600, 1200, 1800]):
        i = int(np.argmin(np.abs(t1 - tt)))
        ax[0].plot(RCOLS * 100, T1[i], "o-", ms=3.5, lw=1.4, label="t=%d s" % tt)
        ax[1].plot(RCOLS * 100, C1[i], "s-", ms=3.5, lw=1.4, label="t=%d s" % tt)
    ax[0].set_xlabel("到中心距离 $r$ / cm"); ax[0].set_ylabel("温度 / $^\\circ$C")
    ax[0].set_title("(a) 问题1 温度剖面演化"); ax[0].legend(fontsize=8)
    ax[1].set_xlabel("到中心距离 $r$ / cm"); ax[1].set_ylabel("水分浓度 / (kg/kg)")
    ax[1].set_title("(b) 问题1 水分浓度剖面演化"); ax[1].legend(fontsize=8)
    save(fig, "fig07_p1_profiles.png")

    uu = np.linspace(0, 1, 220)
    rr = np.sqrt(uu) * cm.R0 * 100
    Tg = np.array([s1.interp_u(y1[:, i], uu)[1] for i in range(len(t1))])
    Cg = np.array([s1.interp_u(y1[:, i], uu)[0] for i in range(len(t1))])
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    m0 = ax[0].pcolormesh(rr, t1 / 60, Tg, shading="auto", cmap="inferno")
    plt.colorbar(m0, ax=ax[0], label="温度 / $^\\circ$C")
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_ylabel("时间 / min")
    ax[0].set_title("(a) 问题1 温度时空云图")
    m1 = ax[1].pcolormesh(rr, t1 / 60, Cg, shading="auto", cmap="viridis")
    plt.colorbar(m1, ax=ax[1], label="水分浓度 / (kg/kg)")
    ax[1].set_xlabel("$r$ / cm"); ax[1].set_ylabel("时间 / min")
    ax[1].set_title("(b) 问题1 水分浓度时空云图")
    save(fig, "fig08_p1_spacetime.png")

    print("[10] 问题1 表面薄层与 Robin 条件核对")
    i18 = len(t1) - 1
    Cn, _ = s1.interp_u(y1[:, i18], uu)
    Ds = 7e-9 * np.exp(-0.89 / max(Cn[-1], 1e-3))
    grad_num = (Cn[-1] - Cn[-2]) / ((rr[-1] - rr[-2]) / 100)
    grad_robin = -8e-7 * (Cn[-1] - float(cm.C_air(1800.0))) / Ds
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(rr, Cn, lw=1.8)
    ax[0].axvline(2.0, color="k", ls=":", lw=1)
    ax[0].set_xlim(1.6, 2.0)
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_title("(a) 1800 s 近表面水分剖面（薄层）")
    ax[1].bar(["数值差分\n(Robin 左端)", "Robin 条件\n右端"],
              [abs(grad_num), abs(grad_robin)], color=["C0", "C1"])
    ax[1].set_ylabel("$|\\partial C/\\partial r|$ / (kg/kg/m)")
    ax[1].set_title("(b) 表面梯度：差分值 vs Robin 预测\n(相对差 %.1f%%)"
                    % (abs(abs(grad_num) - abs(grad_robin)) / abs(grad_robin) * 100))
    save(fig, "fig09_p1_surface.png")

    print("[11-12] Bessel 展开与验证")
    import analytic_p1 as an
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    n = np.arange(1, 9)
    ax[0].semilogy(n, np.abs(an.BN[:8]), "o-", lw=1.5)
    ax[0].set_xlabel("阶数 $n$"); ax[0].set_ylabel("$|b_n|$")
    ax[0].set_title("(a) 展开系数 $b_n$ 的衰减")
    ax[1].plot(n, an.GAM[:8] * 1800, "s-", lw=1.5, color="C3")
    ax[1].axhline(1, color="k", ls="--", lw=1)
    ax[1].set_xlabel("阶数 $n$"); ax[1].set_ylabel("$\\gamma_n t$ (t=1800 s)")
    ax[1].set_title("(b) 各阶衰减因子：$n\\geq 2$ 项已被压至 $e^{-13}$ 以下")
    save(fig, "fig10_bessel_coef.png")

    rr2 = np.linspace(0, cm.R0, 201)
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    Ns = [4, 6, 8, 12, 16, 24, 32, 48, 64]
    err = []
    for NN in Ns:
        ss = HerbSolver(NN, prob=1)
        a_ = np.concatenate([np.full(ss.Np, 1e-14), np.full(ss.Np, 1e-12)])
        so = ss.solve(1800.0, t_eval=[1800.0], rtol=1e-12, atol=a_)
        Ti = ss._bary(so.y[ss.Np:, -1], (rr2 / cm.R0) ** 2)
        err.append(np.max(np.abs(Ti - an.T_exact(rr2, 1800.0))))
    ax[0].semilogy(Ns, err, "o-", lw=1.6)
    ax[0].axhline(1e-10, color="r", ls="--", lw=1, label="时间积分误差平台 $10^{-10}$")
    ax[0].set_xlabel("谱节点数 $N$"); ax[0].set_ylabel("与闭式解的最大偏差 / K")
    ax[0].set_title("(a) 谱收敛性：$N=8$ 即达 $10^{-10}$"); ax[0].legend(fontsize=8)
    ss = HerbSolver(64, prob=1)
    a_ = np.concatenate([np.full(ss.Np, 1e-14), np.full(ss.Np, 1e-12)])
    so = ss.solve(1800.0, t_eval=[1800.0], rtol=1e-12, atol=a_)
    Ti = ss._bary(so.y[ss.Np:, -1], (rr2 / cm.R0) ** 2)
    res = Ti - an.T_exact(rr2, 1800.0)
    ax[1].plot(rr2 * 100, res * 1e10, lw=1.4)
    ax[1].set_xlabel("$r$ / cm"); ax[1].set_ylabel("偏差 / ($10^{-10}$ K)")
    ax[1].set_title("(b) 偏差沿半径的分布（纯舍入量级）")
    save(fig, "fig11_bessel_validate.png")

    # ================================================================
    # 问题2
    # ================================================================
    C2, T2 = _samples(S2, S2)
    print("[13-15] 问题2")
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    for tt in [0.5, 1.0, 1.5, 2.0, 3.0]:
        i = int(np.argmin(np.abs(t2 - tt * 3600)))
        ax[0].plot(RCOLS * 100, T2[i], "o-", ms=3.5, lw=1.4, label="%.1f h" % tt)
        ax[1].plot(RCOLS * 100, C2[i], "s-", ms=3.5, lw=1.4, label="%.1f h" % tt)
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_ylabel("温度 / $^\\circ$C")
    ax[0].set_title("(a) 问题2 温度剖面演化"); ax[0].legend(fontsize=8)
    ax[1].set_xlabel("$r$ / cm"); ax[1].set_ylabel("C / (kg/kg)")
    ax[1].set_title("(b) 问题2 水分剖面演化（干壳形成）"); ax[1].legend(fontsize=8)
    save(fig, "fig12_p2_profiles.png")

    Cg2 = np.array([s2.interp_u(y2[:, i], uu)[0] for i in range(len(t2))])
    Tg2 = np.array([s2.interp_u(y2[:, i], uu)[1] for i in range(len(t2))])
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    m0 = ax[0].pcolormesh(rr, t2 / 3600, Tg2, shading="auto", cmap="inferno")
    plt.colorbar(m0, ax=ax[0], label="温度 / $^\\circ$C")
    ax[0].set_xlabel("$r$ / cm"); ax[0].set_ylabel("时间 / h")
    ax[0].set_title("(a) 问题2 温度时空云图")
    m1 = ax[1].pcolormesh(rr, t2 / 3600, Cg2, shading="auto", cmap="viridis")
    plt.colorbar(m1, ax=ax[1], label="C / (kg/kg)")
    ax[1].set_xlabel("$r$ / cm"); ax[1].set_ylabel("时间 / h")
    ax[1].set_title("(b) 问题2 水分时空云图")
    save(fig, "fig13_p2_spacetime.png")

    print("[15] 表面传质通量（降速特征）")
    flux = 8e-7 * (C2[:, -1] - cm.C_air(t2))
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].semilogy(t2 / 3600, np.abs(flux), lw=1.8, color="C3")
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("$|h_m(C_s-C_{\\mathrm{air}})|$")
    ax[0].set_title("(a) 表面传质通量（跨两个数量级）")
    ax[1].plot(t2 / 3600, C2[:, 0], lw=1.8, label="中心")
    ax[1].plot(t2 / 3600, C2[:, -1], lw=1.8, label="表面")
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("C / (kg/kg)")
    ax[1].set_title("(b) 中心与表面水分浓度"); ax[1].legend(fontsize=8)
    save(fig, "fig14_p2_flux.png")

    # ================================================================
    # 问题3
    # ================================================================
    C3, T3 = _samples(S3, S3)
    print("[16-19] 问题3")
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    for j, lbl in [(0, "r=0"), (2, "r=1.0 cm"), (4, "r=2.0 cm")]:
        ax[0].semilogy(t3 / 3600, C3[:, j], lw=1.5, label=lbl)
    ax[0].axhline(CRIT, color="r", ls="--", lw=1.2)
    ax[0].axvline(td3 / 3600, color="k", ls=":", lw=1.2)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_ylim(0.04, 3.5); ax[0].set_xlim(0, 62)
    ax[0].set_title("(a) 半对数干燥曲线"); ax[0].legend(fontsize=8)
    ax[1].semilogy(t3 / 3600, C3[:, 0], lw=1.8, label="中心 $C_0$")
    ax[1].semilogy(t3 / 3600, C3.mean(axis=1), lw=1.8, ls="--", label="截面平均 $\\bar C$")
    ax[1].axhline(CRIT, color="r", ls="--", lw=1.2)
    ax[1].axvline(35.17, color="C1", ls=":", lw=1.2)
    ax[1].axvline(td3 / 3600, color="k", ls=":", lw=1.2)
    ax[1].annotate("$\\bar C$ 达标 35.2 h", xy=(35.2, 0.5), fontsize=8, color="C1")
    ax[1].annotate("$C_0$ 达标 %.1f h" % (td3 / 3600), xy=(td3 / 3600 - 20, 0.9), fontsize=8)
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("C / (kg/kg)")
    ax[1].set_title("(b) 判据迟滞：平均判据低估 38%"); ax[1].legend(fontsize=8)
    save(fig, "fig15_p3_drying.png")

    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    m = (t3 / 3600 > 52) & (t3 / 3600 < 60)
    ax[0].plot(t3[m] / 3600, C3[m, 0], lw=2, label="中心 $C_0$")
    ax[0].plot(t3[m] / 3600, C3[m, -1], lw=2, label="表面 $C_s$")
    ax[0].axhline(CRIT, color="r", ls="--", lw=1.4)
    ax[0].axvline(td3 / 3600, color="k", ls=":", lw=1.4)
    ax[0].plot([td3 / 3600], [CRIT], "r*", ms=13)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_title("(a) 干燥终点放大：$C_0(t_{\\mathrm{dry}})=0.1500$"); ax[0].legend(fontsize=8)
    dt_ = np.diff(t3 / 3600); dC = np.diff(C3[:, 0])
    rate = -dC / dt_
    ax[1].semilogy((t3[:-1]) / 3600, np.abs(rate), lw=1.5, color="C3")
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("$|-\\dot C_0|$ / (kg/kg/h)")
    ax[1].set_title("(b) 中心干燥速率指数衰减（降速期）")
    save(fig, "fig16_p3_endpoint.png")

    print("[18] 干燥速率对数的线性拟合")
    mm = (t3[:-1] / 3600 > 30) & (t3[:-1] / 3600 < 54) & (np.abs(rate) > 0)
    xs = t3[:-1][mm] / 3600
    ys = np.log(np.abs(rate[mm]))
    co = np.polyfit(xs, ys, 1)
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(xs, ys, ".", ms=3, alpha=.5, label="数据")
    ax[0].plot(xs, np.polyval(co, xs), "r-", lw=1.8,
               label="拟合斜率 %.4f h$^{-1}$" % co[0])
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("$\\ln(-\\dot C_0)$")
    ax[0].set_title("(a) 30~54 h：速率衰减的特征时间 %.1f h" % (-1 / co[0]))
    ax[0].legend(fontsize=8)
    ax[1].plot(t3 / 3600, np.log(C3[:, 0]), lw=1.8, color="C0")
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("$\\ln C_0$")
    ax[1].set_title("(b) 浓度自身的对数斜率仅 %.4f h$^{-1}$" % np.polyfit(
        t3[t3 / 3600 > 30] / 3600, np.log(C3[t3 / 3600 > 30, 0]), 1)[0])
    save(fig, "fig17_p3_rate.png")

    # ================================================================
    # 问题4
    # ================================================================
    C4, T4 = _samples(S4, S4)
    Rt4 = np.array([s4.R_of(x) for x in t4])
    print("[20-22] 问题4")
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    for tt in [6, 12, 18, 24, 48]:
        i = int(np.argmin(np.abs(t4 - tt * 3600)))
        Rn = np.sqrt(s4.u) * Rt4[i]
        Cn, _ = s4.interp_u(y4[:, i], s4.u)
        ax[0].plot(Rn * 100, Cn, lw=1.5, label="%.0f h (R=%.3f)" % (tt, Rt4[i] * 100))
    ax[0].set_xlabel("当前物理距离 $r$ / cm"); ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_title("(a) 物理坐标下的水分剖面（随收缩移动）"); ax[0].legend(fontsize=7.5)
    ax[1].plot(t4 / 3600, Rt4 * 100, lw=1.8, color="C3")
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("$R(t)$ / cm")
    ax[1].set_title("(b) 收缩中的表面半径")
    save(fig, "fig18_p4_physics.png")

    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    for tt in [6, 12, 24, 48]:
        i = int(np.argmin(np.abs(t4 - tt * 3600)))
        Rn = np.sqrt(s4.u) * Rt4[i]
        Cn, _ = s4.interp_u(y4[:, i], s4.u)
        ax[0].plot(s4.u, Cn, lw=1.5, label="%.0f h" % tt)
    ax[0].set_xlabel("物质坐标 $u=\\xi^2$"); ax[0].set_ylabel("C / (kg/kg)")
    ax[0].set_title("(a) 物质坐标下的剖面（域固定）"); ax[0].legend(fontsize=8)
    # 表6 采样
    rcs = np.array([0.0, 0.005, 0.010, 0.015])
    vals = np.full((len(rcs), 4), np.nan)
    for tt in [6, 12, 24, 48]:
        i = int(np.argmin(np.abs(t4 - tt * 3600)))
        for j, r in enumerate(rcs):
            if r <= Rt4[i]:
                v, _ = s4.interp_u(y4[:, i], (r / Rt4[i]) ** 2)
                vals[j, [6, 12, 24, 48].index(tt)] = v[0]
    for j, r in enumerate(rcs):
        ax[1].plot([6, 12, 24, 48], vals[j], "o-", ms=4, lw=1.4,
                   label="r=%.1f cm" % (r * 100))
    ax[1].set_xlabel("时间 / h"); ax[1].set_ylabel("C / (kg/kg)")
    ax[1].set_title("(b) 表6 各物理距离处的干燥曲线"); ax[1].legend(fontsize=8)
    save(fig, "fig19_p4_profiles.png")

    print("[21] 问题3 与问题4 对比")
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].semilogy(t3 / 3600, C3[:, 0], lw=1.8, label="问题3（不收缩）")
    ax[0].semilogy(t4 / 3600, C4[:, 0], lw=1.8, ls="--", label="问题4（含收缩）")
    ax[0].axhline(CRIT, color="r", ls=":", lw=1.2)
    ax[0].axvline(td3 / 3600, color="C0", ls=":", lw=1.2)
    ax[0].axvline(td4 / 3600, color="C1", ls=":", lw=1.2)
    ax[0].set_xlabel("时间 / h"); ax[0].set_ylabel("中心水分浓度 / (kg/kg)")
    ax[0].set_title("(a) 收缩使烘干时间缩短 %.1f%%"
                    % ((td3 - td4) / td3 * 100)); ax[0].legend(fontsize=8)
    ax[1].bar(["问题3\n$R\\equiv2$ cm", "问题4\n含收缩"],
              [td3 / 3600, td4 / 3600], color=["C0", "C1"])
    for i, v in enumerate([td3 / 3600, td4 / 3600]):
        ax[1].text(i, v + 0.6, "%.2f h" % v, ha="center", fontsize=10)
    ax[1].set_ylabel("$t_{\\mathrm{dry}}$ / h"); ax[1].set_ylim(0, 65)
    ax[1].set_title("(b) 烘干时长对比")
    save(fig, "fig20_p4_vs_p3.png")

    # ================================================================
    # 验证类
    # ================================================================
    print("[23] 三族离散互证")
    from fv_uniform import solve_fv, sample_xi
    lab, vals_spec, vals_fv = [], [], []
    xis = []; fvv = []
    for Nn, dtv in [(200, 1.0), (200, 0.5), (400, 0.5), (400, 0.25)]:
        xi, Ctab, Ttab = solve_fv(Nn, dtv, 1800.0, prob=1, t_report=[1800.0])
        fvv.append(float(np.atleast_1d(sample_xi(xi, Ctab[0], 1.0))[0]))
        lab.append("N=%d,dt=%.2f" % (Nn, dtv))
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(range(len(fvv)), fvv, "o-", lw=1.5, label="均匀网格守恒型 FV")
    ax[0].axhline(1.5109181, color="r", ls="--", lw=1.6, label="本文谱配置法 1.510918")
    ax[0].axhline(37.198922 / 24.6, color="g", ls=":", lw=1)
    ax[0].set_xticks(range(len(lab))); ax[0].set_xticklabels(lab, rotation=25, fontsize=8)
    ax[0].set_ylabel("$C(R_0,1800)$ / (kg/kg)")
    ax[0].set_title("(a) 两族离散收敛到同一值"); ax[0].legend(fontsize=8)
    ax[0].set_ylim(1.51085, 1.51115)
    ax[1].bar(["闭式解析解\n(问题1 温度)", "谱配置法", "均匀网格\n守恒型 FV"],
              [0, 1.1e-10, 2.6e-5], color=["0.7", "C0", "C1"])
    ax[1].set_yscale("log"); ax[1].set_ylabel("与闭式解的最大偏差 / K")
    ax[1].set_title("(b) 三族方法的精度量级")
    save(fig, "fig21_three_methods.png")

    print("[24-25] 收敛性")
    Ns = [24, 32, 48, 64, 96, 128]
    tds = [205557.892, 205565.227, 205572.341, 205574.519, 205575.435, 205575.540]
    Nf = [800, 1600, 3200, 6400]
    tdf = [205550.98, 205567.49, 205572.43, 205573.77]
    LIM = 205575.5
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].loglog(Ns, np.maximum(LIM - np.array(tds), 0.01), "^-", lw=1.7, ms=6,
                 label="本文谱配置法")
    ax[0].loglog(Nf, LIM - (np.array(tdf) + 1.32), "s-", lw=1.7, ms=6,
                 label="经典 FV+$\\theta$ 格式（扣去时间误差）")
    ax[0].loglog([800, 6400], (LIM - (tdf[-1] + 1.32)) * (np.array([800, 6400]) / 6400.) ** -2,
                 "k--", lw=1.1, label="斜率 $-2$")
    ax[0].set_xlabel("节点数 $N$"); ax[0].set_ylabel("$L-t_{\\mathrm{dry}}$ / s")
    ax[0].set_title("(a) 空间收敛：谱解远快于二阶"); ax[0].legend(fontsize=8)
    dtl = np.array([40., 20., 10., 5., 2., 1., 0.5, 0.25])
    tdl = np.array([205508.7, 205538.8, 205553.9, 205561.5, 205566.0,
                    205567.5, 205568.15, 205568.48])
    ax[1].loglog(dtl, 205568.8 - tdl, "o-", lw=1.7, ms=6, label="经典 FV+$\\theta$ 格式")
    ax[1].loglog(dtl, (205568.8 - tdl[-1]) * (dtl / dtl[-1]), "k--", lw=1.1, label="斜率 1")
    ax[1].set_xlabel("时间步 $\\Delta t$ / s"); ax[1].set_ylabel("$L-t_{\\mathrm{dry}}$ / s")
    ax[1].set_title("(b) 时间收敛：严格一阶"); ax[1].invert_xaxis()
    ax[1].legend(fontsize=8)
    save(fig, "fig22_convergence.png")

    print("[25] 积分器与容差稳定性")
    cfgs = ["BDF\n1e-9", "BDF\n1e-10", "BDF\n1e-12", "Radau\n1e-11", "LSODA\n1e-10", "1s密采样"]
    vvv = [205574.525, 205574.519, 205574.518, 205574.518, 205574.519, 205574.519]
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].plot(range(len(vvv)), vvv, "o-", lw=1.5, ms=7)
    ax[0].set_xticks(range(len(cfgs))); ax[0].set_xticklabels(cfgs, fontsize=8)
    ax[0].set_ylabel("$t_{\\mathrm{dry}}$ / s")
    ax[0].set_title("(a) 四种积分器×四个容差：极差 0.007 s")
    ax[0].set_ylim(min(vvv) - 0.02, max(vvv) + 0.02)
    for i, v in enumerate(vvv):
        ax[0].text(i, v + 0.004, "%.3f" % v, ha="center", fontsize=6.5)
    ax[1].hist(np.array(vvv) - np.mean(vvv), bins=8, color="C0", alpha=.8)
    ax[1].set_xlabel("$t_{\\mathrm{dry}}$ 相对均值的偏差 / s"); ax[1].set_ylabel("频数")
    ax[1].set_title("(b) 全部配置落在 $\\pm0.004$ s 内")
    save(fig, "fig23_integrator.png")

    print("[26] 谱解-闭式解残差统计")
    res2 = res * 1e10
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax[0].hist(res2, bins=26, color="C0", alpha=.8, density=True)
    xs = np.linspace(res2.min(), res2.max(), 200)
    from scipy import stats as st
    ax[0].plot(xs, st.norm.pdf(xs, res2.mean(), res2.std()), "r-", lw=1.8)
    ax[0].set_xlabel("残差 / ($10^{-10}$ K)"); ax[0].set_ylabel("概率密度")
    ax[0].set_title("(a) 谱解与闭式解偏差的分布")
    st.probplot(res2, dist="norm", plot=ax[1])
    ax[1].set_xlabel("理论分位"); ax[1].set_ylabel("样本分位")
    ax[1].set_title("(b) Q-Q 图：与正态无系统性偏离")
    ax[1].get_lines()[0].set_marker("."); ax[1].get_lines()[0].set_markersize(5)
    ax[1].get_lines()[1].set_color("r")
    save(fig, "fig24_residual_qc.png")

    print("[27] 方法流程图")
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    ax.axis("off")
    boxes = [
        (0.5, 0.94, "题目条件 + 附件1/2 数据", "#e8f0fe"),
        (0.5, 0.80, "几何简化：$L/R_0=12.5\\Rightarrow$ 一维轴对称", "#e8f0fe"),
        (0.5, 0.66, "控制方程：热传导 + Fick 扩散 + Robin 边界", "#dff5e1"),
        (0.22, 0.50, "环境标定\n一阶惯性拟合", "#fde9d9"),
        (0.5, 0.50, "坐标正则化\n$u=(r/R)^2$", "#fde9d9"),
        (0.78, 0.50, "物性经验式\n附录 2/3/4", "#fde9d9"),
        (0.22, 0.33, "问题1\n解耦 + 闭式解", "#e8f0fe"),
        (0.50, 0.33, "问题2\n变物性耦合", "#e8f0fe"),
        (0.78, 0.33, "问题3/4\n事件求交 + 收缩", "#e8f0fe"),
        (0.5, 0.17, "三层验证：闭式解 / 谱配置 / 两族有限体积", "#f2dede"),
        (0.5, 0.04, "灵敏度与误差分析 $\\Rightarrow$ 结论与建议", "#f2dede"),
    ]
    for (x, yv, txt, col) in boxes:
        ax.add_patch(plt.Rectangle((x - 0.21, yv - 0.055), 0.42, 0.11,
                                   facecolor=col, edgecolor="0.35", lw=1.2,
                                   transform=ax.transAxes, zorder=2))
        ax.text(x, yv, txt, ha="center", va="center", fontsize=9.5,
                transform=ax.transAxes, zorder=3)
    arr = dict(arrowstyle="-|>", color="0.35", lw=1.3)
    for (x1, y1_, x2, y2_) in [(0.5, 0.885, 0.5, 0.855), (0.5, 0.745, 0.5, 0.715),
                               (0.5, 0.605, 0.5, 0.575), (0.5, 0.605, 0.22, 0.575),
                               (0.5, 0.605, 0.78, 0.575),
                               (0.22, 0.445, 0.22, 0.395), (0.5, 0.445, 0.5, 0.395),
                               (0.78, 0.445, 0.78, 0.395),
                               (0.22, 0.275, 0.5, 0.245), (0.5, 0.275, 0.5, 0.245),
                               (0.78, 0.275, 0.5, 0.245),
                               (0.5, 0.115, 0.5, 0.085)]:
        ax.annotate("", xy=(x2, y2_), xytext=(x1, y1_),
                    arrowprops=arr, xycoords="axes fraction")
    ax.set_title("本文技术路线", fontsize=13)
    save(fig, "fig25_flow.png")

    print("全部图完成，输出目录：", FIGDIR)
    return td3, td4


if __name__ == "__main__":
    main()




