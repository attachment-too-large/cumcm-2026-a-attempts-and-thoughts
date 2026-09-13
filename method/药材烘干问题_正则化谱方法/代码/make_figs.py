# -*- coding: utf-8 -*-
"""生成 PDF 用的插图 fig1/fig2（中文字体：Microsoft YaHei）。"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import common as cm
from spectral import HerbSolver

# ---- 中文字体 ----
FP = r"C:\Windows\Fonts\msyh.ttc"
fm.fontManager.addfont(FP)
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10.5
plt.rcParams["savefig.dpi"] = 220
plt.rcParams["figure.dpi"] = 220

OUT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "报告", "fig"))
os.makedirs(OUT, exist_ok=True)

RCOLS = np.array([0.0, 0.005, 0.010, 0.015, 0.020])
CRIT = 0.15

# ================== 问题1：剖面演化 ==================
s1 = HerbSolver(64, prob=1)
t1 = np.arange(0.0, 1800.1, 60.0)
at1 = np.concatenate([np.full(s1.Np, 1e-13), np.full(s1.Np, 1e-11)])
sol1 = s1.solve(1800.0, t_eval=t1, rtol=1e-11, atol=at1)
U1 = (RCOLS / cm.R0) ** 2
T1 = np.array([s1.interp_u(sol1.y[:, i], U1)[1] for i in range(len(t1))])
C1 = np.array([s1.interp_u(sol1.y[:, i], U1)[0] for i in range(len(t1))])

fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.5))
rr = RCOLS * 100
for k, tt in enumerate([100, 600, 1200, 1800]):
    i = int(np.argmin(np.abs(t1 - tt)))
    ax[0].plot(rr, T1[i], "o-", ms=3.5, lw=1.4, label="t=%d s" % tt)
    ax[1].plot(rr, C1[i], "s-", ms=3.5, lw=1.4, label="t=%d s" % tt)
ax[0].set_xlabel("到药材中心的距离 r / cm"); ax[0].set_ylabel("温度 T / °C")
ax[0].set_title("(a) 问题1 温度剖面演化")
ax[1].set_xlabel("到药材中心的距离 r / cm"); ax[1].set_ylabel("水分浓度 C / (kg/kg)")
ax[1].set_title("(b) 问题1 水分浓度剖面演化")
for a in ax:
    a.grid(alpha=.3); a.legend(fontsize=8.5)
fig.tight_layout(); fig.savefig(OUT + r"\fig1_p1_profiles.png"); plt.close(fig)

# ================== 问题3：干燥曲线 ==================
s3 = HerbSolver(48, prob=3)
at3 = np.concatenate([np.full(s3.Np, 1e-13), np.full(s3.Np, 1e-11)])


def ev3(t, y):
    return np.max(y[:s3.Np]) - CRIT
ev3.terminal = True; ev3.direction = -1
solE = s3.solve(300000.0, rtol=1e-10, atol=at3, events=ev3)
tdry3 = float(solE.t_events[0][0])

tv = np.linspace(0.0, tdry3 * 1.06, 2000)
sol3 = s3.solve(tdry3 * 1.06, t_eval=list(tv), rtol=1e-10, atol=at3)
U3 = (RCOLS / cm.R0) ** 2
C3 = np.array([s3.interp_u(sol3.y[:, i], U3)[0] for i in range(len(tv))])
hh = tv / 3600.0

fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.5))
for j, r in enumerate(RCOLS):
    ax[0].semilogy(hh, C3[:, j], lw=1.4,
                   label=("r=%.1f cm" % (r * 100)) if j in (0, 2, 4) else None)
ax[0].axhline(CRIT, color="r", ls="--", lw=1.1)
ax[0].axvline(tdry3 / 3600.0, color="k", ls=":", lw=1.1)
ax[0].annotate("判据 0.15", xy=(2, 0.16), color="r", fontsize=9)
ax[0].annotate("$t_{dry}$=%.2f h" % (tdry3 / 3600.0),
               xy=(tdry3 / 3600.0, 1.6), fontsize=9)
ax[0].set_xlabel("时间 t / h"); ax[0].set_ylabel("水分浓度 C / (kg/kg)")
ax[0].set_title("(a) 问题3：各半径处干燥曲线（半对数）")
ax[0].set_xlim(0, 60); ax[0].set_ylim(0.04, 3.2)
ax[0].grid(alpha=.3, which="both"); ax[0].legend(fontsize=8.5, loc="lower left")

zoom = (hh > 50) & (hh < 60)
ax[1].plot(hh[zoom], C3[zoom, 0], lw=1.8, label="中心 C(0,t)")
ax[1].plot(hh[zoom], C3[zoom, -1], lw=1.8, label="表面 C(R,t)")
ax[1].axhline(CRIT, color="r", ls="--", lw=1.2)
ax[1].axvline(tdry3 / 3600.0, color="k", ls=":", lw=1.2)
ax[1].plot([tdry3 / 3600.0], [CRIT], "r*", ms=11)
ax[1].annotate("$t_{dry}$=%.3f h" % (tdry3 / 3600.0),
               xy=(tdry3 / 3600.0 - 8.0, 0.158), fontsize=9)
ax[1].set_xlabel("时间 t / h"); ax[1].set_ylabel("水分浓度 C / (kg/kg)")
ax[1].set_title("(b) 干燥终点附近放大")
ax[1].grid(alpha=.3); ax[1].legend(fontsize=8.5)
fig.tight_layout(); fig.savefig(OUT + r"\fig2_p3_drying.png"); plt.close(fig)

print("fig1/fig2 done")
