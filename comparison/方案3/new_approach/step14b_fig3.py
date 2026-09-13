# -*- coding: utf-8 -*-
"""重画 fig3：把时间误差与空间误差分离，使收敛阶一目了然。"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import os

FP = r"C:\Windows\Fonts\msyh.ttc"
fm.fontManager.addfont(FP)
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10.5
plt.rcParams["savefig.dpi"] = 220
OUT = r"C:\Users\qing1\Desktop\数A\new_approach\fig"
os.makedirs(OUT, exist_ok=True)

# ---- 数据来源：step8 / step10 / step11 报告（paper 格式 FV 与本文谱解）----
LIM = 205575.5            # 连续极限（谱解，四种积分器一致）
L_FV_N1600 = 205568.8     # paper 格式 FV，N=1600，Δt->0

dt = np.array([40., 20., 10., 5., 2., 1., 0.5, 0.25])
td = np.array([205508.7, 205538.8, 205553.9, 205561.5, 205566.0,
               205567.5, 205568.15, 205568.48])
err_dt = L_FV_N1600 - td                       # 纯时间误差

N_fv = np.array([800., 1600., 3200., 6400.])
td_fv = np.array([205550.98, 205567.49, 205572.43, 205573.77])
td_fv_corr = td_fv + 1.32                      # 扣掉 Δt=1s 的一阶时间误差
err_N_fv = LIM - td_fv_corr

N_sp = np.array([24., 32., 48., 64., 96., 128.])
td_sp = np.array([205557.892, 205565.227, 205572.341, 205574.519,
                  205575.435, 205575.540])
err_N_sp = np.maximum(LIM - td_sp, 0.01)   # N=128 已越过参考线，取对数下限

fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.6))

# (a) 时间步收敛
ax[0].loglog(dt, err_dt, "o-", lw=1.7, ms=6, label="paper 格式 FV，N=1600")
ax[0].loglog(dt, err_dt[-1] * (dt / dt[-1]) ** 1.0, "k--", lw=1.1,
             label=r"斜率 1（一阶）")
ax[0].set_xlabel(r"时间步 $\Delta t$ / s")
ax[0].set_ylabel(r"$L_{N=1600}-t_{dry}$  / s")
ax[0].set_title(r"(a) 时间收敛：$\mathcal{O}(\Delta t)$ 一阶")
ax[0].grid(alpha=.3, which="both")
ax[0].legend(fontsize=8.5)
ax[0].invert_xaxis()

# (b) 空间收敛
ax[1].loglog(N_fv, err_N_fv, "s-", lw=1.7, ms=6,
             label=r"paper 格式 FV（扣去时间误差）")
ax[1].loglog(N_sp, err_N_sp, "^-", lw=1.7, ms=6, label="本文谱配置法")
ax[1].loglog(np.array([800., 6400.]),
             err_N_fv[-1] * (np.array([800., 6400.]) / 6400.) ** -2.0,
             "k--", lw=1.1, label=r"斜率 $-2$（二阶）")
ax[1].set_xlabel("节点数 $N$")
ax[1].set_ylabel(r"$L-t_{dry}$  / s")
ax[1].set_title("(b) 空间收敛：谱解远快于二阶")
ax[1].grid(alpha=.3, which="both")
ax[1].legend(fontsize=8.5, loc="lower left")

fig.tight_layout()
fig.savefig(OUT + r"\fig3_tdry_conv.png")
plt.close(fig)

print("err_dt      =", np.array2string(err_dt, precision=2))
print("err_N_fv    =", np.array2string(err_N_fv, precision=3))
print("err_N_sp    =", np.array2string(err_N_sp, precision=6))
print("ok")
