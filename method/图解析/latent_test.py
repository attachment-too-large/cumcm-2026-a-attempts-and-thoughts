# -*- coding: utf-8 -*-
"""第三问加与不加汽化潜热的对照实验。

两种加法：
  甲  表面蒸发冷却：蒸发消耗的潜热出现在温度的表面边界条件里。
      能量平衡是 由内部导到表面的热量 = 对流散失 + 蒸发带走的潜热，即
          -k dT/dr|R = h(Ts - Tair) + L * rho_d * hm * (Cs - Cair)
      换成 u 坐标后，表面通量变成
          g_N = 0.5 * [ h R (Tair - T1) - L rho_d hm R (C1 - Cair) ]
  乙  体积潜热汇：把蒸发当成均匀分布在体内的吸热项，即
          rho cp dT/dt = ... + L * rho_d * dC/dt
      由于 dC/dt 是负的，这一项使温度下降。

对照基准是不加潜热的现行模型。
潜热取值：50 摄氏度下水的汽化潜热约 2.383e6 J/kg。
"""
import sys
import time
import numpy as np

CODEDIR = r"C:\Users\qing1\Desktop\bestway\药材烘干问题_正则化谱方法\代码"
sys.path.insert(0, CODEDIR)
import common as cm
from spectral import HerbSolver

CRIT = 0.15
LAT = 2.383e6          # J/kg，50 摄氏度附近水的汽化潜热


class LatentSurface(HerbSolver):
    """表面蒸发冷却：潜热进入温度边界条件。"""

    L = LAT

    def rhs(self, t, y):
        N, Np = self.N, self.Np
        C = y[:Np]
        T = y[Np:]
        R = self.R_of(t)
        rho, cp, k, D = cm.props(self.prob, C, T)

        gC = self.u * D * (self.Du @ C)
        Cair = self.C_air(t)
        gC[N] = 0.5 * cm.hm * R * (Cair - C[N])
        dC = (4.0 / R ** 2) * (self.Du @ gC)

        gT = self.u * k * (self.Du @ T)
        rho_d = rho / (1.0 + C)                      # 干物质密度
        # 水离开表面的质量通量 = rho_d * hm * (C1 - Cair)
        latent_flux = self.L * rho_d[N] * cm.hm * (C[N] - Cair)
        gT[N] = 0.5 * (cm.h * R * (self.T_air(t) - T[N]) - R * latent_flux)
        dT = (4.0 / (R ** 2 * rho * cp)) * (self.Du @ gT)

        out = np.empty_like(y)
        out[:Np] = dC
        out[Np:] = dT
        return out


class LatentVolumetric(HerbSolver):
    """体积潜热汇：把蒸发当作体内的吸热项。"""

    L = LAT

    def rhs(self, t, y):
        N, Np = self.N, self.Np
        C = y[:Np]
        T = y[Np:]
        R = self.R_of(t)
        rho, cp, k, D = cm.props(self.prob, C, T)

        gC = self.u * D * (self.Du @ C)
        gC[N] = 0.5 * cm.hm * R * (self.C_air(t) - C[N])
        dC = (4.0 / R ** 2) * (self.Du @ gC)

        gT = self.u * k * (self.Du @ T)
        gT[N] = 0.5 * cm.h * R * (self.T_air(t) - T[N])
        dT = (4.0 / (R ** 2 * rho * cp)) * (self.Du @ gT)
        # 体积潜热汇：rho cp dT/dt 里多出 L * rho_d * dC/dt 这一项
        dT = dT + self.L * (rho / (1.0 + C)) * dC / (rho * cp)

        out = np.empty_like(y)
        out[:Np] = dC
        out[Np:] = dT
        return out


def dry(cls, prob=3, shrink=False, rtol=1e-10):
    s = cls(64, prob=prob, shrink=shrink)
    at = np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)])

    def ev(t, y):
        return np.max(y[:s.Np]) - CRIT
    ev.terminal = True
    ev.direction = -1
    sol = s.solve(300000.0, rtol=rtol, atol=at, events=ev)
    if not sol.t_events[0].size:
        return s, sol, float("nan")
    return s, sol, float(sol.t_events[0][0])


print("=" * 74)
print("第三问：加与不加汽化潜热的对照")
print("潜热取 %.3e J/kg" % LAT)
print("=" * 74)

t0 = time.time()
s0, sol0, td0 = dry(HerbSolver)
print("甲 基准，不加潜热            tdry = %10.3f 秒 = %7.4f 小时" % (td0, td0 / 3600))

s1, sol1, td1 = dry(LatentSurface)
print("乙 表面蒸发冷却              tdry = %10.3f 秒 = %7.4f 小时" % (td1, td1 / 3600))

s2, sol2, td2 = dry(LatentVolumetric)
print("丙 体积潜热汇                tdry = %10.3f 秒 = %7.4f 小时" % (td2, td2 / 3600))

print("")
print("相对基准的变化：")
print("  表面蒸发冷却 %+8.3f%%   （%+.1f 小时）" % ((td1 - td0) / td0 * 100, (td1 - td0) / 3600))
print("  体积潜热汇   %+8.3f%%   （%+.1f 小时）" % ((td2 - td0) / td0 * 100, (td2 - td0) / 3600))
print("总用时 %.1f 秒" % (time.time() - t0))

print("")
print("=" * 74)
print("表面温度的对照（摄氏度）")
print("=" * 74)
for name, s, sol in [("不加潜热", s0, sol0), ("表面蒸发冷却", s1, sol1)]:
    tt = np.array([600.0, 1800.0, 3600.0, 7200.0, 21600.0])
    tt = tt[tt <= sol.t[-1]]
    e = s.solve(tt[-1], t_eval=list(tt), rtol=1e-10,
                atol=np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)]))
    Ts = [float(s.interp_u(e.y[:, i], np.array([1.0]))[1][0]) for i in range(len(tt))]
    Tc = [float(s.interp_u(e.y[:, i], np.array([0.0]))[1][0]) for i in range(len(tt))]
    print("%-14s" % name, end="")
    for v in Ts:
        print(" %8.3f" % v, end="")
    print("   （表面）")
    print("%-14s" % "", end="")
    for v in Tc:
        print(" %8.3f" % v, end="")
    print("   （中心）")
print("%-14s" % "时刻 / 秒", end="")
for t in [600.0, 1800.0, 3600.0, 7200.0, 21600.0]:
    print(" %8.0f" % t, end="")
print("")

print("")
print("=" * 74)
print("两项热流的量级对照（早期，单位 W/m2）")
print("=" * 74)
s = HerbSolver(64, prob=3)
e = s.solve(1800.0, t_eval=[1800.0], rtol=1e-11,
            atol=np.concatenate([np.full(s.Np, 1e-13), np.full(s.Np, 1e-11)]))
C = e.y[:s.Np, -1]
T = e.y[s.Np:, -1]
rho, cp, k, D = cm.props(3, C, T)
Tair = float(s.T_air(1800.0))
Cair = float(s.C_air(1800.0))
rho_d = rho / (1.0 + C)
conv = cm.h * (Tair - T[-1])
lat = LAT * rho_d[-1] * cm.hm * (C[-1] - Cair)
print("  一千八百秒时：对流给热 %8.1f，蒸发耗热 %8.1f，两者之比 %.2f" %
      (conv, lat, lat / conv if conv else float("nan")))
print("  同一时刻表面含水率 %.4f，空气含湿量 %.5f" % (C[-1], Cair))
