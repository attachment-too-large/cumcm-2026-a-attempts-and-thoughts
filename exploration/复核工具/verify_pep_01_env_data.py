# -*- coding: utf-8 -*-
"""针对 paper_electronic.pdf（新一版）的独立复核 · 第 1 部分

核验项：
  A. 附件1 环境一阶惯性拟合（论文式(14)）—— 参数、R^2、RMSE，并检查 0-4 h 的系统性偏差
  B. 附件2 半径：表 6 各行取值是否与附件 2 一致
  C. 附录 3/4 物性式的量级复算（论文 5.6 节所引 D 值）
  D. 论文自报"截面平均"的口径定位（0.1232 / 2.2907）
"""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import openpyxl
from scipy.optimize import curve_fit

ROOT = r'C:\Users\qing1\Desktop\A题'
ATT = os.path.join(ROOT, '附件')

# ---------------- A 环境拟合 ----------------
wb = openpyxl.load_workbook(os.path.join(ATT, '附件1.xlsx'), data_only=True)
rows = [r for r in wb['Sheet1'].iter_rows(values_only=True)][1:]
t = np.array([float(r[0]) for r in rows])
Ta = np.array([float(r[1]) for r in rows])
Ca = np.array([float(r[2]) for r in rows])


def f1(tt, Tset, tau):
    return Tset - (Tset - 28.0) * np.exp(-tt / tau)


def g1(tt, Cset, tau):
    return Cset - (Cset - 0.01963) * np.exp(-tt / tau)


pT, _ = curve_fit(f1, t, Ta, p0=[50.2, 1877.0], maxfev=20000)
pC, _ = curve_fit(g1, t, Ca, p0=[0.0509, 2776.0], maxfev=20000)
rT, rC = f1(t, *pT), g1(t, *pC)
print('=' * 74)
print('A. 论文式(14) 环境一阶惯性拟合（初值固定为实测 28 °C / 0.01963）')
print('   温度: Tset = %.4f °C, tau_T = %.1f s, R^2 = %.5f, RMSE = %.4f °C'
      % (pT[0], pT[1], 1 - np.sum((Ta - rT) ** 2) / np.sum((Ta - Ta.mean()) ** 2),
         np.sqrt(np.mean((Ta - rT) ** 2))))
print('   论文: Tset=50.212, tau_T=1877.5, R^2=0.9957, RMSE=0.336')
print('   湿度: Cset = %.6f, tau_C = %.1f s, R^2 = %.4f, RMSE = %.2e'
      % (pC[0], pC[1], 1 - np.sum((Ca - rC) ** 2) / np.sum((Ca - Ca.mean()) ** 2),
         np.sqrt(np.mean((Ca - rC) ** 2))))
print('   论文: Cset=0.050908, tau_C=2776.3, R^2=0.9901, RMSE=8.1e-4')
print()
print('   0-4 h 区间拟合 vs 实测（检查是否系统性偏低）：')
print('   %8s %10s %10s %9s' % ('t/s', '实测', '拟合', '偏差'))
for tt in (60, 300, 600, 900, 1200, 1800, 2700, 3600, 7200, 10800, 14400):
    i = int(np.where(np.isclose(t, tt))[0][0])
    print('   %8d %10.4f %10.4f %+9.4f' % (tt, Ta[i], f1(tt, *pT), f1(tt, *pT) - Ta[i]))
seg = (t >= 3600) & (t <= 14400)
print('   1-4 h（3600-14400 s）拟合偏差: 均值 %+.4f K, 极差 [%+.4f, %+.4f] K'
      % (np.mean(f1(t[seg], *pT) - Ta[seg]), np.min(f1(t[seg], *pT) - Ta[seg]),
         np.max(f1(t[seg], *pT) - Ta[seg])))
print('   全程残差标准差 %.4f K（论文 RMSE 0.336）' % np.sqrt(np.mean((Ta - rT) ** 2)))

# ---------------- B 附件2 半径 ----------------
wb = openpyxl.load_workbook(os.path.join(ATT, '附件2.xlsx'), data_only=True)
rows2 = [r for r in wb['Sheet1'].iter_rows(values_only=True)][1:]
tR = np.array([float(r[0]) for r in rows2])
RR = np.array([float(r[1]) for r in rows2])
print()
print('=' * 74)
print('B. 论文表 6 的 R(t) 列 与 附件 2 对照')
paper_R = {6: 1.374, 12: 1.248, 18: 1.214, 24: 1.204, 30: 1.201, 36: 1.200, 42: 1.200, 48: 1.200}
for h, pr in paper_R.items():
    i = int(np.where(np.isclose(tR, h * 3600))[0][0])
    print('   %2d h: 附件2 = %.3f cm, 论文 = %.3f cm, 差 = %+.4f' % (h, RR[i], pr, pr - RR[i]))
print('   附件2: R(0)=%.3f -> R(end)=%.3f cm; 单调不增: %s' % (RR[0], RR[-1], bool(np.all(np.diff(RR) <= 0))))
i = int(np.argmax(RR <= 1.5))
x1, y1, x2, y2 = tR[i - 1], RR[i - 1], tR[i], RR[i]
print('   R(t)=1.5 cm 交点（分段线性）: t = %.0f s = %.3f h  [论文 9.2/9.4 节写 3.65 h]'
      % (x1 + (1.5 - y1) * (x2 - x1) / (y2 - y1), (x1 + (1.5 - y1) * (x2 - x1) / (y2 - y1)) / 3600))

# ---------------- C 物性式量级 ----------------
print()
print('=' * 74)
print('C. 物性式量级复算（论文 5.6 节引用的 D 值）')
def D3(C, Tc):
    return 2.4e-3 * np.exp(-0.45 / C) * np.exp(-3850.0 / (Tc + 273.15))
def D4(C, Tc):
    return 4.2e-4 * np.exp(-0.30 / C) * np.exp(-3850.0 / (Tc + 273.15))
print('   附录3  C=2.55, T=50 °C : D = %.4e   [论文 1.35e-8]' % D3(2.55, 50))
print('   附录3  C=2.55, T=48.4 °C: D = %.4e   （反推 1.35e-8 对应 48.40 °C）' % D3(2.55, 48.40))
print('   附录3  C=0.15, T=50 °C : D = %.4e   [论文 8.0e-10]' % D3(0.15, 50))
print('   附录3  C=0.15, T=48.4 °C: D = %.4e' % D3(0.15, 48.40))
print('   附录3  C=2.55, T=49.97 °C（3h 烘房温度）: D = %.4e' % D3(2.55, 49.97))
print('   比值 D(2.55)/D(0.15) @50°C = %.1f 倍  [论文"下降17倍"]' % (D3(2.55, 50) / D3(0.15, 50)))
print('   alpha = k/(rho*cp) = %.4e  [论文 1.6886e-7]' % (0.36 / (820 * 2600)))
print('   D(C0)=%.4e  [论文 4.9377e-9]' % D3(2.55, 28.0))
print('   Bi_m = hm*R0/D(C0) = %.4f  [论文 3.2404]' % (8e-7 * 0.02 / D3(2.55, 28.0)))
print('   Bi_m = hm*R0/D(2.55,50°C) = %.4f  [论文 1.19]' % (8e-7 * 0.02 / D3(2.55, 50)))
print('   Le = D(2.55,50)/alpha = %.4f  [论文 0.13]' % (D3(2.55, 50) / (0.36 / (820 * 2600))))
print('   Le = D(C0)/alpha        = %.4f' % (D3(2.55, 28.0) / (0.36 / (820 * 2600))))

# ---------------- D 截面平均口径 ----------------
print()
print('=' * 74)
print('D. 论文自报"截面平均"的口径定位')
print('   论文表5 末行 5 点值 (0/0.5/1/1.5/2 cm): 0.1500, 0.1477, 0.1402, 0.1249, 0.0536')
v5 = np.array([0.1500, 0.1477, 0.1402, 0.1249, 0.0536])
r5 = np.array([0.0, 0.5, 1.0, 1.5, 2.0]) / 100.0
print('   5 点【等权算术平均】  = %.4f      <- 论文报告 0.1232' % v5.mean())
print('   5 点【面积加权 2r/R^2】= %.4f' % np.average(v5, weights=2 * r5 / 0.02 ** 2))
print()
print('   论文表2 t=1800 s 5 点值: 2.5500, 2.5497, 2.5383, 2.3755, 1.5102')
v5b = np.array([2.5500, 2.5497, 2.5383, 2.3755, 1.5102])
print('   5 点【等权算术平均】  = %.4f      <- 论文报告 2.2907' % v5b.mean())
print('   5 点【面积加权】      = %.4f' % np.average(v5b, weights=2 * r5 / 0.02 ** 2))
# 干物质与失水
rho_d0 = 820 / (1 + 2.55)
V = np.pi * 0.02 ** 2 * 0.25
m_d = rho_d0 * V
print('   干物质 = %.4f g [论文 72.57]、初始水 = %.4f g [论文 185.04]' % (m_d * 1e3, m_d * 2.55 * 1e3))
print('   按等权 2.3047 算失水 = %.4f g   [论文 18.71]' % (m_d * (2.55 - v5b.mean()) * 1e3))
print('   按面积加权 2.2534 算失水 = 见第2部分（需完整剖面）')
