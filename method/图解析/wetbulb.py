# -*- coding: utf-8 -*-
"""计算烘房空气状态下的湿球温度，用来判断蒸发冷却能把表面降到多低。

湿球温度是蒸发冷却的物理下限：在恒定干燥条件下，
湿物料表面的温度会趋向湿球温度，不会低于它。

用 ASHRAE 的湿球方程（温度用摄氏度，含湿量 w 用千克水每千克干空气）：
    w = [ (2501 - 2.326 Tw) * ws(Tw) - 1.006 (T - Tw) ]
        / [ 2501 + 1.86 T - 4.186 Tw ]
其中 ws 是饱和含湿量，用 Magnus 公式算饱和蒸汽压：
    psat(T) = 0.61094 * exp(17.625 T / (T + 243.04))   kPa
    ws(T)   = 0.622 * psat / (101.325 - psat)
"""
import numpy as np
from scipy.optimize import brentq


def psat(T):
    return 0.61094 * np.exp(17.625 * T / (T + 243.04))


def ws(T):
    p = psat(T)
    return 0.622 * p / (101.325 - p)


def wetbulb(T, w):
    def f(Tw):
        return ((2501.0 - 2.326 * Tw) * ws(Tw) - 1.006 * (T - Tw)) \
            / (2501.0 + 1.86 * T - 4.186 * Tw) - w
    return brentq(f, -20.0, T - 1e-6, xtol=1e-10)


print("烘房空气状态下的湿球温度")
print("=" * 62)
print("%10s %12s %12s %12s" % ("空气温度", "含湿量", "饱和含湿量", "湿球温度"))
for T, w in [(50.0, 0.0509), (41.70, 0.03455), (35.0, 0.0250), (28.0, 0.01963)]:
    print("%10.2f %12.5f %12.5f %12.2f" % (T, w, ws(T), wetbulb(T, w)))

print("")
print("对照：叠加潜热以后算出的表面温度")
print("  600 秒 时表面 9.94 摄氏度")
print("  1800 秒 时表面 8.37 摄氏度")
print("  3600 秒 时表面 13.41 摄氏度")
print("")
print("结论：这些值远低于对应空气状态的湿球温度，物理上不可能。")
