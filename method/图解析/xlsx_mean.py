# -*- coding: utf-8 -*-
"""直接从交付文件 result3.xlsx 复算截面平均含水率的达标时刻。

文件里第一行是到药材中心的距离（厘米），A 列是时间（秒），
每一行是一个时刻的径向剖面。面积加权平均的定义是
    Cbar = (2/R^2) * 积分 C(r) r dr，从 0 到 R。
"""
import numpy as np
import openpyxl

XL = r"C:\Users\qing1\Desktop\bestway\药材烘干问题_正则化谱方法\结果\result3.xlsx"
wb = openpyxl.load_workbook(XL, data_only=True)
ws = wb.active

rows = list(ws.iter_rows(values_only=True))
head = rows[0]
r_cm = np.array([float(v) for v in head[1:]])
r = r_cm * 1e-2
R = r[-1]

t = []
prof = []
for row in rows[1:]:
    if row[0] is None:
        continue
    vals = [float(v) for v in row[1:]]
    if any(v is None for v in vals):
        continue
    t.append(float(row[0]))
    prof.append(vals)
t = np.array(t)
C = np.array(prof)

print("结果文件：%s" % XL)
print("数据行数 %d，径向列数 %d，半径从 %.2f 到 %.2f 厘米" %
      (len(t), len(r_cm), r_cm[0], r_cm[-1]))
print("时间从 %.0f 到 %.4f 秒" % (t[0], t[-1]))
print("")

# 面积加权平均：梯形积分
mean = np.array([np.trapezoid(C[i] * r, r) * 2.0 / R ** 2 for i in range(len(t))])
# 对照：把 21 个采样点直接算术平均
plain = C.mean(axis=1)

print("第一次低于 0.15 的时刻：")
idx = np.where(mean < 0.15)[0]
if idx.size:
    i = idx[0]
    print("  面积加权平均：%.4f 小时（第 %d 行，t=%.1f 秒）" % (t[i] / 3600.0, i, t[i]))
    if i > 0:
        print("  前一行 t=%.1f 秒 时平均值 %.6f" % (t[i - 1], mean[i - 1]))
idx2 = np.where(plain < 0.15)[0]
if idx2.size:
    print("  直接算术平均：%.4f 小时" % (t[idx2[0]] / 3600.0))

print("")
print("几个时刻的两种平均值：")
print("%-10s %16s %16s" % ("t/h", "面积加权", "算术平均"))
for tt in [6, 12, 24, 30, 36, 42]:
    i = int(np.argmin(np.abs(t - tt * 3600.0)))
    print("%-10d %16.6f %16.6f" % (tt, mean[i], plain[i]))

print("")
print("终点行：t=%.4f 秒，面积加权平均 %.6f，中心 %.6f，表面 %.6f" %
      (t[-1], mean[-1], C[-1, 0], C[-1, -1]))
