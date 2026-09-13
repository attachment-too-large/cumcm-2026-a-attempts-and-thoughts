# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s12_fix_outputs.py
# 作用  : 对已生成的 results/result*.xlsx 做两处"表述层"补正(不改动任何数值),
#         使其与附件 3 模板完全一致:
#           (1) result1~result3 末列表头: 模板为数值 2(即 2 cm, 固定半径下
#               即为表面), 早期版本误写成"药材表面"; result4 模板本身就是
#               "药材表面", 保持不变。
#           (2) result2 的时间列: result2 用 dt=0.1 s 推进, 长时间累加的浮点
#               误差使记录时刻出现 0.9999999999999999 这类值, 早期版本用
#               int() 截断, 把 1 s 写成了 0 s。这里把时间列严格写回
#               1,2,...,10800。
#         数值本身(温度、含水率)一律不修改, 补正后逐列校验。
# 说明  : 生成脚本 scenarios.write_result_xlsx 已同步修正(round 取整 +
#         surface_header 开关), 重新完整运行 s02/s03/s04 会直接得到同样结果;
#         本脚本仅用于免去重复数十小时的求解过程。
# 产出  : 覆盖 results/result1.xlsx ~ result3.xlsx, 并打印校验结果
# 运行  : python s12_fix_outputs.py
# =============================================================================
"""Patch only the label layer of the already-computed result workbooks."""

from __future__ import annotations

import os
import sys

import numpy as np
import openpyxl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scenarios as sc        # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")

# 文件 -> (工作表名, 时间列应为的序列)
PLAN = {
    "result1.xlsx": [("温度", np.arange(1, 1801)), ("水分浓度", np.arange(1, 1801))],
    "result2.xlsx": [("温度", np.arange(1, 10801)), ("水分浓度", np.arange(1, 10801))],
}


def main():
    report = []
    # ---- (1) 末列表头补正 ----
    for fname in ("result1.xlsx", "result2.xlsx", "result3.xlsx"):
        p = os.path.join(RES, fname)
        wb = openpyxl.load_workbook(p)
        for ws in wb.worksheets:
            last = ws.cell(row=1, column=ws.max_column)
            old = last.value
            if old != 2:
                last.value = 2
                report.append(f"{fname}/{ws.title}: 末列表头 '{old}' -> 2")
            # 距离轴列头校验: 0,0.1,...,倒数第二列
            head = [ws.cell(row=1, column=c).value for c in range(2, ws.max_column + 1)]
            expect = list(sc.positions_cm_fixed())[:-1] + [2]
            ok = all(abs(float(a) - float(b)) < 1e-9 for a, b in zip(head, expect))
            report.append(f"{fname}/{ws.title}: 距离轴列头共 {len(head)} 列, "
                          f"与 0~2 cm 每 0.1 cm 一致 = {ok}")
        wb.save(p)

    # ---- (2) result2 时间列补正 ----
    p = os.path.join(RES, "result2.xlsx")
    wb = openpyxl.load_workbook(p)
    for name, seq in PLAN["result2.xlsx"]:
        ws = wb[name]
        n = ws.max_row - 1
        assert n == len(seq), (name, n, len(seq))
        for i, t in enumerate(seq, start=2):
            ws.cell(row=i, column=1).value = int(t)
        report.append(f"result2.xlsx/{name}: 时间列已写回 1..{int(seq[-1])} (共 {n} 行)")
    wb.save(p)

    # ---- 校验 ----
    for fname in ("result1.xlsx", "result2.xlsx", "result3.xlsx", "result4.xlsx"):
        wb = openpyxl.load_workbook(os.path.join(RES, fname), read_only=True)
        for ws in wb.worksheets:
            t = [ws.cell(row=r, column=1).value for r in (2, 3, ws.max_row)]
            report.append(f"{fname}/{ws.title}: {ws.max_row-1} 行 x {ws.max_column} 列, "
                          f"首/次/末时间 = {t}, 末列表头 = "
                          f"{ws.cell(row=1, column=ws.max_column).value}")
        wb.close()

    text = "\n".join(report)
    with open(os.path.join(BASE, "data", "output_files_check.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
