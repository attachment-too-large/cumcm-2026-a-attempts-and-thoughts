# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s14_inverse.py
# 作用  : 只重跑 s05 的 E 组(合成数据回归), 用于在修正二次曲面拟合的数值
#         条件数之后刷新 data/inverse_*.csv, 而不重复 A~D 组。
# 运行  : python s14_inverse.py
# =============================================================================
"""Re-run only the synthetic-data inversion part with the conditioned fit."""

from __future__ import annotations

import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scenarios as sc        # noqa: E402
import s05_verification as v  # noqa: E402


def main():
    lines = []
    push = lines.append
    v.part_e_inverse(push, lines)
    text = "\n".join(lines)
    with io.open(os.path.join(sc.DIR_DATA, "inverse_summary.txt"), "w",
                 encoding="utf-8", newline="\n") as f:
        f.write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
