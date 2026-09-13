# -*- coding: utf-8 -*-
"""控制台 UTF-8 输出修复（Windows 默认 GBK 会把中文打印成乱码）。"""
import sys


def fix_console():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
