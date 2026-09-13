# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: plotstyle.py
# 作用  : 统一论文绘图风格(中文字体、字号、分辨率、色彩), 避免各脚本重复设置。
#         所有出图脚本都必须先调用 setup()。
# 说明  : 中文字体按候选列表探测注册, 保证在离线环境下也能显示中文;
#         Unicode 负号统一关闭, 避免坐标轴出现方框。
# =============================================================================
"""Shared matplotlib style for the paper figures."""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")                      # 无界面后端, 便于批处理
import matplotlib.pyplot as plt            # noqa: E402
from matplotlib import font_manager        # noqa: E402

FONT_CANDIDATES = [
    r"C:\Windows\Fonts\msyh.ttc",          # 微软雅黑
    r"C:\Windows\Fonts\simhei.ttf",        # 黑体
    r"C:\Windows\Fonts\simsun.ttc",        # 宋体
]

PALETTE = {
    "blue": "#1f4e79",
    "red": "#c0392b",
    "green": "#1e8449",
    "orange": "#d68910",
    "purple": "#6c3483",
    "gray": "#566573",
    "cyan": "#117a8b",
}

_CONFIGURED = False


def setup():
    """注册中文字体并设置全局绘图参数(幂等)。"""
    global _CONFIGURED
    if _CONFIGURED:
        return plt.rcParams["font.family"]
    chosen = None
    for path in FONT_CANDIDATES:
        if os.path.exists(path):
            try:
                font_manager.fontManager.addfont(path)
                chosen = font_manager.FontProperties(fname=path).get_name()
                break
            except Exception:
                continue
    if chosen is None:
        chosen = "DejaVu Sans"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [chosen, "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["mathtext.fontset"] = "dejavusans"
    plt.rcParams["figure.dpi"] = 120
    plt.rcParams["savefig.dpi"] = 220
    plt.rcParams["font.size"] = 10.5
    plt.rcParams["axes.titlesize"] = 11.5
    plt.rcParams["axes.labelsize"] = 10.5
    plt.rcParams["legend.fontsize"] = 9.5
    plt.rcParams["xtick.labelsize"] = 9.5
    plt.rcParams["ytick.labelsize"] = 9.5
    plt.rcParams["axes.grid"] = True
    plt.rcParams["grid.alpha"] = 0.28
    plt.rcParams["grid.linestyle"] = "--"
    plt.rcParams["grid.linewidth"] = 0.6
    plt.rcParams["axes.linewidth"] = 0.9
    plt.rcParams["legend.frameon"] = True
    plt.rcParams["legend.framealpha"] = 0.92
    plt.rcParams["savefig.bbox"] = "tight"
    plt.rcParams["figure.autolayout"] = False
    _CONFIGURED = True
    return chosen


def save(fig, path, dpi=220):
    """保存图片并关闭画布, 保证不残留句柄。"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path
