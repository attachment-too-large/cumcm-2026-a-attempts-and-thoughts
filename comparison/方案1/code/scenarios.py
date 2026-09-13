# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: scenarios.py
# 作用  : 统一管理四个问题的"场景"配置与结果落盘格式, 保证各求解脚本使用同一
#         套环境/半径/网格/输出定义, 避免不同脚本各写一份导致口径不一致。
#         内容:
#           (1) 输入读取(附件 1、附件 2)与环境/半径构造函数;
#           (2) 半径插值方式(线性 / PCHIP 保形三次), 用于灵敏度对比;
#           (3) 输出位置定义(0~2 cm 每隔 0.1 cm; 或问题 4 的物理位置);
#           (4) Excel 结果文件模板写入函数(四位小数)。
# 单位  : 与 drying_core 一致(SI), 表头距离用 cm。
# =============================================================================
"""Scenario configuration and result-file writers for the four sub-problems."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

import drying_core as dc

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_IN = os.path.join(BASE, "data", "input")
DIR_DATA = os.path.join(BASE, "data")
DIR_FIG = os.path.join(BASE, "figures")
DIR_RES = os.path.join(BASE, "results")

HEADER_LABEL = "时间\\到药材中心的距离"
SURFACE_LABEL = "药材表面"


def load_env(plateau_rule="tail_mean", plateau_scale_C=1.0, plateau_scale_T=0.0):
    """读取附件 1 并构造环境函数。"""
    df = pd.read_excel(os.path.join(DIR_IN, "attachment1_room.xlsx"))
    df.columns = ["t", "T", "C"]
    return dc.RoomConditions(df["t"].to_numpy(float), df["T"].to_numpy(float),
                             df["C"].to_numpy(float), plateau_rule=plateau_rule,
                             plateau_scale_C=plateau_scale_C,
                             plateau_scale_T=plateau_scale_T)


def load_env_table():
    """返回附件 1 原始表(t, T_air, C_air)。"""
    df = pd.read_excel(os.path.join(DIR_IN, "attachment1_room.xlsx"))
    df.columns = ["t", "T", "C"]
    return df


def load_radius(interp="linear", t_end=None):
    """构造半径函数; interp 为 'linear' 或 'pchip'。

    'linear' 与附件 2 的离散数据严格一致; 'pchip' 为保形三次插值, 用于检验
    半径插值方式对问题 4 结果的影响。
    """
    df = pd.read_excel(os.path.join(DIR_IN, "attachment2_radius.xlsx"))
    df.columns = ["t", "R"]
    t = df["t"].to_numpy(float)
    R = df["R"].to_numpy(float) * 1e-2
    if interp == "linear":
        return dc.RadiusHistory(t, df["R"].to_numpy(float))
    if interp == "pchip":
        f = PchipInterpolator(t, R, extrapolate=False)
        t_last, R_last = float(t[-1]), float(R[-1])

        def radius(t_query):
            if t_query >= t_last:
                return R_last
            return float(f(t_query))

        return radius
    raise ValueError(f"unknown radius interpolation: {interp}")


def fixed_radius():
    """问题 1-3 的常数半径函数。"""
    return dc.constant_radius(dc.R0)


def positions_cm_fixed():
    """问题 1-3 的输出位置(cm): 0, 0.1, ..., 2.0。"""
    return np.round(np.arange(0.0, 2.0 + 1e-9, 0.1), 4)


def positions_xi_fixed():
    """问题 1-3 输出位置对应的归一化坐标 xi = r/R0。"""
    return positions_cm_fixed() / 2.0


def write_result_xlsx(path, times, pos_cm, values, sheet_names,
                      header_label=HEADER_LABEL, surface_last=False,
                      surface_header=False):
    """写出与附件 3 模板同构的结果文件。

    参数
    ----
    times      : (nt,) 时间列 [s]
    pos_cm     : (m,) 距离列 [cm]; 若 surface_last=True, 最后一个位置为药材表面
    values     : (nt, m) 数据(温度或含水率), 已按四位小数取整
    sheet_names: 工作表名列表; 单元素时为单表(问题 3、4), 双元素时两个表
    surface_last   : 最后一列的取值是否为半单元重构的表面值
    surface_header : 最后一列的表头是否写成"药材表面"。
                     附件 3 模板中 result1~result3 的末列表头是数值 2(即 2 cm,
                     对固定半径即为表面), 只有 result4 的模板写"药材表面",
                     因此本参数默认 False, 仅问题 4 传 True。

    时间列用 round 而不是 int 取整: 长时间累加的浮点误差会让 t=1.0 变成
    0.9999999999999999, 直接 int() 截断会把 1 s 标成 0 s。
    """
    from openpyxl import Workbook

    wb = Workbook()
    wb.remove(wb.active)
    header = [header_label] + [int(v) if float(v).is_integer() else float(v)
                               for v in pos_cm]
    if surface_header:
        header[-1] = SURFACE_LABEL
    for name, arr in zip(sheet_names, values):
        ws = wb.create_sheet(title=name)
        ws.append(header)
        for i, t in enumerate(times):
            ws.append([int(round(t))] + [float(f"{v:.4f}") for v in arr[i]])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    wb.save(path)
    return path


def read_position_series(res, pos_cm, surface_last=False):
    """把结果换算到给定物理位置, 返回 (nt, m) 的含水率与温度序列。

    两种记录模式:
      - 采样模式(res.sample_mode=True): 求解时已按固定归一化坐标记录, 直接返回;
      - 全场模式: 逐时刻把单元中心场插值到目标物理半径(问题 4 的 R 随时间变化)。
    若 surface_last=True, 最后一列由半单元重构的表面值填充。
    """
    m = len(pos_cm)
    if getattr(res, "sample_mode", False):
        C_out = np.asarray(res.C, dtype=float).copy()
        T_out = np.asarray(res.T, dtype=float).copy()
        xi_expect = np.asarray(res.record_xi, dtype=float)
        xi_got = np.asarray(pos_cm, dtype=float) * 1e-2 / res.R[0]
        if xi_expect.size != m or not np.allclose(xi_expect, xi_got, atol=1e-12):
            raise ValueError("采样位置与请求位置不一致, 请检查 record_xi")
        if surface_last:
            C_out[:, m - 1] = res.c_surf
            T_out[:, m - 1] = res.t_surf
        return C_out, T_out

    nt = res.C.shape[0]
    C_out = np.full((nt, m), np.nan)
    T_out = np.full((nt, m), np.nan)
    for i in range(nt):
        R = res.R[i]
        tgt = []
        for j, r_cm in enumerate(pos_cm):
            if surface_last and j == m - 1:
                tgt.append(None)          # 表面单独处理
            else:
                tgt.append(min(r_cm * 1e-2 / R, 1.0))
        valid = [k for k in range(m) if tgt[k] is not None]
        if valid:
            idx, w = dc.lagrange_weights_grid(res.n_cells, [tgt[k] for k in valid])
            C_out[i, valid] = dc.apply_weights(res.C[i], idx, w)
            T_out[i, valid] = dc.apply_weights(res.T[i], idx, w)
        if surface_last:
            C_out[i, m - 1] = res.c_surf[i]
            T_out[i, m - 1] = res.t_surf[i]
    return C_out, T_out


def drying_time(res, crit=dc.C_DRY):
    """由 max C 时间序列求干燥时间 [s]; 结果已在求解时插值给出。"""
    if res.t_dry is not None:
        return float(res.t_dry)
    mc = res.step_maxC - crit
    sign = np.sign(mc)
    k = np.where(sign < 0)[0]
    if k.size == 0:
        return float("nan")
    i = int(k[0])
    t0, t1 = res.step_times[i - 1], res.step_times[i]
    m0, m1 = res.step_maxC[i - 1], res.step_maxC[i]
    return float(t0 + (m0 - crit) / (m0 - m1) * (t1 - t0))


def step_dt_schedule(t_switch=86400.0, dt_early=10.0, dt_late=2.0):
    """分段固定时间步(用于长时间问题, 早期步长大、后期步长小)。

    每段步长都能整除 60 s, 保证 60 s 输出点落在整步上。
    """
    return [(t_switch, dt_early), (float("inf"), dt_late)]
