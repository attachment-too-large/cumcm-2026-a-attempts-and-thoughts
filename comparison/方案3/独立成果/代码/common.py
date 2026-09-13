# -*- coding: utf-8 -*-
"""A题 新解法 —— 公共参数、烘房环境、物性经验式。

所有数值均由题目/附件直接给出，或由附件数据独立标定（一阶惯性拟合）得到。
"""
import os
import numpy as np
import openpyxl
from scipy.optimize import curve_fit
from scipy.special import j0, j1
from scipy.optimize import brentq

# ---------------- 数据文件定位（可移植：相对本文件搜索） ----------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_CANDIDATES = [
    os.path.join(_HERE, "..", "题目"),
    os.path.join(_HERE, ".."),
    os.path.join(_HERE, "..", ".."),
    _HERE,
]


def _find_attach(name):
    """在若干候选目录下查找 附件/<name>。"""
    for base in _CANDIDATES:
        p = os.path.abspath(os.path.join(base, "附件", name))
        if os.path.isfile(p):
            return p
    raise FileNotFoundError(
        "找不到 附件/%s；已搜索：%s" % (name, [os.path.abspath(b) for b in _CANDIDATES]))


# ---------------- 几何与传递系数（附录2） ----------------
R0 = 0.02         # m   药材半径
Lc = 0.25         # m   药材长度
h  = 25.0         # W/(m^2*K)   对流换热系数
hm = 8.0e-7       # m/s         对流传质系数

# ---------------- 初值 ----------------
T0 = 28.0         # degC
C0 = 2.55         # kg/kg (干基)

# ---------------- 附件1：烘房环境，一阶惯性拟合 ----------------
def _load_env():
    ws = openpyxl.load_workbook(_find_attach("附件1.xlsx"), data_only=True).active
    rows = list(ws.iter_rows(min_row=2, values_only=True))
    t = np.array([float(r[0]) for r in rows])
    Ta = np.array([float(r[1]) for r in rows])
    Ca = np.array([float(r[2]) for r in rows])
    return t, Ta, Ca


_t_env, _Ta_env, _Ca_env = _load_env()
T_air0 = _Ta_env[0]          # 28.0
C_air0 = _Ca_env[0]          # 0.01963
T_END = float(_t_env[-1])    # 14400 s


def _fT(t, Tset, tau):
    return Tset - (Tset - T_air0) * np.exp(-t / tau)


def _fC(t, Cset, tau):
    return Cset - (Cset - C_air0) * np.exp(-t / tau)


(Tset_fit, tauT_fit), _ = curve_fit(_fT, _t_env, _Ta_env, p0=[50.2, 1800.0])
(Cset_fit, tauC_fit), _ = curve_fit(_fC, _t_env, _Ca_env, p0=[0.0509, 2800.0])

# 标定参数的报告精度：把拟合值按其有效位数保留（T_set 到 0.001 degC、
# tau 到 0.1 s）。该取舍在 T(R,1800s) 上仅造成约 6e-5 degC 的变化，
# 远小于拟合自身的标准误（T_set: 0.030 degC、tau_T: 12.7 s）。
Tset = 50.212
tauT = 1877.5
Cset = 0.050908
tauC = 2776.3
VARIANTS = {
    "reported_precision": dict(Tset=Tset, tauT=tauT, Cset=Cset, tauC=tauC),
    "full_precision_fit": dict(Tset=Tset_fit, tauT=tauT_fit,
                               Cset=Cset_fit, tauC=tauC_fit),
}
# 协方差（用于报告标准误）
_nT, _nC = len(_t_env), len(_t_env)
_sT = np.sqrt(np.sum((_Ta_env - _fT(_t_env, Tset, tauT)) ** 2) / (_nT - 2))
_sC = np.sqrt(np.sum((_Ca_env - _fC(_t_env, Cset, tauC)) ** 2) / (_nC - 2))


# 数值 Jacobian（更稳妥）
def _jac(f, p, t):
    J = np.zeros((len(t), len(p)))
    for i in range(len(p)):
        dp = p.copy(); dh = abs(p[i]) * 1e-7 + 1e-12; dp[i] += dh
        J[:, i] = (f(t, *dp) - f(t, *p)) / dh
    return J
_seT = np.sqrt(np.diag(_sT ** 2 * np.linalg.inv(_jac(_fT, np.array([Tset, tauT]), _t_env).T @
                                                 _jac(_fT, np.array([Tset, tauT]), _t_env))))
_seC = np.sqrt(np.diag(_sC ** 2 * np.linalg.inv(_jac(_fC, np.array([Cset, tauC]), _t_env).T @
                                                 _jac(_fC, np.array([Cset, tauC]), _t_env))))
RMSE_T = _sT
RMSE_C = _sC


def T_air(t):
    """烘房温度：预热平衡段取一阶惯性式，恒温干燥段取设定值。"""
    t = np.asarray(t, dtype=float)
    ramp = Tset - (Tset - T_air0) * np.exp(-t / tauT)
    return np.where(t <= T_END, ramp, Tset)


def C_air(t):
    """烘房水分浓度：同上。"""
    t = np.asarray(t, dtype=float)
    ramp = Cset - (Cset - C_air0) * np.exp(-t / tauC)
    return np.where(t <= T_END, ramp, Cset)


# ---------------- 物性经验式 ----------------
C_FLOOR = 1.0e-3


def props(prob, C, T):
    """返回 (rho, cp, k, D)。

    prob=1 用附录2（常数物性）；prob=2,3 用附录3；prob=4 用附录4。
    C: kg/kg；T: degC（内部转 K 供 Arrhenius 项使用）。
    """
    C = np.asarray(C, dtype=float)
    T = np.asarray(T, dtype=float)
    Cf = np.maximum(C, C_FLOOR)
    TK = T + 273.15
    if prob == 1:
        one = np.ones_like(C)
        rho = 820.0 * one
        cp = 2600.0 * one
        k = 0.36 * one
        D = 7.0e-9 * np.exp(-0.89 / Cf)
    elif prob in (2, 3):
        rho = 650.0 + 128.0 * C
        cp = 1450.0 + 2736.0 * C / (C + 1.0)
        k = 0.21 + 0.38 * C / (C + 1.0)
        D = 2.4e-3 * np.exp(-0.45 / Cf) * np.exp(-3850.0 / TK)
    elif prob == 4:
        rho = 760.0 + 90.0 * C
        cp = 1850.0 + 2150.0 * C / (C + 1.0)
        k = 0.12 + 0.20 * C / (C + 1.0)
        D = 4.2e-4 * np.exp(-0.30 / Cf) * np.exp(-3850.0 / TK)
    else:
        raise ValueError("prob must be 1,2,3,4")
    return rho, cp, k, D


# ---------------- 附件2：半径历史 ----------------
def _load_radius():
    ws = openpyxl.load_workbook(_find_attach("附件2.xlsx"), data_only=True).active
    rows = list(ws.iter_rows(min_row=2, values_only=True))
    t = np.array([float(r[0]) for r in rows])
    R = np.array([float(r[1]) for r in rows]) * 1e-2   # cm -> m
    return t, R


_t_R, _R_hist = _load_radius()


def R_of(t):
    """问题4：半径分段线性插值（t 超出范围则取端点值）。"""
    return np.interp(t, _t_R, _R_hist)


# ---------------- Bessel 特征值 ----------------
def bessel_roots(Bi, n, beta_min=0.0):
    """求 beta*J1(beta) = Bi*J0(beta) 的前 n 个正根（升序）。"""
    from scipy.special import jnp_zeros
    z = jnp_zeros(1, n + 2)                    # J1 的正零点
    br = [0.0] + list(z)                       # 括号端点
    roots = []
    f = lambda b: b * j1(b) - Bi * j0(b)
    # 第一个根在 (0, z0)
    lo = 1e-12
    hi = z[0] - 1e-12
    if f(lo) * f(hi) < 0:
        roots.append(brentq(f, lo, hi, xtol=1e-15, rtol=8.9e-16))
    for i in range(len(z) - 1):
        lo, hi = z[i] + 1e-12, z[i + 1] - 1e-12
        if f(lo) * f(hi) < 0:
            roots.append(brentq(f, lo, hi, xtol=1e-15, rtol=8.9e-16))
        if len(roots) >= n:
            break
    return np.array(roots[:n])
