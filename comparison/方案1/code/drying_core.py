# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: drying_core.py
# 作用  : CUMCM 2026 A 题「药材的烘干问题」的核心物理与数值模块。
#         本模块实现「新模型」的全部要素:
#           (1) 圆柱一维轴对称耦合传热-传质控制方程组;
#           (2) 节点中心有限体积(FVM)离散 + 半单元表面重构(Robin 串联阻力);
#           (3) theta 时间推进(默认 theta=1, 即隐式 Euler) + Picard 迭代,
#               并强制右端项固定在使用旧时间层, 避免"一个时间步内多次推进";
#           (4) 问题 4 的物质坐标收缩域形式 (xi = r/R(t));
#           (5) 干燥终点判据 max C < 0.15 的守恒插值求交。
#
# 单位约定(全部 SI, 仅温度对用户显示用摄氏度):
#   长度 m; 时间 s; 温度 degC(物性公式内部自动换 K);
#   干基含水率 C kg/kg; 密度 kg/m3; 比热 J/(kg K); 导热 k W/(m K);
#   水分扩散系数 D m2/s; 对流换热 h W/(m2 K); 对流传质 h_m m/s.
#
# 关键变量含义:
#   n_cells       径向控制体(单元)数目 N
#   dxi           归一化网格步长, dxi = 1/N
#   xi_c          单元中心归一化坐标 xi_c[i] = (i+0.5)*dxi
#   xi_f          单元界面归一化坐标 xi_f[f] = f*dxi, f=0..N (xi_f[N]=1 为表面)
#   vol           单元归一化体积因子 vol[i] = (xi_f[i+1]^2-xi_f[i]^2)/2
#   g_up          界面导流系数(含面积因子与 1/dxi), 用于三对角装配
#   g_surf        表面串联阻力等效导流系数 g_surf = R*G_eff
#   c_surf        由半单元重构得到的真实表面含水率(不是最后一个单元中心值)
#
# 控制方程(0<=r<=R(t), 一维径向, 轴对称):
#   水分: dC/dt        = (1/r) d/dr ( r D(C,T) dC/dr )
#   温度: rho(C)cp(C) dT/dt = (1/r) d/dr ( r k(C) dT/dr )
#   边界: r=0 对称;  r=R(t): -k dT/dr = h (T_air - Ts),
#                            -D dC/dr = h_m (Cs - C_air)
#   初始: T(r,0)=28 degC, C(r,0)=2.55 kg/kg
#
# 归一化坐标 xi=r/R(t) 下(物质坐标, 材料点 xi 恒定, 纯径向仿射收缩):
#   dC/dt|_xi = (1/R^2) * (1/xi) d/dxi ( xi D dC/dxi )
#   rho cp dT/dt|_xi = (1/R^2) * (1/xi) d/dxi ( xi k dT/dxi )
# 时间导数取固定材料点, 因此不出现"收缩速度 x 梯度"的附加对流项。
# =============================================================================
"""Core physics and finite-volume solver for CUMCM 2026 Problem A.

The module is intentionally self-contained (numpy + scipy only) so that the
same code can be re-run by a reviewer without any project context.
"""

from __future__ import annotations

import time as _time

import numpy as np
from scipy.linalg import solve_banded

# ---------------------------------------------------------------------------
# 题目给定常数(附录 2 与题面)
# ---------------------------------------------------------------------------
R0 = 0.02                 # 药材初始半径 [m] (题面: 2 cm)
L_HERB = 0.25             # 药材长度 [m] (题面: 25 cm, 用于说明端面效应可忽略)
T_INIT_C = 28.0           # 药材初始温度 [degC]
C_INIT = 2.55             # 药材初始干基含水率 [kg/kg]
H_CONV = 25.0             # 对流换热系数 [W/(m2 K)]
HM_CONV = 8.0e-7          # 对流传质系数 [m/s]
KELVIN = 273.15           # 摄氏度到开尔文的换算
C_DRY = 0.15              # 干燥终点判据 [kg/kg]
T_PLATEAU_START = 14400.0  # 附件 1 数据结束时刻 [s]

# exp(-a/C) 在 C->0 时趋于 0；数值上给 C 一个下限保护，避免指数下溢。
# 该下限只影响已经远低于终点判据的极表层单元，对干燥时间影响需在灵敏度分析中量化。
C_FLOOR = 1.0e-4


# ---------------------------------------------------------------------------
# 物性关系(题目附录 2 / 附录 3 / 附录 4)
#   统一签名: props(C, T_C) -> (rho, cp, k, D)
#   T_C 为摄氏温度; 附录 3/4 的 Arrhenius 因子 exp(-3850/T) 中 T 必须用开尔文。
# ---------------------------------------------------------------------------
def props_appendix2(C, T_C):
    """问题 1 物性: 常数密度/比热/导热, 扩散系数仅依赖含水率。"""
    C = np.asarray(C, dtype=float)
    rho = np.full_like(C, 820.0)
    cp = np.full_like(C, 2600.0)
    k = np.full_like(C, 0.36)
    D = 7.0e-9 * np.exp(-0.89 / np.maximum(C, C_FLOOR))
    return rho, cp, k, D


def props_appendix3(C, T_C):
    """问题 2、3 物性: 附录 3 经验公式(密度、比热、导热、扩散系数)。"""
    C = np.asarray(C, dtype=float)
    T_K = np.asarray(T_C, dtype=float) + KELVIN
    rho = 650.0 + 128.0 * C
    cp = 1450.0 + 2736.0 * C / (C + 1.0)
    k = 0.21 + 0.38 * C / (C + 1.0)
    D = 2.4e-3 * np.exp(-0.45 / np.maximum(C, C_FLOOR)) * np.exp(-3850.0 / T_K)
    return rho, cp, k, D


def props_appendix4(C, T_C):
    """问题 4 物性: 附录 4 经验公式。"""
    C = np.asarray(C, dtype=float)
    T_K = np.asarray(T_C, dtype=float) + KELVIN
    rho = 760.0 + 90.0 * C
    cp = 1850.0 + 2150.0 * C / (C + 1.0)
    k = 0.12 + 0.20 * C / (C + 1.0)
    D = 4.2e-4 * np.exp(-0.30 / np.maximum(C, C_FLOOR)) * np.exp(-3850.0 / T_K)
    return rho, cp, k, D


PROPS = {
    "appendix2": props_appendix2,
    "appendix3": props_appendix3,
    "appendix4": props_appendix4,
}


# ---------------------------------------------------------------------------
# 烘房环境(附件 1)与药材半径历史(附件 2)
# ---------------------------------------------------------------------------
class RoomConditions:
    """烘房温湿度环境。

    附件 1 只覆盖 0..14400 s(间隔 60 s)。基准模型对 t 落在数据区间内使用线性
    插值; 对 t 超出末点使用"平台假定", 即恒温干燥阶段环境保持常数。
    plateau_rule 决定平台取值方式:
      'tail_mean'  : 末段 3600 s(最后 60 个点)的算术平均(基准)
      'last_point' : 直接取最后一个数据点
      'tail_mean_1800' / 'tail_mean_7200': 其它平均窗长, 用于灵敏度分析
    """

    def __init__(self, t_arr, T_arr, C_arr, plateau_rule="tail_mean",
                 plateau_scale_C=1.0, plateau_scale_T=0.0):
        self.t = np.asarray(t_arr, dtype=float)
        self.T = np.asarray(T_arr, dtype=float)
        self.C = np.asarray(C_arr, dtype=float)
        self.t_last = float(self.t[-1])
        n_tail = {"tail_mean": 60, "tail_mean_1800": 30,
                  "tail_mean_7200": 120, "last_point": 1}[plateau_rule]
        self.T_plateau = float(np.mean(self.T[-n_tail:]))
        self.C_plateau = float(np.mean(self.C[-n_tail:]))
        # 灵敏度用的比例/偏置扰动(默认不改变基准)
        self.C_plateau *= plateau_scale_C
        self.T_plateau += plateau_scale_T

    def __call__(self, t):
        """返回 (T_air [degC], C_air [kg/kg])。"""
        if t <= self.t_last:
            T = float(np.interp(t, self.t, self.T))
            C = float(np.interp(t, self.t, self.C))
        else:
            T = self.T_plateau
            C = self.C_plateau
        return T, C

    def array(self, t_arr):
        t_arr = np.asarray(t_arr, dtype=float)
        T = np.interp(t_arr, self.t, self.T, left=self.T[0], right=self.T_plateau)
        C = np.interp(t_arr, self.t, self.C, left=self.C[0], right=self.C_plateau)
        return T, C


class RadiusHistory:
    """附件 2 给出的药材半径历史 R(t), 单位由 cm 换算为 m。"""

    def __init__(self, t_arr, r_cm_arr):
        self.t = np.asarray(t_arr, dtype=float)
        self.R = np.asarray(r_cm_arr, dtype=float) * 1.0e-2

    def __call__(self, t):
        if t <= self.t[-1]:
            return float(np.interp(t, self.t, self.R))
        return float(self.R[-1])

    def array(self, t_arr):
        t_arr = np.asarray(t_arr, dtype=float)
        return np.interp(t_arr, self.t, self.R, left=self.R[0], right=self.R[-1])


def constant_radius(R_value=R0):
    """返回常数半径函数, 供问题 1-3 使用。"""
    return lambda t: R_value


# ---------------------------------------------------------------------------
# 网格与插值权重
# ---------------------------------------------------------------------------
def build_grid(n_cells):
    """构造节点中心有限体积网格, 返回几何量的字典。"""
    n = int(n_cells)
    dxi = 1.0 / n
    xi_f = np.linspace(0.0, 1.0, n + 1)          # 界面: 0, dxi, ..., 1
    xi_c = (np.arange(n) + 0.5) * dxi            # 单元中心
    vol = 0.5 * (xi_f[1:] ** 2 - xi_f[:-1] ** 2)  # 单元归一化体积因子
    return {"n": n, "dxi": dxi, "xi_f": xi_f, "xi_c": xi_c, "vol": vol}


def lagrange_weights_grid(n_cells, xi_target):
    """对给定归一化采样位置构造三点二次 Lagrange 插值权重。

    返回 (idx, w): idx 形状 (m,3) 为所用单元中心下标, w 形状 (m,3) 为权重。
    xi_target 可以包含 0 与 1; 对边界点采用最近三点外推, 误差与内部同阶
    (中心处被插值函数为偶函数, 二次外推精度更高)。
    """
    grid = build_grid(n_cells)
    xi_c = grid["xi_c"]
    n = grid["n"]
    idx = np.empty((len(xi_target), 3), dtype=int)
    w = np.empty((len(xi_target), 3), dtype=float)
    for m, xt in enumerate(np.asarray(xi_target, dtype=float)):
        j = int(np.searchsorted(xi_c, xt))
        j0 = min(max(j - 1, 0), n - 3)
        nodes = np.array([j0, j0 + 1, j0 + 2])
        xn = xi_c[nodes]
        # 三点二次 Lagrange 基函数在 xt 处的取值
        wm = np.empty(3)
        for a in range(3):
            num = 1.0
            den = 1.0
            for b in range(3):
                if b != a:
                    num *= (xt - xn[b])
                    den *= (xn[a] - xn[b])
            wm[a] = num / den
        idx[m] = nodes
        w[m] = wm
    return idx, w


def apply_weights(X, idx, w):
    """按预计算权重把单元中心场插值到采样位置。

    X 形状 (n,) 或 (nt,n); 返回 (m,) 或 (nt,m)。
    """
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        return np.einsum("ma,ma->m", w, X[idx])
    return np.einsum("ma,nma->nm", w, X[:, idx])


# ---------------------------------------------------------------------------
# 核心求解器
# ---------------------------------------------------------------------------
class DryingResult:
    """求解结果容器(只保留可复现的少量字段)。"""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def summary(self):
        return {
            "t_end": float(self.times[-1]) if len(self.times) else 0.0,
            "max_C_final": float(self.max_C[-1]) if len(self.max_C) else np.nan,
            "min_C_final": float(self.min_C[-1]) if len(self.min_C) else np.nan,
            "n_steps": int(self.n_steps),
            "picard_mean": float(np.mean(self.picard_iters)),
            "picard_max": int(np.max(self.picard_iters)) if len(self.picard_iters) else 0,
            "wall_seconds": float(self.wall_seconds),
        }


def _assemble_and_solve(cap, theta, dt, g_up, g_surf, outside, x_old, l_old, r_new):
    """装配并求解一个三对角隐式系统(theta 法)。

    参数
    ----
    cap       : 新时间层容量系数数组 (n,), 热方程 rho*cp*vol, 传质方程 vol
    theta     : 时间加权系数 (1.0 = 隐式 Euler, 0.5 = Crank-Nicolson)
    dt        : 时间步长
    g_up      : 界面导流系数, g_up[i] 为单元 i 与 i+1 之间界面(含面积/步长因子)
    g_surf    : 表面等效导流系数(N-1 单元外边界), g_surf = R*G_eff
    outside   : 外部环境值(空气温度或空气含湿量), 标量
    x_old     : 旧时间层单元中心值 (n,), 在 Picard 循环中保持不变
    l_old     : 旧时间层离散算子值 (n,), 使用旧物性、旧几何计算, 同样保持不变
    r_new     : 新时间层半径 [m]

    离散格式(容量系数取新时间层, 与 rho(C)cp(C) 显含 C 的控制方程自洽):
        cap*(x_new-x_old)/dt = theta*L_new(x_new) + (1-theta)*L_old(x_old)
    其中 L_new 的系数由 Picard 当前迭代的 C,T 计算；容量取当前迭代层以保持
    与连续方程的显式容量形式一致。旧时间层算子 l_old 在时间步内冻结。
    """
    n = x_old.size
    diag = np.empty(n)
    upper = np.empty(n - 1)
    lower = np.empty(n - 1)

    # 单元 0: 中心对称, 下界面通量为零
    diag[0] = cap[0] / dt + theta * g_up[0] / r_new ** 2
    upper[0] = -theta * g_up[0] / r_new ** 2
    # 内部单元 1..n-2
    if n > 2:
        diag[1:n - 1] = cap[1:n - 1] / dt + theta * (g_up[1:n - 1] + g_up[0:n - 2]) / r_new ** 2
        upper[1:n - 1] = -theta * g_up[1:n - 1] / r_new ** 2
        lower[0:n - 2] = -theta * g_up[0:n - 2] / r_new ** 2
    # 表面单元 n-1: 外边界为半单元导热 + 表面膜阻力串联
    diag[n - 1] = cap[n - 1] / dt + theta * (g_up[n - 2] + g_surf) / r_new ** 2
    lower[n - 2] = -theta * g_up[n - 2] / r_new ** 2

    # 右端项: 旧时间层取值(时间步内固定) + 表面边界源项
    rhs = (cap / dt) * x_old + (1.0 - theta) * l_old
    rhs[n - 1] += theta * g_surf * outside / r_new ** 2

    ab = np.empty((3, n))
    ab[0, 1:] = upper
    ab[1, :] = diag
    ab[2, :-1] = lower
    return solve_banded((1, 1), ab, rhs)


def _operator(x, g_up, g_surf, outside, r):
    """离散扩散算子 L(x): 返回 (1/R^2)*[净流入], 表面用等效导流系数。"""
    n = x.size
    flux_up = np.zeros(n)
    flux_up[:n - 1] = g_up * (x[1:] - x[:-1])
    flux_dn = np.zeros(n)
    flux_dn[1:] = flux_up[:n - 1]
    out = np.zeros(n)
    out[:n - 1] += flux_up[:n - 1]
    out[1:] -= flux_dn[1:]
    out[n - 1] -= g_surf * (x[n - 1] - outside)
    return out / r ** 2


def simulate(n_cells, t_end, dt_spec, props, env, radius_fn, theta=1.0,
             c_stop=None, record_dt=None, record_xi=None, snapshot_times=None,
             picard_tol=1e-10, picard_max=60, latent_heat=0.0, x0=None,
             h_conv=None, h_mass=None, max_steps=5_000_000):
    """推进耦合传热传质方程。

    参数
    ----
    n_cells    : 径向控制体数目 N
    t_end      : 终止时间 [s]
    dt_spec    : 时间步长; 标量表示全程固定; 也可为 [(t_until, dt), ...] 分段
    props      : 物性函数 f(C, T_C) -> (rho, cp, k, D)
    env        : 环境函数 f(t) -> (T_air, C_air)
    radius_fn  : 半径函数 f(t) -> R(t) [m]
    theta      : 1.0 隐式 Euler; 0.5 Crank-Nicolson
    c_stop     : 若给定, 当 max(C) < c_stop 时提前结束(干燥判据)
    record_dt  : 记录间隔 [s]; None 表示只保留首末
    record_xi  : 归一化采样位置; 给定后 res.T/res.C 为 (nt, len(record_xi)) 的
                 采样序列(可显著降低长时间算例的内存); None 表示记录全场
    snapshot_times : 需要保存全场快照的时刻列表(用于画图), 存入 res.snap_*
    latent_heat: 汽化潜热 [J/kg]; 0 表示基准模型(不含蒸发耗热)
    h_conv/h_mass : 可选的对流换热/传质系数覆盖值(默认取题目附录 2 的 25 与 8e-7)

    返回 DryingResult。
    """
    h_conv = H_CONV if h_conv is None else float(h_conv)
    h_mass = HM_CONV if h_mass is None else float(h_mass)
    grid = build_grid(n_cells)
    n = grid["n"]
    dxi = grid["dxi"]
    vol = grid["vol"]
    xi_c = grid["xi_c"]

    if x0 is None:
        T = np.full(n, T_INIT_C)
        C = np.full(n, C_INIT)
    else:
        T = np.array(x0[0], dtype=float, copy=True)
        C = np.array(x0[1], dtype=float, copy=True)

    # 采样权重: 固定归一化坐标下的三点二次插值(问题 4 的输出位置由后处理换算)
    if record_xi is not None:
        ridx, rw = lagrange_weights_grid(n, np.asarray(record_xi, dtype=float))
        sample_mode = True
    else:
        ridx, rw = None, None
        sample_mode = False

    snap_wanted = np.asarray(snapshot_times, dtype=float) if snapshot_times else None
    snap_T, snap_C, snap_t, snap_R = [], [], [], []

    R = radius_fn(0.0)
    if callable(dt_spec):
        dt_of = dt_spec
    elif np.isscalar(dt_spec):
        dt_of = lambda t: float(dt_spec)  # noqa: E731
    else:
        pts = list(dt_spec)

        def dt_of(t, _pts=pts):
            for t_until, dtv in _pts:
                if t < t_until:
                    return float(dtv)
            return float(_pts[-1][1])

    n_steps = 0
    t = 0.0
    picard_hist = []
    wall0 = _time.perf_counter()

    # 记录容器
    def _pack(field):
        """采样模式下只保留感兴趣的归一化位置, 否则保留全场。"""
        if sample_mode:
            return apply_weights(field, ridx, rw)
        return field.copy()

    rec_t = [0.0]
    rec_T = [_pack(T)]
    rec_C = [_pack(C)]
    rec_R = [R]
    rec_cs = [None]
    rec_ts = [None]
    step_t = [0.0]
    step_maxC = [float(C.max())]
    step_minC = [float(C.min())]
    step_maxT = [float(T.max())]
    step_minT = [float(T.min())]
    step_Tair = [env(0.0)[0]]
    step_Cair = [env(0.0)[1]]
    # 守恒诊断量: 材料坐标下的含水率积分 Z = sum(C_i*vol_i) 与累计表面逸出量
    # F = integral[ g_surf*(C_{N-1}-C_air)/R^2 dt ]。离散质量守恒给出
    #   dZ/dt = -g_surf*(C_{N-1}-C_air)/R^2 ,
    # 因此 Z 的减少量应严格等于 F。注意不能直接用物理总水量
    #   sum(C_i*vol_i)*R^2 做检验: 在问题 4 中 R(t) 随时间变化,
    #   该量的变化还包含"几何收缩"贡献(等价于干物质质量本身在变),
    #   用它检验会把模型假设不自洽的效应误判为数值不守恒。
    water = [float(np.sum(C * vol))]
    cum_flux = [0.0]
    next_rec = record_dt if record_dt else None

    t_dry = None
    dry_snapshot = None

    while t < t_end - 1e-12 and n_steps < max_steps:
        dt = dt_of(t)
        dt = min(dt, t_end - t)
        t_new = t + dt
        R_old = R
        R_new = radius_fn(t_new)
        T_air_new, C_air_new = env(t_new)
        T_air_old, C_air_old = env(t)

        # ---- 旧时间层物性与算子(固定, 不随 Picard 迭代变化) ----
        _, _, k_o, D_o = props(C, T)
        g_up_k_o = _face_conductance(grid, k_o)
        g_up_D_o = _face_conductance(grid, D_o)
        g_surf_k_o = _surface_conductance(grid, k_o[-1], R_old, h_conv)
        g_surf_D_o = _surface_conductance(grid, D_o[-1], R_old, h_mass)
        l_old_T = _operator(T, g_up_k_o, g_surf_k_o, T_air_old, R_old)
        l_old_C = _operator(C, g_up_D_o, g_surf_D_o, C_air_old, R_old)
        # 旧时间层算子(显式部分)固定不变; 容量系数在装配时取新时间层值,
        # 见 _assemble_and_solve 的说明。

        # ---- Picard 迭代 ----
        T_m = T.copy()
        C_m = C.copy()
        for it in range(picard_max):
            rho_m, cp_m, k_m, D_m = props(C_m, T_m)
            if not (np.all(np.isfinite(rho_m)) and np.all(np.isfinite(cp_m))
                    and np.all(np.isfinite(k_m)) and np.all(np.isfinite(D_m))):
                raise FloatingPointError(
                    "non-finite properties: t=%.6g step=%d iter=%d "
                    "T[%.6g,%.6g] C[%.6g,%.6g] R=%.6g"
                    % (t, n_steps, it, np.nanmin(T_m), np.nanmax(T_m),
                       np.nanmin(C_m), np.nanmax(C_m), R_new))
            if not (np.all(np.isfinite(T_m)) and np.all(np.isfinite(C_m))):
                raise FloatingPointError(
                    "non-finite state before step: t=%.6g n_step=%d "
                    "T[%.3g,%.3g] C[%.3g,%.3g] R=%.6g"
                    % (t, n_steps, np.nanmin(T_m), np.nanmax(T_m),
                       np.nanmin(C_m), np.nanmax(C_m), R_new))
            g_up_k = _face_conductance(grid, k_m)
            g_up_D = _face_conductance(grid, D_m)
            g_surf_k = _surface_conductance(grid, k_m[-1], R_new, h_conv)
            g_surf_D = _surface_conductance(grid, D_m[-1], R_new, h_mass)
            cap_T = rho_m * cp_m * vol
            cap_C = vol
            T_new = _assemble_and_solve(cap_T, theta, dt, g_up_k, g_surf_k,
                                        T_air_new, T, l_old_T, R_new)
            C_new = _assemble_and_solve(cap_C, theta, dt, g_up_D, g_surf_D,
                                        C_air_new, C, l_old_C, R_new)
            if latent_heat > 0.0:
                # 可选诊断模型: 由本时间步的含水率变化显式加入蒸发热沉。
                # 该项不属于题面给定的基准模型，因此不与旧时间层算子重复计入。
                # 连续能量方程(单位体积): rho*cp*dT/dt = div(k grad T) + L*rho_d*dC/dt,
                # rho_d = rho/(1+C) 为干物质表观密度; 汽化所需热量由物料自身提供,
                # 故 dC/dt<0 时该项为负(降温)。
                # 归一化离散后: vol*rho*cp*dT/dt = [通量]/R^2 + vol*L*rho_d*dC/dt,
                # 因此每个时间步的温度修正量是
                #     dT = + L*rho_d*(C_new-C_old)/(rho*cp)
                # 分母是 rho*cp，而不是容量系数 cap_T=rho*cp*vol；vol 在源项两侧抵消。
                rho_d = rho_m / (1.0 + C_m)
                T_new = T_new + latent_heat * rho_d * (C_new - C) / (rho_m * cp_m)
            err_T = np.max(np.abs(T_new - T_m)) / max(1.0, float(np.max(np.abs(T_new))))
            err_C = np.max(np.abs(C_new - C_m)) / max(0.01, float(np.max(np.abs(C_new))))
            T_m = T_new
            C_m = C_new
            if max(err_T, err_C) < picard_tol:
                break
        picard_hist.append(it + 1)

        T = T_m
        C = C_m
        t = t_new
        R = R_new
        n_steps += 1

        # 表面重构值(半单元串联阻力), 用当前物性与当前环境
        _, _, k_n, D_n = props(C, T)
        c_surf = _surface_value(grid, C[-1], D_n[-1], R, h_mass, C_air_new)
        t_surf = _surface_value(grid, T[-1], k_n[-1], R, h_conv, T_air_new)

        maxC = float(max(C.max(), c_surf))
        minC = float(min(C.min(), c_surf))
        # 表面逸出(归一化): out_flux = g_surf*(C[-1]-C_air); 除以 R^2 后
        # 即为 dZ/dt 的相反数, 见上面守恒诊断量的说明。
        g_surf_D_now = _surface_conductance(grid, D_n[-1], R, h_mass)
        out_flux = g_surf_D_now * (C[-1] - C_air_new)
        cum_flux.append(cum_flux[-1] + out_flux / R ** 2 * dt)
        water.append(float(np.sum(C * vol)))
        step_t.append(t)
        step_maxC.append(maxC)
        step_minC.append(minC)
        step_maxT.append(float(max(T.max(), t_surf)))
        step_minT.append(float(min(T.min(), t_surf)))
        step_Tair.append(T_air_new)
        step_Cair.append(C_air_new)

        # 干燥终点: max C 首次跌破判据, 在线性插值处求交
        if c_stop is not None and t_dry is None and maxC < c_stop:
            m1, m0 = maxC, step_maxC[-2]
            frac = 0.0 if m0 == m1 else (m0 - c_stop) / (m0 - m1)
            t_dry = step_t[-2] + frac * (t - step_t[-2])
            dry_snapshot = (T.copy(), C.copy(), t, R)

        if next_rec is not None and t >= next_rec - 1e-9:
            # 记录时刻对齐到 record_dt 的整数倍, 避免长时间累加的浮点误差
            # 使 1 s 变成 0.9999999999999999 而被 int() 截断成 0 s
            rec_t.append(round(t / record_dt) * record_dt)
            rec_T.append(_pack(T))
            rec_C.append(_pack(C))
            rec_R.append(R)
            rec_cs.append(c_surf)
            rec_ts.append(t_surf)
            while next_rec <= t + 1e-9:
                next_rec += record_dt

        # 全场快照(仅供画图, 与 record_dt / sample 模式无关)
        if snap_wanted is not None:
            while len(snap_t) < snap_wanted.size and t >= snap_wanted[len(snap_t)] - 1e-9:
                snap_t.append(t)
                snap_T.append(T.copy())
                snap_C.append(C.copy())
                snap_R.append(R)

    wall = _time.perf_counter() - wall0

    # 末态快照始终保留, 与 record_dt 无关
    rho_f, cp_f, k_f, D_f = props(C, T)
    T_air_f, C_air_f = env(t)
    res = DryingResult(
        n_cells=n, t_end=t, times=np.array(rec_t), T=np.array(rec_T),
        C=np.array(rec_C), R=np.array(rec_R),
        c_surf=np.array([np.nan if v is None else v for v in rec_cs]),
        t_surf=np.array([np.nan if v is None else v for v in rec_ts]),
        T_final=T.copy(), C_final=C.copy(), R_final=R,
        t_surf_final=_surface_value(grid, T[-1], k_f[-1], R, h_conv, T_air_f),
        c_surf_final=_surface_value(grid, C[-1], D_f[-1], R, h_mass, C_air_f),
        step_times=np.array(step_t), step_maxC=np.array(step_maxC),
        step_minC=np.array(step_minC), step_Tair=np.array(step_Tair),
        step_maxT=np.array(step_maxT), step_minT=np.array(step_minT),
        step_Cair=np.array(step_Cair), water=np.array(water),
        cum_flux=np.array(cum_flux), t_dry=t_dry, dry_snapshot=dry_snapshot,
        n_steps=n_steps, picard_iters=np.array(picard_hist),
        wall_seconds=wall, theta=theta, dxi=dxi, xi_c=xi_c,
        grid=grid, sample_mode=sample_mode, record_xi=record_xi,
        snap_times=np.array(snap_t), snap_T=np.array(snap_T),
        snap_C=np.array(snap_C), snap_R=np.array(snap_R),
    )
    return res


def _face_conductance(grid, k_cell):
    """界面导流系数: 调和平均有效系数 x 界面半径 / dxi。

    调和平均对应相邻半单元串联阻力，保证正系数下的界面通量守恒。
    """
    k = np.asarray(k_cell, dtype=float)
    km = 2.0 * k[:-1] * k[1:] / (k[:-1] + k[1:])
    return grid["xi_f"][1:-1] * km / grid["dxi"]


def _surface_conductance(grid, k_last, R, h_film):
    """表面半单元等效导流系数 g_surf = R * G_eff。

    a = 2k/dr 为半单元(中心到表面, 厚度 dr/2)导热系数;
    G_eff = h*a/(a+h) 为"半单元导热阻力与表面膜阻力串联"的等效换热系数。
    """
    dr = R * grid["dxi"]
    a = 2.0 * k_last / dr
    g_eff = h_film * a / (a + h_film)
    return R * g_eff


def _surface_value(grid, x_last, k_last, R, h_film, outside):
    """由半单元串联阻力重构真实表面值(不是最后一个单元中心值)。"""
    dr = R * grid["dxi"]
    a = 2.0 * k_last / dr
    return (a * x_last + h_film * outside) / (a + h_film)


def dimensionless_numbers(props, C=1.0, T_C=50.0, R=R0, h=H_CONV, h_m=HM_CONV):
    """返回 (Fo 特征时间, Bi_heat, Bi_mass, 半径扩散时间常数)。

    用于量纲与数量级检查, 不参与求解。
    """
    rho, cp, k, D = props(np.array([C]), np.array([T_C]))
    rho, cp, k, D = float(rho[0]), float(cp[0]), float(k[0]), float(D[0])
    alpha = k / (rho * cp)
    return {
        "rho": rho, "cp": cp, "k": k, "D": D,
        "alpha_thermal": alpha,
        "t_diff_mass": R ** 2 / D,
        "t_diff_heat": R ** 2 / alpha,
        "Bi_heat": h * R / k,
        "Bi_mass": h_m * R / D,
    }
