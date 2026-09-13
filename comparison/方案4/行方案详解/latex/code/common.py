# -*- coding: utf-8 -*-
# ============================================================================
# common.py —— 公共模块：物性经验关系、烘房环境、半径历史、有限体积求解器
#
# 其余脚本均 import 本模块。除注明外一律用 SI 单位；经验式中的温度取绝对
# 温度 K，故代入前先加 273.15。
#
# 题目给定的常量：
#   R0, HERB_LENGTH   药材半径 0.02 m、长度 0.25 m
#   T_INIT_C, C_INIT  初始温度 28 degC、初始干基含水率 2.55 kg/kg
#   H_CONV, HM_CONV   对流换热系数 25 W/(m2 K)、对流传质系数 8e-7 m/s
#
# 求解的方程（xi = r/R(t) 为归一化半径，R(t) 见 RadiusHistory，不收缩时恒为 R0）：
#   水分  dC/dt        = (1/(xi R^2)) d/dxi ( xi D(C,T) dC/dxi )
#   温度  rho cp dT/dt = (1/(xi R^2)) d/dxi ( xi k(C)   dT/dxi )
#   边界  xi=0 对称；xi=1 处 -k dT/dr = h(Ts-Tair)、-D dC/dr = hm(Cs-Cair)
#   收缩问题改用物质坐标后对流项严格为零；R=const 时两式退化为标准圆柱扩散方程。
#
# 数值方法：节点中心有限体积（节点恰好落在 xi=0 与 xi=1 上），theta 法时间推进
#   （theta=0.5 为 Crank-Nicolson，前 4 步改用 theta=1 以抑制初值跳跃引起的振荡），
#   每步解两个三对角方程组。网格与时间步的收敛性见 resolution.py、accuracy.py。
#
# 典型用法：
#   room = RoomConditions()                       # 烘房环境
#   res = simulate(1600, props_problem23, t_end, dt, room,
#                  t_report=times, xi_report=radii)   # 求解并采样
# ============================================================================

import os
import numpy as np
import pandas as pd
from scipy.linalg import solve_banded

# --------------------------------------------------------------------------
# 题目给定的常量
# --------------------------------------------------------------------------
R0 = 0.02            # 药材初始半径                                  [m]
HERB_LENGTH = 0.25   # 药材长度                                      [m]
T_INIT_C = 28.0      # 药材与烘房的初始温度                          [degC]
C_INIT = 2.55        # 药材初始干基含水率                            [kg/kg]
H_CONV = 25.0        # 对流换热系数 h                                [W/(m2 K)]
HM_CONV = 8.0e-7     # 对流传质系数 hm                               [m/s]
T_KELVIN = 273.15    # 摄氏度换成开尔文需加的常数

# 烘房空气状态：预热平衡阶段（0 ~ 14400 s）取附件 1 拟合式，
# 之后烘房稳定在设定值（详见 RoomConditions）。
T_PREHEAT_END = 14400.0   # 预热平衡阶段结束时刻                      [s]
T_SET_C = 50.212          # 恒温干燥阶段设定温度（附件 1 拟合值）     [degC]
C_SET = 0.050908          # 恒温干燥阶段设定含湿量（附件 1 拟合值）   [kg/kg]
TAU_AIR_T = 1877.5        # 温度一阶惯性式的时间常数                  [s]
TAU_AIR_C = 2776.3        # 含湿量一阶惯性式的时间常数                [s]
C_AIR_INIT = 0.01963      # 烘房初始含湿量（附件 1 实测初值）         [kg/kg]

# 含水率下限保护：D 含 exp(-a/C)，C -> 0 时发散，故把 C 截断到 C_FLOOR。
C_FLOOR = 1.0e-3

WORKDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(WORKDIR, "data")
FIG_DIR = os.path.join(WORKDIR, "figures")


def _find_topic_dir():
    """返回题目附件所在目录（即含「附件/附件1.xlsx」的目录）。

    按候选顺序逐个探测，可用环境变量 HERB_TOPIC_DIR 覆盖；全部候选都找不到
    时抛出 FileNotFoundError，避免后续静默读到错误数据。
    """
    env = os.environ.get("HERB_TOPIC_DIR")
    cands = ([env] if env else []) + [
        os.path.join(WORKDIR, "..", "..", "题目", "A题"),
        os.path.join(WORKDIR, "..", "题目", "A题"),
        os.path.join(WORKDIR, "题目", "A题"),
        r"D:\建模国赛\题目\A题",
    ]
    for cand in cands:
        cand = os.path.abspath(cand)
        if os.path.isfile(os.path.join(cand, "附件", "附件1.xlsx")):
            return cand
    raise FileNotFoundError(
        "找不到题目附件目录（应含 附件/附件1.xlsx）。已尝试：%s；"
        "可用环境变量 HERB_TOPIC_DIR 指定。" % cands)


TOPIC_DIR = _find_topic_dir()

for _d in (DATA_DIR, FIG_DIR):
    os.makedirs(_d, exist_ok=True)


# --------------------------------------------------------------------------
# 题目附录 2 / 3 / 4 给出的物性经验关系
# 三个函数签名统一为 f(C, T_K) -> (rho, cp, k, D)，便于求解器统一调用
# --------------------------------------------------------------------------
def props_problem1(C, T_K=None):
    """附录 2（问题 1）：常数物性，仅 D 随含水率变化，D = 7e-9 exp(-0.89/C)。"""
    C = np.maximum(np.asarray(C, dtype=float), C_FLOOR)
    rho = np.full_like(C, 820.0)
    cp = np.full_like(C, 2600.0)
    k = np.full_like(C, 0.36)
    D = 7.0e-9 * np.exp(-0.89 / C)
    return rho, cp, k, D


def props_problem23(C, T_K):
    """附录 3（问题 2、3）：rho、cp、k 随 C 线性/分式变化，D 同时依赖 C 与 T。"""
    C = np.maximum(np.asarray(C, dtype=float), C_FLOOR)
    frac = C / (C + 1.0)
    rho = 650.0 + 128.0 * C
    cp = 1450.0 + 2736.0 * frac
    k = 0.21 + 0.38 * frac
    D = 2.4e-3 * np.exp(-0.45 / C) * np.exp(-3850.0 / T_K)
    return rho, cp, k, D


def props_problem4(C, T_K):
    """附录 4（问题 4）：收缩工况的物性式，形式同附录 3 但系数不同。"""
    C = np.maximum(np.asarray(C, dtype=float), C_FLOOR)
    frac = C / (C + 1.0)
    rho = 760.0 + 90.0 * C
    cp = 1850.0 + 2150.0 * frac
    k = 0.12 + 0.20 * frac
    D = 4.2e-4 * np.exp(-0.30 / C) * np.exp(-3850.0 / T_K)
    return rho, cp, k, D


PROPERTY_SETS = {
    "p1": props_problem1,
    "p23": props_problem23,
    "p4": props_problem4,
}


# --------------------------------------------------------------------------
# 烘房（干燥介质）状态
# --------------------------------------------------------------------------
class RoomConditions:
    """烘房空气的温湿度历程，作为两个场的 Robin 边界条件。

    预热平衡阶段用一阶惯性式拟合附件 1 的实测序列，初值固定为实测值：
        T_air(t) = t_set - (t_set - 28)      * exp(-t / tau_t)
        C_air(t) = c_set - (c_set - 0.01963) * exp(-t / tau_c)
    t >= t_switch 后（hold_after=True）保持恒温干燥阶段的设定值 t_set、c_set。
    """

    def __init__(self, t_switch=T_PREHEAT_END, t_set=T_SET_C, c_set=C_SET,
                 tau_t=TAU_AIR_T, tau_c=TAU_AIR_C, hold_after=True):
        self.t_switch = float(t_switch)
        self.t_set = float(t_set)
        self.c_set = float(c_set)
        self.tau_t = float(tau_t)
        self.tau_c = float(tau_c)
        self.hold_after = bool(hold_after)

    def T_air(self, t):
        t = np.asarray(t, dtype=float)
        val = self.t_set - (self.t_set - T_INIT_C) * np.exp(-t / self.tau_t)
        if self.hold_after:
            val = np.where(t >= self.t_switch, self.t_set, val)
        return val if val.ndim else float(val)

    def C_air(self, t):
        t = np.asarray(t, dtype=float)
        val = self.c_set - (self.c_set - C_AIR_INIT) * np.exp(-t / self.tau_c)
        if self.hold_after:
            val = np.where(t >= self.t_switch, self.c_set, val)
        return val if val.ndim else float(val)


def load_attachment1():
    """读附件 1：烘房温度与水分浓度的实测序列，0 ~ 14400 s 每 60 s 一点。"""
    df = pd.read_excel(os.path.join(TOPIC_DIR, "附件", "附件1.xlsx"))
    df.columns = ["t", "T_air", "C_air"]
    return df


def load_attachment2():
    """读附件 2：药材半径的实测序列，0 ~ 259200 s 每 1800 s 一点。"""
    df = pd.read_excel(os.path.join(TOPIC_DIR, "附件", "附件2.xlsx"))
    df.columns = ["t", "R_cm"]
    return df


class RadiusHistory:
    """药材半径 R(t)，取自附件 2（问题 4 使用；其余问题传 None 即恒为 R0）。

    mode='table'  ：分段线性插值实测表（在每个实测点上精确通过，区间外取端值）；
    mode='smooth' ：单指数拟合 R = rinf + (R0 - rinf) exp(-t/tau)，用于灵敏度对照；
    scale         ：整体缩放因子，用于灵敏度分析中扰动收缩幅度。
    """

    def __init__(self, mode="table", scale=1.0):
        df = load_attachment2()
        self.t_tab = df["t"].to_numpy(float)
        self.R_tab = df["R_cm"].to_numpy(float) * 1.0e-2
        self.mode = mode
        self.scale = float(scale)
        self.R0 = self.R_tab[0] * self.scale
        if mode == "smooth":
            from scipy.optimize import curve_fit

            def model(x, rinf, tau):
                return rinf + (self.R_tab[0] - rinf) * np.exp(-x / tau)

            popt, _ = curve_fit(model, self.t_tab, self.R_tab, p0=[1.2e-2, 1.45e4],
                                maxfev=400000)
            self.rinf, self.tau = float(popt[0]), float(popt[1])
        elif mode != "table":
            raise ValueError("mode must be 'table' or 'smooth'")

    def R(self, t):
        t = np.asarray(t, dtype=float)
        if self.mode == "smooth":
            r = self.rinf + (self.R_tab[0] - self.rinf) * np.exp(-t / self.tau)
            r = np.where(t > self.t_tab[-1], self.R_tab[-1], r)
        else:
            r = np.interp(t, self.t_tab, self.R_tab)
        return (r * self.scale) if r.ndim else float(r * self.scale)

    def dRdt(self, t):
        """用中心差分求 dR/dt [m/s]（步长取 d s，用于能量方程的体积项）。"""
        t = float(t)
        d = 5.0
        return float((self.R(t + d) - self.R(t - d)) / (2.0 * d))


# --------------------------------------------------------------------------
# 有限体积求解器
# --------------------------------------------------------------------------
class DryingSolver:
    """轴对称有限体积求解器，求解域为归一化半径 xi = r/R(t) 的 [0, 1]。

    构造参数：
      n_cell     控制体（= 区间）数目；节点格式下未知量个数为 n_cell+1
      prop       物性函数 f(C, T_K) -> (rho, cp, k, D)
      theta      时间格式：1.0 隐式 Euler，0.5 Crank-Nicolson（默认）
      radius     RadiusHistory 实例；None 表示半径不变（问题 1~3）
      grid       'node' 节点中心（默认，节点落在 xi=0 与 xi=1 上）
                 'cell' 单元中心，作为独立空间离散用于收敛性对照
      volume_term 能量方程是否计入体积变化项（仅问题 4 灵敏度用）
      evap_cooling 表面能量平衡是否扣除蒸发潜热（默认关闭）

    关于 grid 的选择：D(C) = D0 exp(-a/C) 在表面含水率下极小，xi=1 附近的水分
    剖面近乎垂直，故表面值必须直接存放而不是由内部单元外推——单元中心格式
    的表面对流条件只有 O(h) 精度，节点中心格式则精确落在边界上。
    """

    def __init__(self, n_cell, prop, theta=0.5, radius=None, volume_term=False,
                 room=None, n_rannacher=4, grid="node", h_conv=None,
                 hm_conv=None, evap_cooling=False, latent_heat=2.383e6):
        self.n = int(n_cell)
        self.prop = prop
        self.theta = float(theta)
        self.radius = radius
        self.volume_term = bool(volume_term)
        self.room = room if room is not None else RoomConditions()
        self.grid = grid
        self.h_conv = H_CONV if h_conv is None else float(h_conv)
        self.hm_conv = HM_CONV if hm_conv is None else float(hm_conv)
        # 可选的物理修正：把蒸发水带走的潜热从表面能量平衡中扣掉
        self.evap_cooling = bool(evap_cooling)
        self.latent_heat = float(latent_heat)
        # Rannacher 启动：最初 n_rannacher 步用隐式 Euler。初始时刻表面与
        # 环境存在跳跃，Crank-Nicolson 对最高频空间模态几乎不衰减，会产生振荡。
        self.n_rannacher = int(n_rannacher)
        self.n_step_done = 0

        h = 1.0 / self.n
        self.h = h
        if grid == "node":
            # 节点中心格式：n+1 个节点，含轴心 xi=0 与表面 xi=1
            self.xi_c = np.arange(self.n + 1) * h          # storage positions
            self.xi_f = (np.arange(self.n) + 0.5) * h      # interior faces
            vol = np.empty(self.n + 1)
            vol[0] = h * h / 8.0
            vol[1:self.n] = self.xi_c[1:self.n] * h
            vol[self.n] = h / 2.0 - h * h / 8.0
            self.vol = vol
            self.nv = self.n + 1
        elif grid == "cell":
            self.xi_c = (np.arange(self.n) + 0.5) * h
            self.xi_f = np.arange(self.n + 1) * h
            self.vol = self.xi_c * h
            self.nv = self.n
        else:
            raise ValueError("grid must be 'node' or 'cell'")
        # 单元中心格式：d(phi_i)/dt = -M_i (g_{i+1} - g_i)
        self.M = 2.0 / (h ** 2 * (2.0 * np.arange(self.n) + 1.0))
        self.t = 0.0
        self.R = R0 if radius is None else float(radius.R(0.0))
        self.reset()

    # -- 状态量：场变量与物性缓存 -----------------------------------------
    def reset(self, T_init=None, C_init=None):
        self.T = np.full(self.nv, T_INIT_C if T_init is None else T_init)
        self.C = np.full(self.nv, C_INIT if C_init is None else C_init)
        self.t = 0.0
        self.R = R0 if self.radius is None else float(self.radius.R(0.0))
        self.n_step_done = 0
        self._cache()

    def _cache(self):
        rho, cp, k, D = self.prop(self.C, self.T + T_KELVIN)
        self.rho, self.cp, self.k, self.D = rho, cp, k, D
        self.rct = rho * cp

    # -- 离散算子：通量、三对角系数矩阵与右端项组装 -----------------------
    @staticmethod
    def _face_coeff(a):
        """单元中心格式的面系数：内部面取相邻单元的算术平均（长度 n+1，两端取端值）。"""
        af = np.zeros(len(a) + 1)
        af[1:len(a)] = 0.5 * (a[:-1] + a[1:])
        af[0] = a[0]
        af[-1] = a[-1]
        return af

    def _assemble(self, a, cap, coef, env, dt, phi, extra_diag=0.0, theta=None):
        """组装一个场（温度或水分）的 theta 法线性系统。

        连续形式： dphi/dt = 1/(xi R^2) * d/dxi ( xi * a * dphi/dxi )
        其中 a 为扩散系数（水分取 D，温度取 k），cap 为体积热容（水分取 1）。

        表面条件用守恒通量表述：记 g = xi a dphi/dxi，由 a dphi/dxi = -R q_r
        与 q_r|_R = coef (phi_s - phi_env) 得 xi=1 处
            g_n = -R * coef * (phi_s - phi_env)          （coef 为 h 或 hm）

        对第 i 个控制体做有限体积积分：vol_i dphi_i/dt = (g_in - g_out) / R^2，
        即 dphi_i/dt = P_i (g_in - g_out)，P_i = 1/(R^2 vol_i cap_i)。

        返回 scipy.solve_banded 所需的 (1, 1) 带宽矩阵 ab 与右端 rhs。
        """
        th = self.theta if theta is None else theta
        R = self.R

        if self.grid == "node":
            n = self.n
            nv = self.nv
            af = _face_mean(a)                            # n 个内部面系数
            beta = self.xi_f * af / self.h                # 第 j 个面的 beta_j
            P = 1.0 / (R * R * self.vol * cap)

            gf = beta * (phi[1:] - phi[:n])               # g at face j
            g_s = -R * coef * (phi[n] - env)              # g at the surface
            f = np.empty(nv)
            f[0] = P[0] * gf[0]
            f[1:n] = P[1:n] * (gf[1:] - gf[:n - 1])
            f[n] = P[n] * (g_s - gf[n - 1])
            f = f + extra_diag * phi

            dg = np.empty(nv)
            lower = np.empty(nv - 1)
            upper = np.empty(nv - 1)
            dg[0] = -P[0] * beta[0]
            upper[0] = P[0] * beta[0]
            dg[1:n] = -P[1:n] * (beta[:n - 1] + beta[1:])
            lower[:n - 1] = P[1:n] * beta[:n - 1]
            upper[1:] = P[1:n] * beta[1:]
            dg[n] = -P[n] * (beta[n - 1] + R * coef)
            lower[n - 1] = P[n] * beta[n - 1]
            # extra_diag 是与 phi 成正比的一阶项（问题 4 能量方程的体积变化项
            # -2*Rdot/R），必须同时进入显式项 f 与隐式对角 dg，否则按 theta 法
            # 组装时它只按 (1-theta) 起作用：theta=1（Rannacher 前几步）时整项
            # 消失，结果与单元中心布局不一致。cell 分支一直是这样处理的。
            dg = dg + extra_diag
            const = np.zeros(nv)
            const[n] = P[n] * R * coef * env
        else:
            n = self.n
            nv = self.nv
            af = self._face_coeff(a)
            beta = np.zeros(n + 1)
            beta[1:n] = self.xi_f[1:n] * af[1:n] / self.h
            P = self.M / (R * R * cap)

            g = np.zeros(n + 1)
            g[1:n] = beta[1:n] * (phi[1:n] - phi[:n - 1])
            phi_s = 1.5 * phi[n - 1] - 0.5 * phi[n - 2]
            g[n] = -R * coef * (phi_s - env)
            f = P * (g[1:n + 1] - g[:n]) + extra_diag * phi

            w_surf = R * coef
            dg = -P * (beta[:n] + beta[1:n + 1])
            dg[n - 1] = P[n - 1] * (-1.5 * w_surf - beta[n - 1])
            dg = dg + extra_diag
            lower = P[1:] * beta[1:n]
            lower[n - 2] = P[n - 1] * (0.5 * w_surf + beta[n - 1])
            upper = P[:n - 1] * beta[1:n]

            const = np.zeros(n)
            const[n - 1] = P[n - 1] * w_surf * env

        ab = np.zeros((3, nv))
        ab[0, 1:] = -dt * th * upper
        ab[1, :] = 1.0 - dt * th * dg
        ab[2, :-1] = -dt * th * lower

        rhs = phi + dt * (1.0 - th) * f + dt * th * const
        return ab, rhs

    # -- 推进一步：解两个三对角方程组，更新温度场与水分浓度场 ------------
    def step(self, dt):
        """推进一个时间步 dt：物性取上一步的值（冻结系数），两个场隐式求解。"""
        t_new = self.t + dt
        t_mid = self.t + 0.5 * dt
        self.R = R0 if self.radius is None else float(self.radius.R(t_new))
        th = 1.0 if self.n_step_done < self.n_rannacher else self.theta

        T_air = float(self.room.T_air(t_mid))
        C_air = float(self.room.C_air(t_mid))

        abC, rhsC = self._assemble(self.D, np.ones(self.nv), self.hm_conv,
                                   C_air, dt, self.C, theta=th)
        C_new = solve_banded((1, 1), abC, rhsC)

        extra = 0.0
        if self.volume_term and self.radius is not None:
            rdot = float(self.radius.dRdt(t_new))
            extra = -2.0 * rdot / self.R
        abT, rhsT = self._assemble(self.k, self.rct, self.h_conv, T_air, dt,
                                   self.T, extra_diag=extra, theta=th)
        T_new = solve_banded((1, 1), abT, rhsT)

        if self.evap_cooling:
            # 蒸发汇加在最外层控制体上，单位面积质量通量为
            # rho_dry * hm * (C_s - C_air)，乘潜热后即为表面热流
            if self.grid == "node":
                C_s = self.C[-1]
                v_last = self.vol[-1]
                cap_last = self.rct[-1]
            else:
                C_s = 1.5 * self.C[-1] - 0.5 * self.C[-2]
                v_last = self.vol[-1]
                cap_last = self.rct[-1]
            rho_dry = self.rho[-1] / (1.0 + max(self.C[-1], C_FLOOR))
            sink = self.R * self.latent_heat * rho_dry * self.hm_conv * \
                (C_s - C_air)
            T_new[-1] -= dt * sink / (self.R ** 2 * v_last * cap_last)

        self.C = np.maximum(C_new, 0.0)
        self.T = T_new
        self.t = t_new
        self.n_step_done += 1
        self._cache()

    # -- 采样：把节点值插值到指定的输出位置 -------------------------------
    def sample_xi(self, xi):
        """把温度场与水分场插值到给定的归一化半径 xi 上。"""
        return (interp_cells(self.T, self.xi_c, xi),
                interp_cells(self.C, self.xi_c, xi))


def _face_mean(a):
    """节点中心网格内部面的面扩散系数：相邻两节点值的算术平均。

    输入 a 为 n+1 个节点值（xi = 0 ... 1），返回 n 个内部面值。

    面系数怎么取不影响守恒性——守恒来自 g_{i+1/2}-g_{i-1/2} 的通量差分形式，
    与面值本身无关；两种平均也都是空间二阶。差别只在误差常数：本问题 D 随
    含水率呈指数变化（可达几个数量级）但剖面光滑，调和平均在 D 梯度大的地方
    会引入 ∝(D')^2/D 的额外误差项，实测算术平均的误差常数明显更小
    （问题 3、N=1600 时中心含水率误差由 7.7e-6 降到 2.0e-6）。
    若系数本身是分段常数（如多层材料界面），则应当改用调和平均。
    """
    return 0.5 * (a[:-1] + a[1:])


def interp_cells(vals, xi_c, xi_t):
    """把网格上的场值二次重构到任意 xi（输出位置一般不在节点上）。

    分三段处理：
      内部       用相邻三点的 Lagrange 二次插值；
      xi=0 附近  利用偶对称性 phi(xi) = phi(0) + a xi^2；
      xi=1 附近  用最后三点的二次外推（节点中心格式下输出点不会超出 xi=1）。
    """
    vals = np.asarray(vals, dtype=float)
    xi_t = np.atleast_1d(np.asarray(xi_t, dtype=float))
    n = len(vals)
    out = np.empty_like(xi_t)

    # 内部插值：定位落在 [xi_c[k-1], xi_c[k]] 的区间，并限制到 [1, n-2]
    idx = np.searchsorted(xi_c, xi_t)
    k = np.clip(idx, 1, n - 2)

    x0, x1, x2 = xi_c[k - 1], xi_c[k], xi_c[k + 1]
    w0 = (xi_t - x1) * (xi_t - x2) / ((x0 - x1) * (x0 - x2))
    w1 = (xi_t - x0) * (xi_t - x2) / ((x1 - x0) * (x1 - x2))
    w2 = (xi_t - x0) * (xi_t - x1) / ((x2 - x0) * (x2 - x1))
    out[:] = w0 * vals[k - 1] + w1 * vals[k] + w2 * vals[k + 1]

    # 轴心附近：由前两点做偶对称二次式
    a_even = (vals[1] - vals[0]) / (xi_c[1] ** 2 - xi_c[0] ** 2)
    mask0 = xi_t < xi_c[0]
    out[mask0] = vals[0] - a_even * xi_c[0] ** 2 + a_even * xi_t[mask0] ** 2

    # 表面附近及以外：由最后三点做二次外推
    mask1 = xi_t > xi_c[-1]
    if np.any(mask1):
        j = n - 1
        p0, p1, p2 = xi_c[j - 2], xi_c[j - 1], xi_c[j]
        v0, v1, v2 = vals[j - 2], vals[j - 1], vals[j]
        xx = xi_t[mask1]
        c0 = v0 / ((p0 - p1) * (p0 - p2))
        c1 = v1 / ((p1 - p0) * (p1 - p2))
        c2 = v2 / ((p2 - p0) * (p2 - p1))
        out[mask1] = (c0 * (xx - p1) * (xx - p2) + c1 * (xx - p0) * (xx - p2)
                      + c2 * (xx - p0) * (xx - p1))
    return out


# --------------------------------------------------------------------------
# 时间推进驱动：按给定的步长安排推进，并在指定的输出时刻采样
# --------------------------------------------------------------------------
def simulate(n_cell, prop, t_end, dt, room, radius=None, theta=0.5,
             t_report=None, xi_report=None, volume_term=False,
             dt_schedule=None, verbose=False, grid="node", h_conv=None,
             hm_conv=None, evap_cooling=False, latent_heat=2.383e6,
             track_max=False):
    """按给定的步长安排推进求解，并在指定时刻采样输出。

    xi_report 决定在哪些位置取值，支持四种写法：
      None          用网格自身的存储位置；
      一维数组      固定的归一化半径 xi；
      二维数组      形状 (len(t_report), n_out)，每个输出时刻一行；
      可调用对象    f(t) -> 一维数组，在每个输出时刻求值。

    后两种是必要的：输出位置通常按物理距离 r 给定，而求解器工作在物质坐标
    xi = r/R(t) 上，同一个物理半径在每个输出时刻对应不同的 xi。

    track_max=True 时，每个输出时刻额外记录全场的最大含水率 Cmax、其位置
    xi_at_max，以及轴心值 C_axis 与场最小值 C_min。判据要用全场最值而不是
    少量输出点上的最值；剖面仍均匀到舍入精度时 argmax 只是并列的浮点数，
    其位置不携带信息，判断「中心最后干」应比较 Cmax 与 C_axis。

    返回 dict：t 输出时刻，T/C 形状 (len(t_report), n_out)，R 各时刻半径，
    Cmax 与 xi_at_max（仅 track_max），solver 末态，n_step 总步数。
    """
    slv = DryingSolver(n_cell, prop, theta=theta, radius=radius,
                       volume_term=volume_term, room=room, grid=grid,
                       h_conv=h_conv, hm_conv=hm_conv,
                       evap_cooling=evap_cooling, latent_heat=latent_heat)

    if t_report is None:
        t_report = np.arange(0.0, t_end + 0.5 * dt, dt)
    t_report = np.asarray(t_report, dtype=float)
    n_rep = len(t_report)

    # 把 xi_report 统一成可调用对象 xi_of(t) -> 一维数组
    n_out = None
    if xi_report is None:
        xi_fixed = np.asarray(slv.xi_c, dtype=float)
        n_out = len(xi_fixed)
        xi_of = lambda t: xi_fixed                                   # noqa: E731
    elif callable(xi_report):
        probe = np.atleast_1d(np.asarray(xi_report(t_report[0]), dtype=float))
        n_out = len(probe)
        xi_of = lambda t: np.atleast_1d(                            # noqa: E731
            np.asarray(xi_report(t), dtype=float))
    else:
        xi_arr = np.asarray(xi_report, dtype=float)
        if xi_arr.ndim == 2:
            if xi_arr.shape[0] != n_rep:
                raise ValueError("2-D xi_report must have one row per report "
                                 "time (%d), got %d" % (n_rep, xi_arr.shape[0]))
            n_out = xi_arr.shape[1]
            xi_of = lambda t, _a=xi_arr: _a[                          # noqa: E731
                int(np.argmin(np.abs(t_report - t)))]
        else:
            n_out = len(xi_arr)
            xi_of = lambda t: xi_arr                                 # noqa: E731

    steps = [(t_end, dt)] if dt_schedule is None else list(dt_schedule)

    T_out = np.full((n_rep, n_out), np.nan)
    C_out = np.full((n_rep, n_out), np.nan)
    R_out = np.full(n_rep, np.nan)
    Cmax_out = np.full(n_rep, np.nan)
    xi_max_out = np.full(n_rep, np.nan)
    Caxis_out = np.full(n_rep, np.nan)
    Cmin_out = np.full(n_rep, np.nan)

    prev_t = slv.t
    prev_T = slv.T.copy()
    prev_C = slv.C.copy()
    prev_R = slv.R
    k_rep = 0
    n_step = 0
    tol = 1.0e-9

    def flush(t_now):
        nonlocal k_rep, prev_t, prev_T, prev_C, prev_R
        while k_rep < n_rep and t_report[k_rep] <= t_now + tol:
            tr = t_report[k_rep]
            if tr <= prev_t + tol:
                Tf, Cf, Rf = prev_T, prev_C, prev_R
            else:
                w = (tr - prev_t) / (t_now - prev_t)
                Tf = prev_T + w * (slv.T - prev_T)
                Cf = prev_C + w * (slv.C - prev_C)
                Rf = prev_R + w * (slv.R - prev_R)
            xi_t = xi_of(tr)
            T_out[k_rep] = interp_cells(Tf, slv.xi_c, xi_t)
            C_out[k_rep] = interp_cells(Cf, slv.xi_c, xi_t)
            R_out[k_rep] = Rf
            if track_max:
                # 取全场最值（而不是只在采样点上取最值）
                Cmax_out[k_rep] = Cf.max()
                xi_max_out[k_rep] = slv.xi_c[int(np.argmax(Cf))]
                Caxis_out[k_rep] = Cf[0]
                Cmin_out[k_rep] = Cf.min()
            k_rep += 1
        prev_t = slv.t
        prev_T = slv.T.copy()
        prev_C = slv.C.copy()
        prev_R = slv.R

    flush(0.0)
    t_cur = 0.0
    for t_target, h in steps:
        while t_cur < t_target - tol:
            dt_use = min(h, t_target - t_cur)
            slv.step(dt_use)
            t_cur = slv.t
            n_step += 1
            flush(t_cur)
            if verbose and n_step % 50000 == 0:
                print("   t=%.0f s  C_center=%.5f  T_center=%.3f"
                      % (slv.t, slv.C[0], slv.T[0]), flush=True)
    flush(t_cur + 1.0)

    return dict(t=t_report, T=T_out, C=C_out, R=R_out, solver=slv,
                n_step=n_step, Cmax=Cmax_out, xi_at_max=xi_max_out,
                C_axis=Caxis_out, C_min=Cmin_out)


def default_dt_schedule(t_end, switch=7200.0, fine=1.0, coarse=30.0):
    """时间步长安排：热瞬变阶段用细步长，之后放宽（switch 为切换时刻）。"""
    if t_end <= switch:
        return [(t_end, fine)]
    return [(switch, fine), (t_end, coarse)]


# --------------------------------------------------------------------------
# 结果表与论文表格的组装工具
# --------------------------------------------------------------------------
RADII_REPORT_CM = [0.0, 0.5, 1.0, 1.5, 2.0]
RADII_FINE_CM = [round(0.1 * k, 1) for k in range(21)]


def result_frame(times, radii_cm, values, time_header="时间"):
    """按附件 3 模板的版式组装结果表：首列时间，其余列为各半径处的值。"""
    cols = {time_header: np.asarray(times, dtype=float)}
    for j, r in enumerate(radii_cm):
        label = "药材表面" if isinstance(r, str) else ("%.1f" % r)
        cols[label] = np.round(np.asarray(values[:, j], dtype=float), 4)
    return pd.DataFrame(cols)


# --------------------------------------------------------------------------
# 绘图环境设置（中文字体、负号、分辨率）
# --------------------------------------------------------------------------
def setup_matplotlib():
    """设置 matplotlib：中文字体、负号、字号与出图分辨率。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # 把 Microsoft YaHei 放在 SimHei 之前：mathtext 的 default/regular 字体取
    # findfont(该文本的 FontProperties)，即本列表里第一个可用字体；SimHei、
    # SimSun 都没有 U+2212（真减号），对数刻度 10^-3 的负号会被替换成占位符。
    # Microsoft YaHei 与 Noto Sans SC 都含 U+2212，故排在前面，SimHei 作为备选。
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "Noto Sans SC",
                                       "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    # 数学字体本身用 DejaVu Sans（含 U+2212），与上面一起保证负号正常
    plt.rcParams["mathtext.fontset"] = "dejavusans"
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams["font.size"] = 10.5
    plt.rcParams["axes.titlesize"] = 11.5
    plt.rcParams["axes.labelsize"] = 10.5
    plt.rcParams["legend.fontsize"] = 9
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["savefig.dpi"] = 200
    plt.rcParams["savefig.bbox"] = "tight"
    plt.rcParams["axes.grid"] = True
    plt.rcParams["grid.alpha"] = 0.3
    plt.rcParams["grid.linestyle"] = "--"
    return plt


def savefig(fig, name, plt=None):
    """按 200 dpi 把图像保存到插图目录，并在传入 plt 时关闭图像。"""
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, dpi=200)
    if plt is not None:
        plt.close(fig)
    print("   figure ->", path, flush=True)
    return path


# --------------------------------------------------------------------------
# 由模型导出的诊断量（无量纲数与特征时间）
# --------------------------------------------------------------------------
def characteristic_numbers():
    """问题 1 的无量纲数与特征时间（Biot 数、热/质扩散时间）。"""
    alpha = 0.36 / (820.0 * 2600.0)
    bi_heat = H_CONV * R0 / 0.36
    D_ref = 7.0e-9 * np.exp(-0.89 / C_INIT)
    bi_mass = HM_CONV * R0 / D_ref
    return dict(alpha=alpha, bi_heat=bi_heat, D_ref=D_ref, bi_mass=bi_mass,
                tau_heat=R0 ** 2 / alpha, tau_mass=R0 ** 2 / D_ref)
