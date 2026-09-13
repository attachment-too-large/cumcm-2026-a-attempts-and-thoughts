# -*- coding: utf-8 -*-
# ==============================================================================
# common.py —— 公共参数、烘房环境标定与物性经验式（四个问题共用）
# ==============================================================================
# 做什么：集中定义四问共用的几何尺寸、传递系数与初值；把附件 1 的 241 组实测环境
#         数据标定成一阶惯性式 T_air(t)/C_air(t)；按题目附录 2/3/4 给出 rho、cp、
#         k、D 随含水率与温度变化的经验关系。主求解器与各验证脚本都从这里取物性，
#         故论文里写的物性与代码里算的是同一份。
# 输入：附件 1.xlsx（"时间/温度/水分浓度"实测序列，241 行）、附件 2.xlsx（"时间/半径"
#       序列，半径列单位 cm，仅问题 4 使用）；附录 2/3/4 的物性参数直接写在下面。
#       附件由 _find_attach() 在包内候选目录中按文件名搜索，故整个包可整体搬移。
# 输出：不写文件。对外提供 R0/Lc/h/hm/T0/C0 等常量，以及 T_air、C_air、props、
#       R_of、bessel_roots 五个函数。
# 关键变量：
#   R0    药材初始半径                  2.000e-2    [m]
#   Lc    药材长度                      0.25        [m]
#   h     对流换热系数                  25.0        [W/(m2 K)]
#   hm    对流传质系数                  8.0e-7      [m/s]
#   T0    药材初温                      28.0        [degC]
#   C0    药材初含水率（干基）           2.55        [kg/kg]
#   Tset  烘房温度惯性式稳态值          50.212      [degC]
#   tauT  烘房温度惯性时间常数          1877.5      [s]
#   Cset  烘房空气含湿量稳态值          0.050908    [kg/kg]
#   tauC  烘房空气含湿量时间常数        2776.3      [s]
#   T_END 惯性段/恒温段分段时刻 14400 s（附件 1 采样终点 4 h）    [s]
# ==============================================================================
import os
import numpy as np
import openpyxl
from scipy.optimize import curve_fit
from scipy.special import j0, j1
from scipy.optimize import brentq

# ---------------- 数据文件定位（可移植：相对本文件搜索） ----------------
# 为什么这么写：交付包会被解压到任意路径（且可能被整体移动），写死绝对路径必然失效。
# 候选顺序是"包根的题目/ → 包根 → 上一级 → 本目录"，覆盖压缩包原样解压与
# 只取 code/ 单跑两种用法。
_HERE = os.path.dirname(os.path.abspath(__file__))
_CANDIDATES = [
    os.path.join(_HERE, "..", "题目"),
    os.path.join(_HERE, ".."),
    os.path.join(_HERE, "..", ".."),
    _HERE,
]


def _find_attach(name):
    """在候选目录下查找 附件/<name>，返回绝对路径；全部落空则报错并列出搜索过的目录。

    宁可直接抛 FileNotFoundError 也不要返回 None：附件缺失时后续拟合会静默读到
    错误数据，报错比出错数好查。
    """
    for base in _CANDIDATES:
        p = os.path.abspath(os.path.join(base, "附件", name))
        if os.path.isfile(p):
            return p
    raise FileNotFoundError(
        "找不到 附件/%s；已搜索：%s" % (name, [os.path.abspath(b) for b in _CANDIDATES]))


# ---------------- 几何与传递系数（附录 2） ----------------
# 圆柱侧面对流、两端绝热的假设下，只有半径 R0 进入方程，长度 Lc 仅用于报告体积。
R0 = 0.02         # m   药材半径
Lc = 0.25         # m   药材长度
h  = 25.0         # W/(m^2*K)   对流换热系数
hm = 8.0e-7       # m/s         对流传质系数

# ---------------- 初值 ----------------
T0 = 28.0         # degC
C0 = 2.55         # kg/kg (干基)

# ---------------- 附件1：烘房环境，一阶惯性拟合 ----------------
def _load_env():
    """读附件 1 的前三列，返回 (时间 [s], 烘房温度 [degC], 空气含湿量 [kg/kg])。"""
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


# 烘房升温/加湿的物理过程是一阶惯性（送风量有限、围护结构有热容），故用
# phi(t) = phi_set - (phi_set - phi(0))*exp(-t/tau) 拟合，只标定三个量：
# 稳态值、时间常数与初值。初值取实测首点，不作为自由参数。
def _fT(t, Tset, tau):
    """温度惯性式（供 curve_fit 调用）：t [s]，返回烘房温度 [degC]。"""
    return Tset - (Tset - T_air0) * np.exp(-t / tau)


def _fC(t, Cset, tau):
    """含湿量惯性式（供 curve_fit 调用）：t [s]，返回空气含湿量 [kg/kg]。"""
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

# 拟合残差有两种报告口径，看数时别混：论文表 1 与 run_all.py 报的是 RMSE 用
# ÷n 的 0.336 degC / 8.1e-4 kg/kg（把标定值当作对整段序列的整体描述）；
# curve_fit 的 pcov 走的是 ÷(n-2) 的自由度口径，只在需要参数标准误时使用。


# 工艺分段：附件 1 在 t 约 4 h 达到设定值并稳定，故 t>14400 s 直接取设定值。
# 注意该分段在 t=14400 s 处不连续：一阶惯性式在该点的值与设定值尚差
#   T_air: 1.037e-2 degC，C_air: 1.748e-4 kg/kg
# 这是"实测已达稳态、故按稳态取值"这一建模选择的结果，论文 S3.5 已注明量级；
# 若需要严格连续，应把分段点后移或对惯性式做尾部截断，两者都会改动 t_dry 的第 4 位有效数字。
def T_air(t):
    """烘房温度 [degC]：t<=14400 s 取一阶惯性式，其后取稳态设定值。

    入参 t 可为标量或数组；返回同形状数组（后续用 numpy 广播，标量入参也返回 0 维数组）。
    """
    t = np.asarray(t, dtype=float)
    ramp = Tset - (Tset - T_air0) * np.exp(-t / tauT)
    return np.where(t <= T_END, ramp, Tset)


def C_air(t):
    """烘房空气含湿量 [kg/kg]：分段口径与 T_air 完全相同（同一时刻同一个 T_END）。"""
    t = np.asarray(t, dtype=float)
    ramp = Cset - (Cset - C_air0) * np.exp(-t / tauC)
    return np.where(t <= T_END, ramp, Cset)


# ---------------- 物性经验式 ----------------
# 浓度下限：D 含 exp(-a/C)，C→0 时该因子发散。物理上表面含水率也确实不会降到 0
# （烘干判据是 0.15 kg/kg），故把求值用的 C 截断到 1e-3，避免迭代中出现 inf/NaN。
C_FLOOR = 1.0e-3


def props(prob, C, T):
    """返回 (rho, cp, k, D)，四个量与入参 C、T 同形状。

    prob 决定用哪套物性（对应题目三个附录）：
      1 -> 附录 2：常数物性（问题 1 用，温度场与水分场解耦）
      2、3 -> 附录 3：含 C、T 的强非线性关系（问题 2、3 用）
      4 -> 附录 4：另一套系数（问题 4 用，配合实测收缩）

    量纲：C 为干基含水率 [kg/kg]；T 为温度 [degC]，仅 Arrhenius 项内部换算成 K。
    返回：rho [kg/m3]、cp [J/(kg K)]、k [W/(m K)]、D [m2/s]。
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
    """读附件 2，返回 (时间 [s], 半径 [m])；附件以 cm 记录，此处统一乘 1e-2 换成 m。"""
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
def bessel_roots(Bi, n):
    """求 beta*J1(beta) = Bi*J0(beta) 的前 n 个正根（升序）。

    Bi = h*R0/k 为传热 Biot 数；递推关系 J0'=-J1 保证该式正是 Robin 边界条件
    k*dT/dr = -h*T 分离变量后的特征方程。求根用"先用 J1 的零点把实轴分段，
    再在每段内用 brentq 收缩"的办法——不依赖初值猜测，且不会漏根或串根。
    """
    from scipy.special import jnp_zeros
    z = jnp_zeros(1, n + 2)                    # J1 的正零点
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
