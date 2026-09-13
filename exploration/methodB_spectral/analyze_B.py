# -*- coding: utf-8 -*-
"""方法B 独立分析：
 (1) 能量审计：表面对流供热 vs 蒸发潜热需求（复核工作区 §5.5 的 70 倍结论）
 (2) 体积加和律下的收缩一致性：附件2 R(t) <-> 模型 Cbar(t)
 (3) 判据对照：中心(最大) 判据 vs 平均含水率判据
"""
import numpy as np
from hb_core import (Ambient, Radius, load_att1, load_att2, props_q1, props_q23,
                     props_q4, R0, LZ, H_CONV, HM_MASS, T0_C, C0)
from sem_core import SEM1D, Model

RHO_W = 1000.0
L_LAT = 2.4e6


def mean_fields(model, y, t):
    """返回 (Cbar 体积平均, Cs, Ts, rho_dry 体积平均干密度)"""
    S = model.S
    sig, w = S.global_quad()
    T = S.interp(y[:S.ndof], sig)
    C = S.interp(y[S.ndof:], sig)
    vol = w.sum()
    Cbar = float((C * w).sum() / vol)
    rho = model.pr['rho'](C)
    return Cbar, float(C[S.surf_dof]), float(T[S.surf_dof]), rho


def audit(props, radius, tag, t_end, mesh=None, rtol=1e-9, atol=1e-11):
    S = SEM1D(**(mesh or dict(breaks=(0., .6, .9, 1.), n=16)))
    t, T, C = load_att1()
    amb = Ambient(t, T, C)
    m = Model(S, props, amb, radius)
    sol = m.solve(t_end, t_eval=np.linspace(0, t_end, 400), rtol=rtol, atol=atol)
    Area0 = 2 * np.pi * R0 * LZ                     # 初始表面积 (m^2)
    Qconv = 0.0
    mevap = 0.0
    for k in range(len(sol.t) - 1):
        tm = 0.5 * (sol.t[k] + sol.t[k + 1])
        dt = sol.t[k + 1] - sol.t[k]
        ym = sol.sol(tm)
        Tinf, Cinf = amb(tm)
        R = float(radius(tm)) if radius is not None else R0
        A = 2 * np.pi * R * LZ
        Ts = ym[S.surf_dof]
        Cs = ym[S.ndof + S.surf_dof]
        # 单位长度/整根的换算是同一比例，全部按"整根"计
        Qconv += A * H_CONV * (Tinf - Ts) * dt
        mevap += A * HM_MASS * (Cs - Cinf) * dt / 1000.0 * 1000.0   # 体积浓度差 x m/s -> 见下
    # 水的质量流量 [kg/s] = A * km * rho_dry_bulk * (Cs - Cinf)  (rho 抵消口径见文档)
    # 这里按 干基: dm_w/dt = A * km * (Cs-Cinf) * rho_dry(参考)  ; 取 rho_dry 初值
    return sol, m, S, amb


def main():
    t1, T1, C1 = load_att1()
    t2, R2 = load_att2()
    amb = Ambient(t1, T1, C1)
    S = SEM1D(breaks=(0., .6, .9, 1.), n=16)

    print('== (1) 能量审计（单位: 整根药材, 长度 0.25 m）==')
    cases = [('问题1(0-1800s)', props_q1(), None, 1800.0),
             ('问题3(0-57h)', props_q23(), None, 205809.7),
             ('问题4(0-51h)', props_q4(), Radius(fixed=False, t=t2, R=R2), 182950.0)]
    for tag, pr, rad, tend in cases:
        m = Model(S, pr, amb, rad)
        sol = m.solve(tend, t_eval=np.linspace(0, tend, 600), rtol=1e-8, atol=1e-10)
        Qconv = m_evap = 0.0
        for k in range(len(sol.t) - 1):
            tm = 0.5 * (sol.t[k] + sol.t[k + 1]); dt = sol.t[k + 1] - sol.t[k]
            ym = sol.sol(tm)
            Tinf, Cinf = amb(tm)
            R = float(rad(tm)) if rad is not None else R0
            A = 2 * np.pi * R * LZ
            Qconv += A * H_CONV * (Tinf - ym[S.surf_dof]) * dt
            m_evap += A * HM_MASS * (ym[S.ndof + S.surf_dof] - Cinf) * dt
        # 干物质质量（体积加和取干密度初值近似）
        y0 = np.concatenate([np.full(S.ndof, T0_C), np.full(S.ndof, C0)])
        Cb0, _, _, _ = mean_fields(m, y0, 0.0)
        Cb1, Cs1, Ts1, _ = mean_fields(m, sol.y[:, -1], sol.t[-1])
        Wtot = m_evap                                   # kg (干基浓度差 x km x A 积分)
        E_lat = Wtot * L_LAT
        print('  %-16s Q_conv=%.3e J  m_evap=%.4f kg(积分) / %.4f kg(Cbar差)  E_lat=%.3e J  比 %.1f 倍'
              % (tag, Qconv, Wtot, (Cb0 - Cb1), E_lat, E_lat / Qconv))
    return


if __name__ == '__main__':
    main()
