# -*- coding: utf-8 -*-
"""二维轴对称有限体积求解器（端面效应复核用）

物理域: r in [0,R], z in [0,L/2]，z=0 取对称面（等价于整根长 L 的药材）
控制方程（保守形式，长度单位 m）:
    dC/dt      = (1/r) d_r(r D d_r C) + d_z(D d_z C)
    rho cp dT/dt = (1/r) d_r(r k d_r T) + d_z(k d_z T)
离散（单元中心有限体积，面积/体积均按 2*pi 归一）:
    体积 V_ij = Vr_i * dz,  Vr_i = (rf_{i+1}^2 - rf_i^2)/2
    径向导数面面积 Ar_kj = rf_k * dz ；轴向面面积 Az_il = rc_i * dr
    由散度定理  V dC/dt = sum_faces A_f * (D dC/dn)|_face（n 为外法向）:
      右径向面: +Ar_{i+1} * D_f/dr * (C_{i+1}-C_i)
      左径向面: +Ar_i     * D_f/dr * (C_{i-1}-C_i)     (i=0 时面积 Ar_0=0，自动零通量)
      上轴向面: +Az_{j+1} * D_f/dz * (C_{j+1}-C_j)
      下轴向面: +Az_j     * D_f/dz * (C_{j-1}-C_j)     (j=0 时 z=0 为对称面 -> 面积取 0)
    外边界 r=R / z=L/2 用二阶 ghost 消元后的等效 Robin 系数
      beta = h_m / (1 + h_m*h/(2D)),  通量入流 = A*beta*(C_inf - C_cell)
时间: theta 格式，系数取外推中点值（二阶），可做一次 Picard 修正。
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

R0, LZ = 0.02, 0.25
T0C, C0 = 28.0, 2.55
H_CONV, HM_MASS = 25.0, 8.0e-7


def build_grid(Nr, Nz, L=LZ, R=R0):
    dr, dz = R / Nr, (L / 2.0) / Nz
    rf = np.arange(Nr + 1) * dr
    rc = (np.arange(Nr) + 0.5) * dr
    Vr = 0.5 * (rf[1:] ** 2 - rf[:-1] ** 2)                  # (Nr,)
    V = np.broadcast_to(Vr[:, None], (Nr, Nz)).copy() * dz   # (Nr,Nz)
    Ar = np.broadcast_to(rf[:, None], (Nr + 1, Nz)).copy() * dz   # (Nr+1,Nz)
    Az = np.broadcast_to(rc[:, None], (Nr, Nz + 1)).copy() * dr   # (Nr,Nz+1)
    Az[:, 0] = 0.0                                           # z=0 对称面
    return dict(dr=dr, dz=dz, rc=rc, rf=rf, V=V, Ar=Ar, Az=Az, Nr=Nr, Nz=Nz)


def _div_matrix(g, coef, V):
    """组装 L，使  V d(u)/dt = L u - K(u-Yinf)（均不除体积）。
       L 的元素量纲为 m^3/s：a = A_f * coef_f / h 。coef 形状 (Nr,Nz)。"""
    Nr, Nz = g['Nr'], g['Nz']
    dr, dz = g['dr'], g['dz']
    Ar, Az = g['Ar'], g['Az']
    n = Nr * Nz
    idx = np.arange(n).reshape(Nr, Nz)
    rows, cols, vals = [], [], []

    def add(i, j, v):
        rows.append(np.asarray(i).ravel())
        cols.append(np.asarray(j).ravel())
        vals.append(np.broadcast_to(np.asarray(v, float), np.shape(i)).ravel())

    diag = np.zeros((Nr, Nz))
    # 径向内部面 k=1..Nr-1
    kf = 0.5 * (coef[:-1, :] + coef[1:, :])                  # (Nr-1,Nz)
    a = Ar[1:Nr, :] * kf / dr                                # (Nr-1,Nz)  m^3/s
    diag[:-1, :] += a
    diag[1:, :] += a
    add(idx[:-1, :], idx[1:, :], a)
    add(idx[1:, :], idx[:-1, :], a)
    # 轴向内部面 l=1..Nz-1
    kz = 0.5 * (coef[:, :-1] + coef[:, 1:])                  # (Nr,Nz-1)
    b = Az[:, 1:Nz] * kz / dz                                # (Nr,Nz-1)  m^3/s
    diag[:, :-1] += b
    diag[:, 1:] += b
    add(idx[:, :-1], idx[:, 1:], b)
    add(idx[:, 1:], idx[:, :-1], b)
    add(idx, idx, -diag)
    L = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n, n)).tocsr()
    return L


def _robin(g, coef, bnd_coef, side):
    """外边界等效 Robin 面系数 beta（通量入流 = A*beta*(Yinf - Y_cell)）"""
    if side == 'r':
        D = coef[-1, :]
        h = g['dr']
        A = g['Ar'][g['Nr'], :]
        beta = bnd_coef / (1.0 + bnd_coef * h / (2.0 * D))
        return A * beta
    else:
        D = coef[:, -1]
        h = g['dz']
        A = g['Az'][:, g['Nz']]
        beta = bnd_coef / (1.0 + bnd_coef * h / (2.0 * D))
        return A * beta


class Solver2D:
    """props: dict(rho, cp, k, D); amb: callable(t)->(Tinf,Cinf); radius: callable(t)->R"""

    def __init__(self, props, amb, radius=None, Nr=80, Nz=50, L=LZ, R=R0,
                 theta=0.5, bnd_z=True):
        self.pr, self.amb, self.rad = props, amb, radius
        self.g = build_grid(Nr, Nz, L=L, R=R)
        self.theta = theta
        self.bnd_z = bnd_z          # False: z=L/2 设为绝热（用于与一维解对拍）
        self.R = R

    def _rhs_matrices(self, C, T, Tinf, Cinf):
        """给定状态，组装 theta 半步的线性化算子（返回解算 C/T 的矩阵与右端）"""
        g, pr = self.g, self.pr
        Nr, Nz = g['Nr'], g['Nz']
        V = g['V']
        R = self.R
        D = pr['D'](C, T)
        k = pr['k'](C)
        rc = pr['rho'](C) * pr['cp'](C)
        LC = _div_matrix(g, D, V)
        LT = _div_matrix(g, k, V)

        # ---- Robin 外边界：入流 = A*beta*(Yinf - Y_cell) ----
        def assemble(bnd_coef, coef, Yinf):
            kk = np.zeros((Nr, Nz))
            dd = np.zeros((Nr, Nz))
            # r = R 面：所有 j 的最后一个径向单元；K = A*beta  (m^3/s)
            br = g['Ar'][Nr, :] * (bnd_coef / (1.0 + bnd_coef * g['dr'] / (2.0 * coef[-1, :])))
            kk[Nr - 1, :] += br
            dd[Nr - 1, :] += br * Yinf
            if self.bnd_z:
                bz = g['Az'][:, Nz] * (bnd_coef / (1.0 + bnd_coef * g['dz'] / (2.0 * coef[:, -1])))
                kk[:, Nz - 1] += bz
                dd[:, Nz - 1] += bz * Yinf
            return kk.ravel(), dd.ravel()

        kC, dC = assemble(HM_MASS, D, Cinf)
        kT, dT = assemble(H_CONV, k, Tinf)
        return dict(LC=LC, LT=LT, mC=V.ravel(), mT=(V * rc).ravel(),
                    kC=kC, kT=kT, dC=dC, dT=dT, D=D, k=k, rc=rc)

    def solve(self, t_end, dt, t_eval=(), field=False):
        g, pr = self.g, self.pr
        Nr, Nz, n = g['Nr'], g['Nz'], g['Nr'] * g['Nz']
        T = np.full((Nr, Nz), T0C)
        C = np.full((Nr, Nz), C0)
        theta = self.theta
        nt = int(np.ceil(t_end / dt))
        want = sorted(set(t_eval))
        store = {}
        wi = 0
        Tc_prev, Cc_prev = T.copy(), C.copy()
        for step in range(nt):
            t = step * dt
            # 外推中点状态
            if step == 0:
                Tm, Cm = T, C
            else:
                Tm = np.clip(1.5 * T - 0.5 * Tc_prev, 1.0, 1e4)
                Cm = np.clip(1.5 * C - 0.5 * Cc_prev, 1e-6, 1e3)
            Tinf, Cinf = self.amb(t + 0.5 * dt)
            A = self._rhs_matrices(Cm, Tm, Tinf, Cinf)
            # 组装 theta 格式:  (M/dt - th*L - th*K) Y^{n+1}
            #                 = (M/dt + (1-th)L + (1-th)*(-K?)) Y^n + ...
            # 注意 K 为边界吸收项，对 n 与 n+1 取 theta 加权
            Mdt = sp.diags(A['mC'] / dt)
            Kc = sp.diags(A['kC'])
            MC = Mdt - theta * A['LC'] + theta * Kc
            rhsC = (A['mC'] / dt) * C.ravel() \
                + (1 - theta) * (A['LC'] @ C.ravel()) \
                - (1 - theta) * (A['kC'] * C.ravel()) + A['dC']
            Cn = spla.spsolve(MC.tocsc(), rhsC).reshape(Nr, Nz)
            MT = sp.diags(A['mT'] / dt) - theta * A['LT'] + theta * sp.diags(A['kT'])
            rhsT = (A['mT'] / dt) * T.ravel() \
                + (1 - theta) * (A['LT'] @ T.ravel()) \
                - (1 - theta) * (A['kT'] * T.ravel()) + A['dT']
            Tn = spla.spsolve(MT.tocsc(), rhsT).reshape(Nr, Nz)
            Tc_prev, Cc_prev = T, C
            T, C = np.clip(Tn, -50, 200), np.clip(Cn, 0.0, 1e3)
            if not np.all(np.isfinite(T)) or not np.all(np.isfinite(C)):
                raise RuntimeError('2D 求解发散于 t=%.1f s' % t)
            while wi < len(want) and want[wi] <= (step + 1) * dt + 1e-9:
                store[want[wi]] = (float(C[0, 0]), float(T[0, 0]),
                                   float(C[-1, :].mean()), float(T[-1, :].mean()),
                                   C[:, 0].copy() if field else None,
                                   float(self.g['rc'][-1]))
                wi += 1
        return store

    def drying_time(self, dt, t_max=4.0e5, target=0.15):
        """以中心单元 C(0,0) 穿越 target 为判据，逐步推进并线性插值定位"""
        g, pr = self.g, self.pr
        Nr, Nz = g['Nr'], g['Nz']
        T = np.full((Nr, Nz), T0C)
        C = np.full((Nr, Nz), C0)
        Tp, Cp = T.copy(), C.copy()
        nt = int(np.ceil(t_max / dt))
        for step in range(nt):
            t = step * dt
            Tm = T if step == 0 else np.clip(1.5 * T - 0.5 * Tp, 1.0, 1e4)
            Cm = C if step == 0 else np.clip(1.5 * C - 0.5 * Cp, 1e-6, 1e3)
            Tinf, Cinf = self.amb(t + 0.5 * dt)
            A = self._rhs_matrices(Cm, Tm, Tinf, Cinf)
            th = self.theta
            MC = sp.diags(A['mC'] / dt) - th * A['LC'] + th * sp.diags(A['kC'])
            rhsC = (A['mC'] / dt) * C.ravel() + (1 - th) * (A['LC'] @ C.ravel()) \
                - (1 - th) * (A['kC'] * C.ravel()) + A['dC']
            Cn = spla.spsolve(MC.tocsc(), rhsC).reshape(Nr, Nz)
            MT = sp.diags(A['mT'] / dt) - th * A['LT'] + th * sp.diags(A['kT'])
            rhsT = (A['mT'] / dt) * T.ravel() + (1 - th) * (A['LT'] @ T.ravel()) \
                - (1 - th) * (A['kT'] * T.ravel()) + A['dT']
            Tn = spla.spsolve(MT.tocsc(), rhsT).reshape(Nr, Nz)
            Tp, Cp = T, C
            T, C = np.clip(Tn, -50, 200), np.clip(Cn, 0.0, 1e3)
            c0, c1 = Cp[0, 0], C[0, 0]
            if (c0 - 0.15) * (c1 - 0.15) <= 0 and c1 < c0:
                frac = (c0 - 0.15) / max(c0 - c1, 1e-30)
                return (step + frac) * dt, T, C
            if not np.isfinite(c1):
                raise RuntimeError('发散于 t=%.1f' % t)
        raise RuntimeError('未在 %.0f s 内达标' % t_max)


if __name__ == '__main__':
    import sys, os, time
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from hb_core import Ambient, load_att1, props_q23, Radius, load_att2
    t, Ta, Ca = load_att1()
    amb = Ambient(t, Ta, Ca)
    print('=== 2D 有限体积：一维极限自检（z 边界绝热 + 长柱，应等于一维解）===')
    for Nr, Nz in [(40, 8), (80, 8)]:
        s = Solver2D(props_q23(), amb, None, Nr=Nr, Nz=Nz, L=2.0, bnd_z=False)
        t0 = time.time()
        st = s.solve(3600.0, 30.0, t_eval=[600.0, 1800.0, 3600.0])
        print('  Nr=%3d Nz=%2d  C(0,0): 600s=%.6f 1800s=%.6f 3600s=%.6f  (%.1fs)'
              % (Nr, Nz, st[600.0][0], st[1800.0][0], st[3600.0][0], time.time() - t0))
