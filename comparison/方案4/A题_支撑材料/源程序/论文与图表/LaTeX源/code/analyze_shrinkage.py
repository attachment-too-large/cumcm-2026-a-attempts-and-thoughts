# -*- coding: utf-8 -*-
"""附件2 收缩数据 与 模型水分历史 的一致性分析（独立于工作区成果）"""
import numpy as np
import openpyxl
from hb_core import load_att1, load_att2, Ambient, Radius, props_q1, props_q23, props_q4, R0
from sem_core import SEM1D, Model

OUT = 'results'


def load_result_C(path, sheet='水分浓度'):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    t = np.array([r[0] for r in rows[1:]], float)
    M = np.array([[np.nan if v is None else v for v in r[1:]] for r in rows[1:]], float)
    return t, M


def mean_C(t, M, sig=None):
    """体积加权平均 Cbar = 2*int C sigma dsigma"""
    n = M.shape[1]
    s = np.linspace(0, 1, n) if sig is None else sig
    w = s.copy()
    w[1:-1] *= 2
    return (M * w).sum(axis=1) / w.sum() * 1.0


def main():
    t2, R2 = load_att2()
    t1, T1, C1 = load_att1()
    amb = Ambient(t1, T1, C1)
    S = SEM1D((0., 0.5, 0.8, 0.95, 1.), 16)
    sig = np.round(np.arange(0, 20.0001, 1.0) / 10.0, 4) / R0 / 100.0

    print('== 收缩-水分一致性的三种检验 ==')
    # (A) 问题1 参数下前 1800 s 的水分损失
    m = Model(S, props_q1(), amb)
    sol = m.solve(1800.0, t_eval=[0., 300., 600., 900., 1200., 1500., 1800.], rtol=1e-10, atol=1e-12)
    Cb1 = [mean_C(sol.t, np.array([S.interp(sol.y[S.ndof:, k], sig) for k in range(len(sol.t))]))[k]
           for k in range(len(sol.t))]
    print('  Q1 参数: Cbar(1800s)=%.4f  (初始 2.55, 失水 %.4f kg/kg=%.1f%%)'
          % (Cb1[-1], 2.55 - Cb1[-1], (2.55 - Cb1[-1]) / 2.55 * 100))
    # 附件2 在 1800 s 的半径 -> 由体积加和律反推 Cbar
    Rd = R2[-1]
    a = Rd ** 2 / 4.0                      # R^2 = 4[a + (1-a)C/C0]
    def C_from_R(R):
        return (R ** 2 / 4.0 - a) / (1.0 - a) * 2.55
    print('  附件2 半径在 1800 s: R=%.4f cm -> 体积加和律反推 Cbar=%.4f' % (R2[1], C_from_R(R2[1])))
    print('  附件2 半径在 259200 s: R=%.4f -> Cbar=%.4f (完全干燥极限 %.4f)'
          % (Rd, C_from_R(Rd), 0.0))

    # (B) 问题3 参数：把模型计算的 Cbar(t) 映射到体积加和律下的 R(t)，与附件2 比较
    print('  --- 用模型 Cbar(t) 预测 R(t)（体积加和律, 干态半径 1.198 cm）---')
    for tag, pr, tf in [('Q1参数', props_q1(), 14400.), ('Q2/3参数', props_q23(), 14400.)]:
        mm = Model(S, pr, amb)
        s2 = mm.solve(tf, t_eval=t2[t2 <= tf], rtol=1e-9, atol=1e-11)
        Cb = mean_C(s2.t, np.array([S.interp(s2.y[S.ndof:, k], sig) for k in range(len(s2.t))]))
        Rpred = 2.0 * np.sqrt(a + (1 - a) * np.clip(Cb, 0, None) / 2.55)
        Rdat = np.interp(s2.t, t2, R2)
        print('   %s  t/s       Cbar     R模型    R附件2   差' % tag)
        for k in range(0, len(s2.t), max(1, len(s2.t) // 8)):
            print('        %8.0f  %7.4f  %7.4f  %7.4f  %+.4f'
                  % (s2.t[k], Cb[k], Rpred[k], Rdat[k], Rpred[k] - Rdat[k]))


if __name__ == '__main__':
    main()
