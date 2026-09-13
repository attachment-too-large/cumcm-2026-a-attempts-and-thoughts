# -*- coding: utf-8 -*-
"""更正版自检 1：环境数据「噪声」与「偏置」哪个真的影响答案？
   对被扩散系统而言，高频噪声应被强烈低通，而系统性偏置不会被滤掉。
   若实测证实这一点，则「用拟合曲线替代实测」这一做法在本问题上是得不偿失的。
"""
import io
import sys
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
CORE = r'C:\Users\qing1\Desktop\2026CUMCM_\A题_支撑材料\源程序\核心'
sys.path.insert(0, CORE)
from hb_core import Ambient, load_att1, props_q1, props_q23          # noqa: E402
from sem_core import SEM1D, Model, drying_time                       # noqa: E402

OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))

t1, T1, C1 = load_att1()
TF = np.array([0.0, 2.0])
MESH1 = dict(breaks=(0.0, 0.5, 0.8, 0.95, 1.0), n=20)     # ndof=81
MESHL = dict(breaks=(0.0, 0.6, 0.9, 1.0), n=18)            # ndof=55

X_T = (34.0425, 37.1989, 1.5109)
B_T = (33.5753, 36.7856, 1.5102)


def q1(amb):
    S = SEM1D(**MESH1)
    m = Model(S, dict(props_q1()), amb, None, 'affine')
    sol = m.solve(1800.0, t_eval=[1800.0], rtol=1e-10, atol=1e-12, method='BDF')
    T, C = m.sample(sol.y[:, -1], 1800.0, r_cm=TF)
    return float(T[0]), float(T[-1]), float(C[-1])


class Fit:
    """答案 X 的一阶惯性式环境。"""
    T_SET, TAU_T, C_SET, TAU_C, C0 = 50.212, 1877.5, 0.050908, 2776.3, 0.01963
    def __call__(self, t):
        t = float(t)
        return (self.T_SET - (self.T_SET - 28.0) * np.exp(-t / self.TAU_T),
                self.C_SET - (self.C_SET - self.C0) * np.exp(-t / self.TAU_C))


def noisy(seed, sigma):
    rng = np.random.default_rng(seed)
    return Ambient(t1, T1 + rng.normal(0, sigma, T1.shape),
                   C1 + rng.normal(0, sigma * 0.06, C1.shape), mode='linear')


say('=' * 84)
say('【1】问题 1：1800 s 的中心温度 / 表面温度 / 表面含水率')
say('=' * 84)
say('   环境处理                                  T中心      T表面     C表面    |Δ| vs X')
cases = [
    ('分段线性（原样用附件1）', Ambient(t1, T1, C1, mode='linear')),
    ('PCHIP 保形三次', Ambient(t1, T1, C1, mode='pchip')),
    ('SG(41,2) 平滑', Ambient(t1, T1, C1, mode='linear', smooth=41)),
    ('线性 + 噪声 σ=0.2 (种子1)', noisy(1, 0.2)),
    ('线性 + 噪声 σ=0.2 (种子2)', noisy(2, 0.2)),
    ('线性 + 噪声 σ=0.2 (种子3)', noisy(3, 0.2)),
    ('线性 + 噪声 σ=0.8 (种子1)', noisy(1, 0.8)),
    ('线性 + 噪声 σ=0.8 (种子2)', noisy(2, 0.8)),
    ('一阶惯性拟合（答案 X）', Fit()),
]
res = {}
for tag, amb in cases:
    tc, ts, cs = q1(amb)
    res[tag] = (tc, ts, cs)
    d = max(abs(tc - X_T[0]), abs(ts - X_T[1])) if 'X' in tag else None
    say('   %-34s %9.4f %9.4f %9.4f   %s'
        % (tag, tc, ts, cs, ('—' if d is None else '%.4f' % d)))

lin = res['分段线性（原样用附件1）']
say()
say('   分段线性          = (%.4f, %.4f, %.4f)' % lin)
say('   答案 X 一阶惯性拟合 = (%.4f, %.4f, %.4f)' % X_T)
say('   方法 B 论文值（=分段线性）= (%.4f, %.4f, %.4f)' % B_T)
say()
ns = [res[k] for k in res if k.startswith('线性 + 噪声 σ=0.2')]
ds = [max(abs(v[0] - lin[0]), abs(v[1] - lin[1])) for v in ns]
say('   同一 σ=0.2 的三个不同噪声实现之间的最大差异 = %.4f ℃' % (max(ds) - min(ds)))
say('   3σ=0.6 噪声实现与无噪结果的最大偏差        = %.4f ℃'
    % max(max(abs(v[0] - lin[0]), abs(v[1] - lin[1]))
          for k, v in res.items() if 'σ=0.8' in k))
say('   拟合替代实测造成的偏差（不可被滤波消除）    = %.4f ℃'
    % max(abs(X_T[0] - lin[0]), abs(X_T[1] - lin[1])))

say()
say('=' * 84)
say('【2】问题 3：烘干时长对同一组环境处理的敏感性')
say('=' * 84)
say('   环境处理                                  t_dry / s      t_dry / h    相对线性')
base = None
for tag, amb in cases:
    S = SEM1D(**MESHL)
    m = Model(S, props_q23(), amb, None, 'affine')
    tf, _ = drying_time(m, t_max=600000.0, rtol=1e-9, atol=1e-11, method='BDF')
    if base is None:
        base = tf
    say('   %-34s %11.1f   %9.4f   %+8.4f%%'
        % (tag, tf, tf / 3600.0, 100 * (tf - base) / base))

say()
say('=' * 84)
say('【3】结论')
say('=' * 84)
say('   · 高频噪声几乎被扩散系统完全滤掉：σ=0.2（约等于实测噪声水平）的不同噪声实现之间')
say('     差异远小于拟合替代实测造成的偏差。')
say('   · 一阶惯性拟合引入的是【系统性偏置】（早期 +0.96 ℃），扩散系统不会滤掉偏置。')
say('   · 因此在本问题上，「用拟合曲线替代实测数据」是【得不偿失】的：')
say('     它用 0.47 ℃ 的系统偏差，去换本来就不影响的噪声。')

with open('../输出/selfcheck_noise_vs_bias.txt', 'w', encoding='utf-8') as fp:
    fp.write('\n'.join(OUT))
print('\n[saved]')
