# -*- coding: utf-8 -*-
"""把所有中文名的 .py 脚本改成英文名，并同步改掉它们内部写出的输出文件名。"""
import io
import os
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

WS = r'C:\Users\qing1\Desktop\2026CUMCM_'
SCRIPTS = os.path.join(WS, '审阅报告', '脚本')
OUTDIR = os.path.join(WS, '审阅报告', '输出')

# 旧脚本名 -> (新脚本名, 旧输出名, 新输出名)
MAP = [
    ('c01_坐标变换与谱基.py',            'c01_coord_transform_basis.py',        'c01_输出.txt',            'c01_output.txt'),
    ('c02_附件与定解条件.py',            'c02_attachments_and_bc.py',           'c02_附件与定解条件.txt',   'c02_output.txt'),
    ('c03_扩散系数与物性.py',            'c03_diffusivity_and_props.py',        'c03_扩散系数与物性.txt',   'c03_output.txt'),
    ('c03b_数量级说法核验.py',           'c03b_order_of_magnitude_claims.py',   'c03b_数量级说法核验.txt',  'c03b_output.txt'),
    ('c04_弱形式与谱基.py',              'c04_weak_form_and_basis.py',          'c04_弱形式与谱基.txt',     'c04_output.txt'),
    ('c04b_弱形式与谱基.py',             'c04b_weak_form_manufactured.py',      'c04b_弱形式与谱基.txt',    'c04b_output.txt'),
    ('c05_边界层与守恒.py',              'c05_boundary_layer_and_conservation.py', 'c05_边界层与守恒.txt',  'c05_output.txt'),
    ('c05b_问题3复核.py',                'c05b_q3_full_recheck.py',             'c05b_问题3复核.txt',       'c05b_output.txt'),
    ('c05c_守恒恒等式数值检验.py',        'c05c_conservation_identity_numeric.py', 'c05c_守恒恒等式数值检验.txt', 'c05c_output.txt'),
    ('c05d_式10温度式检验.py',           'c05d_eq10_temperature_form.py',       'c05d_式10温度式检验.txt',  'c05d_output.txt'),
    ('c06_结果表一致性.py',              'c06_result_table_consistency.py',     'c06_结果表一致性.txt',     'c06_output.txt'),
    ('c06b_模板与文字核对.py',           'c06b_template_and_text_check.py',     'c06b_模板与文字核对.txt',  'c06b_output.txt'),
    ('c07_问题4与表5表6.py',             'c07_q4_and_tables_5_6.py',            'c07_问题4与表5表6.txt',    'c07_output.txt'),
    ('c08_式13检验.py',                  'c08_eq13_model_check.py',             'c08_式13检验.txt',         'c08_output.txt'),
    ('c09b_式14核验.py',                 'c09b_eq14_model_check.py',            'c09b_式14核验.txt',        'c09b_output.txt'),
    ('c10_能量审计与二维复核.py',         'c10_energy_audit_and_2d.py',          'c10_能量审计与二维复核.txt', 'c10_output.txt'),
    ('c11_二维复核V1V2.py',              'c11_2d_v1_v2_recheck.py',             'c11_二维复核V1V2.txt',     'c11_output.txt'),
    ('c12_灵敏度与表D.py',               'c12_sensitivity_and_tableD.py',       'c12_灵敏度与表D.txt',      'c12_output.txt'),
    ('c13_结构与一致性.py',              'c13_structure_and_consistency.py',    'c13_结构与一致性.txt',     'c13_output.txt'),
    ('d01_比对Q1Q2.py',                  'd01_compare_q1_q2.py',                'd01_比对Q1Q2.txt',         'd01_output.txt'),
    ('d02_不收缩对照.py',                'd02_no_shrinkage_control.py',         'd02_不收缩对照.txt',       'd02_output.txt'),
    ('d03_问题1温度诊断.py',             'd03_q1_temperature_diagnosis.py',     'd03_问题1温度诊断.txt',    'd03_output.txt'),
    ('d04_环境重构影响.py',              'd04_ambient_reconstruction.py',       'd04_环境重构影响.txt',     'd04_output.txt'),
    ('d05_能量口径与收缩分解.py',         'd05_energy_convention_and_shrinkage.py', 'd05_能量口径与收缩分解.txt', 'd05_output.txt'),
    ('e01_环境口径交叉检验.py',           'e01_ambient_cross_test.py',           'e01_环境口径交叉检验.txt',  'e01_output.txt'),
]

done, missing = [], []
for old, new, oldout, newout in MAP:
    src = os.path.join(SCRIPTS, old)
    if not os.path.exists(src):
        missing.append(old); continue
    txt = io.open(src, encoding='utf-8').read()
    # 改内部输出文件名（可能带 ../输出/ 前缀，只替换文件名部分即可）
    txt = txt.replace(oldout, newout)
    io.open(os.path.join(SCRIPTS, new), 'w', encoding='utf-8', newline='\n').write(txt)
    os.remove(src)
    # 改输出文件
    op = os.path.join(OUTDIR, oldout)
    if os.path.exists(op):
        os.rename(op, os.path.join(OUTDIR, newout))
    done.append((old, new))

print('已重命名 %d 个脚本：' % len(done))
for o, n in done:
    print('   %-34s -> %s' % (o, n))
if missing:
    print('未找到（跳过）：', missing)

# 剩余的非 ASCII 名 .py 检查
rest = []
for base in (SCRIPTS, os.path.join(WS, '行方案详解'), os.path.join(WS, '审阅报告')):
    for root, dirs, files in os.walk(base):
        dirs[:] = [d for d in dirs if d not in ('__pycache__',)]
        for f in files:
            if f.endswith('.py') and any(ord(ch) > 127 for ch in f):
                rest.append(os.path.join(root, f))
print('\n仍含中文名的 .py:', rest if rest else '无')

# 输出目录里剩余的含中文名文件
rest2 = [f for f in os.listdir(OUTDIR) if any(ord(ch) > 127 for ch in f)]
print('输出目录仍含中文名的文件:', rest2 if rest2 else '无')
