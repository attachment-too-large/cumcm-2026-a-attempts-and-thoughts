# -*- coding: utf-8 -*-
"""把方法 B（独立解法）的全部产出打包成 zip（含附件，解压即可运行）"""
import io
import os
import shutil
import zipfile

ROOT = r'C:\Users\qing1\Desktop\A题'
STAGE = os.path.join(ROOT, '_pack', 'CUMCM2026A_方法B_独立解法')
ZIP = os.path.join(ROOT, 'CUMCM2026A_方法B_独立解法.zip')

PLAN = [
    ('方法B_独立求解与交叉比对报告.md', '论文与报告'),
    ('方法B_独立求解与交叉比对报告.pdf', '论文与报告'),
    ('论文_中药材烘干_材料坐标谱元法.md', '论文与报告'),
    ('论文_中药材烘干_材料坐标谱元法.pdf', '论文与报告'),
    ('论文_中药材烘干_材料坐标谱元法_TeX版.pdf', '论文与报告'),
    ('paper_tex/paper.tex', '论文与报告/LaTeX源'),
    ('paper_tex/patch_tex.py', '论文与报告/LaTeX源'),
    ('paper_tex/patch_tex2.py', '论文与报告/LaTeX源'),
    ('paper_tex/patch_tex3.py', '论文与报告/LaTeX源'),
    ('paper_tex/patch_tex4.py', '论文与报告/LaTeX源'),
    ('paper_tex/patch_tex5.py', '论文与报告/LaTeX源'),
    ('paper_tex/fix_format.py', '论文与报告/LaTeX源'),
    ('paper_tex/fix_sens_table.py', '论文与报告/LaTeX源'),
    ('paper_tex/make_tables.py', '论文与报告/LaTeX源'),
    ('paper_tools/md2pdf.py', '工具/md2pdf'),
    ('paper_tools/make_figs.py', '工具'),
    ('paper_tools/make_paper.py', '工具'),
    ('附件/附件1.xlsx', '附件'),
    ('附件/附件2.xlsx', '附件'),
    ('附件/附件3/result1.xlsx', '附件/附件3'),
    ('附件/附件3/result2.xlsx', '附件/附件3'),
    ('附件/附件3/result3.xlsx', '附件/附件3'),
    ('附件/附件3/result4.xlsx', '附件/附件3'),
]
DIRS = [('paper_tex/tabs', '论文与报告/LaTeX源/tabs'),
        ('paper_tex/figs', '论文与报告/LaTeX源/figs'),
        ('paper_tex/code', '论文与报告/LaTeX源/code'),
        ('methodB_spectral/results', '结果文件')]
EXTRA = [
    ('cumcm2026-a-herb-drying-main/outputs/tables/tables.json',
     '结果文件/逐格比对/方法A_对照数据_tables.json'),
]
EVIDENCE = ['_cred.txt', '_e2c.txt', '_sens.txt', '_s3.txt', '_fv2d.txt', '_m1fv.txt',
            '_formula_p3_q1.png', '_formula_p4_q23.png', '_formula_p4_q4.png',
            '_out_test_sem.txt']
PDF_EVID = ['_A_statement.txt', '_bbox.txt', '_formula_spans.txt', '_data_full.txt']


def copy(src, dst_dir):
    if not os.path.exists(src):
        print('  MISSING', src)
        return 0
    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(src, dst_dir)
    return 1


def main():
    if os.path.exists(STAGE):
        shutil.rmtree(STAGE)
    os.makedirs(STAGE)
    n = 0
    for s, d in PLAN:
        n += copy(os.path.join(ROOT, s), os.path.join(STAGE, d))
    for s, d in DIRS:
        sp, dp = os.path.join(ROOT, s), os.path.join(STAGE, d)
        if os.path.isdir(sp):
            shutil.copytree(sp, dp, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('__pycache__', '*.aux', '*.log',
                                                          '*.out', '*.toc', '*.html'))
            n += 1
    for s, d in EXTRA:
        n += copy(os.path.join(ROOT, s), os.path.dirname(os.path.join(STAGE, d)))
    # 代码
    code_dir = os.path.join(STAGE, '代码')
    for f in sorted(os.listdir(os.path.join(ROOT, 'methodB_spectral'))):
        if f.endswith('.py'):
            n += copy(os.path.join(ROOT, 'methodB_spectral', f), code_dir)
    # 证据
    for f in EVIDENCE:
        n += copy(os.path.join(ROOT, 'methodB_spectral', f), os.path.join(STAGE, '证据/实验日志'))
    for f in PDF_EVID:
        n += copy(os.path.join(ROOT, f), os.path.join(STAGE, '证据/题面抽取'))
    # 让 compare_AB.py 能在包内找到对照数据（<代码>/results/对照数据/tables/tables.json）
    ref_dst = os.path.join(STAGE, '代码', 'results', '对照数据', 'tables')
    copy(os.path.join(ROOT, 'cumcm2026-a-herb-drying-main', 'outputs', 'tables',
                      'tables.json'), ref_dst)
    io.open(os.path.join(os.path.dirname(ref_dst), '说明.md'), 'w',
            encoding='utf-8').write(
        '# 方法 A 对照数据\n\n'
        '本目录下的 `tables/tables.json` 是**工作区已有实现（方法 A）**的表 1—表 6 结果，\n'
        '仅用于 `compare_AB.py` 做逐格比对，版权归其作者。\n'
        '本包不含方法 A 的源程序；如需重新生成该对照数据，请设置环境变量\n'
        '`CUMCM_WS_OUT` 指向方法 A 工程的 `outputs` 目录。\n')
    print('copied items:', n)

    # README 与附件说明
    io.open(os.path.join(STAGE, 'README.md'), 'w', encoding='utf-8').write(
        io.open(os.path.join(ROOT, 'pack_readme.md'), encoding='utf-8').read())
    io.open(os.path.join(STAGE, '附件', 'README.md'), 'w', encoding='utf-8').write(
        io.open(os.path.join(ROOT, 'pack_attach_readme.md'), encoding='utf-8').read())

    # 文件清单
    lines = ['# 文件清单', '', '| 相对路径 | 大小 / KB |', '|---|---|']
    for dirpath, dirnames, filenames in os.walk(STAGE):
        dirnames[:] = [d for d in dirnames if d != '__pycache__']
        for fn in sorted(filenames):
            if fn == '文件清单.md':
                continue
            fp = os.path.join(dirpath, fn)
            rel = os.path.relpath(fp, STAGE).replace('\\', '/')
            lines.append('| %s | %.1f |' % (rel, os.path.getsize(fp) / 1024))
    io.open(os.path.join(STAGE, '文件清单.md'), 'w', encoding='utf-8').write(
        '\n'.join(lines) + '\n')

    if os.path.exists(ZIP):
        os.remove(ZIP)
    with zipfile.ZipFile(ZIP, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for dirpath, dirnames, filenames in os.walk(STAGE):
            dirnames[:] = [d for d in dirnames if d != '__pycache__']
            for fn in sorted(filenames):
                fp = os.path.join(dirpath, fn)
                z.write(fp, os.path.join(os.path.basename(STAGE),
                                         os.path.relpath(fp, STAGE)))
    print('zip:', ZIP, '%.2f MB' % (os.path.getsize(ZIP) / 1024 / 1024))
    with zipfile.ZipFile(ZIP) as z:
        print('完整性:', 'OK' if z.testzip() is None else 'BAD', '| 条目', len(z.namelist()))


if __name__ == '__main__':
    main()
