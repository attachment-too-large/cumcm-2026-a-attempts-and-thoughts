# -*- coding: utf-8 -*-
"""按竞赛规范（format2026.doc 第九~十一条）整理：
   ① 参赛论文：单个 PDF，第一页为摘要专用页，≤20 MB
   ② 支撑材料：单个 ZIP，含全部可运行源程序 + 中间结果 + 文件列表，≤20 MB
   输出到桌面。
"""
import io
import os
import shutil
import zipfile

import fitz

ROOT = r'C:\Users\qing1\Desktop\A题'
DESK = r'C:\Users\qing1\Desktop'
STAGE = os.path.join(ROOT, '_pack', 'A题_支撑材料')
ZIP = os.path.join(DESK, '2026CUMCM_A题_支撑材料.zip')
PAPER = os.path.join(DESK, '2026CUMCM_A题_参赛论文.pdf')
SRC_PAPER = os.path.join(ROOT, 'paper_tex', 'paper.pdf')
SRC = os.path.join(ROOT, 'methodB_spectral')

CORE = ['hb_core.py', 'sem_core.py', 'run_all.py']
VERIFY = ['selftest.py', 'test_sem.py', 'compare_AB.py', 'credibility.py',
          'sens_extra.py', 's3_plateau.py', 'e2b_surface.py', 'audit_energy.py',
          'analyze_B.py', 'analyze_shrinkage.py',
          'fv2d.py', 'fv2d_run.py', 'm1_fv.py', 'm1_run.py']
TOOLS = [os.path.join(ROOT, 'paper_tools', 'make_figs.py'),
         os.path.join(ROOT, 'paper_tools', 'make_paper.py'),
         os.path.join(ROOT, 'paper_tex', 'make_tables.py'),
         os.path.join(ROOT, 'paper_tools', 'md2pdf.py')]
EVID = ['_cred.txt', '_e2c.txt', '_sens.txt', '_s3.txt', '_fv2d.txt', '_m1fv.txt']
EVID_PNG = ['_formula_p3_q1.png', '_formula_p4_q23.png', '_formula_p4_q4.png']


def cp(s, d, name=None):
    os.makedirs(d, exist_ok=True)
    if not os.path.exists(s):
        print('  MISSING', s)
        return
    shutil.copy2(s, os.path.join(d, name or os.path.basename(s)))


def scrub_pdf(src, dst):
    """复制论文并清除 PDF 元数据中的可能身份信息"""
    doc = fitz.open(src)
    doc.set_metadata({'title': '中药材热风烘干过程的材料坐标谱元建模与求解',
                      'author': '', 'subject': '2026 CUMCM A 题', 'keywords': '',
                      'creator': '', 'producer': '', 'creationDate': '', 'modDate': ''})
    doc.save(dst, garbage=4, deflate=True)
    n = len(doc)
    doc.close()
    return n


def main():
    # ---------- ① 参赛论文 ----------
    n = scrub_pdf(SRC_PAPER, PAPER)
    mb = os.path.getsize(PAPER) / 1048576
    d = fitz.open(PAPER)
    first = d[0].get_text()
    d.close()
    print('参赛论文: %s  (%d 页, %.2f MB)' % (os.path.basename(PAPER), n, mb))
    print('  首页含标题/摘要/关键词:',
          ('摘' in first) and ('关键词' in first))
    print('  满足 ≤20 MB:', mb <= 20)

    # ---------- ② 支撑材料 ----------
    if os.path.exists(STAGE):
        shutil.rmtree(STAGE)
    os.makedirs(STAGE)
    for f in CORE:
        cp(os.path.join(SRC, f), os.path.join(STAGE, '源程序', '核心'))
    for f in VERIFY:
        cp(os.path.join(SRC, f), os.path.join(STAGE, '源程序', '验证'))
    for f in TOOLS:
        cp(f, os.path.join(STAGE, '源程序', '论文与图表'))
    # 论文 LaTeX 源
    shutil.copytree(os.path.join(ROOT, 'paper_tex'), os.path.join(STAGE, '源程序', '论文与图表', 'LaTeX源'),
                    ignore=shutil.ignore_patterns('__pycache__', '*.aux', '*.log', '*.out',
                                                  '*.toc', '*.html', 'paper.pdf'))
    # 结果文件
    for f in sorted(os.listdir(os.path.join(SRC, 'results'))):
        cp(os.path.join(SRC, 'results', f), os.path.join(STAGE, '结果文件'))
    cp(os.path.join(SRC, '_compare_AB.md'), os.path.join(STAGE, '结果文件'))
    cp(os.path.join(ROOT, 'cumcm2026-a-herb-drying-main', 'outputs', 'tables', 'tables.json'),
       os.path.join(STAGE, '结果文件'), '方法A对照数据_tables.json')
    # 中间结果 / 图表
    for f in os.listdir(os.path.join(ROOT, 'paper_tools')):
        if f.endswith('.png'):
            cp(os.path.join(ROOT, 'paper_tools', f), os.path.join(STAGE, '中间结果与图表'))
    for f in EVID_PNG:
        cp(os.path.join(SRC, f), os.path.join(STAGE, '中间结果与图表'), '公式原版面_' + f[9:])
    # 实验日志
    for f in EVID:
        cp(os.path.join(SRC, f), os.path.join(STAGE, '实验日志'))
    for f in ['_A_statement.txt', '_formula_spans.txt', '_bbox.txt', '_data_full.txt']:
        cp(os.path.join(ROOT, f), os.path.join(STAGE, '实验日志', '题面抽取'))
    # 赛题原始数据（同时放到 源程序/附件 下，使代码解压后即可直接运行）
    for f in ['附件1.xlsx', '附件2.xlsx']:
        cp(os.path.join(ROOT, '附件', f), os.path.join(STAGE, '赛题原始数据'))
        cp(os.path.join(ROOT, '附件', f), os.path.join(STAGE, '源程序', '附件'))
    for f in sorted(os.listdir(os.path.join(ROOT, '附件', '附件3'))):
        cp(os.path.join(ROOT, '附件', '附件3', f), os.path.join(STAGE, '赛题原始数据', '附件3'))
        cp(os.path.join(ROOT, '附件', '附件3', f), os.path.join(STAGE, '源程序', '附件', '附件3'))
    cp(os.path.join(ROOT, 'pack_attach_readme.md'), os.path.join(STAGE, '源程序', '附件'),
       'README.md')
    # README
    io.open(os.path.join(STAGE, 'README.md'), 'w', encoding='utf-8').write(
        io.open(os.path.join(ROOT, 'pack_readme.md'), encoding='utf-8').read())

    # 文件列表
    lines = ['# 支撑材料文件列表', '',
             '> 本列表同时列于参赛论文附录 A。', '',
             '| 相对路径 | 大小 / KB |', '|---|---|']
    for dp, dns, fns in os.walk(STAGE):
        dns[:] = [x for x in dns if x != '__pycache__']
        for fn in sorted(fns):
            fp = os.path.join(dp, fn)
            lines.append('| %s | %.1f |' % (
                os.path.relpath(fp, STAGE).replace('\\', '/'), os.path.getsize(fp) / 1024))
    io.open(os.path.join(STAGE, '文件列表.md'), 'w', encoding='utf-8').write(
        '\n'.join(lines) + '\n')

    if os.path.exists(ZIP):
        os.remove(ZIP)
    with zipfile.ZipFile(ZIP, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for dp, dns, fns in os.walk(STAGE):
            dns[:] = [x for x in dns if x != '__pycache__']
            for fn in sorted(fns):
                fp = os.path.join(dp, fn)
                z.write(fp, os.path.join(os.path.basename(STAGE),
                                         os.path.relpath(fp, STAGE)))
    mb = os.path.getsize(ZIP) / 1048576
    with zipfile.ZipFile(ZIP) as z:
        bad = z.testzip()
        cnt = len(z.namelist())
    print('支撑材料: %s  (%d 个条目, %.2f MB, 完整性 %s)'
          % (os.path.basename(ZIP), cnt, mb, 'OK' if bad is None else 'BAD'))
    print('  满足 ≤20 MB:', mb <= 20)


if __name__ == '__main__':
    main()
