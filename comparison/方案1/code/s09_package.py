# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s09_package.py
# 作用  : 交付物打包。把源程序复制为附录用的 .txt 版本、同步插图目录、
#         生成支撑材料压缩包, 并输出交付清单。
#         不做任何删除操作, 只创建/覆盖本工程内部的产物文件。
# 产出  : paper/code_appendix/*.txt, paper/figures/*.png,
#         支撑材料.zip, 交付清单.txt
# 运行  : python s09_package.py
# =============================================================================
"""Package the deliverables: appendix code copies, figures, support archive."""

from __future__ import annotations

import io
import os
import shutil
import sys
import zipfile

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.join(BASE, "code")
PAPER = os.path.join(BASE, "paper")
APP = os.path.join(PAPER, "code_appendix")
FIGS = os.path.join(BASE, "figures")
PAPER_FIGS = os.path.join(PAPER, "figures")

CODE_FILES = [
    ("drying_core.py", "drying_core.txt"),
    ("scenarios.py", "scenarios.txt"),
    ("plotstyle.py", "plotstyle.txt"),
    ("s01_data_analysis.py", "s01_data_analysis.txt"),
    ("s02_solve_p1.py", "s02_solve_p1.txt"),
    ("s03_solve_p23.py", "s03_solve_p23.txt"),
    ("s04_solve_p4.py", "s04_solve_p4.txt"),
    ("s05_verification.py", "s05_verification.txt"),
    ("s06_sensitivity.py", "s06_sensitivity.txt"),
    ("s07_figures.py", "s07_figures.txt"),
    ("s08_make_tables.py", "s08_make_tables.txt"),
    ("s09_package.py", "s09_package.txt"),
    ("s10_fill_numbers.py", "s10_fill_numbers.txt"),
    ("s11_recheck_conservation.py", "s11_recheck_conservation.txt"),
    ("s12_fix_outputs.py", "s12_fix_outputs.txt"),
    ("s13_build_paper.py", "s13_build_paper.txt"),
    ("s14_inverse.py", "s14_inverse.txt"),
    ("s15_inverse_e1.py", "s15_inverse_e1.txt"),
    ("s16_fig_p4.py", "s16_fig_p4.txt"),
    ("s17_compact_paper.py", "s17_compact_paper.txt"),
]

# 附录小节标题(需转义 LaTeX 下划线)
TITLE_MAP = {
    "drying_core": "核心求解模块 drying\\_core.py",
    "scenarios": "场景与结果输出 scenarios.py",
    "plotstyle": "绘图风格 plotstyle.py",
    "s01_data_analysis": "数据探索 s01\\_data\\_analysis.py",
    "s02_solve_p1": "问题 1 求解 s02\\_solve\\_p1.py",
    "s03_solve_p23": "问题 2、3 求解 s03\\_solve\\_p23.py",
    "s04_solve_p4": "问题 4 求解 s04\\_solve\\_p4.py",
    "s05_verification": "模型验证 s05\\_verification.py",
    "s06_sensitivity": "灵敏度与不确定度 s06\\_sensitivity.py",
    "s07_figures": "论文插图 s07\\_figures.py",
    "s08_make_tables": "论文表格生成 s08\\_make\\_tables.py",
    "s09_package": "交付物打包 s09\\_package.py",
    "s10_fill_numbers": "论文数字自动回填 s10\\_fill\\_numbers.py",
    "s11_recheck_conservation": "守恒检验复算 s11\\_recheck\\_conservation.py",
    "s12_fix_outputs": "结果文件标签补正 s12\\_fix\\_outputs.py",
    "s13_build_paper": "论文编译与交付自检 s13\\_build\\_paper.py",
    "s14_inverse": "合成数据反演复算 s14\\_inverse.py",
    "s15_inverse_e1": "反演稳健估计 s15\\_inverse\\_e1.py",
    "s16_fig_p4": "问题 4 插图 s16\\_fig\\_p4.py",
    "s17_compact_paper": "正文页数压缩 s17\\_compact\\_paper.py",
}


def main():
    os.makedirs(APP, exist_ok=True)
    os.makedirs(PAPER_FIGS, exist_ok=True)

    # 1) 源程序副本(附录用; 只改扩展名, 内容逐字节一致)
    n_code = 0
    total_lines = 0
    for src, dst in CODE_FILES:
        sp = os.path.join(CODE, src)
        if not os.path.exists(sp):
            print("[WARN] 缺失源程序:", src)
            continue
        text = io.open(sp, encoding="utf-8-sig").read()
        io.open(os.path.join(APP, dst), "w", encoding="utf-8", newline="\n").write(text)
        n_code += 1
        total_lines += text.count("\n") + 1
    # 1b) 自动生成附录的源程序清单(避免手工维护 \lstinputlisting 列表时漏项)
    lines = ["% 由 s09_package.py 自动生成: 附录源程序清单, 请勿手工编辑",
             "\\section{完整源程序}", "",
             "以下给出建模与求解所使用的\\textbf{全部}源程序，内容与支撑材料中的"
             ".py 文件逐字节一致，可直接运行。",
             "运行顺序：\\texttt{s01} $\\to$ \\texttt{s02} $\\to$ \\texttt{s03} $\\to$ "
             "\\texttt{s04} $\\to$ \\texttt{s05} $\\to$ \\texttt{s06} $\\to$ "
             "\\texttt{s07}（另有 \\texttt{s08}$\\sim$\\texttt{s17} 为出表、回填、"
             "编译、打包等辅助脚本）。", ""]
    for src, dst in CODE_FILES:
        if not os.path.exists(os.path.join(CODE, src)):
            continue
        key = dst[:-4]
        lines += [f"\\subsection{{{TITLE_MAP.get(key, key)}}}",
                  f"\\lstinputlisting[language=Python]{{code_appendix/{dst}}}", ""]
    io.open(os.path.join(BASE, "paper", "sections", "s09_code_listing.tex"),
            "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    print(f"附录源程序: {n_code} 个文件, 共 {total_lines} 行")

    # 2) 插图同步
    n_fig = 0
    for f in sorted(os.listdir(FIGS)):
        if f.lower().endswith(".png"):
            shutil.copy2(os.path.join(FIGS, f), os.path.join(PAPER_FIGS, f))
            n_fig += 1
    print(f"插图同步: {n_fig} 张")

    # 3) 支撑材料压缩包
    zip_path = os.path.join(BASE, "支撑材料.zip")
    items = []
    for sub in ("code", "data", "results", "figures"):
        d = os.path.join(BASE, sub)
        for root, _dirs, files in os.walk(d):
            for f in files:
                if f.endswith(".pyc"):
                    continue
                p = os.path.join(root, f)
                items.append((p, os.path.relpath(p, BASE)))
    for extra in ("docs/自检报告.md", "docs/参考文献来源.md"):
        p = os.path.join(BASE, extra)
        if os.path.exists(p):
            items.append((p, extra))
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for p, arc in items:
            z.write(p, arc)
    size_mb = os.path.getsize(zip_path) / 1024.0 / 1024.0
    print(f"支撑材料: {zip_path} ({len(items)} 个文件, {size_mb:.2f} MB)")

    # 4) 交付清单
    lines = ["交付物清单", "=" * 70, ""]
    for sub in ("code", "data", "docs", "figures", "results", "paper"):
        d = os.path.join(BASE, sub)
        if not os.path.isdir(d):
            continue
        lines.append(f"[{sub}/]")
        for root, _dirs, files in os.walk(d):
            for f in sorted(files):
                if f.endswith((".aux", ".log", ".out", ".toc", ".pyc")):
                    continue
                p = os.path.join(root, f)
                rel = os.path.relpath(p, BASE)
                lines.append(f"  {rel:60s} {os.path.getsize(p):>10d} B")
        lines.append("")
    out = os.path.join(BASE, "交付清单.txt")
    io.open(out, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    print("交付清单:", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
