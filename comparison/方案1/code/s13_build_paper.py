# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s13_build_paper.py
# 作用  : 论文编译与交付自检。
#           (1) 用 xelatex 连续编译两次(解析交叉引用), 生成电子版与纸质版 PDF;
#           (2) 用 pypdf 提取页数与正文文本, 自动核对竞赛格式要求:
#               电子版第一页为摘要页、不含承诺书/编号页、无目录、无身份信息、
#               正文页数不超过 30 页、附录含源程序;
#           (3) 用 pdftoppm 渲染若干页为 PNG, 供人工读图检查排版。
# 产出  : paper/paper_electronic.pdf, paper/paper_print.pdf,
#         paper/build_report.txt, paper/probe/render_*.png
# 运行  : python s13_build_paper.py
# =============================================================================
"""Compile both PDFs and run the format compliance self-check."""

from __future__ import annotations

import io
import os
import re
import subprocess
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(BASE, "paper")
XELATEX = r"D:\miktex\miktex\bin\x64\xelatex.exe"
PDFTOPPM = r"D:\miktex\miktex\bin\x64\pdftoppm.exe"

FORBIDDEN = ["承诺书", "编号专用页", "赛区评阅", "报名参赛队号", "参赛学校",
             "指导教师", "目  录", "目 录", "tableofcontents"]


def compile_tex(tex, passes=2):
    """连续编译 tex 文件, 返回 (成功与否, 日志尾部)。"""
    log_tail = ""
    for _ in range(passes):
        p = subprocess.run([XELATEX, "-interaction=nonstopmode", tex],
                           cwd=PAPER, capture_output=True)
        out = p.stdout.decode("utf-8", errors="replace")
        log_tail = "\n".join(out.splitlines()[-40:])
    pdf = os.path.join(PAPER, os.path.splitext(tex)[0] + ".pdf")
    return os.path.exists(pdf), log_tail


def main():
    lines = []
    push = lines.append
    push("=" * 78)
    push("论文编译与交付自检")
    push("=" * 78)

    results = {}
    for tex in ("paper_electronic.tex", "paper_print.tex"):
        ok, tail = compile_tex(tex, passes=2)
        results[tex] = ok
        push(f"[编译] {tex}: {'成功' if ok else '失败'}")
        errs = [l for l in tail.splitlines() if l.startswith("! ")]
        if errs:
            push("        报错行: " + " | ".join(errs[:5]))

    from pypdf import PdfReader

    for tex, pdf in (("paper_electronic.tex", "paper_electronic.pdf"),
                     ("paper_print.tex", "paper_print.pdf")):
        p = os.path.join(PAPER, pdf)
        if not os.path.exists(p):
            push(f"[检查] {pdf} 不存在, 跳过")
            continue
        try:
            r = PdfReader(p)
        except Exception as e:                       # noqa: BLE001
            push(f"[检查] {pdf} 读取失败: {e}")
            continue
        n = len(r.pages)
        size_mb = os.path.getsize(p) / 1024.0 / 1024.0
        push("")
        push("-" * 78)
        push(f"[{pdf}] 页数 = {n}, 大小 = {size_mb:.2f} MB")
        # 第 1 页文本
        p1 = r.pages[0].extract_text() or ""
        head = p1.strip().splitlines()[:3]
        push(f"  第 1 页开头: {head}")
        push(f"  第 1 页含'摘要' = {'摘要' in p1}; "
             f"含'承诺书' = {'承诺书' in p1}; 含'编号专用页' = {'编号专用页' in p1}")
        # 全文扫描
        full = []
        for pg in r.pages:
            try:
                full.append(pg.extract_text() or "")
            except Exception:                        # noqa: BLE001
                full.append("")
        text = "\n".join(full)
        for kw in FORBIDDEN:
            if kw in text:
                push(f"  [警告] 全文出现敏感词/禁止项: {kw}")
        # 章节完整性
        need = ["问题重述", "问题分析", "模型假设", "符号说明", "模型建立",
                "模型求解", "模型检验", "模型评价", "参考文献", "附录"]
        miss = [k for k in need if k not in text]
        push(f"  必需章节缺失 = {miss if miss else '无'}")
        push(f"  含'drying_core'源程序 = {'drying_core' in text}")
        push(f"  页码/字数: 提取字符数 = {len(text)}")
        # 第二页(电子版应为正文开始)
        if n > 1 and pdf == "paper_electronic.pdf":
            p2 = (r.pages[1].extract_text() or "").strip().splitlines()[:2]
            push(f"  电子版第 2 页开头: {p2}")

    # 渲染若干页供读图检查
    os.makedirs(os.path.join(PAPER, "probe"), exist_ok=True)
    pdf = os.path.join(PAPER, "paper_electronic.pdf")
    if os.path.exists(pdf):
        for page in (1, 2, 3):
            subprocess.run([PDFTOPPM, "-png", "-r", "110", "-f", str(page),
                            "-l", str(page), pdf,
                            os.path.join(PAPER, "probe", f"final_p{page:02d}")],
                           capture_output=True)
        push("")
        push("已渲染 paper/probe/final_p01..p03.png 供读图检查")

    text = "\n".join(lines)
    with io.open(os.path.join(PAPER, "build_report.txt"), "w",
                 encoding="utf-8", newline="\n") as f:
        f.write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
