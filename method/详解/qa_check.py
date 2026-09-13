# -*- coding: utf-8 -*-
"""讲稿格式合规检查。

检查项：
  中文引号、LaTeX 脚注引号、书名号、控制字符、未闭合的美元符号行数、
  以及正文里紧跟数学模式的全角括号。
"""
import io
import os
import re
import sys

W = r"C:\Users\qing1\Desktop\bestway\详解"
FILES = ["part1.tex", "part2.tex", "part3.tex", "part4.tex", "part5.tex",
         "part6.tex", "详解.tex"]

PATS = [
    ("中文引号", re.compile("[\u201c\u201d\u2018\u2019\u300c\u300d\u300e\u300f]")),
    ("LaTeX脚注引号", re.compile("``|''")),
    ("书名号", re.compile("[\u300a\u300b]")),
    ("控制字符", re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f]")),
    ("紧贴数学的全角括号", re.compile(r"\$[\uFF08\uFF09]|[\uFF08\uFF09]\$")),
]

print("%-12s %7s %8s %8s %7s %7s %8s" %
      ("文件", "行数", "中文引号", "脚注引号", "书名号", "控制符", "紧贴括号"))
total = 0
bad = 0
for fn in FILES:
    p = os.path.join(W, fn)
    if not os.path.exists(p):
        print("%-12s (不存在)" % fn)
        continue
    text = io.open(p, encoding="utf-8").read()
    lines = text.count("\n") + 1
    total += lines
    counts = []
    for name, pat in PATS:
        hits = pat.findall(text)
        counts.append(len(hits))
        bad += len(hits)
    print("%-12s %7d %8d %8d %7d %7d %8d" %
          (fn, lines, counts[0], counts[1], counts[2], counts[3], counts[4]))

print("合计行数 %d，违例合计 %d" % (total, bad))

# 逐行报告紧贴括号的位置，便于定点修
print("")
print("紧贴数学模式的全角括号明细：")
found = 0
for fn in FILES:
    p = os.path.join(W, fn)
    if not os.path.exists(p):
        continue
    for k, line in enumerate(io.open(p, encoding="utf-8"), 1):
        if PATS[4][1].search(line):
            print("  %s:%d  %s" % (fn, k, line.strip()[:90]))
            found += 1
print("合计 %d 处" % found)
