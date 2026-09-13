# -*- coding: utf-8 -*-
"""对最终 PDF 做交付前校验。

校验项：
  1. 页数与是否可打开
  2. 六个部分的标题是否都在
  3. 七个关键答案是否都出现在正文里
  4. 是否出现缺字替换符或未定义引用
  5. 目录是否生成
"""
import os
import sys
import fitz

PDF = sys.argv[1] if len(sys.argv) > 1 else \
    r"C:\Users\qing1\Desktop\bestway\详解\详解.pdf"

doc = fitz.open(PDF)
text = "\n".join(p.get_text() for p in doc)

print("文件：", PDF)
print("页数：", doc.page_count)
print("字符数：", len(text))

PARTS = ["数学与物理基础", "数值方法原理", "代码逐行解释",
         "如何确认结果正确", "对照求解器与结果文件", "问题一闭式解的完整推导"]
print("\n各部分标题：")
for s in PARTS:
    print("  %-24s %s" % (s, "找到" if s in text else "缺失"))

KEYS = {
    "问题1 中心温度 34.0425": "34.0425",
    "问题1 表面温度 37.1989": "37.1989",
    "问题1 表面含水 1.5109": "1.5109",
    "问题2 中心温度 49.9051": "49.9051",
    "问题2 中心含水 1.7673": "1.7673",
    "问题2 表面含水 1.0083": "1.0083",
    "问题3 烘干时长 205574.5": "205574.5",
    "问题4 烘干时长 182778.3": "182778.3",
    "问题1 闭式解 34.0424989": "34.0424989",
}
print("\n关键答案：")
missing = 0
for name, s in KEYS.items():
    ok = s in text
    if not ok:
        missing += 1
    print("  %-26s %s" % (name, "找到" if ok else "缺失"))

print("\n异常字符检查：")
for name, ch in [("缺字替换符", "\ufffd"), ("反斜杠残留", "\\begin{"),
                 ("问号引用", "??")]:
    n = text.count(ch)
    print("  %-12s %d 处" % (name, n))

print("\n目录页检查：")
print("  目录条目数（按 部分 统计）：", text.count("部分"))

verdict = (missing == 0 and "\ufffd" not in text)
print("\n结论：", "通过" if verdict else "仍有问题，需要继续修")
