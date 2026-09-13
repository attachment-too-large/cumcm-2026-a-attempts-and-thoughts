# -*- coding: utf-8 -*-
"""校验图解析文档：页数、七张图是否都排入、关键读数是否齐全。"""
import fitz

PDF = r"C:\Users\qing1\Desktop\bestway\图解析\图解析.pdf"
doc = fitz.open(PDF)
t = "".join(p.get_text() for p in doc)

print("页数：", doc.page_count)
print("字符数：", len(t))

print("")
print("插图检查（每页图片数）：")
n = 0
for i, p in enumerate(doc):
    k = len(p.get_images(full=True))
    if k:
        print("  第 %d 页：%d 张" % (i + 1, k))
        n += k
print("  图片总数：", n, "（应为 7）")

print("")
print("七张图的小节标题：")
for s in ["图 1：控制体与两种径向离散", "图 2：问题一温度与水分剖面",
          "图 3：坐标变换的几何意义", "图 4：切比雪夫节点的来源",
          "图 5：问题三干燥曲线与终点放大", "图 6：通量形式的谱配置",
          "图 7：两种离散方式的收敛对比"]:
    print("  %-28s %s" % (s, "找到" if s in t else "缺失"))

print("")
print("关键读数检查：")
for k in ["7.59", "2.88", "205574.519", "57.1040", "35.9", "35.845",
          "26.121", "1.5109", "37.1989", "34.0425"]:
    print("  %-12s %s" % (k, "在" if k in t else "缺失"))

print("")
print("异常检查：缺字替换符 %d 处，问号引用 %d 处" %
      (t.count("\ufffd"), t.count("??")))
