# -*- coding: utf-8 -*-
"""检查七张插图是否都排进了 PDF，以及它们落在第几页。"""
import fitz

PDF = r"C:\Users\qing1\Desktop\bestway\详解\详解.pdf"
doc = fitz.open(PDF)

markers = {
    "控制体几何": "圆柱环壳控制体与两种径向离散的对照",
    "剖面演化": "问题 1 的温度剖面与水分剖面演化",
    "坐标变换": "坐标变换 u 等于半径比值的平方",
    "切比雪夫节点": "切比雪夫节点的来源",
    "干燥曲线": "问题 3 的干燥曲线",
    "通量形式": "通量形式的谱配置",
    "收敛对比": "两种离散方式的空间收敛",
}

print("页数：", doc.page_count)
for i, page in enumerate(doc):
    txt = page.get_text().replace("\n", "")
    for name, key in markers.items():
        if key.replace(" ", "") in txt.replace(" ", ""):
            print("  %-12s 第 %d 页" % (name, i + 1))

print("")
print("各页图片数量统计：")
cnt = 0
for i, page in enumerate(doc):
    n = len(page.get_images(full=True))
    if n:
        print("  第 %d 页：%d 张" % (i + 1, n))
        cnt += n
print("图片总数：", cnt)
