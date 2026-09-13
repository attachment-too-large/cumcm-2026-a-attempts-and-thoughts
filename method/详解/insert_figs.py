# -*- coding: utf-8 -*-
"""把教学插图插入讲稿正文的指定位置。

每个插入点用一段唯一的正文片段作为锚点，脚本先检查锚点唯一，
再在锚点之前插入 figure 环境。若锚点不唯一或找不到，脚本报错退出，
不会写出半成品。
"""
import io
import os
import sys

W = r"C:\Users\qing1\Desktop\bestway\详解"

# (文件, 锚点, 图文件名, 图标题, 标签, 宽度)
PLAN = [
    ("part1.tex", "\\subsection{能量守恒，一步不跳}", "fig_cv.png",
     "圆柱环壳控制体与两种径向离散的对照。(a) 取半径 r 到 r+dr 的一层薄壳作控制体，"
     "侧面积是 2 pi r L；(b) 常规 r 坐标下的节点中心有限体积，第一个控制体只有 "
     "半格宽，其体积测度是 h 平方除以 8。",
     "fig:cv", 0.96),

    ("part1.tex", "\\section{边界条件与初始条件}", "fig_p1.png",
     "问题 1 的温度剖面与水分剖面演化。(a) 温度在 30 分钟内已经拉开明显梯度；"
     "(b) 水分在同样时间内几乎只动表层。两张图的差别来自热扩散率与水分扩散系数的量级差别。",
     "fig:p1", 0.96),

    ("part2.tex", "\\subsection{第一步：求出两个导数之间的换算关系}", "fig_map.png",
     "坐标变换 u 等于半径比值的平方。(a) 等距的半径映射成不等距的 u；"
     "(b) 同一批节点在两种坐标下的位置，上半行是 u 坐标，下半行是物理半径。",
     "fig:map", 0.96),

    ("part2.tex", "\\subsection{插值多项式的唯一性}", "fig_cheb.png",
     "切比雪夫节点的来源。(a) 在单位圆上按等角距取点，再垂直投影到直径上；"
     "(b) 投影结果是两端密、中间疏的分布，两端间距约为 N 平方分之一量级。",
     "fig:cheb", 0.96),

    ("part2.tex", "\\subsection{判据的数学表述}", "fig_p3.png",
     "问题 3 的干燥曲线。(a) 各半径处水分浓度随时间下降，横轴到 60 小时；"
     "(b) 终点附近放大，截面平均值先达标，中心最后达标，"
     "所以判据必须取全场最大值而不是平均值。",
     "fig:p3", 0.96),

    ("part3.tex", "\\subsection{右端项函数：核心中的核心}", "fig_flux.png",
     "通量形式的谱配置。内部节点的通量由谱微分矩阵算出，"
     "表面节点的通量直接取自 Robin 边界条件的物理值，"
     "轴心节点的通量因为 u 等于零而自动为零。",
     "fig:flux", 0.96),

    ("part4.tex", "\\subsection{空间方向的加密}", "fig_conv.png",
     "两种离散方式的空间收敛。(a) 谱配置法：节点参数 N 从 4 增到 6 时误差从 "
     "7.6e-6 掉到 2.9e-10，此后落在时间积分容差平台上不再下降；"
     "(b) 有限体积对照格式只达到二阶收敛。",
     "fig:conv", 0.96),
]


def figblock(name, caption, label, width):
    return (
        "\\begin{figure}[H]\n"
        "\\centering\n"
        "\\includegraphics[width=%.2f\\textwidth]{fig/%s}\n"
        "\\caption{%s}\n"
        "\\label{%s}\n"
        "\\end{figure}\n\n" % (width, name, caption, label)
    )


def main():
    cache = {}
    for fname, anchor, name, caption, label, width in PLAN:
        path = os.path.join(W, fname)
        if path not in cache:
            with io.open(path, encoding="utf-8") as fh:
                cache[path] = fh.read()
        text = cache[path]
        n = text.count(anchor)
        if n != 1:
            print("FAIL %s : 锚点出现 %d 次 : %s" % (fname, n, anchor))
            return 1
        if "fig/" + name in text:
            print("SKIP %s : 已插入 %s" % (fname, name))
            continue
        text = text.replace(anchor, figblock(name, caption, label, width) + anchor, 1)
        cache[path] = text
        print("OK   %s <- %s" % (fname, name))

    for path, text in cache.items():
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
    print("所有插图已写入")
    return 0


if __name__ == "__main__":
    sys.exit(main())
