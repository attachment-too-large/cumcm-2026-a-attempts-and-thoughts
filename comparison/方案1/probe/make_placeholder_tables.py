# -*- coding: utf-8 -*-
# 临时脚本: 为尚未生成的真实表格写入占位文件, 以便提前做 LaTeX 编译冒烟测试。
# s08_make_tables.py 会覆盖真实表格; 占位表一旦被覆盖即自然消失。
import io
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(BASE, "paper", "tables")
NAMES = ["tab_p1_T", "tab_p1_C", "tab_p2_T", "tab_p2_C", "tab_p3_C", "tab_p4_C",
         "tab_convergence_p1", "tab_convergence_p3", "tab_convergence_p4",
         "tab_analytic", "tab_conservation", "tab_scheme",
         "tab_sensitivity_oat", "tab_plateau", "tab_mc_contrib", "tab_inverse"]
DUMMY = "\n".join([
    "\\begin{table}[H]", "\\centering",
    "\\caption{占位表（编译冒烟测试用）}", "\\label{tab:placeholder}",
    "\\small", "\\begin{tabular}{@{}lrr@{}}", "\\toprule",
    "项目 & A & B \\\\", "\\midrule", "占位 & 0.0000 & 0.0000 \\\\",
    "\\bottomrule", "\\end{tabular}", "\\end{table}",
])
os.makedirs(TAB, exist_ok=True)
for n in NAMES:
    p = os.path.join(TAB, n + ".tex")
    io.open(p, "w", encoding="utf-8", newline="\n").write(DUMMY + "\n")
print("wrote", len(NAMES), "placeholder tables")
