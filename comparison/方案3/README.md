# A题「药材的烘干问题」—— 新解法打包说明

本包是「与 paper 不同的新解法」的完整交付物：推导报告、全部代码、逐项验证报告、
结果文件与原始题目资料。

---

## 1. 先看什么

| 顺序 | 文件 | 说明 |
|---|---|---|
| ① | `report/report.pdf` | **主报告（20 页）**：模型推导、新解法论证、验证与归因 |
| ② | `新解法说明.md` | 速览版：方法对比表、答案对照表、自查记录 |
| ③ | `题目与附件/联想截图_20260912035339.png` | 对照基准（目标答案图） |
| ④ | `题目与附件/附件/附件3/result1.xlsx`…`result4.xlsx` | 本文生成的四个结果文件 |

---

## 2. 核心结论

| 项目 | 目标答案 | 本文结果 | 判定 |
|---|---|---|---|
| 问题1 中心温度 | 34.0425 | 34.042499 | **一致** |
| 问题1 表面温度 | 37.1989 | 37.198922 | **一致** |
| 问题1 表面含水 | 1.5109 | 1.510918 | **一致** |
| 问题2 中心温度 | 49.9051 | 49.905089 | **一致** |
| 问题2 中心含水 | 1.7673 | 1.767275 | **一致** |
| 问题2 表面含水 | 1.0083 | 1.008288 | **一致** |
| 问题3 t_dry | 205 559 s | 205 575.5 s | 差 0.008% |
| 问题4 t_dry | 182 764.7 s | 182 778.3 s | 差 0.007% |

7 项中 5 项完全一致；问题 1 的表 1（35 个值）与参考解完全相同。

问题 3、4 的差异已完整归因：**原样复现参考解的有限体积格式可将其表 7 序列逐值重现**
（问题 3 四档误差 ≤0.2 s；问题 4 的 Richardson 值 182764.72 对 182764.7），
证明模型完全一致；差异来自参考解长时间段 **Δt = 10 s 的一阶时间误差**，
其 Richardson 外推未收敛到连续极限。详见 `report/report.pdf` 第 9 章。

---

## 3. 新解法一句话

**坐标正则化 `u = (r/R)²` + Chebyshev 谱配置 + 自适应高阶隐式积分 + Bessel–Duhamel 闭式解**

四处与参考解不同：

1. **坐标**：不用 `ξ = r/R`，取 `u = ξ²`，使
   `(1/r)∂_r(rD∂_r) = (4/R²)∂_u(uD∂_u)` 在 `u=0` **完全正则**，轴心无需任何特殊处理；
2. **空间**：Chebyshev–Gauss–Lobatto 谱配置，**N=8 即达 10⁻¹⁰**；
3. **时间**：变阶变步长隐式 BDF（物性每步重估），非定步长冻结系数 θ 法；
4. **问题1 温度场**：准静态分解给出闭式级数解，**根本不离散**，2 项即达机器精度。

---

## 4. 目录结构

```
.
├── README.md                      本文件
├── report/
│   ├── report.pdf                 ★ 主报告（20 页，含 3 幅插图）
│   ├── report.tex                 LaTeX 源码
│   ├── report_extract.txt         PDF 文字提取（便于检索）
│   └── fig/                       报告插图（PNG）
├── 新解法说明.md                    速览版说明
├── code/                          全部求解与验证代码
│   ├── common.py                  参数、烘房环境、物性经验式、附件读取
│   ├── spectral.py                ★ 新解法核心：正则化谱配置 + 自适应 BDF
│   ├── ref_fv.py                  独立参照：均匀 r 网格守恒型有限体积
│   ├── paper_fv.py                参考解格式复现：节点中心 FV + θ 法
│   ├── step1 …step14*.py          13 个分步脚本（见下表）
│   ├── step*.txt                  每个脚本的原始输出报告
│   └── fig/                       插图生成结果
├── 题目与附件/
│   ├── A题.pdf                     题目原文
│   ├── 联想截图_20260912035339.png 目标答案图
│   └── 附件/
│       ├── 附件1.xlsx              烘房温湿度序列
│       ├── 附件2.xlsx              半径收缩序列
│       └── 附件3/result1–4.xlsx    ★ 本文生成的结果文件
├── 参考论文/
│   └── paper_electronic.pdf       参考论文（对照用）
└── 资料提取/
    ├── _A_ti.txt                   题目 PDF 文字提取
    ├── _paper.txt                  论文 PDF 文字提取（含其源代码）
    ├── _data12.txt                 附件1/2 数据转储
    └── _templates.txt              结果模板结构
```

---

## 5. 分步脚本与原始报告

按执行顺序编号，每个脚本对应一个同名 `.txt` 报告。运行环境：Python 3.9，
`numpy / scipy / pandas / openpyxl / matplotlib / pymupdf`；报告 PDF 另需 MiKTeX（xelatex）。

| 脚本 | 作用 | 报告 |
|---|---|---|
| `step1_fit.py` | 独立复核附件 1 的一阶惯性拟合 | `step1_fit.txt` |
| `step2_p1_analytic.py` | 问题1 温度场 Bessel–Duhamel 闭式解 | `step2_p1_analytic.txt` |
| `step2b`（记录在 txt） | 常数口径敏感性（解释 8×10⁻⁵ 差异） | `step2b_const_check.txt` |
| `step3_spectral_conv.py` | 谱收敛性 vs 解析解 | `step3_spectral_conv.txt` |
| `step4_p1_full.py` | 问题1 完整求解 + 表1/表2 | `step4_p1_full.txt` |
| `step5_ref_check.py` | 有限体积参照解（问题1） | `step5_ref_check.txt` |
| `step5b_fine.py` | 高精度逐点对照（临界第 4 位小数） | `step5b_fine.txt` |
| `step6_p2.py` | 问题2 求解 + 表3/表4 | `step6_p2.txt` |
| `step7_p3.py` | 问题3 t_dry + 表5 | `step7_p3.txt` |
| `step8_paper_repro.py` | **复现参考解表 7（问题3）** | `step8_paper_repro.txt` |
| `step9_p4.py` | 问题4 求解 + 表6 | `step9_p4.txt` |
| `step10_verify_tdry.py` | t_dry 数值可信度剖析 | `step10_verify_tdry.txt` |
| `step11_p4_repro.py` | **复现参考解表 7（问题4）** | `step11_p4_repro.txt` |
| `step12_make_results.py` | 生成 result1–4.xlsx | `step12_results.txt` |
| `step13`（记录在 txt） | 结果文件回读校验 | `step13_verify_results.txt` |
| `step14_make_figs.py` / `step14b_fig3.py` | 报告插图 | `fig/*.png` |

依赖顺序：`common.py` 被所有脚本引用；`spectral.py` / `ref_fv.py` / `paper_fv.py` 为三个求解器。
运行方式：`cd code && python stepN_xxx.py`（脚本内使用绝对路径指向题目附件，
如需移植请修改 `common.py` 顶部的 `_XLS` / `_XLS2` 两个常量）。

---

## 6. 报告 PDF 的编译

```bash
cd report
xelatex report.tex
xelatex report.tex      # 第二遍生成目录与交叉引用
```
需 TeX 发行版支持 `ctex` 宏包与 `xelatex`（本例用 MiKTeX）。

---

## 7. 自查与勘误

推导与实现过程中主动排查并修正了 10 项问题，完整清单见
`report/report.pdf` 第 10 章；其中最关键的 3 项：

* **建模**：收缩问题**不能**在 `r` 坐标下直接写扩散方程，必须补对流项
  `−Ṙ r ∂_r C / R`（已给出证明与 D=0 极限检验）；
* **口径**：参考解程序把拟合参数四舍五入为 `50.212 / 1877.5`，
  该 1×10⁻⁴ 的变化足以翻转第 4 位小数；
* **结论**：t_dry 的 16 s 差异不是模型差异（已用格式复现证明），
  而是参考解的一阶时间误差。
