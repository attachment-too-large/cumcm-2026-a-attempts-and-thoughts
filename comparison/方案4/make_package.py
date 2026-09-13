# -*- coding: utf-8 -*-
"""把全部内容（含更正版）打进一个 zip。"""
import io
import os
import shutil
import sys
import zipfile

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

WS = r'C:\Users\qing1\Desktop\2026CUMCM_'
RUN = r'C:\Users\qing1\Desktop\对比\行'
STAGE = os.path.join(WS, '_stage', '2026CUMCM_A题_完整交付包')
ZIP = os.path.join(WS, '2026CUMCM_A题_完整交付包.zip')

SKIP_EXT = {'.aux', '.log', '.out', '.toc', '.pyc', '.html', '.synctex.gz', '.fls', '.fdb_latexmk'}
SKIP_DIRS = {'__pycache__', '.idea', '.git', '_stage', '_打包暂存'}


def copy_tree(src, dst, skip_ext=SKIP_EXT):
    n = 0
    for root, dirs, files in os.walk(src):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        rel = os.path.relpath(root, src)
        outdir = os.path.join(dst, rel) if rel != '.' else dst
        os.makedirs(outdir, exist_ok=True)
        for f in files:
            if os.path.splitext(f)[1].lower() in skip_ext:
                continue
            shutil.copy2(os.path.join(root, f), os.path.join(outdir, f))
            n += 1
    return n


if os.path.exists(STAGE):
    shutil.rmtree(STAGE)
os.makedirs(STAGE)

README = r'''# 2026 CUMCM A 题「药材的烘干问题」完整交付包

包含：**交付报告 4 份**、**两套独立生产代码与全部结果**、**赛题原始数据**、**全部核查脚本与原始输出**。

---

## 目录

| 目录 | 内容 |
|---|---|
| `A_交付报告/` | 四份 PDF 报告及其 LaTeX / Markdown 源 |
| `A0_更正版解答与自查报告/` | **15 页**：按审阅结论做出的 11 项更正、正确性验证、与答案 X 的偏离归因与自我反思 |
| `A1_逐步推导详解/` | **37 页**：从零开始、不跳步的完整推导论文（XeLaTeX + cumcmthesis） |
| `A2_两方案对比检验报告/` | **8 页**：两套独立方案的对比检验（2×2 环境互换交叉实验） |
| `A3_论文逐式审阅报告/` | **39 页**：对方法 B 论文 14 个编号公式的逐条审计 + 26 项问题 + 19 条优化措施 |
| `B_生产代码/` | 两套完整可运行的生产代码与结果 |
| `B1_行方案/` | 节点中心有限体积 + 一阶惯性环境（code / data / figures / paper / 修改说明） |
| `B2_方法B方案/` | 分级谱元弱形式 + 分段线性环境（源程序 / 结果文件 / 实验日志 / 中间结果与图表 / 赛题数据） |
| `C_赛题原始数据/` | 附件 1、附件 2、附件 3（4 个结果模板） |
| `D_核查脚本与原始输出/` | 34 个 Python 核查脚本 + 原始输出（**文件名全部为英文**） |
| `E_论文编译校验/` | 用「行」的 LaTeX 源重编，验证可复现性（93 页 / 4098 KB / 0 错误） |

> 所有 `.py` 脚本文件名为纯英文，避免跨平台与命令行传参的编码问题；报告与数据保留中文名便于识别。

---

## 核心结论速览

### 1. 四个问题的最终答案

| 问题 | 答案 |
|---|---|
| 问题 1 | 1800 s：中心 **33.5753 ℃**、表面 **36.7856 ℃**；中心含水 2.5500、表面 **1.5102** |
| 问题 2 | 3 h：中心 **49.8495 ℃**；中心含水 1.7662、表面 1.0081 |
| 问题 3 | 烘干时长 **57.17 h**，含环境测量噪声 ±0.3 h；模型形式区间 **53.1—57.2 h** |
| 问题 4 | 烘干时长 **50.82 h**；**收缩使时长缩短 60.6%**（不收缩对照 129.09 h） |

### 2. 与答案 X 的偏离全部由同一原因解释

**2×2 环境互换交叉实验**：

| | 行方案环境（一阶惯性式） | 方法B环境（分段线性） |
|---|---|---|
| 行方案 FV N=1600 | 34.0425 / 37.1989 / 1.5109 | 33.5753 / 36.7856 / 1.5102 |
| 方法B 谱元 ndof=81 | 34.0425 / 37.1989 / 1.5109 | 33.5753 / 36.7856 / 1.5102 |

同一环境下两种完全不同的数值方法四位小数一致。

**偏离判定（分项）**：

| 项目 | 偏差 | 与之比较的不确定度带 | 判定 |
|---|---|---|---|
| 问题 1 温度 | 0.47 ℃ | 环境噪声带 0.036 ℃ | **显著**（13 倍） |
| 问题 3/4 时长 | 0.10%~0.12% | 环境噪声带 ±0.44% | **不显著**（1/4） |
| 问题 4 收缩效应 | 0 | — | **完全一致** |

### 3. 噪声-偏置判定性实验（本轮最关键的自检）

| 扰动 | 对问题 1 温度的影响 | 对问题 3 时长的影响 |
|---|---|---|
| 环境噪声（同一 σ=0.2，3 个实现） | 0.036 ℃ 散布 | **±0.25 h（±0.44%）** |
| 拟合替代实测（系统偏置） | 0.467 ℃ | −0.11% |

→ **两个问题的最优环境处理不同**：问题 1/2 用实测插值（噪声被低通滤掉，偏置滤不掉）；
问题 3/4 两种口径都可用，但**必须把 ±0.25 h 的噪声不确定度报出来**（它比口径差大 4 倍）。

### 4. 逐式审计结果（方法 B 论文）

14 个编号公式：**9 个完全正确**、**3 个有推导缺口但结论可用**、**2 个有实质错误或错配**。
216 个表格格子与两条烘干时长逐位复现无误；26 项问题分级；19 条优化措施。

---

## 运行环境

- Python 3.9 + numpy 2.0.2 + scipy 1.13.1 + openpyxl + pandas + PyMuPDF
- MiKTeX 25.12（XeLaTeX）

## 复现方式

```powershell
# 重编报告
cd "A_交付报告\A0_更正版解答与自查报告\LaTeX源"; xelatex corrected_report.tex   # 两遍

# 复现「行」论文
cd E_论文编译校验; xelatex paper_electronic.tex                                   # 两遍

# 跑核查脚本
cd D_核查脚本与原始输出\脚本
python s01_noise_vs_bias.py       # 噪声 vs 偏置判定性实验（约 20 分钟）
python s02_result_reproduce.py    # 端到端重算 result1-4（约 25 分钟）
python c04b_weak_form_manufactured.py   # 制造解验证（约 5 秒）
```

## 脚本命名对照

| 前缀 | 含义 |
|---|---|
| `c01`–`c13` | 对方法 B 论文的逐式核查 |
| `d01`–`d05` | 与「行」方案的比对与差异溯源 |
| `e01` | 环境口径 2×2 交叉检验（决定性实验） |
| `s01`–`s02` | 更正版自检（噪声-偏置实验、端到端重算） |
| `md2tex` / `make_latex` / `rebuild_compare` / `build_corrected` | 报告构建 |
'''
io.open(os.path.join(STAGE, 'README.md'), 'w', encoding='utf-8', newline='\n').write(README)

X = os.path.join(WS, '行方案详解')
H = os.path.join(WS, '更正版')

# ---------------- A ----------------
d = os.path.join(STAGE, 'A_交付报告'); os.makedirs(d)

a0 = os.path.join(d, 'A0_更正版解答与自查报告'); os.makedirs(a0)
shutil.copy2(os.path.join(H, '输出', '更正版解答与自查报告.pdf'), a0)
shutil.copy2(os.path.join(H, 'README.md'), os.path.join(a0, 'README.md'))
os.makedirs(os.path.join(a0, 'LaTeX源'))
for f in ('corrected_report.tex', 'body_correct.tex', 'body_correct.md', 'body_patch.tex'):
    p = os.path.join(H, 'latex', f)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(a0, 'LaTeX源', f))
os.makedirs(os.path.join(a0, 'Markdown源'))
shutil.copy2(os.path.join(H, '报告', '00_更正版解答与自查报告.md'), os.path.join(a0, 'Markdown源'))
os.makedirs(os.path.join(a0, '代码更正补丁'))
shutil.copy2(os.path.join(H, '补丁', '代码更正补丁.md'), os.path.join(a0, '代码更正补丁'))

a1 = os.path.join(d, 'A1_逐步推导详解'); os.makedirs(a1)
for f in ('中药材烘干_行方案_逐步推导详解_LaTeX版.pdf', '中药材烘干_行方案_逐步推导详解.pdf'):
    shutil.copy2(os.path.join(X, '输出', f), a1)
os.makedirs(os.path.join(a1, 'LaTeX源'))
for f in ('中药材烘干_行方案_逐步推导详解_LaTeX版.tex', 'body_detail.tex', 'body_detail.md', 'cumcmthesis.cls'):
    p = os.path.join(X, 'latex', f)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(a1, 'LaTeX源', f))
os.makedirs(os.path.join(a1, 'Markdown源'))
for f in ('00_导读与方程推导.md', '01_收缩与数值方法.md', '02_环境求解与验证.md'):
    shutil.copy2(os.path.join(X, '报告', f), os.path.join(a1, 'Markdown源', f))

a2 = os.path.join(d, 'A2_两方案对比检验报告'); os.makedirs(a2)
for f in ('中药材烘干_两方案对比检验报告_LaTeX版.pdf', '中药材烘干_两方案对比检验报告.pdf'):
    shutil.copy2(os.path.join(X, '输出', f), a2)
os.makedirs(os.path.join(a2, 'LaTeX源'))
for f in ('中药材烘干_两方案对比检验报告_LaTeX版.tex', 'body_compare.tex', 'body_compare.md'):
    p = os.path.join(X, 'latex', f)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(a2, 'LaTeX源', f))
shutil.copy2(os.path.join(X, '报告', '10_对比报告.md'), os.path.join(a2, 'Markdown源.md'))
shutil.copy2(os.path.join(X, 'README与自检报告.md'), os.path.join(a2, '自检报告.md'))

a3 = os.path.join(d, 'A3_论文逐式审阅报告'); os.makedirs(a3)
shutil.copy2(os.path.join(WS, '审阅报告', '输出', '2026CUMCM_A题论文独立审阅报告.pdf'), a3)
shutil.copy2(os.path.join(WS, '审阅报告', 'README.md'), os.path.join(a3, 'README.md'))
copy_tree(os.path.join(WS, '审阅报告', '报告'), os.path.join(a3, '报告源文件'))

# ---------------- B ----------------
d = os.path.join(STAGE, 'B_生产代码'); os.makedirs(d)
n1 = copy_tree(os.path.join(RUN, 'code'), os.path.join(d, 'B1_行方案', 'code'))
n2 = copy_tree(os.path.join(RUN, 'data'), os.path.join(d, 'B1_行方案', 'data'))
n3 = copy_tree(os.path.join(RUN, 'paper'), os.path.join(d, 'B1_行方案', 'paper'))
n4 = copy_tree(os.path.join(RUN, '修改说明'), os.path.join(d, 'B1_行方案', '修改说明'))
n5 = copy_tree(os.path.join(RUN, 'figures'), os.path.join(d, 'B1_行方案', 'figures'))
n6 = copy_tree(os.path.join(WS, 'A题_支撑材料'), os.path.join(d, 'B2_方法B方案'))

# ---------------- C ----------------
d = os.path.join(STAGE, 'C_赛题原始数据'); os.makedirs(d)
n7 = copy_tree(os.path.join(WS, 'A题_支撑材料', '赛题原始数据'), d)

# ---------------- D ----------------
d = os.path.join(STAGE, 'D_核查脚本与原始输出'); os.makedirs(d)
nd1 = copy_tree(os.path.join(WS, '审阅报告', '脚本'), os.path.join(d, '脚本'))
for src in (X, H, os.path.join(H, '脚本')):
    if not os.path.isdir(src):
        continue
    for f in sorted(os.listdir(src)):
        if f.endswith('.py'):
            shutil.copy2(os.path.join(src, f), os.path.join(d, '脚本', f)); nd1 += 1
for f in ('rename_scripts.py', 'make_package.py'):
    p = os.path.join(WS, f)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(d, '脚本', f)); nd1 += 1
os.makedirs(os.path.join(d, '原始输出'))
nd2 = 0
for base in (os.path.join(WS, '审阅报告', '输出'), os.path.join(X, '输出'),
             os.path.join(X, 'latex'), os.path.join(H, '输出')):
    if not os.path.isdir(base):
        continue
    for f in sorted(os.listdir(base)):
        if os.path.splitext(f)[1].lower() in ('.txt', '.json', '.csv'):
            shutil.copy2(os.path.join(base, f), os.path.join(d, '原始输出', f)); nd2 += 1

# ---------------- E ----------------
d = os.path.join(STAGE, 'E_论文编译校验'); os.makedirs(d)
src = os.path.join(X, 'latex', '行论文编译校验')
ne = 0
for f in sorted(os.listdir(src)):
    p = os.path.join(src, f)
    if os.path.isfile(p) and os.path.splitext(f)[1].lower() in ('.tex', '.cls', '.pdf', '.bib'):
        shutil.copy2(p, d); ne += 1
shutil.copy2(os.path.join(X, 'latex', 'cumcmthesis.cls'), d)

# ---------------- zip ----------------
if os.path.exists(ZIP):
    os.remove(ZIP)
total = 0
with zipfile.ZipFile(ZIP, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for root, dirs, files in os.walk(STAGE):
        dirs[:] = [x for x in dirs if x not in SKIP_DIRS]
        for f in files:
            p = os.path.join(root, f)
            arc = os.path.join('2026CUMCM_A题_完整交付包', os.path.relpath(p, STAGE))
            z.write(p, arc)
            total += os.path.getsize(p)

print('B1行方案 code=%d data=%d paper=%d 修改说明=%d figures=%d' % (n1, n2, n3, n4, n5))
print('B2方法B方案 %d ；C赛题数据 %d ；E编译校验 %d' % (n6, n7, ne))
print('D 脚本 %d 个、原始输出 %d 份' % (nd1, nd2))
print()
print('ZIP: %s' % ZIP)
print('  未压缩 %.1f MB  ->  压缩后 %.1f MB' % (total / 1e6, os.path.getsize(ZIP) / 1e6))
with zipfile.ZipFile(ZIP) as z:
    names = z.namelist()
bad = [x for x in names if x.lower().endswith(('.html', '.aux', '.log', '.out', '.pyc'))]
cn = [x for x in names if x.endswith('.py') and any(ord(c) > 127 for c in os.path.basename(x))]
print('  条目 %d ；临时文件 %s ；中文名 .py %s'
      % (len(names), bad if bad else '无', cn if cn else '无'))
import collections
c = collections.Counter(x.split('/')[1] for x in names if x.count('/') >= 1)
for k in sorted(c):
    print('    %-28s %d 项' % (k, c[k]))
