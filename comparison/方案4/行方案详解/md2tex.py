# -*- coding: utf-8 -*-
"""Markdown -> LaTeX 转换器（修正版）
要点：
  1) 先按数学段切分，文本段才做转义与 Unicode 替换（旧版把 $2\\times3$ 换成 $2$\\times$3$，是 90 处报错的根因）
  2) 支持跨行的 $$...$$ 独立公式块
  3) 支持 markdown 表格 -> booktabs tabular
"""
import io
import os
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

SPECIAL = {'&': r'\&', '%': r'\%', '#': r'\#', '_': r'\_',
           '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}',
           '^': r'\textasciicircum{}'}
UNI = {'→': r'$\rightarrow$', '←': r'$\leftarrow$', '×': r'$\times$',
       '≥': r'$\ge$', '≤': r'$\le$', '≈': r'$\approx$', '±': r'$\pm$',
       '⇒': r'$\Rightarrow$', '⇔': r'$\Leftrightarrow$', '√': r'$\surd$',
       'α': r'$\alpha$', 'β': r'$\beta$', 'δ': r'$\delta$', 'ρ': r'$\rho$',
       'σ': r'$\sigma$', 'ξ': r'$\xi$', 'τ': r'$\tau$', 'θ': r'$\theta$',
       'μ': r'$\mu$', 'Δ': r'$\Delta$', '·': r'$\cdot$'}


def fmt_text(s):
    """只处理数学段之外的文本。"""
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if c == '\\' and i + 1 < len(s):
            out.append(s[i:i + 2]); i += 2; continue
        if c in UNI:
            out.append(UNI[c]); i += 1; continue
        out.append(SPECIAL.get(c, c)); i += 1
    s = ''.join(out)
    s = re.sub(r'\*\*(.+?)\*\*', r'\\textbf{\1}', s, flags=re.S)
    s = re.sub(r'(?<!\*)\*(?!\s)([^*\n]+?)(?<!\s)\*(?!\*)', r'\\emph{\1}', s)
    s = re.sub(r'`([^`\n]+)`', r'\\texttt{\1}', s)
    return s


def inline(s):
    parts = re.split(r'(\$\$.+?\$\$|\$[^$\n]+?\$)', s, flags=re.S)
    return ''.join(p if p.startswith('$') else fmt_text(p) for p in parts)


def convert_table(rows):
    head = [c.strip() for c in rows[0].strip().strip('|').split('|')]
    ncol = len(head)
    out = [r'\begin{center}', r'\small',
           r'\begin{tabular}{%s}' % ('l' * ncol), r'\toprule',
           ' & '.join(inline(c) for c in head) + r' \\', r'\midrule']
    for r in rows[2:]:
        cells = [c.strip() for c in r.strip().strip('|').split('|')]
        cells = (cells + [''] * ncol)[:ncol]
        out.append(' & '.join(inline(c) for c in cells) + r' \\')
    out += [r'\bottomrule', r'\end{tabular}', r'\end{center}', '']
    return out


def convert(md):
    lines = md.split('\n')
    out, i = [], 0
    while i < len(lines):
        raw = lines[i]
        s = raw.strip()

        if re.match(r'^</?(div|span|br|p)\b', s):
            i += 1; continue
        # 代码块
        if s.startswith('```'):
            i += 1
            buf = []
            while i < len(lines) and not lines[i].strip().startswith('```'):
                buf.append(lines[i]); i += 1
            i += 1
            out += [r'\begin{lstlisting}', *buf, r'\end{lstlisting}', '']
            continue
        # 表格
        if s.startswith('|') and i + 1 < len(lines) and re.match(r'^\|[\s:\-|]+\|$', lines[i + 1].strip()):
            rows = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                rows.append(lines[i]); i += 1
            out += convert_table(rows)
            continue
        # 跨行的 $$ ... $$ 独立公式
        if s.startswith('$$') and not (s.endswith('$$') and len(s) > 4):
            buf = [s[2:]]
            i += 1
            while i < len(lines) and '$$' not in lines[i]:
                buf.append(lines[i]); i += 1
            if i < len(lines):
                buf.append(lines[i].split('$$')[0]); i += 1
            body = '\n'.join(buf).strip()
            env = 'equation' if '\\tag{' in body else 'equation*'
            out += [r'\begin{%s}' % env, body, r'\end{%s}' % env, '']
            continue
        if s.startswith('$$') and s.endswith('$$') and len(s) > 4:
            body = s[2:-2].strip()
            env = 'equation' if '\\tag{' in body else 'equation*'
            out += [r'\begin{%s}' % env, body, r'\end{%s}' % env, '']
            i += 1; continue
        # 已有的 LaTeX 环境（\[ \]、equation、align）
        if s.startswith(r'\[') or s.startswith(r'\begin{equation}') or s.startswith(r'\begin{align}'):
            buf = [raw]
            endpat = r'\]' if s.startswith(r'\[') else r'\\end\{(equation|align)\}'
            while not re.search(endpat, '\n'.join(buf)) and i + 1 < len(lines):
                i += 1; buf.append(lines[i])
            out += buf + ['']
            i += 1; continue
        # 标题
        m = re.match(r'^(#{1,4})\s+(.*)$', s)
        if m:
            lvl, txt = len(m.group(1)), m.group(2)
            cmd = {1: 'section', 2: 'subsection', 3: 'subsubsection', 4: 'paragraph'}[lvl]
            txt = re.sub(r'^\d+(\.\d+)*\s*', '', txt)
            out += [r'\%s{%s}' % (cmd, inline(txt)), '']
            i += 1; continue
        if re.match(r'^-{3,}$', s):
            out += [r'\vspace{1mm}\hrule\vspace{1mm}', '']
            i += 1; continue
        if s.startswith('>'):
            buf = []
            while i < len(lines) and lines[i].strip().startswith('>'):
                buf.append(lines[i].strip().lstrip('>').strip()); i += 1
            out += [r'\begin{quote}', inline(' '.join(buf)), r'\end{quote}', '']
            continue
        if re.match(r'^[-*]\s+', s) or re.match(r'^\d+[.)]\s+', s):
            env = 'itemize' if re.match(r'^[-*]\s+', s) else 'enumerate'
            out.append(r'\begin{%s}' % env)
            while i < len(lines) and (re.match(r'^\s*[-*]\s+', lines[i]) or
                                      re.match(r'^\s*\d+[.)]\s+', lines[i])):
                item = re.sub(r'^\s*(?:[-*]|\d+[.)])\s+', '', lines[i])
                out.append(r'\item ' + inline(item.strip()))
                i += 1
            out += [r'\end{%s}' % env, '']
            continue
        if not s:
            out.append(''); i += 1; continue
        out += [inline(s), '']
        i += 1
    return '\n'.join(out)


if __name__ == '__main__':
    md = io.open(sys.argv[1], encoding='utf-8').read()
    io.open(sys.argv[2], 'w', encoding='utf-8', newline='\n').write(convert(md))
    print('转换完成:', os.path.basename(sys.argv[1]))
