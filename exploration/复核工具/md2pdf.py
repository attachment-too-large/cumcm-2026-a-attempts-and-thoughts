# -*- coding: utf-8 -*-
"""把本目录的说明文档（Markdown）渲染为可打印 HTML，再用 Chrome/Edge 打印为 PDF。

用法: python _md2pdf.py <input.md> <output.pdf> [标题]
"""
import sys, io, os, re, html, subprocess, tempfile

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

CSS = r"""
@page { size: A4; margin: 18mm 16mm 16mm 16mm; }
* { box-sizing: border-box; }
body {
  font-family: "Microsoft YaHei", "Noto Sans SC", "SimSun", sans-serif;
  font-size: 10.5pt; line-height: 1.75; color: #1a1a1a;
  margin: 0; padding: 0; -webkit-print-color-adjust: exact; print-color-adjust: exact;
}
h1 { font-size: 19pt; text-align: center; margin: 0 0 4mm 0; padding-bottom: 3mm;
     border-bottom: 2.5px solid #1f4e79; color: #12395c; letter-spacing: .5px; }
h2 { font-size: 14pt; color: #12395c; margin: 8mm 0 3mm 0; padding: 1.6mm 0 1.6mm 3mm;
     border-left: 4px solid #1f4e79; background: #eef4fa; page-break-after: avoid; }
h3 { font-size: 12pt; color: #1f4e79; margin: 5.5mm 0 2mm 0; page-break-after: avoid; }
h4 { font-size: 11pt; color: #333; margin: 4mm 0 1.5mm 0; page-break-after: avoid; }
p { margin: 1.8mm 0; }
strong { color: #0d2f4d; }
ul, ol { margin: 1.5mm 0 1.5mm 0; padding-left: 7mm; }
li { margin: .9mm 0; }
code { font-family: Consolas, "Courier New", monospace; font-size: 9.4pt;
       background: #f2f4f7; padding: .3mm 1mm; border-radius: 2px; color: #b0004a; }
table { border-collapse: collapse; width: 100%; margin: 2.5mm 0 3.5mm 0;
        font-size: 9.3pt; page-break-inside: avoid; }
th { background: #1f4e79; color: #fff; font-weight: 600; padding: 1.5mm 1.6mm;
     border: 1px solid #7f9db9; text-align: left; }
td { padding: 1.3mm 1.6mm; border: 1px solid #b9c9d9; vertical-align: top; }
tbody tr:nth-child(even) td { background: #f6f9fc; }
blockquote { margin: 2.5mm 0; padding: 2mm 3mm; background: #fffbf0;
             border-left: 4px solid #e0a83c; color: #5a4a20; font-size: 9.8pt; }
blockquote p { margin: 1mm 0; }
hr { border: none; border-top: 1px dashed #b9c9d9; margin: 5mm 0; }
.math { font-family: "Cambria Math", "Times New Roman", serif; font-size: 10.5pt; }
.mathblock { text-align: center; margin: 3mm 0; font-family: "Cambria Math", "Times New Roman", serif;
             font-size: 11pt; page-break-inside: avoid; }
sub, sup { font-size: 7.5pt; }
.meta { text-align: center; color: #5a6b7a; font-size: 9.2pt; margin-bottom: 5mm; }
"""


def esc(s):
    return html.escape(s, quote=False)


GREEK = [
    (r'\\le\b', '≤'), (r'\\ge\b', '≥'),
    (r'\\partial', '∂'), (r'\\nabla', '∇'), (r'\\infty', '∞'), (r'\\propto', '∝'),
    (r'\\approx', '≈'), (r'\\times', '×'), (r'\\cdot', '·'), (r'\\quad', '  '),
    (r'\\qquad', '    '), (r'\\to', '→'), (r'\\rightarrow', '→'), (r'\\le\b', '≤'),
    (r'\\ge\b', '≥'), (r'\\neq', '≠'), (r'\\pm', '±'), (r'\\circ', '°'),
    (r'\\in\b', '∈'), (r'\\forall', '∀'), (r'\\sim', '~'), (r'\\ll', '≪'), (r'\\gg', '≫'),
    (r'\\sum', 'Σ'), (r'\\int_0\^R', '∫₀^R'), (r'\\int\b', '∫'), (r'\\ne\b', ' ≠ '), (r'\\neq', ' ≠ '),
    (r'\\equiv', '≡'), (r'\\Rightarrow', '⇒'), (r'\\int_0\^1', '∫₀¹'),
    (r'\\alpha', 'α'), (r'\\beta', 'β'), (r'\\gamma', 'γ'), (r'\\delta', 'δ'),
    (r'\\Delta', 'Δ'), (r'\\epsilon', 'ε'), (r'\\theta', 'θ'), (r'\\lambda', 'λ'),
    (r'\\mu', 'μ'), (r'\\nu', 'ν'), (r'\\xi', 'ξ'), (r'\\pi', 'π'), (r'\\rho', 'ρ'),
    (r'\\sigma', 'σ'), (r'\\tau', 'τ'), (r'\\phi', 'φ'), (r'\\Phi', 'Φ'),
    (r'\\omega', 'ω'), (r'\\Omega', 'Ω'),
    (r'\\mathrm\{([^{}]*)\}', r'\1'), (r'\\text\{([^{}]*)\}', r'\1'),
    (r'\\ln\b', ' ln '), (r'\\log\b', ' log '), (r'\\exp\b', ' exp '),
    (r'\\left', ''), (r'\\right', ''),
]


def latex_to_text(s):
    """把 LaTeX 片段转成可读的 Unicode 文本（不依赖 MathJax/KaTeX）。"""
    frac = re.compile(r'\\d?frac\{([^{}]*)\}\{([^{}]*)\}')
    for _ in range(6):
        new = frac.sub(lambda m: '(%s)/(%s)' % (m.group(1).strip(), m.group(2).strip()), s)
        if new == s:
            break
        s = new
    s = re.sub(r'\\sqrt\{([^{}]*)\}', r'√(\1)', s)
    s = re.sub(r'\\bar\{([^{}]*)\}', r'\1̄', s)
    s = re.sub(r'\\hat\{([^{}]*)\}', r'\1̂', s)
    s = re.sub(r'\\tilde\{([^{}]*)\}', r'\1̃', s)
    for pat, rep in GREEK:
        s = re.sub(pat, rep, s)
    # 负指数/分数指数用 ^(...) 表达，避免歧义
    s = re.sub(r'\^\s*\{([^{}]*)\}', lambda m: expo(m.group(1)), s)
    s = re.sub(r'\^\s*(-?[0-9])', lambda m: expo(m.group(1)), s)
    s = re.sub(r'_\s*\{([^{}]*)\}', lambda m: subv(m.group(1)), s)
    s = re.sub(r'_\s*([A-Za-z0-9])', lambda m: subv(m.group(1)), s)
    s = s.replace('\\,', ' ').replace('\\;', ' ').replace('\\!', '')
    s = s.replace('\\ ', ' ').replace('{', '').replace('}', '')
    s = re.sub(r'\s{2,}', ' ', s)
    s = sub2uni(s)
    return s.strip()


def _plain(t):
    return t.strip().lstrip('\\')


def expo(t):
    t = t.strip()
    if t.startswith('\\'):
        t = t[1:]
    if not t:
        return ''
    if re.fullmatch(r'[0-9]+', t):
        return t.translate(str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹'))
    if re.fullmatch(r'[+-]', t):
        return t.translate(str.maketrans('+-', '⁺⁻'))
    if re.fullmatch(r'-[0-9]+', t):
        return t.translate(str.maketrans('0123456789-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻'))
    return '^' + ('(%s)' % t if len(t) > 1 else t)


def subv(t):
    t = _plain(t)
    if not t:
        return ''
    table = str.maketrans('0123456789+-=()aeioxhklmnpst∞',
                          '₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑᵢₒₓₕₖₗₘₙₚₛₜ∞')
    if all(ch in '0123456789+-=()aeioxhklmnpst∞' for ch in t):
        return t.translate(table)
    return '_{%s}' % t


def sub2uni(s):
    """把散落在正文里的裸下标（如 T_∞、R_0）转成 Unicode 下标。"""
    def rep(m):
        v = subv(m.group(1))
        return v if not v.startswith('_') else m.group(0)
    return re.sub(r'(?<=[A-Za-z])_\{?([A-Za-z0-9∞]+)\}?', rep, s)


def degfix(s):
    """把 $...$ 里残留的 ^\\circ 之类收尾。"""
    s = re.sub(r'\^\s*(?:\\circ|°|o|O)\s*C', '°C', s)
    return s


def inline(s):
    """行内标记：行内代码、行内公式、加粗。"""
    stash = []

    def keep(m):
        stash.append(m.group(1))
        return '\x00%d\x00' % (len(stash) - 1)

    s = re.sub(r'`([^`]+)`', keep, s)
    # 行内公式 $...$（含跨行）
    s = re.sub(r'\$([^$]+?)\$', lambda m: sub2uni(degfix(latex_to_text(m.group(1)))), s, flags=re.S)
    s = sub2uni(s)
    s = re.sub(r'\\d?frac\{[^{}]*\}\{[^{}]*\}',
               lambda m: sub2uni(latex_to_text(m.group(0))), s)
    s = re.sub(r'\\(?:alpha|beta|gamma|delta|Delta|epsilon|theta|lambda|mu|nu|xi|pi|rho|'
               r'sigma|tau|phi|Phi|omega|Omega|partial|infty|propto|approx|times|cdot|circ|in|'
               r'sum|int|ne|neq|equiv|forall|sim|ll|gg|le|ge|to|pm|neq)\b',
               lambda m: latex_to_text(m.group(0)).replace('\\int', '∫').replace('\\sum', 'Σ'), s)
    s = esc(s)
    s = re.sub(r'\^\s*°C', '°C', s).replace(' ^°C', ' °C').replace('^ C', ' C')
    s = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', s)
    s = re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)', r'<em>\1</em>', s)
    for i, code in enumerate(stash):
        s = s.replace('\x00%d\x00' % i, '<code>%s</code>' % esc(code))
    return s


def md_to_html(text):
    # 先做公式预处理（避免数学内容被二次转换）：
    #   1) 保护行内代码
    stash = []

    def keep(m):
        stash.append(m.group(1))
        return '\x00%d\x00' % (len(stash) - 1)

    text = re.sub(r'`([^`]+)`', keep, text)
    #   2) 块级公式 $$...$$
    text = re.sub(r'\$\$(.+?)\$\$',
                  lambda m: '\n\nMATHBLOCK\x01' + degfix(latex_to_text(
                      ' '.join(m.group(1).split()))) + '\x02\n\n',
                  text, flags=re.S)
    #   3) 行内公式 $...$
    text = re.sub(r'\$([^$\n]+)\$',
                  lambda m: degfix(latex_to_text(m.group(1))), text)
    #   4) 还原行内代码
    for i, code in enumerate(stash):
        text = text.replace('\x00%d\x00' % i, '`%s`' % code)

    lines = text.split('\n')
    out = []
    i = 0
    n = len(lines)
    while i < n:
        ln = lines[i]
        if not ln.strip():
            i += 1
            continue
        # 块级公式（预处理后的占位）
        if 'MATHBLOCK\x01' in ln:
            buf = []
            while i < n and '\x02' not in lines[i]:
                buf.append(lines[i].replace('MATHBLOCK\x01', ''))
                i += 1
            if i < n:
                buf.append(lines[i].replace('MATHBLOCK\x01', '').split('\x02')[0])
                i += 1
            out.append('<div class="mathblock">%s</div>' % esc(sub2uni(' '.join(buf)).strip()))
            continue
        # 分隔线
        if re.match(r'^\s*---+\s*$', ln):
            out.append('<hr>')
            i += 1
            continue
        # 标题
        m = re.match(r'^(#{1,4})\s+(.*)$', ln)
        if m:
            lvl = len(m.group(1))
            out.append('<h%d>%s</h%d>' % (lvl, inline(m.group(2)), lvl))
            i += 1
            continue
        # 公式块 $$...$$（已在预处理中转换，此处仅兜底）
        if ln.strip().startswith('$$'):
            buf = []
            body = ln.strip()[2:]
            if body.endswith('$$') and len(body) > 2:
                buf.append(body[:-2])
                i += 1
            else:
                i += 1
                while i < n and '$$' not in lines[i]:
                    buf.append(lines[i]); i += 1
                if i < n:
                    buf.append(lines[i].split('$$')[0]); i += 1
            out.append('<div class="mathblock">%s</div>'
                       % esc(sub2uni(degfix(latex_to_text(' '.join(x.strip() for x in buf))))))
            continue
        # 表格
        if ln.lstrip().startswith('|'):
            tbl = []
            while i < n and lines[i].lstrip().startswith('|'):
                tbl.append(lines[i].strip())
                i += 1
            rows = []
            for r in tbl:
                cells = [c.strip() for c in r.strip('|').split('|')]
                if all(re.match(r'^:?-{2,}:?$', c) for c in cells if c):
                    continue
                rows.append(cells)
            if rows:
                out.append('<table>')
                head, body = rows[0], rows[1:]
                out.append('<thead><tr>' + ''.join('<th>%s</th>' % inline(c) for c in head) + '</tr></thead>')
                if body:
                    out.append('<tbody>')
                    for r in body:
                        out.append('<tr>' + ''.join('<td>%s</td>' % inline(c) for c in r) + '</tr>')
                    out.append('</tbody>')
                out.append('</table>')
            continue
        # 引用
        if ln.lstrip().startswith('>'):
            buf = []
            while i < n and lines[i].lstrip().startswith('>'):
                buf.append(re.sub(r'^\s*>\s?', '', lines[i]))
                i += 1
            out.append('<blockquote>%s</blockquote>' % '<br>'.join(inline(x) for x in buf))
            continue
        # 列表
        if re.match(r'^\s*([-*+]|\d+\.)\s+', ln):
            ordered = bool(re.match(r'^\s*\d+\.\s+', ln))
            tag = 'ol' if ordered else 'ul'
            items = []
            while i < n and re.match(r'^\s*([-*+]|\d+\.)\s+', lines[i]):
                items.append(re.sub(r'^\s*([-*+]|\d+\.)\s+', '', lines[i]))
                i += 1
            out.append('<%s>%s</%s>' % (tag, ''.join('<li>%s</li>' % inline(x) for x in items), tag))
            continue
        # 普通段落
        buf = [ln]
        i += 1
        while i < n and lines[i].strip() and not re.match(r'^(#{1,4}\s|\s*[-*+]\s|\s*\d+\.\s|\s*>|\s*\||\s*---+\s*$|\s*\$\$)', lines[i]):
            buf.append(lines[i]); i += 1
        out.append('<p>%s</p>' % inline(' '.join(x.strip() for x in buf)))
    return '\n'.join(out)


def find_browser():
    cands = [
        r'C:\Program Files\Google\Chrome\Application\chrome.exe',
        r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe',
        r'C:\Program Files\Microsoft\Edge\Application\msedge.exe',
        r'C:\Program Files (x86)\Google\Chrome\Application\chrome.exe',
    ]
    for c in cands:
        if os.path.exists(c):
            return c
    raise SystemExit('未找到 Chrome/Edge')


def main():
    src, dst = sys.argv[1], sys.argv[2]
    title = sys.argv[3] if len(sys.argv) > 3 else os.path.splitext(os.path.basename(src))[0]
    text = open(src, encoding='utf-8').read()
    body = md_to_html(text)
    doc = ('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">'
           '<title>%s</title><style>%s</style></head><body>%s</body></html>'
           % (esc(title), CSS, body))
    htmlpath = os.path.splitext(dst)[0] + '.html'
    with open(htmlpath, 'w', encoding='utf-8') as f:
        f.write(doc)
    chrome = find_browser()
    with tempfile.TemporaryDirectory() as prof:
        cmd = [chrome, '--headless=new', '--disable-gpu', '--no-sandbox',
               '--no-pdf-header-footer', '--run-all-compositor-stages-before-draw',
               '--virtual-time-budget=8000',
               '--user-data-dir=' + prof,
               '--print-to-pdf=' + os.path.abspath(dst),
               'file:///' + os.path.abspath(htmlpath).replace('\\', '/')]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if r.returncode != 0:
            print('chrome stderr:', r.stderr[-2000:])
    print('HTML ->', htmlpath)
    print('PDF  ->', dst, os.path.getsize(dst), 'bytes' if os.path.exists(dst) else 'MISSING')


if __name__ == '__main__':
    main()
