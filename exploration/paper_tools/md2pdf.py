# -*- coding: utf-8 -*-
"""md2pdf：Markdown(含 LaTeX 公式) -> HTML -> 无头 Edge 打印 -> PDF(带页码)

用法:
    python md2pdf.py <input.md> <output.pdf> [--title 标题] [--style paper|report]

依赖: markdown_it, PyMuPDF, Microsoft Edge/Chrome
"""
import io
import os
import re
import sys
import html as _html
import subprocess
import tempfile

import fitz
from markdown_it import MarkdownIt

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

EDGE = [r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe',
        r'C:\Program Files\Microsoft\Edge\Application\msedge.exe',
        r'C:\Program Files\Google\Chrome\Application\chrome.exe',
        r'C:\Program Files (x86)\Google\Chrome\Application\chrome.exe']

MATHJAX = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js'

CSS_COMMON = r"""
@page { size: A4; margin: 20mm 18mm 20mm 18mm; }
* { box-sizing: border-box; }
body { margin: 0; padding: 0; color: #111; -webkit-print-color-adjust: exact;
       print-color-adjust: exact; font-family: "SimSun", "Noto Serif SC", serif; }
h1 { font-family: "SimHei","Microsoft YaHei",sans-serif; font-size: 16pt; text-align: center;
     margin: 0 0 3mm 0; line-height: 1.5; }
h2 { font-family: "SimHei","Microsoft YaHei",sans-serif; font-size: 13pt; color: #000;
     margin: 6.5mm 0 2.5mm 0; page-break-after: avoid; }
h3 { font-family: "SimHei","Microsoft YaHei",sans-serif; font-size: 12pt; color: #000;
     margin: 4.5mm 0 2mm 0; page-break-after: avoid; }
h4 { font-family: "SimHei","Microsoft YaHei",sans-serif; font-size: 11pt; color: #222;
     margin: 3.5mm 0 1.5mm 0; page-break-after: avoid; }
p  { margin: 1.6mm 0; text-align: justify; line-height: 1.62; }
ul, ol { margin: 1.4mm 0 1.4mm 0; padding-left: 8mm; }
li { margin: .8mm 0; line-height: 1.55; }
strong { font-family: "SimHei","Microsoft YaHei",sans-serif; font-weight: 600; }
em { font-style: italic; }
hr { border: none; border-top: 1px solid #999; margin: 5mm 0; }
img { max-width: 100%; }
mjx-container { overflow: visible !important; }
"""

CSS_PAPER = CSS_COMMON + r"""
body { font-size: 10.5pt; }
.abstract { margin: 4mm 0 3mm 0; }
.abstract h2 { font-size: 12pt; margin: 0 0 2mm 0; }
.keywords { margin-top: 2.5mm; }
table { border-collapse: collapse; width: 100%; margin: 2mm auto 4mm auto;
        font-size: 8.4pt; page-break-inside: avoid; font-family: "SimSun", serif; }
th { border-top: 1.4pt solid #000; border-bottom: .7pt solid #000; padding: 1.3mm 1mm;
     font-family: "SimHei","Microsoft YaHei",sans-serif; font-weight: 600; text-align: center; }
td { padding: 1.1mm 1mm; text-align: center; }
tbody tr:last-child td { border-bottom: 1.4pt solid #000; }
caption, .tcap { font-family: "SimHei","Microsoft YaHei",sans-serif; font-size: 9.5pt;
                 text-align: center; margin: 3mm 0 1mm 0; }
.caption { font-family: "SimHei","Microsoft YaHei",sans-serif; font-size: 9.5pt; text-align: center;
           margin: 1mm 0 3mm 0; page-break-after: avoid; }
blockquote { margin: 2mm 0; padding: 1.5mm 3mm; background: #f7f7f7;
             border-left: 3px solid #888; font-size: 9.5pt; }
code { font-family: Consolas, monospace; font-size: 9pt; background: #f2f2f2; padding: 0 .6mm; }
pre { background: #f6f6f6; border: .5pt solid #ccc; padding: 2mm; font-size: 8.4pt;
      white-space: pre-wrap; page-break-inside: avoid; }
pre code { background: none; }
.fig { text-align: center; margin: 3.5mm 0; page-break-inside: avoid; }
.fig img { max-width: 90%; }
img { max-width: 92%; }
"""

CSS_REPORT = CSS_COMMON + r"""
body { font-size: 10pt; font-family: "Microsoft YaHei","SimSun",sans-serif; }
h1 { border-bottom: 2.5px solid #1f4e79; color: #12395c; padding-bottom: 2.5mm; }
h2 { color: #12395c; border-left: 4px solid #1f4e79; background: #eef4fa;
     padding: 1.4mm 0 1.4mm 3mm; }
h3 { color: #1f4e79; }
table { border-collapse: collapse; width: 100%; margin: 2mm 0 3.5mm 0; font-size: 8.8pt; }
th { background: #1f4e79; color: #fff; font-weight: 600; padding: 1.4mm 1.5mm;
     border: 1px solid #7f9db9; text-align: left; }
td { padding: 1.2mm 1.5mm; border: 1px solid #b9c9d9; vertical-align: top; }
tbody tr:nth-child(even) td { background: #f6f9fc; }
blockquote { margin: 2mm 0; padding: 1.6mm 3mm; background: #fffbf0;
             border-left: 4px solid #e0a83c; color: #5a4a20; font-size: 9.5pt; }
code { font-family: Consolas, monospace; font-size: 8.8pt; background: #f2f4f7; color: #b0004a;
       padding: .2mm .8mm; border-radius: 2px; }
pre { background: #f6f8fa; border: 1px solid #d6dde5; border-radius: 3px; padding: 2mm;
      font-size: 8.2pt; white-space: pre-wrap; page-break-inside: avoid; }
"""

TEMPLATE = """<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>@@TITLE@@</title>
<script>
window.MathJax = {
  tex: { inlineMath: [['$','$']], displayMath: [['$$','$$']],
          packages: {'[+]': ['ams','boldsymbol','noerrors']} },
  svg: { fontCache: 'local', scale: 0.98 },
  options: { enableMenu: false },
  startup: { typeset: true }
};
</script>
<script src="@@MATHJAX@@" id="MathJax-script"></script>
<style>@@CSS@@</style></head>
<body>
@@BODY@@
</body></html>
"""


def protect_math(md):
    """把 $$..$$ / $..$ 换成占位符，避免 markdown 解析破坏公式"""
    store = []

    def keep(m):
        store.append(m.group(0))
        return '@@MATH%d@@' % (len(store) - 1)

    md = re.sub(r'\$\$.+?\$\$', keep, md, flags=re.S)
    md = re.sub(r'(?<!\\)\$(?!\$).+?(?<!\\)\$', keep, md, flags=re.S)
    return md, store


def restore_math(html, store):
    def put(m):
        s = store[int(m.group(1))]
        if s.startswith('$$'):
            return '<div style="text-align:center;margin:3mm 0;">%s</div>' % s
        return s
    return re.sub(r'@@MATH(\d+)@@', put, html)


def md_to_html(md, base_dir='.'):
    md_p, store = protect_math(md)
    parser = MarkdownIt('commonmark', {'html': True, 'linkify': False})
    parser.enable(['table', 'strikethrough'])
    body = parser.render(md_p)
    body = restore_math(body, store)

    # 图片路径：相对 Markdown 所在目录解析为绝对 file:// 路径
    def fix_src(m):
        src = m.group(1)
        if re.match(r'^(https?:|data:|file:)', src):
            return m.group(0)
        p = os.path.abspath(os.path.join(base_dir, src))
        if not os.path.exists(p):
            raise RuntimeError('插图不存在: %s' % p)
        return 'src="file:///%s"' % p.replace('\\', '/')
    body = re.sub(r'src="([^"]+)"', fix_src, body)

    # 图片 + 紧随其后的题注 -> 不可分割的图块
    body = re.sub(r'<p>(<img[^>]*/?>)</p>\s*<p><strong>(图\s*\d[^<]*)</strong></p>',
                  r'<div class="fig">\1<div class="caption">\2</div></div>', body)
    # 表格题注：单独一行的 **表 x ...** -> caption
    body = re.sub(r'<p><strong>((?:表|图)\s*\d[^<]*)</strong></p>',
                  r'<div class="caption">\1</div>', body)
    return body


def find_browser():
    for p in EDGE:
        if os.path.exists(p):
            return p
    raise RuntimeError('未找到 Edge/Chrome')


def print_pdf(html_path, pdf_path, wait_ms=30000):
    exe = find_browser()
    html_path = os.path.abspath(html_path)
    pdf_path = os.path.abspath(pdf_path)
    cmd = [exe, '--headless=new', '--disable-gpu', '--no-sandbox', '--no-pdf-header-footer',
           '--run-all-compositor-stages-before-draw',
           '--virtual-time-budget=%d' % wait_ms,
           '--print-to-pdf-no-header', '--print-to-pdf=%s' % pdf_path,
           'file:///' + html_path.replace('\\', '/')]
    r = subprocess.run(cmd, capture_output=True, timeout=600)
    if not os.path.exists(pdf_path):
        raise RuntimeError('打印失败: %s' % (r.stderr[-800:].decode('utf-8', 'replace')))


def stamp_page_numbers(pdf_path, prefix='', skip_first=False):
    """用 PyMuPDF 在页脚居中盖页码"""
    doc = fitz.open(pdf_path)
    n = len(doc)
    for i, page in enumerate(doc):
        if skip_first and i == 0:
            continue
        txt = '第 %d 页 / 共 %d 页' % (i + 1, n)
        w = page.rect.width
        page.insert_text((w / 2 - 22, page.rect.height - 26), txt,
                         fontsize=9, fontname='china-s', color=(0.25, 0.25, 0.25))
    doc.save(pdf_path + '.tmp', garbage=3, deflate=True)
    doc.close()
    os.replace(pdf_path + '.tmp', pdf_path)
    return n


def convert(md_path, pdf_path, title=None, style='paper', page_numbers=True):
    raw = open(md_path, 'rb').read()
    if raw.count(0) > 0:
        raise RuntimeError('输入 Markdown 已损坏（含 %d 个 NUL 字节）: %s' % (raw.count(0), md_path))
    md = raw.decode('utf-8')
    if len(md) < 50:
        raise RuntimeError('输入 Markdown 内容异常（仅 %d 字符）' % len(md))
    title = title or os.path.splitext(os.path.basename(md_path))[0]
    css = CSS_PAPER if style == 'paper' else CSS_REPORT
    body = md_to_html(md, os.path.dirname(os.path.abspath(md_path)))
    html = (TEMPLATE.replace('@@TITLE@@', _html.escape(title))
            .replace('@@MATHJAX@@', MATHJAX)
            .replace('@@CSS@@', css)
            .replace('@@BODY@@', body))
    tmp = os.path.splitext(pdf_path)[0] + '.html'
    with open(tmp, 'w', encoding='utf-8', newline='\n') as f:
        f.write(html)
    print_pdf(tmp, pdf_path)
    n = stamp_page_numbers(pdf_path) if page_numbers else len(fitz.open(pdf_path))
    print('%s -> %s  (%d 页, %.0f KB)' % (os.path.basename(md_path), os.path.basename(pdf_path),
                                         n, os.path.getsize(pdf_path) / 1024))
    return n


if __name__ == '__main__':
    args = sys.argv[1:]
    style = 'paper'
    if '--style' in args:
        i = args.index('--style')
        style = args[i + 1]
        del args[i:i + 2]
    title = None
    if '--title' in args:
        i = args.index('--title')
        title = args[i + 1]
        del args[i:i + 2]
    convert(args[0], args[1], title, style)
