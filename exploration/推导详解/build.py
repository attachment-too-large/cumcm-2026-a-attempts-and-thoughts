# -*- coding: utf-8 -*-
"""Markdown(含 LaTeX) -> HTML(KaTeX 预渲染) -> Chrome 打印 -> PDF(自动目录 + 页脚页码)。

用法:  python build.py 详解.md 输出.pdf "标题"

两遍编译:
  pass1  目录页码填 "?"，打印 PDF，用 PyMuPDF 按字号定位各标题真实页码
  pass2  把真实页码写回目录（列宽固定，版面不位移），再打印并补页脚页码
"""
import sys, io, os, re, html, json, base64, shutil, subprocess, tempfile

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

HERE = os.path.dirname(os.path.abspath(__file__))
KATEX_DIR = os.path.join(HERE, 'katex')
KATEX_JS = os.path.join(KATEX_DIR, 'katex.min.js')
RENDER_JS = os.path.join(HERE, 'render_math.js')
BUILD_DIR = os.path.join(HERE, '_build')
TOC_MARK = '<!--TOC-->'

CSS = r"""
@page { size: A4; margin: 17mm 15mm 16mm 15mm; }
* { box-sizing: border-box; }
html { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
body { font-family: "Microsoft YaHei", "MS YaHei", "SimSun", sans-serif;
       font-size: 10.4pt; line-height: 1.82; color: #16191d; margin: 0; padding: 0; }
p { margin: 1.6mm 0; text-align: justify; }
h1 { font-size: 20pt; text-align: center; margin: 0 0 3mm 0; color: #10375c; letter-spacing: .5px; }
h2 { font-size: 14.5pt; color: #fff; background: #1f4e79; margin: 7mm 0 3mm 0;
     padding: 2mm 3mm; border-radius: 1.5mm; page-break-after: avoid; }
h3 { font-size: 12.4pt; color: #10375c; margin: 5mm 0 2mm 0; padding-left: 2.5mm;
     border-left: 4px solid #2e75b6; page-break-after: avoid; }
h4 { font-size: 11pt; color: #1f4e79; margin: 4mm 0 1.5mm 0; page-break-after: avoid; }
h5 { font-size: 10.4pt; color: #333; margin: 3mm 0 1mm 0; page-break-after: avoid; }
strong { color: #0b2f4f; }
em { color: #7a4b00; font-style: normal; background: #fff8e6; padding: 0 .6mm; }
ul, ol { margin: 1.4mm 0; padding-left: 6.5mm; }
li { margin: .7mm 0; text-align: justify; }
li > ul, li > ol { margin: .5mm 0; }
code { font-family: Consolas, "Courier New", monospace; font-size: 9.3pt;
       background: #f1f3f6; padding: .2mm .8mm; border-radius: 1px; color: #a8064a; }
pre { background: #f7f8fa; border: 1px solid #dde3ea; border-left: 3px solid #2e75b6;
      padding: 2mm 3mm; margin: 2.5mm 0; page-break-inside: avoid; }
pre code { background: none; color: #24292f; font-size: 8.8pt; line-height: 1.5;
           white-space: pre-wrap; word-break: break-all; }
table { border-collapse: collapse; width: 100%; margin: 2.5mm 0 3.5mm 0;
        font-size: 9.2pt; page-break-inside: avoid; }
th { background: #1f4e79; color: #fff; font-weight: 600; padding: 1.4mm 1.5mm;
     border: 1px solid #7f9db9; text-align: left; }
td { padding: 1.2mm 1.5mm; border: 1px solid #b9c9d9; vertical-align: top; }
tbody tr:nth-child(even) td { background: #f6f9fc; }
blockquote { margin: 2.5mm 0; padding: 1.8mm 3mm; background: #fffbf0;
             border-left: 4px solid #e0a83c; color: #5a4a20; font-size: 9.9pt; }
blockquote p { margin: .8mm 0; }
hr { border: none; border-top: 1px dashed #b9c9d9; margin: 4.5mm 0; }

.katex { font-size: 1.045em; }
.katex-display { margin: 2.6mm 0 !important; }
.katex-display > .katex { font-size: 1.1em; }
.mblock { position: relative; page-break-inside: avoid; margin: 1mm 0; }
.eqno { position: absolute; right: 0; top: 50%; transform: translateY(-50%);
        font-size: 9.6pt; color: #10375c; font-family: Consolas, monospace; }

.box { margin: 2.8mm 0; padding: 2mm 3mm; border-radius: 1.2mm;
       page-break-inside: avoid; font-size: 9.9pt; }
.box p:first-child { margin-top: 0; }
.box p:last-child { margin-bottom: 0; }
.box .bt { font-weight: 700; display: block; margin-bottom: .8mm; }
.box.key  { background: #eaf5ec; border-left: 4px solid #2e8b57; }
.box.key .bt { color: #1c6b3e; }
.box.warn { background: #fdecec; border-left: 4px solid #c0392b; }
.box.warn .bt { color: #a5281b; }
.box.tip  { background: #eef4fa; border-left: 4px solid #2e75b6; }
.box.tip .bt { color: #1a5276; }
.box.ask  { background: #f4f0fa; border-left: 4px solid #6b4fa8; }
.box.ask .bt { color: #513a86; }

.cover { text-align: center; padding-top: 40mm; }
.cover .t1 { font-size: 27pt; font-weight: 700; color: #10375c; letter-spacing: 2px; line-height: 1.5; }
.cover .t2 { font-size: 14.5pt; color: #1f4e79; margin-top: 7mm; letter-spacing: 1px; }
.cover .t3 { font-size: 11pt; color: #5a6b7a; margin-top: 3mm; }
.cover .rule { width: 62mm; height: 2.5px; background: #1f4e79; margin: 8mm auto; }
.cover .meta { font-size: 10.5pt; color: #34455a; line-height: 2.3; margin-top: 14mm; }
.cover .foot { margin-top: 34mm; font-size: 9.5pt; color: #7b8794; line-height: 1.9; }
.pagebreak { page-break-after: always; }

.toc h2 { background: #10375c; }
.toc table { font-size: 9.9pt; margin-top: 3mm; }
.toc td { border: none; padding: .85mm 1mm; }
.toc tr:nth-child(even) td { background: none; }
.toc tr.l1 td { font-weight: 700; color: #10375c; padding-top: 2mm; }
.toc tr.l2 td:nth-child(1) { padding-left: 8mm; color: #26333f; font-weight: 400; }
.toc td.dots { border-bottom: 1px dotted #a9b7c6; width: auto; }
.toc td.pg { text-align: right; width: 11mm; font-family: Consolas, monospace; color: #10375c; }

figure { margin: 3mm 0; text-align: center; page-break-inside: avoid; }
figure img { max-width: 100%; }
figcaption { font-size: 9.2pt; color: #495a6b; margin-top: 1.2mm; text-align: center; }
"""

_CSS_CACHE = os.path.join(BUILD_DIR, 'katex_inline.css')


def katex_css_inline():
    """把 katex.min.css 内联，并把 woff2 字体转成 base64，避免任何相对路径问题。"""
    if os.path.exists(_CSS_CACHE):
        return open(_CSS_CACHE, encoding='utf-8').read()
    css = open(os.path.join(KATEX_DIR, 'katex.min.css'), encoding='utf-8').read()
    fonts_dir = os.path.join(KATEX_DIR, 'fonts')

    def font_uri(name):
        p = os.path.join(fonts_dir, name)
        with open(p, 'rb') as f:
            return 'data:font/woff2;base64,' + base64.b64encode(f.read()).decode('ascii')

    # 只保留 woff2（Chrome 全支持），丢掉 woff/ttf 备选，缩小体积
    def fix_src(m):
        block = m.group(0)
        w2 = re.search(r'url\(fonts/([^)]+\.woff2)\)', block)
        if not w2:
            return block
        uri = font_uri(w2.group(1))
        return re.sub(r'src:[^;}]*?(?=[;}])', 'src:url(%s) format("woff2")' % uri, block, count=1)

    css = re.sub(r'@font-face\{[^}]*\}', fix_src, css)
    css = re.sub(r'url\(fonts/[^)]+\.(woff|ttf)\)\s*format\("(woff|truetype)"\),?', '', css)
    os.makedirs(BUILD_DIR, exist_ok=True)
    open(_CSS_CACHE, 'w', encoding='utf-8').write(css)
    return css

# ------------------------------------------------------------------ 行内/块级
def esc(s):
    return html.escape(s, quote=False)


def img_data_uri(path):
    with open(path, 'rb') as f:
        b = f.read()
    ext = os.path.splitext(path)[1].lower().lstrip('.')
    mime = {'png': 'image/png', 'jpg': 'image/jpeg', 'jpeg': 'image/jpeg',
            'svg': 'image/svg+xml', 'gif': 'image/gif'}[ext]
    return 'data:%s;base64,%s' % (mime, base64.b64encode(b).decode('ascii'))


class Stash(object):
    def __init__(self):
        self.math, self.code = [], []

    def add_math(self, tex, display):
        self.math.append({'tex': tex, 'display': display})
        return '\x01M%d\x02' % (len(self.math) - 1)

    def add_code(self, txt):
        self.code.append(txt)
        return '\x00C%d\x00' % (len(self.code) - 1)


def restore_code(s, st):
    return re.sub(r'\x00C(\d+)\x00', lambda m: '<code>%s</code>' % esc(st.code[int(m.group(1))]), s)


def inline(s, st):
    s = re.sub(r'`([^`]+)`', lambda m: st.add_code(m.group(1)), s)
    s = re.sub(r'\$([^$\n]+?)\$', lambda m: st.add_math(m.group(1).strip(), False), s)
    s = esc(s)
    s = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', s)
    s = re.sub(r'(?<!\*)\*([^*\n]+)\*(?!\*)', r'<em>\1</em>', s)
    return restore_code(s, st)


def stash_fences(text, st):
    return re.sub(r'```([^\n`]*)\n(.*?)```', lambda m: st.add_code(m.group(2)), text, flags=re.S)


def parse(md, st):
    lines, out, i, n = md.split('\n'), [], 0, len(md.split('\n'))
    while i < n:
        ln = lines[i]
        s = ln.strip()
        if not s:
            i += 1
            continue
        if s == '\\newpage':
            out.append('<div class="pagebreak"></div>'); i += 1; continue
        if s == TOC_MARK:
            out.append(TOC_MARK); i += 1; continue
        m = re.match(r'^:{2,3}\s*(\w+)\s*(.*)$', s)
        if m:
            kind, title = m.group(1), m.group(2).strip()
            i += 1
            buf = []
            while i < n and lines[i].strip() not in ('::', ':::'):
                buf.append(lines[i]); i += 1
            i += 1
            if kind == 'rawhtml':
                out.append('\n'.join(buf))
                continue
            t = '<span class="bt">%s</span>' % inline(title, st) if title else ''
            out.append('<div class="box %s">%s%s</div>' % (kind, t, parse('\n'.join(buf), st)))
            continue
        if re.match(r'^-{3,}$', s):
            out.append('<hr>'); i += 1; continue
        if s.startswith('$$'):
            buf = []
            body = s[2:]
            if body.endswith('$$') and len(body) > 2:
                buf.append(body[:-2]); i += 1
            else:
                if body:
                    buf.append(body)
                i += 1
                while i < n and '$$' not in lines[i]:
                    buf.append(lines[i]); i += 1
                if i < n:
                    buf.append(lines[i].split('$$')[0]); i += 1
            tex = ' '.join(x.strip() for x in buf).strip()
            num, mm = '', re.search(r'\\tag\{([^{}]*)\}\s*$', tex)
            if mm:
                num, tex = mm.group(1), tex[:mm.start()].strip()
            ph = st.add_math(tex, True)
            out.append('<div class="mblock">%s%s</div>'
                       % (ph, '<div class="eqno">(%s)</div>' % esc(num) if num else ''))
            continue
        m = re.match(r'^\x00C(\d+)\x00$', s)
        if m:
            out.append('<pre><code>%s</code></pre>' % esc(st.code[int(m.group(1))])); i += 1; continue
        m = re.match(r'^(#{1,5})\s+(.*)$', ln)
        if m:
            lvl = len(m.group(1))
            out.append('<h%d>%s</h%d>' % (lvl, inline(m.group(2), st), lvl)); i += 1; continue
        if s.startswith('|'):
            tbl = []
            while i < n and lines[i].strip().startswith('|'):
                tbl.append(lines[i].strip()); i += 1
            rows = []
            for r in tbl:
                cells = split_row(r)
                if all(re.match(r'^:?-{2,}:?$', c) for c in cells if c):
                    continue
                rows.append(cells)
            if rows:
                head, body = rows[0], rows[1:]
                out.append('<table><thead><tr>' +
                           ''.join('<th>%s</th>' % inline(c, st) for c in head) + '</tr></thead>')
                if body:
                    out.append('<tbody>')
                    for r in body:
                        out.append('<tr>' + ''.join('<td>%s</td>' % inline(c, st) for c in r) + '</tr>')
                    out.append('</tbody>')
                out.append('</table>')
            continue
        if s.startswith('>'):
            buf = []
            while i < n and lines[i].strip().startswith('>'):
                buf.append(re.sub(r'^\s*>\s?', '', lines[i])); i += 1
            out.append('<blockquote>%s</blockquote>' % parse('\n'.join(buf), st)); continue
        m = re.match(r'^!\[([^\]]*)\]\(([^)]+)\)$', s)
        if m:
            cap, src = m.group(1), m.group(2)
            p = os.path.normpath(src if os.path.isabs(src) else os.path.join(HERE, src))
            uri = img_data_uri(p) if os.path.exists(p) else ''
            if not uri:
                print('!! 缺图:', p)
            out.append('<figure><img src="%s"><figcaption>%s</figcaption></figure>'
                       % (uri, inline(cap, st)))
            i += 1; continue
        if re.match(r'^(\s*)([-*+]|\d+\.)\s+', ln):
            i = parse_list(lines, i, out, st); continue
        buf = [ln]; i += 1
        while i < n and lines[i].strip() and not re.match(
                r'^(#{1,5}\s|\s*([-*+]|\d+\.)\s|\s*>|\s*\||-{3,}$|\$\$|:{2,3}|!\[|\\newpage|' + TOC_MARK + r')',
                lines[i]):
            buf.append(lines[i]); i += 1
        out.append('<p>%s</p>' % inline(' '.join(x.strip() for x in buf), st))
    return '\n'.join(out)


_TBL_MATH = re.compile(r'\$[^$\n]*\$')


def split_row(r):
    """按 | 切分表格行，但先保护 $...$ 里的竖线（\\big| 、\\right| 、|x| 等）。"""
    stash = []

    def keep(m):
        stash.append(m.group(0))
        return '\x03%d\x03' % (len(stash) - 1)

    body = _TBL_MATH.sub(keep, r.strip())
    if body.startswith('|'):
        body = body[1:]
    if body.endswith('|'):
        body = body[:-1]
    out = []
    for c in body.split('|'):
        c = c.strip()
        c = re.sub(r'\x03(\d+)\x03', lambda m: stash[int(m.group(1))], c)
        out.append(c)
    return out


def parse_list(lines, i, out, st):
    n = len(lines)
    base = len(re.match(r'^(\s*)', lines[i]).group(1))
    ordered = bool(re.match(r'^\s*\d+\.\s+', lines[i]))
    tag = 'ol' if ordered else 'ul'
    items = []
    while i < n:
        m = re.match(r'^(\s*)([-*+]|\d+\.)\s+(.*)$', lines[i])
        if m and len(m.group(1)) == base:
            items.append([m.group(3)]); i += 1
        elif items and lines[i].strip() and len(re.match(r'^(\s*)', lines[i]).group(1)) > base:
            items[-1].append(lines[i]); i += 1
        else:
            break
    out.append('<%s>' % tag)
    for it in items:
        body = inline(it[0], st) + (parse('\n'.join(it[1:]), st) if len(it) > 1 else '')
        out.append('<li>%s</li>' % body)
    out.append('</%s>' % tag)
    return i


def scan_headings(md):
    """收录 H2/H3 标题；frags 用于在 PDF 中定位页码，raw 用于目录显示（保留公式）。"""
    ents = []
    for m in re.finditer(r'^(#{2,3})\s+(.*)$', md, flags=re.M):
        raw = m.group(2).strip()
        frags = [re.sub(r'\s+', '', s) for s in re.split(r'\$[^$]*\$', raw)]
        frags = [s for s in frags if s]
        ents.append({'id': len(ents), 'level': len(m.group(1)), 'raw': raw,
                     'text': re.sub(r'\$[^$]*\$', '', raw).strip(),
                     'frags': frags})
    return ents


def toc_html(ents, pages, st):
    rows = ['<div class="toc">', '<h2>目　录</h2>', '<table>']
    for e in ents:
        pg = pages.get(str(e['id']), pages.get(e['id'], '?')) if pages else '?'
        cls = 'l1' if e['level'] == 2 else 'l2'
        rows.append('<tr class="%s"><td>%s</td><td class="dots"></td><td class="pg">%s</td></tr>'
                    % (cls, inline(e['raw'], st), esc(str(pg))))
    rows += ['</table>', '</div>']
    return '\n'.join(rows)


def build_html(md, title, ents, pages):
    st = Stash()
    body = parse(stash_fences(md, st), st)
    body = body.replace(TOC_MARK, toc_html(ents, pages, st))
    doc = ('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8"><title>%s</title>'
           '<style>%s</style><style>%s</style></head>'
           '<body>%s</body></html>' % (esc(title), katex_css_inline(), CSS, body))
    return doc, st.math


def chrome_print(html_path, pdf_path):
    from urllib.parse import quote
    cands = [r'C:\Program Files\Google\Chrome\Application\chrome.exe',
             r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe',
             r'C:\Program Files\Microsoft\Edge\Application\msedge.exe']
    exe = next(p for p in cands if os.path.exists(p))
    prof = tempfile.mkdtemp(prefix='dshpdf_')
    url = 'file:///' + quote(os.path.abspath(html_path).replace('\\', '/'))
    cmd = [exe, '--headless=new', '--disable-gpu', '--no-sandbox', '--no-first-run',
           '--disable-extensions', '--font-render-hinting=none',
           '--no-pdf-header-footer', '--run-all-compositor-stages-before-draw',
           '--virtual-time-budget=20000', '--user-data-dir=' + prof,
           '--print-to-pdf=' + os.path.abspath(pdf_path), url]
    r = subprocess.run(cmd, capture_output=True, timeout=600)
    out = (r.stdout or b'').decode('utf-8', 'replace')
    err = (r.stderr or b'').decode('utf-8', 'replace')
    shutil.rmtree(prof, ignore_errors=True)
    if not os.path.exists(pdf_path):
        print('chrome rc=%s\nURL=%s\nOUT:%s\nERR:%s' % (r.returncode, url, out[-1500:], err[-1500:]))
        raise SystemExit('打印 PDF 失败')


def scan_pages(pdf_path, ents):
    """按字号定位标题所在页码（1 基）。

    只认 H2/H3 的大号字，避开目录行（目录行是 9.9pt）。
    含行内公式的标题在文本提取时会断成多行，故把整页大号字拼成一个串再按片段顺序匹配。
    """
    import fitz
    doc = fitz.open(pdf_path)
    big = []
    for pg in doc:
        chunks = []
        for blk in pg.get_text('dict')['blocks']:
            for ln in blk.get('lines', []):
                sz = max([sp['size'] for sp in ln['spans']] or [0])
                if sz >= 11.9:
                    chunks.append(''.join(sp['text'] for sp in ln['spans']))
        big.append(re.sub(r'\s+', '', ''.join(chunks)))
    res, cursor = {}, 0
    for e in ents:
        frags = [re.sub(r'\s+', '', f) for f in (e.get('frags') or [e['text']])]
        frags = [f for f in frags if f]
        found = ''
        for pi in range(cursor, len(big)):
            c = big[pi]
            pos, ok = 0, bool(frags)
            for fr in frags:
                k = c.find(fr, pos)
                if k < 0:
                    ok = False
                    break
                pos = k + len(fr)
            if ok:
                found, cursor = pi + 1, pi
                break
        res[e['id']] = found
    doc.close()
    return res


def stamp_pages(pdf_path, out_path):
    import fitz
    doc = fitz.open(pdf_path)
    n = doc.page_count
    for k, pg in enumerate(doc):
        if k == 0:
            continue
        r, txt = pg.rect, '第 %d 页 / 共 %d 页' % (k + 1, n)
        pg.draw_line(fitz.Point(42, r.height - 33), fitz.Point(r.width - 42, r.height - 33),
                     color=(0.80, 0.85, 0.90), width=0.5)
        pg.insert_text((r.width / 2 - 32, r.height - 22), txt,
                       fontname='china-s', fontsize=8.6, color=(0.42, 0.47, 0.53))
    doc.save(out_path, garbage=3, deflate=True)
    doc.close()


def main():
    src, dst = sys.argv[1], sys.argv[2]
    title = sys.argv[3] if len(sys.argv) > 3 else '推导详解'
    md = open(src, encoding='utf-8').read()
    ents = scan_headings(md)
    print('标题条目:', len(ents))
    os.makedirs(BUILD_DIR, exist_ok=True)

    pages = None
    for pas in (1, 2):
        doc, math_items = build_html(md, title, ents, pages)
        hp = os.path.join(BUILD_DIR, 'doc%d.html' % pas)
        open(hp, 'w', encoding='utf-8').write(doc)
        json.dump(math_items, open(os.path.join(BUILD_DIR, 'math%d.json' % pas), 'w', encoding='utf-8'),
                  ensure_ascii=False)
        final = os.path.join(BUILD_DIR, 'final%d.html' % pas)
        r = subprocess.run(['node', RENDER_JS, hp, os.path.join(BUILD_DIR, 'math%d.json' % pas),
                            final, KATEX_JS], capture_output=True, timeout=600)
        out = (r.stdout or b'').decode('utf-8', 'replace').strip()
        err = (r.stderr or b'').decode('utf-8', 'replace').strip()
        print('[pass%d] %s %s' % (pas, out, err[:2000]))
        tmp = os.path.join(BUILD_DIR, 'pass%d.pdf' % pas)
        chrome_print(final, tmp)
        if pas == 1:
            pages = scan_pages(tmp, ents)
            hits = sum(1 for v in pages.values() if v)
            print('  定位页码 %d/%d' % (hits, len(ents)))
            for e in ents:
                if not pages.get(e['id']):
                    print('   未定位:', e['text'][:40])
            json.dump(pages, open(os.path.join(BUILD_DIR, 'pages.json'), 'w', encoding='utf-8'),
                      ensure_ascii=False)

    stamp_pages(os.path.join(BUILD_DIR, 'pass2.pdf'), dst)
    import fitz
    d = fitz.open(dst)
    print('PDF ->', dst, '页数 =', d.page_count, os.path.getsize(dst), 'bytes')
    d.close()


if __name__ == '__main__':
    main()
