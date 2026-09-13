# -*- coding: utf-8 -*-
"""把 Markdown 里被反引号包住的反斜杠命令改成普通文字，避免 LaTeX 的 \\texttt{\\cite} 炸栈。"""
import io
import os
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__))

REPL = {
    '`\\cite`': 'cite 命令',
    '`\\lstinputlisting{../code/*.py}`': 'lstinputlisting 引用',
    '`\\lstinputlisting`': 'lstinputlisting',
    '`\\cite` 出现': 'cite 命令出现',
}

for base in ('报告', '补丁'):
    d = os.path.join(HERE, base)
    if not os.path.isdir(d):
        continue
    for f in sorted(os.listdir(d)):
        if not f.endswith('.md'):
            continue
        p = os.path.join(d, f)
        s = io.open(p, encoding='utf-8').read()
        n0 = len(re.findall(r'`\\[A-Za-z]+', s))
        for a, b in REPL.items():
            s = s.replace(a, b)
        # 兜底：任何残留的 `\xxx` 一律去掉反引号
        s = re.sub(r'`(\\[A-Za-z]+)`', r'\1', s)
        n1 = len(re.findall(r'`\\[A-Za-z]+', s))
        io.open(p, 'w', encoding='utf-8', newline='\n').write(s)
        print('%-40s 原有 %d 处，剩余 %d 处' % (f, n0, n1))
