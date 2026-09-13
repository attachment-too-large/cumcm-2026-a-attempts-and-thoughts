# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import fitz
d = fitz.open(sys.argv[1])
print('pages:', d.page_count)
for i in range(d.page_count):
    print('=' * 20, 'PAGE', i + 1, '=' * 20)
    print(d[i].get_text('text'))
