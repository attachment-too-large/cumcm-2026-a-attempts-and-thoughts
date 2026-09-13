# -*- coding: utf-8 -*-
"""把 md/ 下的分章 Markdown 按序合并，再调用 build.py 生成 PDF。"""
import os, sys, subprocess, glob, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__))
MD = os.path.join(HERE, 'md')

ORDER = ['01_封面目录预备.md', '02_几何与传热.md', '03_传质与边界.md',
         '04_物质坐标.md', '05_物性与无量纲.md', '06_数值格式.md',
         '06b_独立数值验证.md', '07_FAQ附录.md']


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        HERE, '药材烘干模型_推导逐步骤详解.pdf')
    parts = []
    for name in ORDER:
        p = os.path.join(MD, name)
        if not os.path.exists(p):
            print('  [跳过] %s 不存在' % name)
            continue
        parts.append(open(p, encoding='utf-8').read().rstrip())
        print('  [合并] %-24s %6d 字符' % (name, len(parts[-1])))
    text = '\n\n'.join(parts) + '\n'
    dst_md = os.path.join(HERE, '_merged.md')
    open(dst_md, 'w', encoding='utf-8').write(text)
    print('合并后 %d 字符 -> %s' % (len(text), dst_md))
    r = subprocess.run([sys.executable, os.path.join(HERE, 'build.py'), dst_md, out,
                        '药材烘干模型 推导逐步骤详解'], cwd=HERE)
    raise SystemExit(r.returncode)


if __name__ == '__main__':
    main()
