# -*- coding: utf-8 -*-
# ==============================================================================
# console_utf8.py —— 把控制台输出切到 UTF-8（Windows 默认 GBK 会把中文打成乱码）
# ==============================================================================
# 做什么：只做一件事——把当前进程的 stdout/stderr 重新配置成 UTF-8。Windows 控制台
#         默认按 GBK 解码，脚本 print 出来的中文进度信息会变成乱码，而一键复跑跑到
#         哪一步、哪一节报错，全靠这些中文行判断。
# 输入：不读数据文件、不依赖本项目任何模块，只调用 sys.stdout/stderr.reconfigure()。
# 输出：不写文件；唯一效果是改掉本进程这两个标准流的字符编码。
# 运行方式：不单独运行。在调用方脚本的最前面（任何 print 之前）写两行：
#             from console_utf8 import fix_console
#             fix_console()
#           它不参与任何数值计算，放在哪个脚本前面都不影响结果。
# 关键变量：无模块级常量；两个标准流的 errors 一律取 "replace"（原因见 fix_console）。
# 文件本身保持 UTF-8 无 BOM、行尾 LF、行尾不留空格，以便任何平台按纯文本读取。
# ==============================================================================
import sys


def fix_console():
    """把 stdout/stderr 切到 UTF-8；失败时静默跳过，不阻断主流程。

    参数与返回：均为空（就地重配两个标准流）。
    errors="replace" 是关键：即使下游终端不接受 UTF-8，也只是把个别字显示成问号，
    不会让整段复跑因为一次 print 抛 UnicodeEncodeError 而中断。
    """
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        # 需要 Python 3.7+ 才有的 reconfigure，且个别环境的重定向流不支持；
        # 这种情况退回平台默认编码即可，不影响任何计算结果。
        pass
