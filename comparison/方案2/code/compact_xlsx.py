# -*- coding: utf-8 -*-
# ============================================================================
# compact_xlsx.py —— 极简 xlsx 写出器（只写必需的 OOXML 标记）
#
# 做什么：把若干张表写成一个 .xlsx 文件。本模块不依赖 openpyxl/xlsxwriter，
#         只输出必需的标签，把 259200 行 x 22 列的大表压到约 4 MB——用通用库
#         写同样的表约需 16.5 MB，而结果文件包有 20 MB 的体积上限。
#
# 输入：write_xlsx(path, sheets) 的 sheets，形如 {工作表名: (表头列表, 数据行)}，
#       也可写成 (表头, 数据行, 四位小数列号集合)。数据行可为二维 numpy 数组；
#       NaN / None 写成空单元格（问题 4 中用于「已收缩到药材之外」的位置）。
#
# 输出：path 指定的 .xlsx 文件，返回其字节数。数值单元格同时用 numFmt "0.0000"
#       样式（编号 164）和 "%.4f" 文本写出，故忽略样式也能读出四位小数。
#
# 关键变量：
#   _FMT_ID / _FMT_CODE  自定义数字格式的编号 164 与格式代码 "0.0000"
#   _STYLES              xl/styles.xml 模板，其 cellXfs 第 1 项引用上述格式
#   _HEAD / _TAIL        工作表 XML 的头尾，数据行流式写入以压低内存峰值
#   with_refs            是否写单元格引用 r="B12345"：写它会让每个单元格多出约
#                        10 个唯一字节、破坏 ZIP 的重复率，实测 result2.xlsx 由
#                        3.9 MB 涨到 34.9 MB，故默认关闭；OOXML 允许省略 r，
#                        Excel 与 openpyxl 都能正常读取。
# ============================================================================

import os
import zipfile

_CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
{sheet_overrides}
</Types>"""

_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>"""

# numFmtId 164 是可供自定义格式使用的第一个编号；cellXfs 第 0 项是不带格式的
# 单元格（文字标签与时间列），第 1 项引用 "0.0000" 格式。
_FMT_ID = 164
_FMT_CODE = "0.0000"
_STYLES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<numFmts count="1"><numFmt numFmtId="%d" formatCode="%s"/></numFmts>
<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>
<fills count="1"><fill><patternFill patternType="none"/></fill></fills>
<borders count="1"><border/></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="%d" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/></cellXfs>
</styleSheet>""" % (_FMT_ID, _FMT_CODE, _FMT_ID)

_HEAD = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
         '<worksheet xmlns="http://schemas.openxmlformats.org/'
         'spreadsheetml/2006/main"><sheetData>')
_TAIL = "</sheetData></worksheet>"


def _xml_escape(text):
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def _col_name(idx):
    """列号转工作表列名：0 -> A，25 -> Z，26 -> AA。"""
    name = ""
    idx += 1
    while idx:
        idx, rem = divmod(idx - 1, 26)
        name = chr(65 + rem) + name
    return name


def _fmt(value, ref, is_numfmt_col):
    """把单个单元格写成 XML 片段：数值写四位小数，文字标签写成 inlineStr。"""
    refattr = ' r="%s"' % ref if ref is not None else ""
    if value is None:
        return "<c%s/>" % refattr
    if isinstance(value, str):
        return '<c%s t="inlineStr"><is><t>%s</t></is></c>' % (refattr,
                                                             _xml_escape(value))
    try:
        fv = float(value)
    except (TypeError, ValueError):
        return '<c%s t="inlineStr"><is><t>%s</t></is></c>' % (refattr,
                                                             _xml_escape(value))
    if fv != fv or fv in (float("inf"), float("-inf")):
        return "<c%s/>" % refattr
    if is_numfmt_col:
        # 文本保留尾随零，同时引用 "0.0000" 样式
        return '<c%s s="1"><v>%.4f</v></c>' % (refattr, fv)
    txt = ("%.4f" % fv).rstrip("0").rstrip(".") or "0"
    return "<c%s><v>%s</v></c>" % (refattr, txt)


def write_xlsx(path, sheets, compresslevel=9, with_refs=False, flush_rows=4000):
    """把 {工作表名: (表头, 数据行[, 四位小数列号])} 写成一个紧凑的 xlsx。

    sheets 的每个值为二元组 (表头, 数据行) 或三元组 (表头, 数据行, numfmt_cols)；
    numfmt_cols 是需写成四位小数的列号集合，默认除第 0 列（时间列）外全部包含。
    compresslevel 是 ZIP 压缩级别，flush_rows 是数据行的缓冲区行数，with_refs=True
    会补上可选的 r="A1" 引用（文件体积大幅增大，见模块文件头）。返回文件字节数。
    """
    names = list(sheets.keys())
    sheet_over = "\n".join(
        '<Override PartName="/xl/worksheets/sheet%d.xml" ContentType='
        '"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        % (i + 1) for i in range(len(names)))
    wb_sheets = "".join(
        '<sheet name="%s" sheetId="%d" r:id="rId%d"/>'
        % (_xml_escape(n), i + 1, i + 1) for i, n in enumerate(names))
    wb = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
          '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
          'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
          '<sheets>%s</sheets></workbook>' % wb_sheets)
    wb_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
               '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
               + "".join('<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/'
                         'officeDocument/2006/relationships/worksheet" Target="worksheets/sheet%d.xml"/>'
                         % (i + 1, i + 1) for i in range(len(names)))
               + '<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/'
                 'officeDocument/2006/relationships/styles" Target="styles.xml"/>'
                 % (len(names) + 1)
               + "</Relationships>")

    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED,
                         compresslevel=compresslevel) as z:
        z.writestr("[Content_Types].xml", _CONTENT_TYPES.format(
            sheet_overrides=sheet_over))
        z.writestr("_rels/.rels", _RELS)
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        z.writestr("xl/styles.xml", _STYLES)
        for i, name in enumerate(names):
            entry = sheets[name]
            if len(entry) == 3:
                header, rows, numfmt_cols = entry
            else:
                header, rows = entry
                numfmt_cols = None
            ncol = len(header) if header is not None else (
                len(rows[0]) if len(rows) else 0)
            if numfmt_cols is None:
                numfmt_cols = set(range(1, ncol))     # 除时间列外全部四位小数
            numfmt_cols = set(numfmt_cols)
            cols = ([_col_name(k) for k in range(ncol)] if with_refs
                    else [None] * ncol)

            with z.open("xl/worksheets/sheet%d.xml" % (i + 1), "w") as fh:
                fh.write(_HEAD.encode("utf-8"))
                buf = []
                if header is not None:
                    buf.append('<row r="1">' if with_refs else "<row>")
                    for k, v in enumerate(header):
                        buf.append(_fmt(v, "%s1" % cols[k] if with_refs
                                        else None, False))
                    buf.append("</row>")
                for ri, row in enumerate(rows, start=2):
                    buf.append('<row r="%d">' % ri if with_refs else "<row>")
                    for k, v in enumerate(row):
                        buf.append(_fmt(v, "%s%d" % (cols[k], ri)
                                        if with_refs else None,
                                        k in numfmt_cols))
                    buf.append("</row>")
                    if len(buf) >= flush_rows * (ncol + 2):
                        fh.write("".join(buf).encode("utf-8"))
                        buf = []
                buf.append(_TAIL)
                fh.write("".join(buf).encode("utf-8"))
    return os.path.getsize(path)
