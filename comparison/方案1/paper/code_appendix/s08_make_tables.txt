# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s08_make_tables.py
# 作用  : 把 data/ 下的机器可读结果转换为论文使用的 LaTeX 表格, 保证论文中的
#         每一个数字都直接来自实算输出, 避免手工抄写造成不一致。
# 产出  : paper/tables/tab_*.tex
# 运行  : python s08_make_tables.py
# =============================================================================
"""Convert machine-readable results into LaTeX tables for the paper."""

from __future__ import annotations

import io
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scenarios as sc        # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_TAB = os.path.join(BASE, "paper", "tables")


def save_tex(name, text):
    os.makedirs(DIR_TAB, exist_ok=True)
    path = os.path.join(DIR_TAB, name)
    with io.open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text + "\n")
    return path


def esc(s):
    """把数据文件里的普通字符串转义成 LaTeX 安全文本。

    必须处理 LaTeX 特殊字符, 尤其是:
      % -> \\%   否则会注释掉该行余下内容(包括 & 与 \\\\)，引发
                 "Extra alignment tab" / "Misplaced \\cr" 连锁报错;
      _ -> \\_   否则文本模式下报 "Missing $ inserted"。
    参数名形如 h_m+1%、tail_mean_1800s 时这两个坑必然触发。
    已经是数学模式的字段不要经过本函数。
    """
    out = str(s)
    for a, b in (("\\", "\\textbackslash{}"), ("&", "\\&"), ("%", "\\%"),
                 ("_", "\\_"), ("#", "\\#"), ("$", "\\$"),
                 ("{", "\\{"), ("}", "\\}")):
        out = out.replace(a, b)
    return out


def frame_table(caption, label, body, colspec, size="\\footnotesize", landscape=False):
    """把表格主体包进 table 环境。"""
    return "\n".join([
        "\\begin{table}[H]",
        "\\centering",
        f"\\caption{{{caption}}}",
        f"\\label{{{label}}}",
        size,
        f"\\begin{{tabular}}{{{colspec}}}",
        "\\toprule",
        body,
        "\\bottomrule",
        "\\end{tabular}",
        "\\end{table}",
    ])


def df_to_rows(df, fmt="%.4f", index_header=None, index_fmt="%s"):
    """把 DataFrame 转成 booktabs 行, 支持多级表头。"""
    cols = list(df.columns)
    header = " & ".join(([index_header] if index_header else []) + [str(c) for c in cols])
    lines = [header + " \\\\", "\\midrule"]
    for idx, row in df.iterrows():
        vals = []
        for v in row.to_numpy():
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                vals.append(fmt % v)
            else:
                vals.append(str(v))
        head = (index_fmt % idx) if index_header else None
        lines.append(" & ".join(([head] if head else []) + vals) + " \\\\")
    return "\n".join(lines)


def simple_table(csv_name, tex_name, caption, label, fmt="%.4f",
                 index_header=None, index_fmt="%s", colspec=None, size="\\small",
                 skip=0, rename=None, header_replace=None):
    df = pd.read_csv(os.path.join(sc.DIR_DATA, csv_name))
    if skip:
        df = df.iloc[skip:]
    if rename:
        df = df.rename(columns=rename)
    if index_header:
        df = df.set_index(df.columns[0])
    body = df_to_rows(df, fmt=fmt, index_header=index_header, index_fmt=index_fmt)
    if header_replace:
        for a, b in header_replace.items():
            body = body.replace(a, b)
    if colspec is None:
        ncol = df.shape[1] + (1 if index_header else 0)
        colspec = "@{}l" + "r" * (ncol - 1) + "@{}"
    save_tex(tex_name, frame_table(caption, label, body, colspec, size=size))
    return tex_name


def main():
    made = []

    # ------------------------------------------------ 四问结果表(表 1~表 6)
    made.append(simple_table("p1_table1_T.csv", "tab_p1_T.tex",
                             "表 1\\quad 30 分钟内药材的温度（单位：$^\\circ$C）",
                             "tab:p1T", index_header="时间/s",
                             index_fmt="%s"))
    made.append(simple_table("p1_table2_C.csv", "tab_p1_C.tex",
                             "表 2\\quad 30 分钟内药材的水分浓度（单位：kg/kg）",
                             "tab:p1C", index_header="时间/s", index_fmt="%s"))
    made.append(simple_table("p2_table3_T.csv", "tab_p2_T.tex",
                             "表 3\\quad 3 小时内药材的温度（单位：$^\\circ$C）",
                             "tab:p2T", index_header="时间/h", index_fmt="%s"))
    made.append(simple_table("p2_table4_C.csv", "tab_p2_C.tex",
                             "表 4\\quad 3 小时内药材的水分浓度（单位：kg/kg）",
                             "tab:p2C", index_header="时间/h", index_fmt="%s"))

    t5 = pd.read_csv(os.path.join(sc.DIR_DATA, "p3_table5_C.csv"))
    body = df_to_rows(t5.set_index(t5.columns[0]), index_header="时间/h",
                      index_fmt="%s")
    save_tex("tab_p3_C.tex", frame_table(
        "表 5\\quad 药材烘干过程的水分浓度（单位：kg/kg）", "tab:p3C", body,
        "@{}lrrrrr@{}"))
    made.append("tab_p3_C.tex")

    t6 = pd.read_csv(os.path.join(sc.DIR_DATA, "p4_table6_C.csv"))
    body = df_to_rows(t6.set_index(t6.columns[0]), index_header="时间/h",
                      index_fmt="%s")
    save_tex("tab_p4_C.tex", frame_table(
        "表 6\\quad 药材烘干过程的水分浓度（单位：kg/kg）", "tab:p4C", body,
        "@{}lrrrr@{}"))
    made.append("tab_p4_C.tex")

    # -------------------------------------------------------- 收敛性汇总表
    rows = []
    m = pd.read_csv(os.path.join(sc.DIR_DATA, "p1_mesh_convergence.csv"))
    for _, r in m.iterrows():
        rows.append(("问题 1 网格", f"$N$={int(r['N'])}",
                     f"{r['max|dC|']:.2e}", f"{r['max|dT|']:.2e}"))
    d = pd.read_csv(os.path.join(sc.DIR_DATA, "p1_dt_convergence.csv"))
    for _, r in d.iterrows():
        rows.append(("问题 1 时间步", f"$\\Delta t$={r['dt/s']:g} s", "--",
                     f"$T_c$={r['T_center']:.6f}\\,$^\\circ$C"))
    body_rows = ["指标 & 水平 & 含水率偏差 & 温度偏差 \\\\", "\\midrule"]
    for a, b, c, e in rows:
        body_rows.append(f"{a} & {b} & {c} & {e} \\\\")
    save_tex("tab_convergence_p1.tex", frame_table(
        "问题 1 的网格与时间步收敛证据（相对最细网格/最小步长的最大偏差）",
        "tab:conv1", "\n".join(body_rows), "@{}llrr@{}"))

    p23 = pd.read_csv(os.path.join(sc.DIR_DATA, "p23_convergence.csv"))
    body_rows = ["类别 & 参数 & 干燥时间 $t_3$ / h \\\\", "\\midrule"]
    for _, r in p23.iterrows():
        if r["kind"] == "dt":
            body_rows.append(f"时间步 & $N$=1600, $\\Delta t$={r['dt/s']:g} s & "
                             f"{r['t_dry/h']:.4f} \\\\")
        else:
            body_rows.append(f"网格 & $N$={int(r['dt/s'])}, $\\Delta t$=10 s & "
                             f"{r['t_dry/h']:.4f} \\\\")
    save_tex("tab_convergence_p3.tex", frame_table(
        "问题 3 干燥时间的双重收敛证据", "tab:conv3", "\n".join(body_rows),
        "@{}llr@{}"))

    p4 = pd.read_csv(os.path.join(sc.DIR_DATA, "p4_convergence.csv"))
    body_rows = ["类别 & 参数 & 干燥时间 $t_4$ / h \\\\", "\\midrule"]
    for _, r in p4.iterrows():
        if r["kind"] == "dt":
            body_rows.append(f"时间步 & $N$=1600, $\\Delta t$={r['value']:g} s & "
                             f"{r['t_dry/h']:.4f} \\\\")
        else:
            body_rows.append(f"网格 & $N$={int(r['value'])}, $\\Delta t$=10 s & "
                             f"{r['t_dry/h']:.4f} \\\\")
    save_tex("tab_convergence_p4.tex", frame_table(
        "问题 4 干燥时间的双重收敛证据", "tab:conv4", "\n".join(body_rows),
        "@{}llr@{}"))

    # ---------------------------------------------------------- 验证类表格
    va = pd.read_csv(os.path.join(sc.DIR_DATA, "verification_analytic.csv"))
    body_rows = ["场 & $N$ & 全部采样点最大误差 & 内部点最大误差 & 中心点误差 & "
                 "表面点误差 & 末单元中心当表面值的假误差 \\\\", "\\midrule"]
    for _, r in va.iterrows():
        nm = "传质" if r["field"] == "mass" else "传热"
        body_rows.append(f"{nm} & {int(r['N'])} & {r['max_abs_err']:.3e} & "
                         f"{r['max_err_interior']:.3e} & {r['err_center_last']:.3e} & "
                         f"{r['err_surface_last']:.3e} & "
                         f"{r['err_lastcell_as_surface']:.3e} \\\\")
    save_tex("tab_analytic.tex", frame_table(
        "常物性算例下数值解与 Bessel 级数解析解的比较", "tab:analytic",
        "\n".join(body_rows), "@{}lrrrrrr@{}"))

    vc = pd.read_csv(os.path.join(sc.DIR_DATA, "verification_conservation.csv"))
    body_rows = ["算例 & 水量收支相对闭合误差 & $\\max C$ 单调不增 & "
                 "$C\\ge0$ & 温度不越界 \\\\", "\\midrule"]
    for _, r in vc.iterrows():
        body_rows.append(f"{esc(r['case'])} & {r['rel_closure_err']:.2e} & "
                         f"{'是' if r['maxC_monotone'] else '否'} & "
                         f"{'是' if r['minC_positive'] else '否'} & "
                         f"{'是' if r['T_within_air'] else '否'} \\\\")
    save_tex("tab_conservation.tex", frame_table(
        "守恒性与物理界检验", "tab:cons", "\n".join(body_rows), "@{}lrrrr@{}"))

    if os.path.exists(os.path.join(sc.DIR_DATA, "verification_scheme.csv")):
        vs = pd.read_csv(os.path.join(sc.DIR_DATA, "verification_scheme.csv"))
        body_rows = ["时间格式 & $\\Delta t$ / s & 中心温度 / $^\\circ$C & "
                     "表面含水率 / (kg/kg) \\\\", "\\midrule"]
        for _, r in vs.iterrows():
            body_rows.append(f"{r['method']} & {r['dt/s']:g} & "
                             f"{r['T_center']:.6f} & {r['C_surface']:.6f} \\\\")
        save_tex("tab_scheme.tex", frame_table(
            "两种时间积分格式的互证（问题 1，$N=800$）", "tab:scheme",
            "\n".join(body_rows), "@{}lrrr@{}"))

    if os.path.exists(os.path.join(sc.DIR_DATA, "inverse_estimate.csv")):
        e = pd.read_csv(os.path.join(sc.DIR_DATA, "inverse_estimate.csv")).iloc[0]
        body_rows = [
            "项目 & 真值 & 估计值 & 标准差 & 95\\% 置信区间 \\\\", "\\midrule",
            f"$h_m$ / ($10^{{-7}}$ m/s)，单参数反演(E1，$D$ 固定) & "
            f"{e['hm_true']*1e7:.4f} & {e['hm_fit']*1e7:.4f} & "
            f"{e['hm_sd']*1e7:.4f} & "
            f"[{e['hm_lo95']*1e7:.4f}, {e['hm_hi95']*1e7:.4f}] \\\\",
            "\\midrule",
            "项目 & 观测窗口 & $\\mathrm{Fo}=Dt/R^2$ & 参数相关系数 & \\\\",
            "\\midrule",
            f"$h_m$ 与 $D$ 尺度 & 30 min (E2) & $\\approx0.023$ & "
            f"{e['corr']:+.3f} & \\\\",
            f"$h_m$ 与 $D$ 尺度 & 24 h (E3) & $\\approx1.1$ & "
            f"{e['corr_long']:+.3f} & \\\\",
            "\\midrule",
            f"\\multicolumn{{5}}{{l}}{{残差标准差 = {e['resid_sd']:.4e} "
            f"(注入噪声 {e['noise_sd']:.4e}, 比值 {e['resid_sd']/e['noise_sd']:.3f}); "
            f"残差一阶自相关 = {e['resid_ac1']:+.3f}}} \\\\",
        ]
        save_tex("tab_inverse.tex", frame_table(
            "合成数据回归：$h_m$ 的反演与两参数可辨识性", "tab:inv",
            "\n".join(body_rows), "@{}lrrrr@{}"))

    if os.path.exists(os.path.join(sc.DIR_DATA, "sensitivity_oat.csv")):
        oat = pd.read_csv(os.path.join(sc.DIR_DATA, "sensitivity_oat.csv"))
        body_rows = ["问题 & 扰动 & $\\Delta t_{\\mathrm{dry}}$ / h & 相对变化 / \\% \\\\",
                     "\\midrule"]
        for _, r in oat.iterrows():
            if r["param"] == "base":
                continue
            body_rows.append(f"{r['problem']} & {esc(r['param'])} & "
                             f"{r['delta_h']:+.4f} & {r['rel_pct']:+.3f} \\\\")
        save_tex("tab_sensitivity_oat.tex", frame_table(
            "关键参数单因素扰动对干燥时间的影响", "tab:oat",
            "\n".join(body_rows), "@{}llrr@{}"))

    if os.path.exists(os.path.join(sc.DIR_DATA, "sensitivity_plateau.csv")):
        pl = pd.read_csv(os.path.join(sc.DIR_DATA, "sensitivity_plateau.csv"))
        body_rows = ["平台口径 & $T_{\\mathrm{air}}$ / $^\\circ$C & "
                     "$C_{\\mathrm{air}}$ & $t_3$ / h & $t_4$ / h \\\\", "\\midrule"]
        piv = pl.pivot_table(index=["variant", "T_plateau", "C_plateau"],
                             columns="problem", values="t_dry_h").reset_index()
        for _, r in piv.iterrows():
            body_rows.append(f"{esc(r['variant'])} & {r['T_plateau']:.4f} & "
                             f"{r['C_plateau']:.6f} & {r['P3']:.4f} & {r['P4']:.4f} \\\\")
        save_tex("tab_plateau.tex", frame_table(
            "恒温干燥阶段平台边界口径的影响", "tab:plat",
            "\n".join(body_rows), "@{}lrrrr@{}"))

    if os.path.exists(os.path.join(sc.DIR_DATA, "sensitivity_mc_contrib.csv")):
        ct = pd.read_csv(os.path.join(sc.DIR_DATA, "sensitivity_mc_contrib.csv"))
        piv = ct.pivot_table(index="param", columns="problem",
                             values="share_pct").reset_index()
        names = {"h": "$h$", "h_m": "$h_m$", "D_scale": "$D$ 乘性尺度",
                 "T_plateau": "$T_{\\mathrm{air}}$ 平台值",
                 "C_plateau": "$C_{\\mathrm{air}}$ 平台值"}
        body_rows = ["输入参数 & 对 $t_3$ 方差的贡献 / \\% & "
                     "对 $t_4$ 方差的贡献 / \\% \\\\", "\\midrule"]
        for _, r in piv.iterrows():
            body_rows.append(f"{esc(names.get(r['param'], r['param']))} & "
                             f"{r['P3']:.1f} & {r['P4']:.1f} \\\\")
        save_tex("tab_mc_contrib.tex", frame_table(
            "Monte Carlo 不确定度传播的方差贡献分解", "tab:mc",
            "\n".join(body_rows), "@{}lrr@{}"))

    print("生成 LaTeX 表格 %d 个:" % len(made))
    for m in sorted(set(made)):
        print("  paper/tables/" + m)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
