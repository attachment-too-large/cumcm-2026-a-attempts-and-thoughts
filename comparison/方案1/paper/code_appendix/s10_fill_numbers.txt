# -*- coding: utf-8 -*-
# =============================================================================
# 文件名: s10_fill_numbers.py
# 作用  : 把论文源文件中的占位符 %%XXX%% 替换为实算结果中的真实数字,
#         保证论文正文与结果文件完全一致, 避免手工抄写错误。
#         所有数值均从 data/ 下的 csv/txt 读取, 不在此文件中硬编码。
# 产出  : 直接改写 paper/sections/*.tex (可用 --dry-run 只打印不改写)
# 运行  : python s10_fill_numbers.py [--dry-run]
# =============================================================================
"""Fill the paper's numeric placeholders from the actual computation outputs."""

from __future__ import annotations

import io
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scenarios as sc        # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SECT = os.path.join(BASE, "paper", "sections")


def read_csv(name):
    p = os.path.join(sc.DIR_DATA, name)
    if not os.path.exists(p):
        return None
    return pd.read_csv(p)


def build_map():
    """构造 占位符 -> 字符串 的映射。缺失数据时给出显式占位, 不编造。

    占位符定界符用 @@NAME@@ 而不是 %%NAME%%: 后者含 LaTeX 的注释字符 %,
    会把该行余下的内容(包括闭合花括号)全部注释掉, 引发 "Runaway argument"。
    """
    m = {}
    miss = []

    def put(key, value):
        m[key] = value

    # ---- 问题 3 / 问题 4 的干燥时间 ----
    for prob, tag in (("p23", "3"), ("p4", "4")):
        txt_path = os.path.join(sc.DIR_DATA, f"{prob}_summary.txt")
        if not os.path.exists(txt_path):
            miss.append(f"{prob}_summary.txt")
            continue
        txt = io.open(txt_path, encoding="utf-8").read()
        import re
        r = re.search(r"t[34]\s*=\s*([0-9.]+)\s*s\s*=\s*([0-9.]+)\s*h", txt)
        if r:
            put(f"T{tag}S", f"{float(r.group(1)):.1f}")
            put(f"T{tag}H", f"{float(r.group(2)):.2f}")
            put(f"T{tag}DAY", f"{float(r.group(2))/24:.2f}")
        r2 = re.search(r"已写出\s+results/result" + tag + r"\.xlsx:\s*(\d+)\s*行", txt)
        if r2:
            put(f"N{tag}ROWS", r2.group(1))
    if "T3H" in m and "T4H" in m:
        put("DT43", f"{float(m['T3H'])-float(m['T4H']):.2f}")

    # 结果文件行数(直接数 csv, 比解析文本更稳)
    for f, key in (("p3_table5_C.csv", None),):
        pass

    # ---- 解析解对照 ----
    va = read_csv("verification_analytic.csv")
    if va is not None:
        top = va[va["N"] == va["N"].max()]
        put("ERRANA", f"{top['max_abs_err'].max():.2e}")
        put("ERRFAKE", f"{top['err_lastcell_as_surface'].max():.2e}")
    else:
        miss.append("verification_analytic.csv")

    # ---- 网格与时间步收敛 ----
    mc = read_csv("p1_mesh_convergence.csv")
    if mc is not None:
        row = mc[mc["N"] == 1600]
        if len(row):
            put("MESHC", f"{row['max|dC|'].iloc[0]:.2e}")
            put("MESHT", f"{row['max|dT|'].iloc[0]:.2e}")
    else:
        miss.append("p1_mesh_convergence.csv")

    ri = read_csv("p1_richardson.csv")
    dc_tab = read_csv("p1_dt_convergence.csv")
    if dc_tab is not None:
        # 直接从原始收敛序列重算观测阶与 Richardson 外推,
        # 观测阶 p = ln(e_prev/e_last)/ln(dt 差分之比), 不假设步长折半。
        v = dc_tab["T_center"].to_numpy(float)
        d = dc_tab["dt/s"].to_numpy(float)
        e_last = abs(v[-2] - v[-1])
        e_prev = abs(v[-3] - v[-2])
        r = (d[-3] - d[-2]) / (d[-2] - d[-1])
        p = np.log(e_prev / e_last) / np.log(r)
        ext = v[-1] - e_last / (r ** p - 1.0)
        put("ORDER1", f"{p:.2f}")
        put("RICH", f"{ext:.6f}")
        put("RICHBIAS", f"{abs(v[-1]-ext):.2e}")
    elif ri is not None:
        row = ri[ri["Unnamed: 0"] == "T_center"]
        if len(row):
            put("ORDER1", f"{float(row['observed_order'].iloc[0]):.2f}")
            put("RICH", f"{float(row['extrapolated'].iloc[0]):.6f}")
            put("RICHBIAS", f"{float(row['bias_at_base'].iloc[0]):.2e}")
    else:
        miss.append("p1_dt_convergence.csv")

    # ---- 两种时间格式互证 ----
    vs = read_csv("verification_scheme.csv")
    if vs is not None:
        be = vs[vs["method"] == "BE"]["T_center"].to_numpy(float)
        cn = vs[vs["method"] == "CN"]["T_center"].to_numpy(float)
        put("SCHEME", f"{abs(be.min() - cn.min()):.2e} $^\\circ$C")
    else:
        miss.append("verification_scheme.csv")

    # ---- 合成数据反演 ----
    ie = read_csv("inverse_estimate.csv")
    if ie is not None:
        e = ie.iloc[0]
        put("NOISE", f"{float(e['noise_sd']):.4f}")
        put("CORR", f"{float(e['corr']):+.3f}")
        put("CORRLONG", f"{float(e['corr_long']):+.3f}")
        put("HMFIT", f"{float(e['hm_fit'])*1e7:.4f}")
        put("HMSD", f"{float(e['hm_sd'])*1e7:.4f}")
        put("HMLO", f"{float(e['hm_lo95'])*1e7:.4f}")
        put("HMHI", f"{float(e['hm_hi95'])*1e7:.4f}")
        put("HMIN", "落在区间内" if float(e['hm_lo95']) <= float(e['hm_true']) <= float(e['hm_hi95'])
            else "不落在区间内")
        put("RESIDSD", f"{float(e['resid_sd']):.4e}")
        put("RESIDRATIO", f"{float(e['resid_sd'])/float(e['noise_sd']):.3f}")
        put("AC1", f"{float(e['resid_ac1']):+.3f}")
    else:
        miss.append("inverse_estimate.csv")

    # ---- 干物质守恒 ----
    dm = read_csv("drymass_consistency.csv")
    if dm is not None:
        m0, m1 = float(dm["m_dry"].iloc[0]), float(dm["m_dry"].iloc[-1])
        put("MDRY0", f"{m0:.6f} kg/m")
        put("MDRY1", f"{m1:.6f} kg/m")
        put("MDRYDEV", f"{100.0*(m1-m0)/m0:+.2f}\\%")
    else:
        miss.append("drymass_consistency.csv")

    # ---- 平台口径灵敏度 ----
    pl = read_csv("sensitivity_plateau.csv")
    if pl is not None:
        for prob, tag in (("P3", "3"), ("P4", "4")):
            sub = pl[pl["problem"] == prob]
            lo, hi = float(sub["t_dry_h"].min()), float(sub["t_dry_h"].max())
            put(f"PLATRANGE{tag}", f"{hi-lo:.3f}")
            base = float(sub[sub["variant"] == "tail_mean"]["t_dry_h"].iloc[0])
            put(f"PLATRANGE{tag}P", f"{100.0*(hi-lo)/base:.2f}")
            if prob == "P3":
                fit = float(sub[sub["variant"] == "exp_fit_asymptote"]["t_dry_h"].iloc[0])
                put("PLATT3FIT", f"{fit:.2f}")
                put("PLATT3BASE", f"{base:.2f}")
                put("PLATDIFF", f"{'慢' if fit > base else '快'} {abs(fit-base):.2f}")
    else:
        miss.append("sensitivity_plateau.csv")

    # ---- 空气含湿量口径 ----
    ca = read_csv("sensitivity_cair.csv")
    if ca is not None:
        sub = ca[(ca["problem"] == "P3") & (ca["t_dry_h"].notna())].sort_values("C_air_scale")
        if len(sub):
            put("CAIRLO", f"{float(sub['t_dry_h'].iloc[0]):.2f}")
            put("CAIRHI", f"{float(sub['t_dry_h'].iloc[-1]):.2f}")
        else:
            put("CAIRLO", "n/a")
            put("CAIRHI", "n/a")
        un = ca[~ca["t_dry_h"].notna()]
        put("CAIRUNREACH", f"{float(un['C_plateau'].iloc[0]):.3f}"
            if len(un) else "无")
    else:
        miss.append("sensitivity_cair.csv")

    # ---- 噪声注入 ----
    nz = read_csv("sensitivity_noise.csv")
    if nz is not None:
        for prob, tag in (("P3", "3"), ("P4", "4")):
            sub = nz[nz["problem"] == prob]
            put(f"NOISET{tag}", f"{float(sub['t_dry_h'].std(ddof=1)):.4f}")
    else:
        miss.append("sensitivity_noise.csv")

    # ---- Monte Carlo ----
    mcm = read_csv("sensitivity_montecarlo.csv")
    if mcm is not None:
        for prob, tag in (("P3", "3"), ("P4", "4")):
            sub = mcm[mcm["problem"] == prob]["t_dry_h"].to_numpy(float)
            put(f"MCMEAN{tag}", f"{sub.mean():.2f}")
            put(f"MCSD{tag}", f"{sub.std(ddof=1):.3f}")
            put(f"MCLO{tag}", f"{np.percentile(sub, 2.5):.2f}")
            put(f"MCHI{tag}", f"{np.percentile(sub, 97.5):.2f}")
    else:
        miss.append("sensitivity_montecarlo.csv")

    ctc = read_csv("sensitivity_mc_contrib.csv")
    if ctc is not None:
        names = {"h": "$h$（对流换热系数）", "h_m": "$h_m$（对流传质系数）",
                 "D_scale": "$D$ 的乘性尺度", "T_plateau": "环境温度平台值",
                 "C_plateau": "环境含湿量平台值"}
        sub = ctc[ctc["problem"] == "P3"]
        top = sub.loc[sub["share_pct"].idxmax()]
        put("MCDOM", names.get(top["param"], top["param"]))
    else:
        miss.append("sensitivity_mc_contrib.csv")

    return m, miss


def main():
    dry = "--dry-run" in sys.argv
    values, missing = build_map()
    files = [f for f in sorted(os.listdir(SECT)) if f.endswith(".tex")]
    replaced_total = 0
    leftover = {}
    for f in files:
        p = os.path.join(SECT, f)
        s = io.open(p, encoding="utf-8").read()
        n = 0
        for k, v in values.items():
            tok = "@@" + k + "@@"
            if tok in s:
                n += s.count(tok)
                s = s.replace(tok, v)
        left = [t for t in set(__import__("re").findall(r"@@[A-Z0-9]+@@", s))]
        if left:
            leftover[f] = sorted(left)
        if n and not dry:
            io.open(p, "w", encoding="utf-8", newline="\n").write(s)
        replaced_total += n
    print(f"替换占位符 {replaced_total} 处, 涉及 {len(files)} 个文件")
    if missing:
        print("[缺少数据源]", ", ".join(missing))
    if leftover:
        print("[仍有未填充的占位符]")
        for f, toks in leftover.items():
            print("  ", f, toks)
    else:
        print("[OK] 无残留占位符")
    if dry:
        print("(dry-run, 未改写文件)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
