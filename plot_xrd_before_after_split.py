#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制 XRD 0 点补正的前后对比，并可同时导出 CSV：
- 叠加图：<stem>_before_after.png
- 可选：分别输出：<stem>_before.png, <stem>_after.png
- 额外对比图：差分 <stem>_compare_diff.png；比值 <stem>_compare_ratio.png
- 另存 CSV（--save-csv）：
    <stem>_before.csv         [two_theta_deg, intensity]
    <stem>_after.csv          [two_theta_deg, intensity, delta_deg_used]
    <stem>_compare.csv        [two_theta_deg, intensity_before, intensity_after, diff, ratio]

自动读取 δ（优先级）：
1) 命令行 --delta
2) 命令行 --delta-report <REPORT.txt>（解析 '2THETA ORIGIN = ...'）
3) 输入 CSV 含列 delta_deg_used / delta / origin / twotheta_origin（取中位数）
4) 自动查找同名报告（<stem>.out.txt / <stem>.report.txt / <stem>_report.txt）
5) --delta-default（默认 0.0）

公式：2θ_after = 2θ_before - δ（单位：度）
"""
import argparse, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

NUM_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")

def read_text_guess_encoding(path: Path) -> str:
    for enc in ("utf-8","utf-8-sig","cp932","shift_jis","latin1"):
        try:
            return path.read_text(encoding=enc)
        except Exception:
            continue
    return path.read_bytes().decode("latin1", errors="ignore")

def parse_two_cols_from_text(txt: str):
    rows = []
    for ln in txt.splitlines():
        s = ln.replace(",", " ").replace("\t", " ").strip()
        if not s: continue
        ms = list(NUM_RE.finditer(s))
        if len(ms) >= 2:
            try:
                x = float(ms[0].group(0)); y = float(ms[1].group(0))
                rows.append((x, y))
            except: pass
    if not rows:
        raise ValueError("未在文本中解析到两列数值。")
    arr = np.asarray(rows, float)
    x, y = arr[:,0], arr[:,1]
    def is_theta(col): return np.all(np.diff(col) >= 0) and (-10 <= col.min() <= 190) and (col.max() <= 190)
    if is_theta(x): return x, y
    if is_theta(y): return y, x
    return x, y

def load_csv_with_guess(path: Path):
    for enc in ("utf-8","utf-8-sig","cp932","shift_jis","latin1"):
        try:
            return pd.read_csv(path, encoding=enc)
        except Exception:
            continue
    return pd.read_csv(path, encoding="latin1", on_bad_lines="skip")

def _normalize(s: str) -> str:
    return re.sub(r'[^a-z0-9]+','', s.lower())

def _get_series_by_spec(df: pd.DataFrame, spec):
    cols = list(df.columns)
    if isinstance(spec, int) or (isinstance(spec, str) and str(spec).isdigit()):
        idx = int(spec)
        if not (0 <= idx < len(cols)):
            raise KeyError(f"列序号越界: {idx}, 总列数={len(cols)}")
        return pd.to_numeric(df.iloc[:, idx], errors="coerce")
    low_map = {c.lower(): c for c in cols}
    if isinstance(spec, str) and spec.lower() in low_map:
        return pd.to_numeric(df[low_map[spec.lower()]], errors="coerce")
    norm = _normalize(str(spec))
    norm_map = {_normalize(c): c for c in cols}
    if norm in norm_map:
        return pd.to_numeric(df[norm_map[norm]], errors="coerce")
    raise KeyError(f"找不到列 {spec}. 可用列: {cols}")

def load_xy(path: Path, theta_spec=None, intensity_spec=None):
    ext = path.suffix.lower()
    if ext in {".ras",".txt",".dat",".xy"}:
        txt = read_text_guess_encoding(path)
        tt, ii = parse_two_cols_from_text(txt)
        order = np.argsort(tt, kind="mergesort")
        return tt[order], ii[order]
    df = load_csv_with_guess(path)
    if theta_spec is None or intensity_spec is None:
        low = {c.lower(): c for c in df.columns}
        def pick(names):
            for n in names:
                if n in low: return low[n]
            return None
        theta_col = pick(["two_theta_deg","two_theta","2theta","twotheta","theta","tth"])
        intensity_col = pick(["intensity","counts","y","i"])
        if theta_col is None or intensity_col is None:
            df_num = df.apply(pd.to_numeric, errors="coerce")
            num_cols = [c for c in df_num.columns if pd.api.types.is_numeric_dtype(df_num[c])]
            if len(num_cols) >= 2:
                theta_col = theta_col or num_cols[0]
                intensity_col = intensity_col or num_cols[1]
                df = df_num
            else:
                raise ValueError(f"自动识别失败。可用列：{list(df.columns)}。请用 --theta-col / --intensity-col 指定。")
        tt = pd.to_numeric(df[theta_col], errors="coerce")
        ii = pd.to_numeric(df[intensity_col], errors="coerce")
    else:
        tt = _get_series_by_spec(df, theta_spec)
        ii = _get_series_by_spec(df, intensity_spec)
    mask = np.isfinite(tt) & np.isfinite(ii)
    tt = tt[mask].to_numpy(); ii = ii[mask].to_numpy()
    order = np.argsort(tt, kind="mergesort")
    return tt[order], ii[order]

def estimate_step(tt):
    d = np.diff(tt)
    d = d[np.isfinite(d) & (d > 0)]
    if d.size == 0: return 0.02
    return float(np.median(d))

def common_grid(tt_a, tt_b, step=None):
    lo = max(float(np.min(tt_a)), float(np.min(tt_b)))
    hi = min(float(np.max(tt_a)), float(np.max(tt_b)))
    if step is None: step = estimate_step(tt_a)
    if hi - lo <= step: step = max(step, 1e-3)
    return np.arange(lo, hi + 1e-9, step, dtype=float)

def interp_to_grid(tt, I, grid):
    order = np.argsort(tt, kind="mergesort")
    return np.interp(grid, tt[order], I[order], left=np.nan, right=np.nan)

def parse_delta_from_report(report_path: Path):
    txt = read_text_guess_encoding(report_path)
    # 典型行： "2THETA ORIGIN =   -0.30000"
    m = re.search(r"2\s*THETA\s*ORIGIN\s*=\s*([+-]?\d+(?:\.\d*)?)", txt, flags=re.IGNORECASE)
    if m:
        return float(m.group(1))
    # 备选：若报告里有 “delta” 字样
    m = re.search(r"\bdelta\b\s*[:=]\s*([+-]?\d+(?:\.\d*)?)", txt, flags=re.IGNORECASE)
    if m:
        return float(m.group(1))
    raise ValueError(f"无法从报告解析 δ：{report_path}")

def try_delta_from_csv_meta(csv_path: Path):
    try:
        df = load_csv_with_guess(csv_path)
    except Exception:
        return None
    low = {c.lower(): c for c in df.columns}
    for key in ("delta_deg_used", "delta", "origin", "twotheta_origin"):
        if key in low:
            vals = pd.to_numeric(df[low[key]], errors="coerce").dropna().to_numpy()
            if vals.size:
                return float(np.median(vals))
    return None

def try_delta_from_sibling_reports(stem_base: Path):
    candidates = [
        stem_base.with_suffix(".out.txt"),
        stem_base.with_suffix(".report.txt"),
        stem_base.with_name(stem_base.name + "_report.txt"),
    ]
    for p in candidates:
        if p.exists():
            try:
                return parse_delta_from_report(p)
            except Exception:
                pass
    return None

def save_overlay(two_theta_raw, intensity, two_theta_after, delta, out_png, xlim=None, dpi=150):
    plt.figure(figsize=(8,5))
    plt.plot(two_theta_raw, intensity, label="Before (raw)")
    plt.plot(two_theta_after, intensity, label=f"After (δ = {delta:+.3f}°)")
    plt.xlabel("2θ (deg)"); plt.ylabel("Intensity (a.u.)")
    plt.title("XRD Zero-Point Correction: Before vs After")
    plt.legend()
    if xlim: plt.xlim(xlim[0], xlim[1])
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)

def save_single(two_theta, intensity, title, out_png, xlim=None, dpi=150):
    plt.figure(figsize=(8,5))
    plt.plot(two_theta, intensity)
    plt.xlabel("2θ (deg)"); plt.ylabel("Intensity (a.u.)")
    plt.title(title)
    if xlim: plt.xlim(xlim[0], xlim[1])
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)

def save_compare_diff_ratio(two_theta_raw, intensity, two_theta_after, delta, base, xlim=None, dpi=150, modes=("diff","ratio"), save_csv=False):
    g = common_grid(two_theta_raw, two_theta_after)
    I_before = interp_to_grid(two_theta_raw, intensity, g)
    I_after  = interp_to_grid(two_theta_after, intensity, g)
    m = np.isfinite(I_before) & np.isfinite(I_after)
    g = g[m]; I_before = I_before[m]; I_after = I_after[m]
    if "diff" in modes:
        plt.figure(figsize=(8,5))
        plt.plot(g, I_after - I_before)
        plt.xlabel("2θ (deg)"); plt.ylabel("ΔIntensity (after - before)")
        plt.title(f"Difference on Common Grid (δ = {delta:+.3f}°)")
        if xlim: plt.xlim(xlim[0], xlim[1])
        plt.tight_layout()
        plt.savefig(f"{base}_compare_diff.png", dpi=dpi)
    if "ratio" in modes:
        eps = max(1e-12, float(np.nanmax(I_before)) * 1e-12)
        denom = np.where(np.abs(I_before) < eps, np.nan, I_before)
        ratio = I_after / denom
        plt.figure(figsize=(8,5))
        plt.plot(g, ratio)
        plt.xlabel("2θ (deg)"); plt.ylabel("Intensity Ratio (after / before)")
        plt.title(f"Ratio on Common Grid (δ = {delta:+.3f}°)")
        if xlim: plt.xlim(xlim[0], xlim[1])
        plt.tight_layout()
        plt.savefig(f"{base}_compare_ratio.png", dpi=dpi)
    if save_csv:
        dfc = pd.DataFrame({
            "two_theta_deg": g,
            "intensity_before": I_before,
            "intensity_after": I_after,
            "diff_after_minus_before": I_after - I_before,
            "ratio_after_over_before": I_after / np.where(np.abs(I_before) < 1e-300, np.nan, I_before),
            "delta_deg_used": np.full_like(g, delta, dtype=float),
        })
        dfc.to_csv(f"{base}_compare.csv", index=False)

def save_csv_before_after(two_theta_raw, intensity, two_theta_after, delta, base):
    df_b = pd.DataFrame({"two_theta_deg": two_theta_raw, "intensity": intensity})
    df_a = pd.DataFrame({"two_theta_deg": two_theta_after, "intensity": intensity, "delta_deg_used": delta})
    df_b.to_csv(f"{base}_before.csv", index=False)
    df_a.to_csv(f"{base}_after.csv", index=False)

def decide_delta(p_input: Path, cli_delta, delta_report, delta_default):
    # 1) --delta
    if cli_delta is not None:
        return float(cli_delta), "from --delta"
    # 2) --delta-report
    if delta_report:
        rp = Path(delta_report)
        if not rp.exists():
            raise FileNotFoundError(f"指定的报告不存在：{rp}")
        return parse_delta_from_report(rp), f"from report: {rp.name}"
    # 3) CSV meta
    if p_input.suffix.lower() == ".csv":
        dv = try_delta_from_csv_meta(p_input)
        if dv is not None:
            return dv, f"from CSV column in {p_input.name}"
    # 4) sibling reports
    dv = try_delta_from_sibling_reports(p_input.with_suffix(""))
    if dv is not None:
        return dv, "from sibling report"
    # 5) default
    return float(delta_default), f"default {delta_default}"

def main():
    ap = argparse.ArgumentParser(description="XRD 0点补正对比图：可输出叠加、分别、差分/比值图，并导出 CSV。")
    ap.add_argument("-i","--input", required=True, help="输入文件：.csv/.txt/.dat/.xy/.ras")
    ap.add_argument("-d","--delta", type=float, default=None, help="零点 δ（度）")
    ap.add_argument("--delta-report", help="从 urlap 风格报告中解析 δ（优先行：'2THETA ORIGIN = ...'）")
    ap.add_argument("--delta-default", type=float, default=0.0, help="当无法解析时使用的默认 δ（度），默认 0.0）")
    ap.add_argument("-o","--out", default=None, help="叠加图输出 PNG（默认：<stem>_before_after.png）")
    ap.add_argument("--theta-col", help="2θ列（CSV 时可用：列名或序号）")
    ap.add_argument("--intensity-col", help="强度列（CSV 时可用：列名或序号）")
    ap.add_argument("--xlim", nargs=2, type=float, metavar=("XMIN","XMAX"), help="x 轴范围（度）")
    ap.add_argument("--dpi", type=int, default=150, help="DPI（默认 150）")
    ap.add_argument("--separate", action="store_true", help="另外保存分别的前/后两张图")
    ap.add_argument("--separate-only", action="store_true", help="只保存分别的前/后两张图，不保存叠加图")
    ap.add_argument("--compare", choices=["overlay","diff","ratio","both"], default="overlay",
                    help="额外对比图类型（默认 overlay 叠加；both=同时输出 diff 与 ratio）")
    ap.add_argument("--save-csv", action="store_true", help="同时导出 CSV（before/after/compare）")
    args = ap.parse_args()

    p = Path(args.input)
    two_theta_raw, intensity = load_xy(p, args.theta_col, args.intensity_col)

    delta, source = decide_delta(p, args.delta, args.delta_report, args.delta_default)
    print(f"[delta] {delta:+.6f} deg ({source})")

    two_theta_after = two_theta_raw - delta

    overlay_png = args.out or f"{p.stem}_before_after.png"
    base_for_others = Path(args.out).stem if (args.out and args.out.lower().endswith(".png")) else p.stem

    if args.compare in ("overlay", "both") and not args.separate_only:
        save_overlay(two_theta_raw, intensity, two_theta_after, delta, overlay_png, args.xlim, args.dpi)
        print(f"Saved overlay figure -> {overlay_png}")

    if args.separate or args.separate_only:
        before_png = f"{base_for_others}_before.png"
        after_png  = f"{base_for_others}_after.png"
        save_single(two_theta_raw, intensity, "XRD: Before (raw)", before_png, args.xlim, args.dpi)
        save_single(two_theta_after, intensity, f"XRD: After (δ = {delta:+.3f}°)", after_png, args.xlim, args.dpi)
        print(f"Saved separate figures -> {before_png}, {after_png}")

    if args.save_csv:
        save_csv_before_after(two_theta_raw, intensity, two_theta_after, delta, base_for_others)
        print(f"Saved CSV -> {base_for_others}_before.csv, {base_for_others}_after.csv")

    if args.compare in ("diff","ratio","both"):
        modes = ("diff","ratio") if args.compare == "both" else (args.compare,)
        save_compare_diff_ratio(two_theta_raw, intensity, two_theta_after, delta, base_for_others, args.xlim, args.dpi, modes, save_csv=args.save_csv)
        if "diff" in modes:
            print(f"Saved difference figure -> {base_for_others}_compare_diff.png")
        if "ratio" in modes:
            print(f"Saved ratio figure -> {base_for_others}_compare_ratio.png")
        if args.save_csv:
            print(f"Saved CSV -> {base_for_others}_compare.csv")

if __name__ == "__main__":
    main()
