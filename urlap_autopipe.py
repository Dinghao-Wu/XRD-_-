#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
urlap_autopipe.py  —  一键管线：用 urlap_report_v3.py 找到 δ，然后出图 & 导出 CSV

需要同目录存在：
  - urlap_report_v3.py
  - plot_xrd_before_after_split.py

基本用法：
  # 扫描 δ（例如 -0.5..0.5, step=0.1），以不确定度乘积最小（uncprod）为准；
  # 然后根据最佳 δ 对谱线做补正并绘图；导出 CSV；x 轴范围 5–90°；输出前/后/差分/比值图。
  python urlap_autopipe.py \
    --peaks LATTICE9.TXT \
    --spectrum before.csv \
    --mode scan --delta-range -0.5 0.5 --step 0.1 --select-by uncprod \
    --xlim 5 90 --compare both --separate --save-csv \
    --out-prefix sample1

  # 固定 δ（若输入 RAW 的第 3 行已写 δ 也行）
  python urlap_autopipe.py \
    --peaks TETRA_INPUT.TXT \
    --spectrum before.csv \
    --mode fixed \
    --xlim 5 90 --separate --compare both --save-csv
"""
import argparse, sys, shlex, subprocess
from pathlib import Path

def require_exists(p: Path, what: str):
    if not p.exists():
        print(f"[ERROR] 找不到 {what}: {p}", file=sys.stderr)
        sys.exit(2)

def main():
    ap = argparse.ArgumentParser(description="管线：urlap 找 δ + 出图/导出 CSV")
    ap.add_argument("--peaks", required=True, help="urlap RAW 峰文件（含晶系、λ、δ 或扫描参数）")
    ap.add_argument("--spectrum", required=True, help="谱线文件（.csv/.txt/.dat/.xy/.ras）")
    ap.add_argument("--mode", choices=["fixed","scan"], default="scan", help="urlap 模式（默认 scan）")
    ap.add_argument("--wavelength", type=float, help="覆盖 λ（可选）")
    ap.add_argument("--delta", type=float, help="固定模式下覆盖 δ（可选）")
    ap.add_argument("--delta-range", nargs=2, type=float, metavar=("MIN","MAX"), default=[-0.5, 0.5])
    ap.add_argument("--step", type=float, default=0.1, help="扫描步长（默认 0.1°）")
    ap.add_argument("--select-by", choices=["rms","uncprod"], default="uncprod")
    ap.add_argument("--out-prefix", default=None, help="输出基名（不带后缀）；默认取谱线文件名 stem")
    # 绘图与导出参数（原样传给 plot 脚本）
    ap.add_argument("--xlim", nargs=2, type=float, metavar=("XMIN","XMAX"))
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--compare", choices=["overlay","diff","ratio","both"], default="overlay")
    ap.add_argument("--separate", action="store_true")
    ap.add_argument("--separate-only", action="store_true")
    ap.add_argument("--save-csv", action="store_true")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    urlap = here / "urlap_report_v3.py"
    plot  = here / "plot_xrd_before_after_split.py"
    require_exists(urlap, "urlap_report_v3.py")
    require_exists(plot,  "plot_xrd_before_after_split.py")

    peaks = Path(args.peaks).resolve()
    spec  = Path(args.spectrum).resolve()
    require_exists(peaks, "峰文件")
    require_exists(spec,  "谱线文件")

    # 1) 运行 urlap_report_v3.py 生成报表，并拿到 δ（plot 脚本会直接从报表解析 δ）
    out_prefix = args.out_prefix or spec.stem
    report_path = peaks.with_suffix(peaks.suffix + ".out.txt")  # 默认命名逻辑
    # 允许按基名定制报表名：<out_prefix>_urlap.out.txt
    report_path = spec.with_name(f"{out_prefix}_urlap.out.txt")

    cmd_parts = [sys.executable, str(urlap), "-i", str(peaks), "--mode", args.mode, "--out", str(report_path)]
    if args.wavelength is not None:
        cmd_parts += ["--wavelength", str(args.wavelength)]
    if args.mode == "fixed":
        if args.delta is not None:
            cmd_parts += ["--delta", str(args.delta)]
    else:
        cmd_parts += ["--delta-range", str(args.delta_range[0]), str(args.delta_range[1]), "--step", str(args.step), "--select-by", args.select_by]

    print("[RUN] urlap:", shlex.join(cmd_parts))
    r1 = subprocess.run(cmd_parts, capture_output=True, text=True)
    if r1.returncode != 0:
        print(r1.stdout)
        print(r1.stderr, file=sys.stderr)
        sys.exit(r1.returncode)
    # urlap_report_v3.py 会打印输出路径；我们直接使用 report_path 即可
    print(f"[OK] 报表 -> {report_path}")

    # 2) 调用绘图脚本，delta 从报表自动解析（--delta-report）
    out_png = spec.with_name(f"{out_prefix}_before_after.png")
    cmd2 = [sys.executable, str(plot), "-i", str(spec), "--delta-report", str(report_path),
            "-o", str(out_png), "--compare", args.compare, "--dpi", str(args.dpi)]
    if args.separate: cmd2.append("--separate")
    if args.separate_only: cmd2.append("--separate-only")
    if args.save_csv: cmd2.append("--save-csv")
    if args.xlim: cmd2 += ["--xlim", str(args.xlim[0]), str(args.xlim[1])]

    print("[RUN] plot :", shlex.join(cmd2))
    r2 = subprocess.run(cmd2, capture_output=True, text=True)
    print(r2.stdout)
    if r2.returncode != 0:
        print(r2.stderr, file=sys.stderr)
        sys.exit(r2.returncode)

    print("\n[PIPE DONE] 主要输出：")
    base = out_png if out_png.name.endswith(".png") else spec.with_name(f"{out_prefix}_before_after.png")
    print(" -", base)
    if args.separate or args.separate_only:
        print(" -", spec.with_name(f"{out_prefix}_before.png"))
        print(" -", spec.with_name(f"{out_prefix}_after.png"))
    if args.compare in ("diff","ratio","both"):
        if args.compare in ("diff","both"):
            print(" -", spec.with_name(f"{out_prefix}_compare_diff.png"))
        if args.compare in ("ratio","both"):
            print(" -", spec.with_name(f"{out_prefix}_compare_ratio.png"))
    if args.save_csv:
        print(" -", spec.with_name(f"{out_prefix}_before.csv"))
        print(" -", spec.with_name(f"{out_prefix}_after.csv"))
        print(" -", spec.with_name(f"{out_prefix}_compare.csv"))
    print(" -", report_path)

if __name__ == "__main__":
    main()
