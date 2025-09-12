
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xrd_to_dataset.py  (default: numeric-only patterns CSV)
------------------------------------------------------
Changes vs earlier version:
- ALWAYS write `calcd_patterns.csv` as a pure numeric matrix:
    * no 'file' column
    * no header row (angle labels)
    * no index column
- `targets.csv` keeps the mapping: file -> label (row index)

Supported inputs:
- ASCII two-column (.txt/.csv/.xy/.xye/.dat), auto delimiter & header skip
- XRDML (.xrdml)
- Vendor formats via xylib (.ras/.uxd/.raw/.cpi/...) if `pip install xylib-py`

Usage
    python xrd_to_dataset.py ./mydata --xmin 10 --xmax 80 --step 0.02

Deps
    pip install numpy pandas
    # optional for vendor formats:
    pip install xylib-py
"""
from __future__ import annotations

import os, re, sys, glob, math, argparse
from typing import List, Tuple, Optional
import numpy as np
import pandas as pd

XYLIB_OK=False
try:
    import xylib  # type: ignore
    XYLIB_OK=True
except Exception:
    XYLIB_OK=False

COMMENT_PREFIXES=("#",";","!","*","//")

def is_number(s:str)->bool:
    try: float(s); return True
    except Exception: return False

def read_ascii_two_columns(path: str):
    xs=[]; ys=[]
    with open(path,"r",encoding="utf-8",errors="ignore") as f:
        for line in f:
            s=line.strip()
            if not s: continue
            if s.startswith(COMMENT_PREFIXES): continue
            parts=re.split(r"[,\s;]+", s)
            if len(parts)<2: continue
            if not (is_number(parts[0]) and is_number(parts[1])): continue
            xs.append(float(parts[0])); ys.append(float(parts[1]))
    if not xs: raise ValueError("no numeric two-column rows")
    x=np.asarray(xs,float); y=np.asarray(ys,float)
    o=np.argsort(x); x=x[o]; y=y[o]
    uniq_x, idx=np.unique(x, return_index=True)
    return x[idx], y[idx]

def read_with_xylib(path: str):
    if not XYLIB_OK: raise ValueError("xylib not installed")
    lib=xylib.load_file(path)
    if lib is None or lib.get_block_count()<1: raise ValueError("xylib: no block")
    b=lib.get_block(0); n=b.get_point_count()
    xs=np.zeros(n,float); ys=np.zeros(n,float)
    for i in range(n):
        p=b.get_point(i); xs[i]=p.x; ys[i]=p.y
    return xs, ys

def read_xrdml(path: str):
    import xml.etree.ElementTree as ET
    tree=ET.parse(path); root=tree.getroot()
    dp=None
    for e in root.iter():
        if e.tag.endswith("dataPoints"): dp=e; break
    if dp is None: raise ValueError("XRDML: dataPoints missing")
    intens=None; xpos=None
    for c in dp.iter():
        if c.tag.endswith("intensities"):
            toks=((c.text or "").replace("\n"," ").split())
            intens=np.array([float(t) for t in toks], float)
        if c.tag.endswith("positions"):
            toks=((c.text or "").replace("\n"," ").split())
            if toks: xpos=np.array([float(t) for t in toks], float)
    if intens is None: raise ValueError("XRDML: intensities missing")
    if xpos is None:
        start=end=None
        for c in dp.iter():
            if c.tag.endswith("startPosition"): start=float(c.text)
            if c.tag.endswith("endPosition"): end=float(c.text)
        if start is not None and end is not None:
            n=len(intens); step=(end-start)/(n-1) if n>1 else 0.0
            xpos=start+np.arange(n)*step
    if xpos is None: raise ValueError("XRDML: positions missing")
    return np.asarray(xpos,float), np.asarray(intens,float)

def build_grid(xmin: float, xmax: float, step: float)->np.ndarray:
    if xmax<=xmin: raise ValueError("--xmax must be > --xmin")
    if step<=0: raise ValueError("--step must be > 0")
    n=int(math.floor((xmax-xmin)/step))+1
    return xmin+np.arange(n)*step

def resample_to_grid(x,y,grid): return np.interp(grid,x,y,left=0.0,right=0.0)

def l2_normalize(v, eps: float=1e-12):
    n=float(np.linalg.norm(v))
    if not np.isfinite(n) or n<=eps: return None
    return v/n

ASCII_EXTS=(".txt",".csv",".xy",".xye",".dat")
XYLIB_EXTS=(".ras",".uxd",".raw",".rd",".cpi",".udf",".xdd")
XRDML_EXTS=(".xrdml",)

def discover_files(input_dir: str)->List[str]:
    pats=[]
    for ext in ASCII_EXTS+XYLIB_EXTS+XRDML_EXTS:
        pats += [f"*{ext}", f"*{ext.upper()}"]
    files=[]
    for p in pats: files += glob.glob(os.path.join(input_dir,p))
    if not files:
        files=[os.path.join(input_dir,f) for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir,f))]
    return sorted(set(files))

def read_any_xrd(path: str):
    ext=os.path.splitext(path)[1].lower()
    try: return read_ascii_two_columns(path)
    except Exception: pass
    if ext in XRDML_EXTS:
        try: return read_xrdml(path)
        except Exception: pass
    if XYLIB_OK:
        try: return read_with_xylib(path)
        except Exception: pass
    if XYLIB_OK and ext not in XYLIB_EXTS:
        try: return read_with_xylib(path)
        except Exception: pass
    raise ValueError("unrecognized or unreadable XRD format")

def main(argv: Optional[List[str]] = None)->int:
    ap=argparse.ArgumentParser(description="XRD -> dataset (common grid + L2 normalize)")
    ap.add_argument("input_dir")
    ap.add_argument("--xmin", type=float, default=10.0)
    ap.add_argument("--xmax", type=float, default=80.0)
    ap.add_argument("--step", type=float, default=0.02)
    args=ap.parse_args(argv)

    if not os.path.isdir(args.input_dir):
        sys.stderr.write(f"[ERR] Not a directory: {args.input_dir}\n"); return 1
    try: grid=build_grid(args.xmin, args.xmax, args.step)
    except Exception as e: sys.stderr.write(f"[ERR] Grid error: {e}\n"); return 1

    files=discover_files(args.input_dir)
    if not files:
        sys.stderr.write("[ERR] No candidate XRD files.\n"); return 1

    rows=[]; labels=[]; skipped=[]
    for fp in files:
        base=os.path.basename(fp)
        try:
            x,y=read_any_xrd(fp)
            yi=resample_to_grid(x,y,grid)
            yi=l2_normalize(yi)
            if yi is None: raise ValueError("norm==0")
        except Exception as e:
            skipped.append((base, str(e))); continue
        rows.append(yi); labels.append(base)

    if not rows:
        sys.stderr.write("[ERR] No valid patterns.\n"); return 1
    X=np.vstack(rows)

    # patterns CSV: numeric-only (no header, no index, no file column)
    patterns_csv=os.path.join(os.getcwd(),"calcd_patterns.csv")
    np.savetxt(patterns_csv, X, delimiter=",")

    # targets.csv keeps mapping: file -> label (row index)
    tgt_df=pd.DataFrame({"file": labels, "label": list(range(len(labels))) })
    targets_csv=os.path.join(os.getcwd(), "targets.csv")
    tgt_df.to_csv(targets_csv, index=False)

    print(f"[OK] files={len(labels)} skipped={len(skipped)} gridpts={len(grid)}")
    if skipped:
        for n,r in skipped[:10]: print(f"   [skip] {n}: {r}")
        if len(skipped)>10: print(f"   ... and {len(skipped)-10} more")
    print(f"[OK] calcd_patterns -> {patterns_csv}  (numeric-only)")
    print(f"[OK] targets -> {targets_csv}")
    if not XYLIB_OK:
        print("[NOTE] xylib-py not installed; vendor formats (.ras/.uxd/.raw/...) may be skipped.")
    return 0

if __name__=="__main__":
    sys.exit(main())
