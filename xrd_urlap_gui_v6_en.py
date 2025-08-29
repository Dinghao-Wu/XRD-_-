# -*- coding: utf-8 -*-
"""
XRD URLAP GUI v6-en
- Full-system URLAP (ported from your original script)
- CIF/Materials Project -> auto crystal system detection -> reference peaks
- PDXL-style Peak List table: No / 2θ_obs / d_obs / Intensity% / h k l / 2θ_calc / |Δ2θ| / Phase / Use
- Nearest-neighbour matching: match CIF (hkl, 2θ_calc) to Peak List 2θ_obs (with tolerance, unique assignment)
- δ: Scan first, then Fixed fit to generate URLAP report
- Color selection; double-click table row toggles "Use"

Run: python xrd_urlap_gui_v6_en.py
Deps: numpy, matplotlib, tkinter; (optional) pymatgen (for CIF/MP)
"""

import os, re, math, sys, csv
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict

import numpy as np

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# Optional: pymatgen
try:
    from pymatgen.core import Structure
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.ext.matproj import MPRester
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    PYMATGEN_OK = True
except Exception:
    PYMATGEN_OK = False

# -----------------------------
# Utilities
# -----------------------------

def _try_float(t: str) -> Optional[float]:
    if t is None:
        return None
    s = str(t).strip()
    if not s:
        return None
    if ("," in s) and ("." not in s):
        s = s.replace(",", ".")
    try:
        return float(s)
    except Exception:
        return None

def read_xy_file(path: str) -> Tuple[np.ndarray, np.ndarray]:
    xs: List[float] = []
    ys: List[float] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            parts = re.split(r"[\s,;\t]+", s)
            if len(parts) < 2:
                continue
            x = _try_float(parts[0]); y = _try_float(parts[1])
            if (x is None) or (y is None):
                continue
            xs.append(x); ys.append(y)
    if not xs:
        raise ValueError("No 2θ-Intensity columns detected in the file.")
    x = np.asarray(xs, float)
    y = np.asarray(ys, float)
    if np.max(y) > 0:
        y = y / np.max(y) * 100.0
    return x, y

def median_step(x: np.ndarray) -> float:
    if x.size < 2:
        return 0.02
    dx = np.diff(x)
    dx = dx[dx > 0]
    if dx.size == 0:
        return 0.02
    return float(np.median(dx))

def moving_average(y: np.ndarray, win_pts: int) -> np.ndarray:
    win_pts = max(1, int(win_pts))
    if win_pts == 1:
        return y.copy()
    kernel = np.ones(win_pts, dtype=float) / win_pts
    yy = np.convolve(y, kernel, mode="same")
    return yy

def parabolic_refine(x: np.ndarray, y: np.ndarray, i: int) -> float:
    if i <= 0 or i >= len(x)-1:
        return float(x[i])
    x1, y1 = float(x[i-1]), float(y[i-1])
    x2, y2 = float(x[i]),   float(y[i])
    x3, y3 = float(x[i+1]), float(y[i+1])
    denom = (x1-x2)*(x1-x3)*(x2-x3)
    if abs(denom) < 1e-18:
        return float(x2)
    a = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / denom
    b = (x3*x3*(y1-y2) + x2*x2*(y3-y1) + x1*x1*(y2-y3)) / denom
    if abs(a) < 1e-18:
        return float(x2)
    xv = -b/(2*a)
    if min(x1,x3) <= xv <= max(x1,x3):
        return float(xv)
    return float(x2)

def find_peaks_simple(x: np.ndarray, y: np.ndarray, win_deg: float, min_height_pct: float, min_sep_deg: float, refine: bool) -> List[Tuple[float, float]]:
    """Return [(two_obs, intensity_pct), ...] on smoothed signal, scale intensity to 0-100%."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    y = y - np.min(y)  # baseline shift
    dx = median_step(x)
    win_pts = max(1, int(max(win_deg,0.0)/max(dx,1e-6) + 0.5))
    if win_pts % 2 == 0:
        win_pts += 1
    ys = moving_average(y, win_pts)

    # local maxima
    cond = (ys[1:-1] > ys[:-2]) & (ys[1:-1] >= ys[2:])
    cand_idx = np.where(cond)[0] + 1
    if cand_idx.size == 0:
        return []

    # intensity threshold
    ymax = float(np.max(ys))
    if ymax <= 0:
        return []
    thr = float(min_height_pct)/100.0 * ymax
    cand_idx = cand_idx[ys[cand_idx] >= thr]
    if cand_idx.size == 0:
        return []

    # non-maximum suppression by spacing
    order = np.argsort(ys[cand_idx])[::-1]
    chosen: List[int] = []
    sep = max(min_sep_deg, 0.0)
    for j in order:
        i = int(cand_idx[j])
        if len(chosen) == 0:
            chosen.append(i); continue
        if all(abs(x[i]-x[k]) >= sep for k in chosen):
            chosen.append(i)
    chosen.sort()

    peaks: List[Tuple[float, float]] = []
    for i in chosen:
        xi = parabolic_refine(x, ys, i) if refine else float(x[i])
        yi = float(ys[i]) / ymax * 100.0
        peaks.append((xi, yi))
    return peaks

def d_from_two(two_deg: float, wavelength: float, delta_deg: float) -> float:
    theta = math.radians((two_deg + delta_deg)/2.0)
    s = math.sin(theta)
    if s <= 0:
        return float("inf")
    return wavelength / (2.0*s)

# -----------------------------
# URLAP (multi-crystal systems; from your original file)
# -----------------------------

class URLAP:
    @staticmethod
    def parse_raw(path: str):
        lines = [ln.rstrip("\n") for ln in open(path, "r", encoding="utf-8", errors="ignore")]
        if not lines:
            raise ValueError("RAW 文件为空")
        title = lines[0].strip()
        system_code = None
        if len(lines) >= 2:
            parts = lines[1].split()
            if parts:
                try:
                    system_code = int(float(parts[0]))
                except Exception:
                    pass
        wav = None
        delta = None
        if len(lines) >= 3:
            parts = lines[2].split()
            if len(parts) >= 1:
                try:
                    wav = float(parts[0])
                except Exception:
                    pass
            if len(parts) >= 2:
                try:
                    delta = float(parts[1])
                except Exception:
                    pass
        peaks = []
        for ln in lines[3:]:
            if not ln.strip():
                continue
            parts = ln.split()
            if len(parts) < 4:
                continue
            try:
                h = int(float(parts[0])); k = int(float(parts[1])); l = int(float(parts[2]))
                if h >= 1000:
                    break
                two = float(parts[3])
                peaks.append([h, k, l, two])
            except Exception:
                continue
        if not peaks:
            raise ValueError("RAW 中没有解析到 (h k l 2theta) 峰列表")
        return title, system_code, wav, delta, np.asarray(peaks, float)

    @staticmethod
    def obs_to_Q(two_deg, wavelength):
        theta = np.radians(0.5 * two_deg)
        return (2.0 * np.sin(theta) / wavelength) ** 2  # equals 1/d^2

    @staticmethod
    def two_from_dinv2(dinv2, wavelength):
        d = 1.0 / np.sqrt(np.maximum(dinv2, 1e-300))
        arg = np.clip(wavelength / (2.0 * d), -1.0, 1.0)
        th = np.degrees(np.arcsin(arg))
        return 2.0 * th

    @staticmethod
    def solve_ls(X, y):
        XtX = X.T @ X
        XtY = X.T @ y
        try:
            m = np.linalg.solve(XtX, XtY)
            pinv = np.linalg.inv(XtX)
        except np.linalg.LinAlgError:
            m = np.linalg.lstsq(X, y, rcond=None)[0]
            pinv = np.linalg.pinv(XtX)
        resid = y - X @ m
        p = X.shape[1]
        dof = max(len(y) - p, 1)
        RSS = float(np.sum(resid ** 2))
        sigma2 = RSS / dof
        cov_m = sigma2 * pinv
        return m, cov_m, RSS, dof

    # ---- system models ----
    @staticmethod
    def design_triclinic(hkl):
        h = hkl[:, 0]; k = hkl[:, 1]; l = hkl[:, 2]
        X = np.column_stack([h*h, k*k, l*l, 2*h*k, 2*h*l, 2*k*l])
        return X

    @staticmethod
    def extract_triclinic(m):
        G11, G22, G33, G12, G13, G23 = m.tolist()
        Gstar = np.array([[G11, G12, G13], [G12, G22, G23], [G13, G23, G33]], float)
        g = np.linalg.inv(Gstar)
        a = math.sqrt(g[0, 0]); b = math.sqrt(g[1, 1]); c = math.sqrt(g[2, 2])
        cos_alpha = g[1, 2] / (b * c); cos_beta = g[0, 2] / (a * c); cos_gamma = g[0, 1] / (a * b)
        cos_alpha = float(np.clip(cos_alpha, -1.0, 1.0))
        cos_beta  = float(np.clip(cos_beta , -1.0, 1.0))
        cos_gamma = float(np.clip(cos_gamma, -1.0, 1.0))
        alpha = math.degrees(math.acos(cos_alpha))
        beta  = math.degrees(math.acos(cos_beta))
        gamma = math.degrees(math.acos(cos_gamma))
        return a, b, c, alpha, beta, gamma, Gstar

    @staticmethod
    def dinv2_triclinic(hkl, Gstar):
        H = hkl.astype(float)
        vals = np.einsum('...i,ij,...j->...', H, Gstar, H)
        return vals

    @staticmethod
    def design_monoclinic_b(hkl):
        h = hkl[:, 0].astype(float); k = hkl[:, 1].astype(float); l = hkl[:, 2].astype(float)
        X = np.column_stack([h*h, k*k, l*l, h*l])
        return X

    @staticmethod
    def extract_monoclinic_b(m):
        m1, m2, m3, m4 = m.tolist()
        cos2 = (m4*m4) / (4.0 * max(m1, 1e-300) * max(m3, 1e-300))
        cos2 = float(np.clip(cos2, 0.0, 1.0))
        sin2 = max(1.0 - cos2, 1e-12)
        s = math.sqrt(sin2)
        a = 1.0 / math.sqrt(max(m1, 1e-300) * sin2)
        b = 1.0 / math.sqrt(max(m2, 1e-300) * sin2)
        c = 1.0 / math.sqrt(max(m3, 1e-300) * sin2)
        cosb = -0.5 * m4 * a * c * sin2
        cosb = float(np.clip(cosb, -1.0, 1.0))
        beta = math.degrees(math.acos(cosb))
        return a, b, c, 90.0, beta, 90.0

    @staticmethod
    def dinv2_monoclinic_b(hkl, a, b, c, beta_deg):
        beta = math.radians(beta_deg)
        s2 = (math.sin(beta))**2
        cb = math.cos(beta)
        h = hkl[:,0].astype(float); k = hkl[:,1].astype(float); l = hkl[:,2].astype(float)
        return (h*h)/(a*a*s2) + (k*k)/(b*b*s2) + (l*l)/(c*c*s2) - (2*h*l*cb)/(a*c*s2)

    @staticmethod
    def design_monoclinic_c(hkl):
        h = hkl[:, 0].astype(float); k = hkl[:, 1].astype(float); l = hkl[:, 2].astype(float)
        X = np.column_stack([h*h, k*k, l*l, h*k])
        return X

    @staticmethod
    def extract_monoclinic_c(m):
        m1, m2, m3, m4 = m.tolist()
        cos2 = (m4*m4) / (4.0 * max(m1, 1e-300) * max(m2, 1e-300))
        cos2 = float(np.clip(cos2, 0.0, 1.0))
        sin2 = max(1.0 - cos2, 1e-12)
        a = 1.0 / math.sqrt(max(m1, 1e-300) * sin2)
        b = 1.0 / math.sqrt(max(m2, 1e-300) * sin2)
        c = 1.0 / math.sqrt(max(m3, 1e-300) * sin2)
        cosg = -0.5 * m4 * a * b * sin2
        cosg = float(np.clip(cosg, -1.0, 1.0))
        gamma = math.degrees(math.acos(cosg))
        return a, b, c, 90.0, 90.0, gamma

    @staticmethod
    def dinv2_monoclinic_c(hkl, a, b, c, gamma_deg):
        g = math.radians(gamma_deg)
        s2 = (math.sin(g)) ** 2
        cg = math.cos(g)
        h = hkl[:, 0].astype(float); k = hkl[:, 1].astype(float); l = hkl[:, 2].astype(float)
        return (h*h)/(a*a*s2) + (k*k)/(b*b*s2) + (l*l)/(c*c*s2) - (2*h*k*cg)/(a*b*s2)

    @staticmethod
    def design_orthorhombic(hkl):
        h = hkl[:, 0]; k = hkl[:, 1]; l = hkl[:, 2]
        X = np.column_stack([h*h, k*k, l*l])
        return X

    @staticmethod
    def extract_orthorhombic(m):
        ma, mb, mc = m.tolist()
        a = 1.0 / math.sqrt(max(ma, 1e-300))
        b = 1.0 / math.sqrt(max(mb, 1e-300))
        c = 1.0 / math.sqrt(max(mc, 1e-300))
        return a, b, c, 90.0, 90.0, 90.0

    @staticmethod
    def dinv2_orthorhombic(hkl, a, b, c):
        h = hkl[:, 0].astype(float); k = hkl[:, 1].astype(float); l = hkl[:, 2].astype(float)
        return (h*h)/(a*a) + (k*k)/(b*b) + (l*l)/(c*c)

    @staticmethod
    def design_tetragonal(hkl):
        h = hkl[:, 0]; k = hkl[:, 1]; l = hkl[:, 2]
        Sxy = (h*h + k*k).astype(float); Sz = (l*l).astype(float)
        X = np.column_stack([Sxy, Sz])
        return X

    @staticmethod
    def extract_tetragonal(m):
        ma, mc = m.tolist()
        a = 1.0 / math.sqrt(max(ma, 1e-300))
        c = 1.0 / math.sqrt(max(mc, 1e-300))
        return a, a, c, 90.0, 90.0, 90.0

    @staticmethod
    def dinv2_tetragonal(hkl, a, c):
        h = hkl[:, 0].astype(float); k = hkl[:, 1].astype(float); l = hkl[:, 2].astype(float)
        return (h*h + k*k)/(a*a) + (l*l)/(c*c)

    @staticmethod
    def design_hex(hkl):
        h = hkl[:, 0].astype(float); k = hkl[:, 1].astype(float); l = hkl[:, 2].astype(float)
        H = (4.0/3.0) * (h*h + h*k + k*k)
        L = l*l
        X = np.column_stack([H, L])
        return X

    @staticmethod
    def extract_hex(m):
        ma, mc = m.tolist()
        a = 1.0 / math.sqrt(max(ma, 1e-300))
        c = 1.0 / math.sqrt(max(mc, 1e-300))
        return a, a, c, 90.0, 90.0, 120.0

    @staticmethod
    def dinv2_hex(hkl, a, c):
        h = hkl[:, 0].astype(float); k = hkl[:, 1].astype(float); l = hkl[:, 2].astype(float)
        return (4.0/3.0) * (h*h + h*k + k*k)/(a*a) + (l*l)/(c*c)

    @staticmethod
    def design_cubic(hkl):
        h = hkl[:, 0]; k = hkl[:, 1]; l = hkl[:, 2]
        S = (h*h + k*k + l*l).astype(float)
        X = S.reshape(-1, 1)
        return X

    @staticmethod
    def extract_cubic(m):
        ma = float(m[0])
        a = 1.0 / math.sqrt(max(ma, 1e-300))
        return a, a, a, 90.0, 90.0, 90.0

    @staticmethod
    def dinv2_cubic(hkl, a):
        h = hkl[:, 0].astype(float); k = hkl[:, 1].astype(float); l = hkl[:, 2].astype(float)
        return (h*h + k*k + l*l) / (a*a)

    @staticmethod
    def system_design_and_extract(code):
        if code == 1:
            return URLAP.design_triclinic, URLAP.extract_triclinic, URLAP.dinv2_triclinic, "TRICLINIC"
        if code == 2:
            return URLAP.design_monoclinic_b, URLAP.extract_monoclinic_b, URLAP.dinv2_monoclinic_b, "MONOCLINIC(b)"
        if code == 3:
            return URLAP.design_monoclinic_c, URLAP.extract_monoclinic_c, URLAP.dinv2_monoclinic_c, "MONOCLINIC(c)"
        if code == 4:
            return URLAP.design_orthorhombic, URLAP.extract_orthorhombic, URLAP.dinv2_orthorhombic, "ORTHORHOMBIC"
        if code == 5:
            return URLAP.design_tetragonal, URLAP.extract_tetragonal, URLAP.dinv2_tetragonal, "TETRAGONAL"
        if code == 6:
            return URLAP.design_cubic, URLAP.extract_cubic, URLAP.dinv2_cubic, "CUBIC"
        if code == 7:
            return URLAP.design_hex, URLAP.extract_hex, URLAP.dinv2_hex, "TRIGONAL"
        if code == 8:
            return URLAP.design_hex, URLAP.extract_hex, URLAP.dinv2_hex, "HEXAGONAL"
        raise NotImplementedError("System code must be 1..8")

    @staticmethod
    def numeric_jacobian(func, m, eps=1e-6):
        m = np.asarray(m, float)
        base = np.asarray(func(m), float)
        J = np.zeros((base.shape[0], m.size), float)
        for j in range(m.size):
            mj = m.copy()
            step = eps * max(1.0, abs(mj[j]))
            mj[j] += step
            f2 = np.asarray(func(mj), float)
            J[:, j] = (f2 - base) / step
        return J, base

    @staticmethod
    def compute_from_system(peaks, wavelength, delta, code):
        hkl = peaks[:, :3].astype(int)
        obs_corr = peaks[:, 3] - delta
        Q = URLAP.obs_to_Q(obs_corr, wavelength)
        design, extract, dinv2_fun, label = URLAP.system_design_and_extract(code)
        X = design(hkl)
        m, cov_m, RSS, dof = URLAP.solve_ls(X, Q)
        # derive cell
        a, b, c, alpha, beta, gamma, *extra = extract(m)

        # uncertainties for (a,b,c) by numerical Jacobian
        def map_to_abc(mm):
            a2, b2, c2, *_ = extract(mm)
            return np.array([a2, b2, c2], float)
        J, abc = URLAP.numeric_jacobian(map_to_abc, m, eps=1e-6)
        cov_abc = J @ cov_m @ J.T
        da, db, dc = [float(math.sqrt(max(cov_abc[i, i], 0.0))) for i in range(3)]

        # build table
        if code == 1:
            Gstar = extra[0]
            dinv2 = URLAP.dinv2_triclinic(hkl, Gstar)
        elif code == 2:
            dinv2 = URLAP.dinv2_monoclinic_b(hkl, a, b, c, beta)
        elif code == 3:
            dinv2 = URLAP.dinv2_monoclinic_c(hkl, a, b, c, gamma)
        elif code == 4:
            dinv2 = URLAP.dinv2_orthorhombic(hkl, a, b, c)
        elif code == 5:
            dinv2 = URLAP.dinv2_tetragonal(hkl, a, c)
        elif code == 6:
            dinv2 = URLAP.dinv2_cubic(hkl, a)
        else:  # 7 or 8
            dinv2 = URLAP.dinv2_hex(hkl, a, c)
        cal = URLAP.two_from_dinv2(dinv2, wavelength)
        resid_deg = obs_corr - cal
        d_hkl = 1.0 / np.sqrt(dinv2 + 1e-300)
        Sdisp = (hkl[:, 0]**2 + hkl[:, 1]**2 + hkl[:, 2]**2).astype(float)

        info = {
            "label": label, "a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma,
            "da": da, "db": db, "dc": dc, "cov_m": cov_m, "m": m,
            "table": np.column_stack([hkl, obs_corr, cal, resid_deg, d_hkl, Sdisp])
        }
        return info

    @staticmethod
    def unc_product(info, code):
        if code == 6:  # cubic
            return info["da"] ** 3
        if code in (5, 7, 8):  # tetragonal/trigonal/hex
            return (info["da"] ** 2) * info["dc"]
        return info["da"] * info["db"] * info["dc"]

    @staticmethod
    def rms_residual(info):
        resid = info["table"][:, 5]  # OBS - CAL
        return float(np.sqrt(np.mean(resid ** 2)))

    @staticmethod
    def render_report(title, wavelength, delta, code, info):
        label = info["label"]
        a, b, c = info["a"], info["b"], info["c"]
        da, db, dc = info["da"], info["db"], info["dc"]
        alpha, beta, gamma = info["alpha"], info["beta"], info["gamma"]
        table = info["table"]
        lines = []
        lines.append(f" TITLE : {title:<62}")
        lines.append("")
        lines.append(f" LATTICE SYSTEM     {label}")
        lines.append("")
        lines.append(f"   WAVE LENGTH =    {wavelength:0.5f}     2THETA ORIGIN =   {delta:0.5f}")
        lines.append("")
        lines.append("")
        lines.append("        H    K    L   2THETA_OBS   2THETA_CAL    CAL-OBS")
        lines.append("")
        for row in table:
            h, k, l, obs_show, cal, resid, d_hkl, S = row
            lines.append(f"{int(h):9d}{int(k):5d}{int(l):5d}{obs_show:12.7f}{cal:13.7f}{(obs_show-cal):12.7f}")
        lines.append("")
        lines.append("")
        lines.append(" DIRECT CELL CONSTANT ")
        lines.append("")
        lines.append("    A             DA           B            DB           C            DC       ")
        lines.append(f"{a:9.7f}    {da:0.7f}    {b:9.7f}    {db:0.7f}    {c:9.7f}    {dc:0.7f}")
        lines.append("")
        lines.append("  ALPHA        DALPHA       BETA         DBETA       GAMMA        DGAMMA       ")
        lines.append(f"{alpha:9.6f}     {0.000000:0.6f}  {beta:9.6f}     {0.000000:0.6f}  {gamma:9.6f}     {0.000000:0.6f}")
        return "\n".join(lines)

    @staticmethod
    def select_delta(peaks, wavelength, code, dmin, dmax, step, criterion):
        deltas = np.arange(dmin, dmax + 1e-15, step, dtype=float)
        best = None
        for d in deltas:
            info = URLAP.compute_from_system(peaks, wavelength, d, code)
            score = URLAP.rms_residual(info) if criterion == "rms" else URLAP.unc_product(info, code)
            if (best is None) or (score < best["score"]):
                best = {"delta": d, "info": info, "score": score}
        return best



# -----------------------------
# GUI
# -----------------------------

@dataclass
class PeakRow:
    two_obs: float
    inten: float
    d_obs: float
    h: Optional[int] = None
    k: Optional[int] = None
    l: Optional[int] = None
    two_calc: Optional[float] = None
    d2: Optional[float] = None
    phase: str = ""
    use: bool = False

class URLAPGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("XRD URLAP v6-en: Peak List + Matching + URLAP")
        self.geometry("1280x820")

        self._xy: Optional[Tuple[np.ndarray, np.ndarray]] = None
        self._calc_peaks: List[Tuple[float, Tuple[int,int,int]]] = []  # (2θ_calc, (h,k,l))
        self._auto_code: Optional[int] = None
        self._auto_name: str = ""
        self._struct = None

        # state
        self._peaklist: List[PeakRow] = []
        self._phase_name: str = "Phase"

        self._make_widgets()

    # ---------- UI ----------
    def _make_widgets(self):
        nb = ttk.Notebook(self); nb.pack(fill="both", expand=True)

        # Tab 1: Data + Ref Peaks + Plot
        tab_main = ttk.Frame(nb); nb.add(tab_main, text="Data / Ref Peaks / Plot")
        left = ttk.Frame(tab_main); left.pack(side=tk.LEFT, fill="y", padx=8, pady=8)
        right = ttk.Frame(tab_main); right.pack(side=tk.RIGHT, fill="both", expand=True, padx=8, pady=8)

        ttk.Label(left, text="XRD 2-column file (2θ, Intensity)").pack(anchor="w")
        self.e_xy = ttk.Entry(left, width=40); self.e_xy.pack(anchor="w")
        ttk.Button(left, text="Browse…", command=self._pick_xy).pack(anchor="w", pady=(2,6))

        ttk.Label(left, text="Wavelength (Å)").pack(anchor="w")
        self.e_wav = ttk.Entry(left, width=12); self.e_wav.insert(0, "1.5406"); self.e_wav.pack(anchor="w", pady=(0,6))

        ttk.Label(left, text="Sample / Title").pack(anchor="w")
        self.e_title = ttk.Entry(left, width=40); self.e_title.pack(anchor="w", pady=(0,6))

        ttk.Separator(left, orient="horizontal").pack(fill="x", pady=6)

        ttk.Label(left, text="Phase source").pack(anchor="w")
        self.phase_src = tk.StringVar(value="cif")
        ttk.Radiobutton(left, text="CIF file", variable=self.phase_src, value="cif").pack(anchor="w")
        ttk.Radiobutton(left, text="Materials Project", variable=self.phase_src, value="mp").pack(anchor="w")

        self.e_cif = ttk.Entry(left, width=40); self.e_cif.pack(anchor="w")
        ttk.Button(left, text="Browse CIF…", command=self._pick_cif).pack(anchor="w", pady=(2,6))

        ttk.Label(left, text="MP API Key").pack(anchor="w")
        self.e_mpkey = ttk.Entry(left, width=40); self.e_mpkey.pack(anchor="w")
        ttk.Label(left, text="mp-id or formula").pack(anchor="w")
        self.e_mpid = ttk.Entry(left, width=40); self.e_mpid.pack(anchor="w", pady=(0,6))

        ttk.Label(left, text="Phase name (for table)").pack(anchor="w")
        self.e_phase = ttk.Entry(left, width=40); self.e_phase.insert(0, "Phase"); self.e_phase.pack(anchor="w", pady=(0,6))

        ttk.Button(left, text="Compute reference peaks + auto-detect crystal system", command=self._calc_reference).pack(anchor="w", pady=(6,6))
        self.lbl_sys = ttk.Label(left, text="Crystal system (auto): --"); self.lbl_sys.pack(anchor="w", pady=(0,6))

        # colors
        colors = ["blue","red","green","black","gray","orange","purple","brown","pink","cyan"]
        ttk.Label(left, text="Colors: measured / reference / matched").pack(anchor="w")
        self.combo_c_meas = ttk.Combobox(left, values=colors, width=10); self.combo_c_meas.set("blue"); self.combo_c_meas.pack(anchor="w")
        self.combo_c_ref  = ttk.Combobox(left, values=colors, width=10); self.combo_c_ref.set("red"); self.combo_c_ref.pack(anchor="w")
        self.combo_c_match= ttk.Combobox(left, values=colors, width=10); self.combo_c_match.set("green"); self.combo_c_match.pack(anchor="w")

        # δ scan/fixed
        ttk.Separator(left, orient="horizontal").pack(fill="x", pady=6)
        ttk.Label(left, text="δ scan range (min,max,step)").pack(anchor="w")
        self.e_d_range = ttk.Entry(left, width=18); self.e_d_range.insert(0, "-0.5,0.5,0.1"); self.e_d_range.pack(anchor="w")
        ttk.Label(left, text="Objective").pack(anchor="w")
        self.combo_crit = ttk.Combobox(left, values=["rms", "uncprod"], width=10); self.combo_crit.set("rms"); self.combo_crit.pack(anchor="w")
        ttk.Button(left, text="SCAN: grid search δ & set best", command=self._scan_delta).pack(anchor="w", pady=(6,6))

        ttk.Label(left, text="Δ (fixed)").pack(anchor="w")
        self.e_delta = ttk.Entry(left, width=12); self.e_delta.insert(0, "0.0"); self.e_delta.pack(anchor="w", pady=(0,6))
        ttk.Button(left, text="FIXED: fit with fixed δ (report)", command=self._run_urlap).pack(anchor="w", pady=(6,8))

        # plot
        fig = Figure(figsize=(8.8,6.2), dpi=100)
        ax = fig.add_subplot(111)
        ax.set_xlabel("2θ (°)"); ax.set_ylabel("Relative Intensity (norm to 100)")
        self._fig = fig; self._ax = ax
        self._canvas = FigureCanvasTkAgg(fig, master=right)
        self._canvas.get_tk_widget().pack(fill="both", expand=True)

        # Tab 2: Peak List
        tab_peaks = ttk.Frame(nb); nb.add(tab_peaks, text="Peak List / Matching")
        top = ttk.Frame(tab_peaks); top.pack(side=tk.TOP, fill="x", padx=8, pady=6)
        frm_table = ttk.Frame(tab_peaks); frm_table.pack(side=tk.TOP, fill="both", expand=True, padx=8, pady=(0,8))

        ttk.Label(top, text="Smoothing window (°)").pack(side=tk.LEFT)
        self.e_smooth = ttk.Entry(top, width=8); self.e_smooth.insert(0, "0.15"); self.e_smooth.pack(side=tk.LEFT, padx=(2,8))

        ttk.Label(top, text="Min height (%)").pack(side=tk.LEFT)
        self.e_minH = ttk.Entry(top, width=8); self.e_minH.insert(0, "1.0"); self.e_minH.pack(side=tk.LEFT, padx=(2,8))

        ttk.Label(top, text="Min spacing (°)").pack(side=tk.LEFT)
        self.e_minSep = ttk.Entry(top, width=8); self.e_minSep.insert(0, "0.10"); self.e_minSep.pack(side=tk.LEFT, padx=(2,8))

        self.var_refine = tk.BooleanVar(value=True)
        ttk.Checkbutton(top, text="Parabolic refine", variable=self.var_refine).pack(side=tk.LEFT, padx=(4,10))

        ttk.Button(top, text="Build Peak List", command=self._gen_peaklist).pack(side=tk.LEFT, padx=(4,10))

        ttk.Separator(top, orient="vertical").pack(side=tk.LEFT, fill="y", padx=8)

        ttk.Label(top, text="Match tolerance (°)").pack(side=tk.LEFT)
        self.e_tol = ttk.Entry(top, width=8); self.e_tol.insert(0, "0.15"); self.e_tol.pack(side=tk.LEFT, padx=(2,8))

        ttk.Button(top, text="Match to CIF peaks (nearest)", command=self._do_match).pack(side=tk.LEFT, padx=(4,10))

        ttk.Button(top, text="Export Peak List…", command=self._export_peaklist).pack(side=tk.LEFT, padx=(8,0))

        ttk.Button(top, text="Select all", command=lambda: self._set_use_for_all(True)).pack(side=tk.RIGHT, padx=(8,4))
        ttk.Button(top, text="Select none", command=lambda: self._set_use_for_all(False)).pack(side=tk.RIGHT)

        # table
        cols = ("no","two","d","inten","h","k","l","two_calc","d2","phase","use")
        heads = {
            "no":"No.","two":"2θ (deg)","d":"d (Å)","inten":"Intensity (%)",
            "h":"h","k":"k","l":"l","two_calc":"2θ calc","d2":"|Δ2θ|",
            "phase":"Phase","use":"Use"
        }
        self.tv = ttk.Treeview(frm_table, columns=cols, show="headings", height=16)
        vsb = ttk.Scrollbar(frm_table, orient="vertical", command=self.tv.yview)
        self.tv.configure(yscroll=vsb.set)
        self.tv.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        frm_table.grid_rowconfigure(0, weight=1)
        frm_table.grid_columnconfigure(0, weight=1)
        for c in cols:
            self.tv.heading(c, text=heads[c])
        widths = {"no":50,"two":90,"d":90,"inten":100,"h":40,"k":40,"l":40,"two_calc":90,"d2":70,"phase":140,"use":60}
        for c in cols:
            self.tv.column(c, width=widths.get(c,80), anchor=tk.CENTER if c in ("no","h","k","l","use") else tk.W)
        self.tv.bind("<Double-1>", self._on_tv_double)

        # Tab 3: Report
        tab_rep = ttk.Frame(nb); nb.add(tab_rep, text="Report")
        self.txt_report = tk.Text(tab_rep, wrap="none")
        self.txt_report.pack(fill="both", expand=True, padx=6, pady=6)

    # ---------- File/Ref peaks ----------
    def _pick_xy(self):
        p = filedialog.askopenfilename(title="Choose XRD 2-column file", filetypes=[("Text/CSV","*.txt *.csv"), ("All","*.*")])
        if p:
            self.e_xy.delete(0, tk.END); self.e_xy.insert(0, p)
            try:
                self._xy = read_xy_file(p)
                self._redraw_plot()
            except Exception as e:
                messagebox.showerror("Read failed", str(e))

    def _pick_cif(self):
        p = filedialog.askopenfilename(title="Choose CIF file", filetypes=[("CIF","*.cif"), ("All","*.*")])
        if p:
            self.e_cif.delete(0, tk.END); self.e_cif.insert(0, p)

    def _calc_reference(self):
        if not PYMATGEN_OK:
            messagebox.showerror("Missing dependency", "Please install 'pymatgen' to read CIF/MP and compute reference peaks."); return
        if self._xy is None:
            messagebox.showwarning("Info", "Please load measured XRD first."); return

        try:
            wav = float(self.e_wav.get())
        except Exception:
            wav = 1.5406
        self._phase_name = (self.e_phase.get() or "Phase").strip()

        struct = None
        if self.phase_src.get() == "cif":
            path = self.e_cif.get().strip()
            if not path:
                messagebox.showerror("Error", "Please choose a CIF file."); return
            try:
                struct = Structure.from_file(path)
                # default phase name from file name
                if not self._phase_name or self._phase_name == "Phase":
                    base = os.path.basename(path)
                    self._phase_name = os.path.splitext(base)[0]
                    self.e_phase.delete(0, tk.END); self.e_phase.insert(0, self._phase_name)
            except Exception as e:
                messagebox.showerror("CIF read failed", str(e)); return
        else:
            key = self.e_mpkey.get().strip()
            q = self.e_mpid.get().strip()
            if not key or not q:
                messagebox.showerror("Error", "Please input MP API Key and mp-id/formula."); return
            try:
                with MPRester(key) as mpr:
                    if q.lower().startswith("mp-"):
                        docs = mpr.summary.search(material_ids=[q], fields=["structure"], num_chunks=1, chunk_size=1)
                        struct = docs[0].structure if docs else None
                    else:
                        docs = mpr.summary.search(formula=q, fields=["structure"], num_chunks=1, chunk_size=1)
                        struct = docs[0].structure if docs else None
                if not self._phase_name or self._phase_name == "Phase":
                    self._phase_name = q
                    self.e_phase.delete(0, tk.END); self.e_phase.insert(0, self._phase_name)
            except Exception as e:
                messagebox.showerror("Materials Project error", str(e)); return
            if struct is None:
                messagebox.showerror("Not found", "Material not found on Materials Project."); return

        # auto-detect crystal system
        code, name = self._infer_system(struct)
        self._auto_code = code; self._auto_name = name; self._struct = struct
        self.lbl_sys.config(text=f"Crystal system (auto): {name} (code={code})")

        # compute reference peaks
        try:
            calc = XRDCalculator(wavelength=wav)
            patt = calc.get_pattern(struct)
            self._calc_peaks = []
            for two, hkls in zip(patt.x, patt.hkls):
                if not hkls: continue
                h,k,l = [int(v) for v in hkls[0]["hkl"]]
                self._calc_peaks.append((float(two), (h,k,l)))
        except Exception as e:
            messagebox.showerror("XRD calculation failed", str(e)); return

        self._redraw_plot()

    def _infer_system(self, struct):
        try:
            sga = SpacegroupAnalyzer(struct, symprec=1e-3)
            sys_name = sga.get_crystal_system()
        except Exception:
            sys_name = None
        a = struct.lattice.a; b = struct.lattice.b; c = struct.lattice.c
        alpha = struct.lattice.alpha; beta = struct.lattice.beta; gamma = struct.lattice.gamma
        def near90(x): return abs(x-90.0) < 0.5
        if sys_name is None:
            if (not near90(alpha)) or (not near90(beta)) or (not near90(gamma)):
                non90 = sum([not near90(alpha), not near90(beta), not near90(gamma)])
                sys_name = "monoclinic" if non90 == 1 else "triclinic"
            else:
                if abs(a-b)<1e-6 and abs(b-c)<1e-6: sys_name="cubic"
                elif abs(a-b)<1e-6 or abs(b-c)<1e-6 or abs(a-c)<1e-6: sys_name="tetragonal"
                else: sys_name="orthorhombic"
        m = {"triclinic":1,"monoclinic":2,"orthorhombic":4,"tetragonal":5,"cubic":6,"trigonal":7,"hexagonal":8}
        if sys_name=="monoclinic":
            return 2, "MONOCLINIC(b)"
        return m.get(sys_name,4), sys_name.upper()

    # ---------- Peak List ----------
    def _gen_peaklist(self):
        if self._xy is None:
            messagebox.showwarning("Info", "Please load measured XRD first."); return
        try:
            wav = float(self.e_wav.get())
        except Exception:
            wav = 1.5406
        try:
            win = float(self.e_smooth.get()); minH = float(self.e_minH.get()); sep = float(self.e_minSep.get())
        except Exception:
            messagebox.showerror("Invalid format", "Check smoothing/min height/min spacing."); return
        peaks = find_peaks_simple(self._xy[0], self._xy[1], win, minH, sep, self.var_refine.get())
        self._peaklist = []
        for idx, (two, inten) in enumerate(peaks, start=1):
            d = d_from_two(two, wav, 0.0)  # δ handled during fitting; display-only here
            self._peaklist.append(PeakRow(two_obs=two, inten=float(inten), d_obs=float(d)))
        self._refresh_table()
        self._redraw_plot()

    def _do_match(self):
        if not self._calc_peaks:
            messagebox.showwarning("Info", "Please compute reference peaks first (CIF/MP)."); return
        if not self._peaklist:
            messagebox.showwarning("Info", "Please build the Peak List first."); return
        try:
            tol = float(self.e_tol.get())
        except Exception:
            tol = 0.15
        # nearest neighbour + unique assignment
        obs = np.array([r.two_obs for r in self._peaklist], float)
        used = np.zeros(obs.shape[0], dtype=bool)
        for two_calc, (h,k,l) in self._calc_peaks:
            i = int(np.argmin(np.abs(obs - two_calc)))
            if used[i]:
                continue
            if abs(obs[i]-two_calc) <= tol:
                used[i] = True
                r = self._peaklist[i]
                r.h, r.k, r.l = int(h), int(k), int(l)
                r.two_calc = float(two_calc)
                r.d2 = abs(float(obs[i]-two_calc))
                r.phase = self._phase_name
                r.use = True
        self._refresh_table()
        self._redraw_plot()

    def _export_peaklist(self):
        if not self._peaklist:
            messagebox.showwarning("No data", "Please build the Peak List first."); return
        p = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV","*.csv"), ("All","*.*")], title="Export Peak List as…")
        if not p:
            return
        with open(p, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["No","2theta_obs_deg","d_obs_A","Intensity_%","h","k","l","2theta_calc_deg","abs_d2theta_deg","phase","use"])
            for i, r in enumerate(self._peaklist, start=1):
                w.writerow([i, r.two_obs, r.d_obs, r.inten, r.h, r.k, r.l, r.two_calc, r.d2, r.phase, int(bool(r.use))])
        messagebox.showinfo("Exported", os.path.basename(p))

    def _set_use_for_all(self, flag: bool):
        for r in self._peaklist:
            r.use = bool(flag)
        self._refresh_table()
        self._redraw_plot()

    def _on_tv_double(self, event):
        iid = self.tv.focus()
        if not iid:
            return
        idx = int(self.tv.set(iid, "no")) - 1
        if 0 <= idx < len(self._peaklist):
            self._peaklist[idx].use = not self._peaklist[idx].use
            self._refresh_table_row(idx)
            self._redraw_plot()

    def _refresh_table(self):
        self.tv.delete(*self.tv.get_children())
        for i, r in enumerate(self._peaklist, start=1):
            self.tv.insert("", "end", values=(
                i,
                f"{r.two_obs:.5f}",
                f"{r.d_obs:.5f}",
                f"{r.inten:.1f}",
                "" if r.h is None else r.h,
                "" if r.k is None else r.k,
                "" if r.l is None else r.l,
                "" if r.two_calc is None else f"{r.two_calc:.5f}",
                "" if r.d2 is None else f"{r.d2:.5f}",
                r.phase or "",
                "✓" if r.use else ""
            ))

    def _refresh_table_row(self, idx: int):
        iid = self.tv.get_children()[idx]
        r = self._peaklist[idx]
        self.tv.item(iid, values=(
            idx+1,
            f"{r.two_obs:.5f}",
            f"{r.d_obs:.5f}",
            f"{r.inten:.1f}",
            "" if r.h is None else r.h,
            "" if r.k is None else r.k,
            "" if r.l is None else r.l,
            "" if r.two_calc is None else f"{r.two_calc:.5f}",
            "" if r.d2 is None else f"{r.d2:.5f}",
            r.phase or "",
            "✓" if r.use else ""
        ))

    # ---------- Fitting ----------
    def _get_peaks_for_fit(self) -> Optional[np.ndarray]:
        rows = [(r.h, r.k, r.l, r.two_obs) for r in self._peaklist if r.use and (r.h is not None)]
        if rows:
            try:
                arr = np.asarray(rows, float)
                if arr.ndim == 2 and arr.shape[1] == 4:
                    return arr
            except Exception:
                pass
        return None

    def _scan_delta(self):
        peaks = self._get_peaks_for_fit()
        if peaks is None:
            messagebox.showerror("No usable peaks", "Build & match the Peak List, then select rows to use."); return
        if self._auto_code is None:
            messagebox.showerror("No crystal system", "Compute reference peaks to auto-detect the crystal system first."); return
        try:
            wav = float(self.e_wav.get())
        except Exception:
            wav = 1.5406
        try:
            dmin, dmax, step = [float(v) for v in (self.e_d_range.get().strip() or "-0.5,0.5,0.1").split(",")]
        except Exception:
            messagebox.showerror("Invalid format", "δ range must be 'min,max,step'."); return
        crit = self.combo_crit.get().strip() or "rms"
        code = self._auto_code
        try:
            best = URLAP.select_delta(peaks, wav, code, dmin, dmax, step, crit)
        except Exception as e:
            messagebox.showerror("Scan failed", str(e)); return
        if not best:
            messagebox.showwarning("No better δ", "No better δ found."); return
        d_best = float(best["delta"])
        self.e_delta.delete(0, tk.END); self.e_delta.insert(0, f"{d_best:.6f}")
        messagebox.showinfo("Scan done", f"Best δ = {d_best:.6f} (objective={crit})")

    def _run_urlap(self):
        peaks = self._get_peaks_for_fit()
        if peaks is None:
            messagebox.showerror("No usable peaks", "Build & match the Peak List, then select rows to use."); return
        if self._auto_code is None:
            messagebox.showerror("No crystal system", "Compute reference peaks to auto-detect the crystal system first."); return
        try:
            wav = float(self.e_wav.get())
        except Exception:
            wav = 1.5406
        try:
            delta = float(self.e_delta.get())
        except Exception:
            delta = 0.0
        code = self._auto_code
        title = self.e_title.get().strip() or "Sample"
        try:
            info = URLAP.compute_from_system(peaks, wav, delta, code)
            report = URLAP.render_report(title, wav, delta, code, info)
        except Exception as e:
            messagebox.showerror("URLAP failed", str(e)); return
        self.txt_report.delete("1.0", tk.END)
        self.txt_report.insert(tk.END, report)

    # ---------- Plot ----------
    def _redraw_plot(self):
        ax = self._ax; ax.clear()
        ax.set_xlabel("2θ (°)"); ax.set_ylabel("Relative Intensity (norm to 100)")
        c_meas = self.combo_c_meas.get().strip() or None
        c_ref  = self.combo_c_ref.get().strip() or None
        c_m    = self.combo_c_match.get().strip() or None

        if self._xy is not None:
            x,y = self._xy
            ax.plot(x, y, lw=1.0, color=c_meas if c_meas else None)

        # reference peaks
        if self._calc_peaks:
            ymax = 100.0
            for two, (h,k,l) in self._calc_peaks:
                ax.vlines(two, 0, ymax*0.2, linestyles="dashed", colors=c_ref if c_ref else None)
                ax.text(two, ymax*0.22, f"{h}{k}{l}", rotation=90, fontsize=7, ha="center", va="bottom", color=c_ref if c_ref else None)

        # matched & selected peaks
        for r in self._peaklist:
            if r.use and (r.h is not None):
                ax.plot([r.two_obs], [90.0], marker="s", ms=5, color=c_m if c_m else None)
                ax.text(r.two_obs, 92.0, f"({r.h}{r.k}{r.l})", fontsize=8, ha="center", color=c_m if c_m else None)

        self._canvas.draw()

if __name__ == "__main__":
    app = URLAPGUI()
    app.mainloop()
