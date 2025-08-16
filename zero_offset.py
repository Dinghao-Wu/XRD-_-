
# urlap_report_v3.py
# -*- coding: utf-8 -*-
"""
Lattice codes:
 1: triclinic
 2: monoclinic(b)  (unique b, beta != 90°)
 3: monoclinic(c)  (unique c, gamma != 90°)
 4: orthorhombic
 5: tetragonal     (a=b!=c)
 6: cubic          (a=b=c)
 7: trigonal       (hex setting; same d-spacing as hexagonal)
 8: hexagonal

Conventions:
- delta in degrees, shown as "2THETA ORIGIN".
- Shown 2THETA_OBS = OBS_raw - delta.
- 2THETA_CAL computed from lattice only.
- "CAL-OBS" prints OBS_shown - CAL (matches your previous outputs).

Fitting per system:
- Linear regression in Q-space: Q = (2 sin(theta)/lambda)^2 = 1/d^2.
- Each system maps to a linear model Q = X @ m (through origin). Solve m by least squares.
- From m (metric parameters) derive (a,b,c,angles) and 2theta calc.
- Uncertainties: from cov(m) = sigma^2 (X'X)^-1 with sigma^2 = RSS/(n-p),
  propagate to std(a), std(b), std(c) via numerical Jacobian J of (a,b,c) wrt m:
  Var([a,b,c]) ≈ J cov(m) J^T ⇒ take diagonal sqrt for da,db,dc.
- Scan mode: delta in [min,max] with step; choose delta by criterion:
  "rms": minimize RMS(CAL - OBS_shown)
  "uncprod": minimize product of (da,db,dc) with degeneracy rules:
     cubic: da^3; tetragonal/trigonal/hex: da^2*dc; orthorhombic: da*db*dc;
     monoclinic: da*db*dc; triclinic: da*db*dc.

CLI:
  python urlap_report_v3.py -i INPUT.TXT --mode fixed
  python urlap_report_v3.py -i INPUT.TXT --mode scan --delta-range -0.5 0.5 --step 0.1 --select-by uncprod
"""

import argparse, sys, math
import numpy as np

def parse_raw(path: str):
    lines = [ln.rstrip('\n') for ln in open(path, 'r', encoding='utf-8', errors='ignore')]
    title = lines[0].strip()
    system_code = None
    if len(lines) >= 2:
        parts = lines[1].split()
        if parts:
            try: system_code = int(float(parts[0]))
            except: pass
    wav = None; delta = None
    if len(lines) >= 3:
        parts = lines[2].split()
        if len(parts) >= 1:
            try: wav = float(parts[0])
            except: pass
        if len(parts) >= 2:
            try: delta = float(parts[1])
            except: pass
    peaks = []
    for ln in lines[3:]:
        if not ln.strip(): continue
        parts = ln.split()
        if len(parts) < 4: continue
        try:
            h = int(float(parts[0])); k = int(float(parts[1])); l = int(float(parts[2]))
            if h >= 1000: break
            two = float(parts[3])
            peaks.append([h,k,l,two])
        except:
            continue
    if not peaks:
        raise ValueError("No peaks parsed from RAW")
    return title, system_code, wav, delta, np.asarray(peaks, float)

def obs_to_Q(two_deg, wavelength):
    theta = np.radians(0.5*two_deg)
    return (2.0*np.sin(theta)/wavelength)**2  # equals 1/d^2

def two_from_dinv2(dinv2, wavelength):
    # dinv2 = 1/d^2
    d = 1.0/np.sqrt(np.maximum(dinv2, 1e-300))
    arg = np.clip(wavelength/(2.0*d), -1.0, 1.0)
    th = np.degrees(np.arcsin(arg))
    return 2.0*th

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
    RSS = float(np.sum(resid**2))
    sigma2 = RSS / dof
    cov_m = sigma2 * pinv
    return m, cov_m, RSS, dof

# System-specific linear designs and parameter extraction

def design_triclinic(hkl):
    h = hkl[:,0]; k = hkl[:,1]; l = hkl[:,2]
    X = np.column_stack([h*h, k*k, l*l, 2*h*k, 2*h*l, 2*k*l])
    names = ["G11","G22","G33","G12","G13","G23"]
    return X, names

def extract_triclinic(m):
    G11,G22,G33,G12,G13,G23 = m.tolist()
    Gstar = np.array([[G11,G12,G13],[G12,G22,G23],[G13,G23,G33]], float)
    g = np.linalg.inv(Gstar)
    a = math.sqrt(g[0,0]); b = math.sqrt(g[1,1]); c = math.sqrt(g[2,2])
    cos_alpha = g[1,2]/(b*c); cos_beta = g[0,2]/(a*c); cos_gamma = g[0,1]/(a*b)
    cos_alpha = float(np.clip(cos_alpha, -1.0, 1.0))
    cos_beta  = float(np.clip(cos_beta , -1.0, 1.0))
    cos_gamma = float(np.clip(cos_gamma, -1.0, 1.0))
    alpha = math.degrees(math.acos(cos_alpha))
    beta  = math.degrees(math.acos(cos_beta))
    gamma = math.degrees(math.acos(cos_gamma))
    return a,b,c,alpha,beta,gamma,Gstar

def dinv2_triclinic(hkl, Gstar):
    # 1/d^2 = [h k l] G* [h k l]^T
    H = hkl.astype(float)
    vals = np.einsum('...i,ij,...j->...', H, Gstar, H)
    return vals

def design_monoclinic_b(hkl):
    # unique b, beta != 90°, formula:
    # 1/d^2 = (h^2/a^2 + k^2/b^2 + l^2/c^2 - 2 h l cosβ / (a c)) / sin^2 β
    # Define linear parameters:
    # m1 = 1/(a^2 s^2), m2 = 1/(b^2 s^2), m3 = 1/(c^2 s^2), m4 = -2 cosβ/(a c s^2)
    h = hkl[:,0].astype(float); k = hkl[:,1].astype(float); l = hkl[:,2].astype(float)
    X = np.column_stack([h*h, k*k, l*l, h*l])
    names = ["m1","m2","m3","m4"]
    return X, names

def extract_monoclinic_b(m):
    m1,m2,m3,m4 = m.tolist()
    # Compute sin^2 β from relation E^2/(4 m1 m3) = cos^2 β
    cos2 = (m4*m4) / (4.0 * max(m1,1e-300) * max(m3,1e-300))
    cos2 = float(np.clip(cos2, 0.0, 1.0))
    sin2 = max(1.0 - cos2, 1e-12)
    s = math.sqrt(sin2)
    a = 1.0 / math.sqrt(max(m1,1e-300) * sin2)
    b = 1.0 / math.sqrt(max(m2,1e-300) * sin2)
    c = 1.0 / math.sqrt(max(m3,1e-300) * sin2)
    cosb = -0.5 * m4 * a * c * sin2  # from m4 = -2 cosβ/(a c s^2) → cosβ = -m4 a c s^2 / 2
    cosb = float(np.clip(cosb, -1.0, 1.0))
    beta = math.degrees(math.acos(cosb))
    return a,b,c, 90.0, beta, 90.0

def dinv2_monoclinic_b(hkl, a,b,c, beta_deg):
    beta = math.radians(beta_deg)
    s2 = (math.sin(beta))**2
    cb = math.cos(beta)
    h = hkl[:,0].astype(float); k = hkl[:,1].astype(float); l = hkl[:,2].astype(float)
    return (h*h)/(a*a*s2) + (k*k)/(b*b*s2) + (l*l)/(c*c*s2) - (2*h*l*cb)/(a*c*s2)

def design_monoclinic_c(hkl):
    # unique c, gamma != 90°:
    # 1/d^2 = (h^2/a^2 + k^2/b^2 + l^2/c^2 - 2 h k cosγ / (a b)) / sin^2 γ
    # m1 = 1/(a^2 s^2), m2 = 1/(b^2 s^2), m3 = 1/(c^2 s^2), m4 = -2 cosγ/(a b s^2)
    h = hkl[:,0].astype(float); k = hkl[:,1].astype(float); l = hkl[:,2].astype(float)
    X = np.column_stack([h*h, k*k, l*l, h*k])
    names = ["m1","m2","m3","m4"]
    return X, names

def extract_monoclinic_c(m):
    m1,m2,m3,m4 = m.tolist()
    cos2 = (m4*m4) / (4.0 * max(m1,1e-300) * max(m2,1e-300))
    cos2 = float(np.clip(cos2, 0.0, 1.0))
    sin2 = max(1.0 - cos2, 1e-12)
    a = 1.0 / math.sqrt(max(m1,1e-300) * sin2)
    b = 1.0 / math.sqrt(max(m2,1e-300) * sin2)
    c = 1.0 / math.sqrt(max(m3,1e-300) * sin2)
    cosg = -0.5 * m4 * a * b * sin2
    cosg = float(np.clip(cosg, -1.0, 1.0))
    gamma = math.degrees(math.acos(cosg))
    return a,b,c, 90.0, 90.0, gamma

def dinv2_monoclinic_c(hkl, a,b,c, gamma_deg):
    g = math.radians(gamma_deg)
    s2 = (math.sin(g))**2
    cg = math.cos(g)
    h = hkl[:,0].astype(float); k = hkl[:,1].astype(float); l = hkl[:,2].astype(float)
    return (h*h)/(a*a*s2) + (k*k)/(b*b*s2) + (l*l)/(c*c*s2) - (2*h*k*cg)/(a*b*s2)

def design_orthorhombic(hkl):
    h = hkl[:,0]; k = hkl[:,1]; l = hkl[:,2]
    X = np.column_stack([h*h, k*k, l*l])
    names = ["1/a^2","1/b^2","1/c^2"]
    return X, names

def extract_orthorhombic(m):
    ma, mb, mc = m.tolist()
    a = 1.0/math.sqrt(max(ma,1e-300))
    b = 1.0/math.sqrt(max(mb,1e-300))
    c = 1.0/math.sqrt(max(mc,1e-300))
    return a,b,c, 90.0, 90.0, 90.0

def dinv2_orthorhombic(hkl, a,b,c):
    h = hkl[:,0].astype(float); k = hkl[:,1].astype(float); l = hkl[:,2].astype(float)
    return (h*h)/(a*a) + (k*k)/(b*b) + (l*l)/(c*c)

def design_tetragonal(hkl):
    h = hkl[:,0]; k = hkl[:,1]; l = hkl[:,2]
    Sxy = (h*h + k*k).astype(float); Sz = (l*l).astype(float)
    X = np.column_stack([Sxy, Sz])
    names = ["1/a^2","1/c^2"]
    return X, names

def extract_tetragonal(m):
    ma, mc = m.tolist()
    a = 1.0/math.sqrt(max(ma,1e-300))
    c = 1.0/math.sqrt(max(mc,1e-300))
    return a,a,c, 90.0, 90.0, 90.0

def dinv2_tetragonal(hkl, a,c):
    h = hkl[:,0].astype(float); k = hkl[:,1].astype(float); l = hkl[:,2].astype(float)
    return (h*h + k*k)/(a*a) + (l*l)/(c*c)

def design_hex(hkl):
    h = hkl[:,0].astype(float); k = hkl[:,1].astype(float); l = hkl[:,2].astype(float)
    H = (4.0/3.0)*(h*h + h*k + k*k)
    L = l*l
    X = np.column_stack([H, L])
    names = ["1/a^2","1/c^2"]
    return X, names

def extract_hex(m):
    ma, mc = m.tolist()
    a = 1.0/math.sqrt(max(ma,1e-300))
    c = 1.0/math.sqrt(max(mc,1e-300))
    return a,a,c, 90.0, 90.0, 120.0

def dinv2_hex(hkl, a,c):
    h = hkl[:,0].astype(float); k = hkl[:,1].astype(float); l = hkl[:,2].astype(float)
    return (4.0/3.0)*(h*h + h*k + k*k)/(a*a) + (l*l)/(c*c)

def design_cubic(hkl):
    h = hkl[:,0]; k = hkl[:,1]; l = hkl[:,2]
    S = (h*h + k*k + l*l).astype(float)
    X = S.reshape(-1,1)
    names = ["1/a^2"]
    return X, names

def extract_cubic(m):
    ma = float(m[0])
    a = 1.0/math.sqrt(max(ma,1e-300))
    return a,a,a, 90.0, 90.0, 90.0

def dinv2_cubic(hkl, a):
    h = hkl[:,0].astype(float); k = hkl[:,1].astype(float); l = hkl[:,2].astype(float)
    return (h*h + k*k + l*l)/(a*a)

def system_design_and_extract(code):
    if code == 1:
        return design_triclinic, extract_triclinic, dinv2_triclinic, "TRICLINIC"
    if code == 2:
        return design_monoclinic_b, extract_monoclinic_b, dinv2_monoclinic_b, "MONOCLINIC(b)"
    if code == 3:
        return design_monoclinic_c, extract_monoclinic_c, dinv2_monoclinic_c, "MONOCLINIC(c)"
    if code == 4:
        return design_orthorhombic, extract_orthorhombic, dinv2_orthorhombic, "ORTHORHOMBIC"
    if code == 5:
        return design_tetragonal, extract_tetragonal, dinv2_tetragonal, "TETRAGONAL"
    if code == 6:
        return design_cubic, extract_cubic, dinv2_cubic, "CUBIC"
    if code == 7:
        return design_hex, extract_hex, dinv2_hex, "TRIGONAL"
    if code == 8:
        return design_hex, extract_hex, dinv2_hex, "HEXAGONAL"
    raise NotImplementedError("System code must be 1..8")

def numeric_jacobian(func, m, eps=1e-6):
    m = np.asarray(m, float)
    base = func(m)
    base = np.asarray(base, float)
    J = np.zeros((base.shape[0], m.size), float)
    for j in range(m.size):
        mj = m.copy()
        step = eps * max(1.0, abs(mj[j]))
        mj[j] += step
        f2 = np.asarray(func(mj), float)
        J[:, j] = (f2 - base) / step
    return J, base

def compute_from_system(peaks, wavelength, delta, code):
    hkl = peaks[:, :3].astype(int)
    obs_corr = peaks[:,3] - delta
    Q = obs_to_Q(obs_corr, wavelength)

    design, extract, dinv2_fun, label = system_design_and_extract(code)
    X, names = design(hkl)
    m, cov_m, RSS, dof = solve_ls(X, Q)

    # derive cell
    a,b,c, alpha, beta, gamma, *extra = extract(m)
    # uncertainties via numerical Jacobian of (a,b,c) wrt m
    def map_to_abc(mm):
        a2,b2,c2, *_ = extract(mm)
        return np.array([a2,b2,c2], float)
    J, abc = numeric_jacobian(map_to_abc, m, eps=1e-6)
    cov_abc = J @ cov_m @ J.T
    da, db, dc = [float(math.sqrt(max(cov_abc[i,i], 0.0))) for i in range(3)]

    # Build table
    if code == 1:
        Gstar = extra[0]
        dinv2 = dinv2_triclinic(hkl, Gstar)
    elif code == 2:
        dinv2 = dinv2_monoclinic_b(hkl, a,b,c, beta)
    elif code == 3:
        dinv2 = dinv2_monoclinic_c(hkl, a,b,c, gamma)
    elif code == 4:
        dinv2 = dinv2_orthorhombic(hkl, a,b,c)
    elif code == 5:
        dinv2 = dinv2_tetragonal(hkl, a,c)
    elif code == 6:
        dinv2 = dinv2_cubic(hkl, a)
    else:  # 7 or 8
        dinv2 = dinv2_hex(hkl, a,c)
    cal = two_from_dinv2(dinv2, wavelength)
    resid_deg = obs_corr - cal
    d_hkl = 1.0/np.sqrt(dinv2 + 1e-300)
    Sdisp = (hkl[:,0]**2 + hkl[:,1]**2 + hkl[:,2]**2).astype(float)

    info = {
        "label": label, "a":a, "b":b, "c":c, "alpha":alpha, "beta":beta, "gamma":gamma,
        "da":da, "db":db, "dc":dc, "cov_m": cov_m, "m": m,
        "table": np.column_stack([hkl, obs_corr, cal, resid_deg, d_hkl, Sdisp])
    }
    return info

def unc_product(info, code):
    # define uncertainty product per system
    if code == 6:  # cubic
        return info["da"]**3
    if code in (5,7,8):  # tetragonal/trigonal/hex
        return (info["da"]**2) * info["dc"]
    # others
    return info["da"] * info["db"] * info["dc"]

def rms_residual(info):
    resid = info["table"][:,5]  # OBS - CAL
    return float(np.sqrt(np.mean(resid**2)))

def render_report_generic(title, wavelength, delta, code, info):
    label = info["label"]
    a,b,c = info["a"], info["b"], info["c"]
    da,db,dc = info["da"], info["db"], info["dc"]
    alpha,beta,gamma = info["alpha"], info["beta"], info["gamma"]
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
        h,k,l, obs_show, cal, resid, d_hkl, S = row
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

def select_delta(peaks, wavelength, code, dmin, dmax, step, criterion):
    deltas = np.arange(dmin, dmax + 1e-15, step, dtype=float)
    best = None
    for d in deltas:
        info = compute_from_system(peaks, wavelength, d, code)
        score = rms_residual(info) if criterion=="rms" else unc_product(info, code)
        if (best is None) or (score < best["score"]):
            best = {"delta": d, "info": info, "score": score}
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i","--input", required=True)
    ap.add_argument("--mode", choices=["fixed","scan"], default="fixed")
    ap.add_argument("--wavelength", type=float, default=None)
    ap.add_argument("--delta", type=float, default=None)
    ap.add_argument("--delta-range", nargs=2, type=float, default=[-0.5, 0.5])
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--select-by", choices=["rms","uncprod"], default="uncprod")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    title, system_code, w_in, d_in, peaks = parse_raw(args.input)
    if system_code is None:
        raise ValueError("Lattice system code (1..8) must be provided on line 2 (first integer).")
    wavelength = args.wavelength if args.wavelength is not None else (w_in if w_in is not None else 1.5402)
    if args.mode == "fixed":
        delta = args.delta if args.delta is not None else (d_in if d_in is not None else 0.0)
        info = compute_from_system(peaks, wavelength, delta, system_code)
        report = render_report_generic(title, wavelength, delta, system_code, info)
    else:
        best = select_delta(peaks, wavelength, system_code, args.delta_range[0], args.delta_range[1], args.step, args.select_by)
        delta = best["delta"]; info = best["info"]
        report = render_report_generic(title, wavelength, delta, system_code, info)
    out_path = args.out or (args.input + ".out.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(report)
    print(out_path)

if __name__ == "__main__":
    main()
