#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10_linewidth_TP_sweep.py  (robust file loader)
Temperature/Pressure sweep on a single H2O HITRAN line (Voigt),
checks Doppler ~ T^0.5, collisional width ~ T^{-n_air}, and γ ∝ P.

Usage:
  python scripts/10_linewidth_TP_sweep.py ^
    --hitran data\raw\H2O_lines_7150-7220.csv ^
    --auto-strongest ^
    --T 250,300,350,400,500,700 ^
    --P 0.2,0.5,1.0 ^
    --instr-fwhm 0.10 ^
    --win 1.0 ^
    --out results\MOL-01_H2O_TP_Voigt
"""

import argparse, os, math, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import wofz
from scipy.optimize import curve_fit

# ---------------- Voigt + helpers ----------------

def voigt_profile(nu, nu0, sigma_g, gamma_l):
    sg = max(float(sigma_g), 1e-15)
    x = (nu - nu0) / (sg * np.sqrt(2.0))
    y = gamma_l / (sg * np.sqrt(2.0))
    z = x + 1j*y
    return np.real(wofz(z)) / (sg * np.sqrt(2*np.pi))

def fwhm_gaussian_from_sigma(sigma):
    return 2.0*np.sqrt(2.0*np.log(2.0)) * sigma

def sigma_from_fwhm_gaussian(fwhm):
    return fwhm / (2.0*np.sqrt(2.0*np.log(2.0)))

def fwhm_lorentz_from_gamma(gamma):
    return 2.0*gamma

def doppler_sigma_from_formula(nu0_cm1, T, m_kg):
    kB = 1.380649e-23; c = 299792458.0
    fwhm = nu0_cm1 * math.sqrt(8*kB*T*math.log(2)/(m_kg*c*c))
    return sigma_from_fwhm_gaussian(fwhm)

def model_absorbance(nu, amp, center, sigma_g, gamma_l, instr_sigma):
    sigma_eff = np.sqrt(sigma_g**2 + instr_sigma**2)
    return amp * voigt_profile(nu, center, sigma_eff, gamma_l)

# ---------------- Robust HITRAN file loader ----------------

LIKELY_SEPS = [",", r"\s+"]

def _try_load(path):
    # First try CSV with header row
    try:
        df = pd.read_csv(path)
        return df, "csv_header"
    except Exception:
        pass
    # Try whitespace with header
    try:
        df = pd.read_csv(path, sep=r"\s+")
        return df, "ws_header"
    except Exception:
        pass
    # Try headerless CSV
    try:
        df = pd.read_csv(path, header=None)
        return df, "csv_noheader"
    except Exception:
        pass
    # Try headerless whitespace
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    return df, "ws_noheader"

def load_hitran_minimal(path):
    df, mode = _try_load(path)
    cols = [str(c).strip().lower() for c in df.columns]

    # If the first "column names" are numeric, we likely have NO header.
    def _is_number(s):
        try:
            float(str(s).replace("e","E").replace("+","").replace("-","-"))
            return True
        except Exception:
            return False

    numeric_header = all(_is_number(c) for c in cols)

    # Map with header if we can
    name_map = {}
    if not numeric_header:
        # normalize some common variants
        alias = { "nu": None, "wavenumber": "nu", "wn": "nu",
                  "gamma_air": None, "γair": "gamma_air", "gair": "gamma_air",
                  "n_air": None, "nair": "n_air",
                  "sw": None, "s": "sw" }
        for c in df.columns:
            lc = str(c).strip().lower()
            if lc in alias and alias[lc] is not None:
                name_map[c] = alias[lc]
            elif lc in ["nu","gamma_air","n_air","sw"]:
                name_map[c] = lc
        df = df.rename(columns=name_map)
        missing = [k for k in ["nu","gamma_air","n_air"] if k not in df.columns]
        if missing:
            # fall through to headerless mapping
            numeric_header = True

    if numeric_header:
        # Assume first columns are: 0=nu, 1=sw, 2=gamma_air, 3=n_air
        # This matches your sample line where first 4 tokens were:
        # 7150.06732, 4.42e-24, 0.0827, 0.71, ...
        df = df.rename(columns={df.columns[0]:"nu", df.columns[1]:"sw",
                                df.columns[2]:"gamma_air", df.columns[3]:"n_air"})

    # Keep only needed + helpful
    keep = [c for c in ["nu","sw","gamma_air","n_air"] if c in df.columns]
    df = df[keep].copy()
    # Force numeric
    for c in keep:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["nu","gamma_air","n_air"])
    if "sw" not in df.columns:
        # add dummy sw to allow --auto-strongest
        df["sw"] = 1.0
    return df

# ---------------- Line selection ----------------

def pick_line(df, center=None, auto_strongest=False, isolation_margin=0.2):
    if center is not None:
        idx = (df["nu"] - center).abs().idxmin()
        chosen = df.loc[idx]
    elif auto_strongest:
        idx = df["sw"].astype(float).idxmax()
        chosen = df.loc[idx]
    else:
        raise ValueError("Provide --center or --auto-strongest")
    nu0 = float(chosen["nu"])
    neighbors = df["nu"].sort_values().values
    dists = np.sort(np.abs(neighbors - nu0))
    nearest = dists[1] if len(dists) > 1 else np.inf
    return chosen, bool(nearest >= isolation_margin), float(nearest)

# ---------------- Main ----------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hitran", required=True, help="HITRAN file (.csv/.out, headerless ok)")
    ap.add_argument("--center", type=float, default=None, help="Target line center (cm^-1). If omitted, use --auto-strongest.")
    ap.add_argument("--auto-strongest", action="store_true", help="Pick strongest line in the file (uses 'sw' if present)")
    ap.add_argument("--T", default="250,300,350,400,500,700", help="Temperatures in K (comma-separated)")
    ap.add_argument("--P", default="1.0", help="Pressures in atm (comma-separated)")
    ap.add_argument("--Tref", type=float, default=296.0)
    ap.add_argument("--Pref", type=float, default=1.0)
    ap.add_argument("--instr-fwhm", type=float, default=0.20)
    ap.add_argument("--win", type=float, default=1.0)
    ap.add_argument("--points", type=int, default=2001)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    df = load_hitran_minimal(args.hitran)

    chosen, isolated, nearest = pick_line(
    df,
    center=args.center,
    auto_strongest=args.auto_strongest,   # <-- correct
    isolation_margin=0.2
)
    nu0 = float(chosen["nu"])
    gamma_ref = float(chosen["gamma_air"])
    n_exp = float(chosen["n_air"])

    # constants
    m_H2O = 18.01528 * 1.66053906660e-27
    T_list = np.array([float(t) for t in args.T.split(",")])
    P_list = np.array([float(p) for p in args.P.split(",")])
    instr_sigma = sigma_from_fwhm_gaussian(args.instr_fwhm)
    nu_grid = np.linspace(nu0 - args.win, nu0 + args.win, args.points)
    rng = np.random.default_rng(1234)
    noise_sigma = 1e-4
    amp_true = 0.02

    rows = []
    for P in P_list:
        for T in T_list:
            sigma_g_true = doppler_sigma_from_formula(nu0, T, m_H2O)
            gamma_l_true = gamma_ref * (args.Tref/T)**n_exp * (P/args.Pref)

            sigma_eff = np.sqrt(sigma_g_true**2 + instr_sigma**2)
            A = amp_true * voigt_profile(nu_grid, nu0, sigma_eff, gamma_l_true)
            A_obs = A + rng.normal(0.0, noise_sigma, size=nu_grid.size)

            def model(nu, amp, center, sigma_g, gamma_l):
                return model_absorbance(nu, amp, center, sigma_g, gamma_l, instr_sigma)

            p0 = [amp_true*0.9, nu0, max(1e-6, sigma_g_true*1.1), max(1e-6, gamma_l_true*0.9)]
            bounds = ([0.0, nu0-0.1, 0.0, 0.0], [1.0, nu0+0.1, 0.5, 1.0])
            popt, pcov = curve_fit(model, nu_grid, A_obs, p0=p0, bounds=bounds, maxfev=20000)

            amp_fit, center_fit, sigma_g_fit, gamma_l_fit = popt
            resid = A_obs - model(nu_grid, *popt)
            chi2 = np.sum((resid/noise_sigma)**2); ndf = len(nu_grid)-len(popt)
            rows.append({
                "nu0_cm^-1": nu0, "isolated": bool(isolated), "nearest_line_distance_cm^-1": float(nearest),
                "T_K": float(T), "P_atm": float(P),
                "gamma_ref_cm^-1@Tref": gamma_ref, "n_exp": n_exp,
                "instr_fwhm_cm^-1": args.instr_fwhm,
                "sigma_g_fit": float(sigma_g_fit),
                "gamma_l_fit": float(gamma_l_fit),
                "fwhm_g_fit": float(fwhm_gaussian_from_sigma(sigma_g_fit)),
                "fwhm_l_fit": float(fwhm_lorentz_from_gamma(gamma_l_fit)),
                "chi2_ndf": float(chi2/ndf if ndf>0 else np.nan)
            })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(os.path.join(args.out, "TP_fit_params.csv"), index=False)

    # exponent checks per pressure
    summaries = []
    for P in P_list:
        d = out_df[out_df["P_atm"]==P].sort_values("T_K")
        slope_d, _ = np.polyfit(np.log(d["T_K"]), np.log(d["fwhm_g_fit"]), 1)   # expect ~ +0.5
        slope_l, _ = np.polyfit(np.log(d["T_K"]), np.log(d["gamma_l_fit"]), 1) # expect ~ -n_exp
        summaries.append({
            "P_atm": float(P),
            "doppler_exp_fit": float(slope_d), "expected_doppler": 0.5,
            "lorentz_exp_fit": float(slope_l), "expected_lorentz": float(-n_exp)
        })
    pd.DataFrame(summaries).to_csv(os.path.join(args.out, "TP_exponents_by_pressure.csv"), index=False)

    # pressure proportionality at median T
# ---- Pressure proportionality at ~Tmid (robust) ----
# use the actually present temperature closest to the median
    T_values = np.sort(out_df["T_K"].unique().astype(float))
    Tmid = float(np.median(T_values))
    T_near = float(T_values[np.argmin(np.abs(T_values - Tmid))])

    dT = out_df[np.isclose(out_df["T_K"].astype(float), T_near)].sort_values("P_atm")

    if len(dT) >= 2:
        slopeP, interceptP = np.polyfit(dT["P_atm"].astype(float), dT["gamma_l_fit"].astype(float), 1)
    else:
        slopeP, interceptP = float("nan"), float("nan")

    expected_slope = gamma_ref * (args.Tref / T_near) ** n_exp


    # plots (simple)
    plt.figure(figsize=(6,5))
    for P in P_list:
        d = out_df[out_df["P_atm"]==P].sort_values("T_K")
        plt.loglog(d["T_K"], d["fwhm_g_fit"], marker="o", label=f"P={P} atm")
    plt.xlabel("T (K)"); plt.ylabel("Gaussian FWHM (cm$^{-1}$)"); plt.title("Doppler width vs T"); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.out, "doppler_vs_T.png"), dpi=150); plt.close()

    plt.figure(figsize=(6,5))
    for P in P_list:
        d = out_df[out_df["P_atm"]==P].sort_values("T_K")
        plt.loglog(d["T_K"], d["gamma_l_fit"], marker="s", label=f"P={P} atm")
    plt.xlabel("T (K)"); plt.ylabel("Lorentz HWHM γ (cm$^{-1}$)"); plt.title("Collisional width vs T"); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.out, "lorentz_vs_T.png"), dpi=150); plt.close()

    plt.figure(figsize=(6,5))
    plt.plot(dT["P_atm"], dT["gamma_l_fit"], marker="^")
    plt.xlabel("P (atm)"); plt.ylabel(f"γ(T≈{int(Tmid)}K) (cm$^{{-1}}$)"); plt.title("Collisional width vs Pressure"); plt.tight_layout()
    plt.savefig(os.path.join(args.out, "lorentz_vs_P_Tmid.png"), dpi=150); plt.close()

    # summary
    summary = {
        "nu0_cm^-1": nu0, "isolated_line": True, "nearest_line_distance_cm^-1": nearest,
        "gamma_ref_cm^-1@Tref": gamma_ref, "n_exp_from_HITRAN": n_exp,
        "exponents_by_pressure": summaries,
        "pressure_slope_fit_gamma_vs_P_at_Tmid": float(slopeP),
        "expected_pressure_slope_at_Tmid": float(expected_slope)
    }
    with open(os.path.join(args.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

if __name__ == "__main__":
    main()
