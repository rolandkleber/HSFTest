#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- config ----------
ROOT   = Path(__file__).resolve().parents[1]
INFILE = ROOT / "data" / "raw" / "O2_Bogumil(2003)_293K_235-389nm.txt"
OUTDIR = ROOT / "results"; FIGDIR = OUTDIR / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

# Use the band in the filename
CROP_MIN_NM, CROP_MAX_NM = 121.0, 193.0

# Gas/path (only needed for absolute alpha; KK shape works even unscaled)
P_atm, T_K, mix = 1.0, 370.0, 1.0  # 370 K per dataset
# ---------------------------

kB = 1.380649e-23
c  = 2.99792458e8
atm_to_Pa = 101325.0

def load_sigma_txt(path: Path) -> pd.DataFrame:
    import re
    rows=[]
    with open(path,"r",encoding="utf-8",errors="ignore") as f:
        for ln in f:
            s=ln.strip()
            if not s or s[0] in "#;": 
                continue
            parts=re.split(r"[,\s]+", s)
            if len(parts)<2: 
                continue
            try:
                wl=float(parts[0]); sig=float(parts[1])
            except:
                continue
            rows.append((wl, sig))
    if not rows:
        raise SystemExit("No numeric rows parsed.")

    df=pd.DataFrame(rows, columns=["wl_nm","sigma_cm2"])

    # basic sanity: positive & reasonable range for O2 UV/Vis (in cm^2/molecule)
    df=df[(df["sigma_cm2"]>0) & (df["sigma_cm2"]<1e-18)]

    # average any duplicate wavelengths
    df=(df.groupby("wl_nm", as_index=False)
          .agg(sigma_cm2=("sigma_cm2","mean"))
          .sort_values("wl_nm")
          .reset_index(drop=True))
    return df


def kk_n_from_kappa(omega, kappa):
    """
    Returns n_minus_1(omega) via:
      n(ω) - 1 = (2/π) P ∫ Ω κ(Ω) / (Ω^2 - ω^2) * Ω dΩ
    Non-uniform ΔΩ handled by local half-interval weights.
    """
    omg = omega.astype(float)
    kap = kappa.astype(float)
    n = len(omg)

    # local integration weights (half-intervals)
    dΩ = np.zeros(n, dtype=float)
    dΩ[1:-1] = 0.5*(omg[2:] - omg[:-2])
    dΩ[0]    = (omg[1]-omg[0])
    dΩ[-1]   = (omg[-1]-omg[-2])

    out = np.zeros(n, dtype=float)
    for k in range(n):
        denom = (omg**2 - omg[k]**2)
        mask  = denom != 0.0
        out[k] = (2/np.pi) * np.sum( (omg[mask]*kap[mask]) * dΩ[mask] / denom[mask] )
    # remove tiny offset so that median ~ 0 (finite-window baseline)
    out -= np.median(out)
    return out

def to_omega(df: pd.DataFrame):
    wl_m  = df["wl_nm"].values * 1e-9
    omega = 2*np.pi * c / wl_m
    return omega, df["sigma_cm2"].values

def sigma_to_alpha(sig_cm2, P_atm=1.0, T_K=370.0, mix=1.0):
    kB = 1.380649e-23
    atm_to_Pa = 101325.0
    N = (P_atm*atm_to_Pa) / (kB*T_K) * mix      # [1/m^3]
    sig_m2 = sig_cm2 * 1e-4                     # cm^2 -> m^2
    alpha_m = N * sig_m2                        # [1/m]
    return alpha_m                              # [1/m] (keep SI here)

def alpha_to_kappa(alpha_m, omega):
    c = 2.99792458e8
    return alpha_m * c / (2.0*omega)            # dimensionless

def kk_pv_dispersion(omega, alpha):
    """ Principal-value discrete KK: (1/pi) ∑_{j≠k} α(ω_j) Δω_j / (ω_j-ω_k).
        Handles non-uniform grids by using local Δω. """
    omg = omega.astype(float)
    a   = alpha.astype(float)
    n   = len(omg)
    disp = np.zeros(n, dtype=float)
    # construct Δω_j as half-intervals
    dω = np.zeros(n, dtype=float)
    dω[1:-1] = 0.5*(omg[2:] - omg[:-2])
    dω[0]    = omg[1]-omg[0]
    dω[-1]   = omg[-1]-omg[-2]
    for k in range(n):
        diff = omg - omg[k]
        mask = diff != 0.0
        disp[k] = (1/np.pi) * np.sum( a[mask] * dω[mask] / diff[mask] )
    return disp

def main():
    # --- load σ(λ) and convert to ω, α, κ ---
    df = load_sigma_txt(INFILE)                      # wl_nm, sigma_cm2
    omega_raw, sigma = to_omega(df)                  # ω [rad/s], σ [cm^2/molecule]
    alpha_m = sigma_to_alpha(sigma, P_atm=P_atm, T_K=T_K, mix=mix)  # [1/m]
    kappa_raw = alpha_to_kappa(alpha_m, omega_raw)   # dimensionless

    # --- put κ on a uniform ω grid (needed for Hilbert FFT) ---
    npts = max(2048, 8*len(kappa_raw))
    omg_u = np.linspace(omega_raw.min(), omega_raw.max(), npts)
    kappa_u = np.interp(omg_u, omega_raw, kappa_raw)

    # light linear detrend to reduce DC leakage (on κ, not α)
    x = np.arange(npts, dtype=float)
    kappa_u = kappa_u - np.polyval(np.polyfit(x, kappa_u, 1), x)

    # --- strong taper + zero padding, then Hilbert ---
    from scipy.signal import windows, hilbert
    taper = windows.tukey(npts, alpha=0.6)           # strong edges
    kap_t = kappa_u * taper

    pad = 8 * npts                                   # generous zero-pad
    kap_pad = np.r_[kap_t, np.zeros(pad)]
    analytic = hilbert(kap_pad)                      # analytic signal
    disp_pad = np.imag(analytic)                     # Hilbert twin of κ
    disp = disp_pad[:npts] / (taper + 1e-12)         # undo taper for display

    # --- normalize shapes; allow sign flip; best scale ---
    a_shape = (kappa_u - np.mean(kappa_u)) / (np.max(np.abs(kappa_u)) + 1e-12)
    d_shape = (disp     - np.mean(disp))  / (np.max(np.abs(disp))  + 1e-12)
    if np.corrcoef(a_shape, d_shape)[0, 1] < 0:
        d_shape = -d_shape
    s_opt = float(np.dot(a_shape, d_shape) / (np.dot(d_shape, d_shape) + 1e-12))
    resid = a_shape - s_opt * d_shape
    rms   = float(np.sqrt(np.mean(resid**2)))
    corr  = float(np.corrcoef(a_shape, s_opt*d_shape)[0, 1])

    # --- plots ---
    # absorption (for context) — plot raw α (no detrend) in 1/cm
    plt.figure(figsize=(7,3.5))
    plt.plot(omega_raw, alpha_m*1e-2, lw=1)
    plt.xlabel("Angular frequency ω [rad/s]")
    plt.ylabel("Absorption α [1/cm]")
    plt.title(f"Measured absorption ({INFILE.name}, T≈{T_K:.0f} K)")
    plt.tight_layout(); plt.savefig(FIGDIR/"kk_absorption.png", dpi=200); plt.close()

    # κ vs KK twin on uniform grid
    plt.figure(figsize=(7,3.5))
    plt.plot(omg_u, a_shape,       label="κ(ω) (norm.)", lw=1)
    plt.plot(omg_u, s_opt*d_shape, label="KK twin (norm., scaled)", lw=1)
    plt.xlabel("ω [rad/s]"); plt.ylabel("normalized")
    plt.title("KK: extinction κ vs. dispersive twin (Hilbert, tapered)")
    plt.legend(); plt.tight_layout()
    plt.savefig(FIGDIR/"kk_dispersion.png", dpi=200); plt.close()

    # --- summary ---
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with open(OUTDIR/"kk_summary.txt","w",encoding="utf-8") as f:
        f.write(f"File: {INFILE.name}\n")
        f.write(f"Uniform grid points: {npts}\n")
        f.write(f"Optimal scale s ≈ {s_opt:.3f}\n")
        f.write(f"RMS mismatch ≈ {rms:.3e}\n")
        f.write(f"|corr(κ, KK)| ≈ {abs(corr):.3f}\n")
        f.write("Tip: crop to a single hump (e.g., 135–185 nm) if edges still dominate.\n")

    print("[OK] wrote figures →", FIGDIR, "and summary →", OUTDIR)


if __name__ == "__main__":
    main()
