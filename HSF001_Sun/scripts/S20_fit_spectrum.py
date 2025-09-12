#!/usr/bin/env python3
"""
Per-day GHz spectrum fits for RSTN noontime flux (P0: Sun-function).
- Fits log10 S_nu ~ alpha * log10(nu) + const across 2695–15400 MHz (GHz regime)
- Computes residual ratio at 2.695, 4.995, 8.800, 15.400 GHz
- Ranks "non-power-law" days by largest |residual-1|
- Writes a candidates CSV with quick B[G] estimates for s=2,3

Usage:
  python src/fit_spectrum.py data/rstn/2018.txt -o results/candidates/hsf_candidates_2018.csv
"""

from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress

from S10_rstn_parse import load_year, RSTN_FREQS

GHz_BANDS = [2695, 4995, 8800, 15400]
GHz_COLS  = [f"sfu_{f}MHz" for f in GHz_BANDS]
ALL_COLS  = [f"sfu_{f}MHz" for f in RSTN_FREQS]

def fit_day(row: pd.Series):
    fx, fy = [], []
    for f in GHz_BANDS:
        v = row[f"sfu_{f}MHz"]
        if pd.notna(v) and v > 0:
            fx.append(f); fy.append(v)
    if len(fx) < 3:
        return pd.Series({
            "alpha": np.nan, "r2": np.nan,
            "resid_2695": np.nan, "resid_4995": np.nan,
            "resid_8800": np.nan, "resid_15400": np.nan,
            "max_abs_resid": np.nan, "which_band": np.nan
        })

    X = np.log10(np.asarray(fx, float))
    Y = np.log10(np.asarray(fy, float))
    slope, intercept, r, p, stderr = linregress(X, Y)
    model = 10**(slope*np.log10(np.asarray(GHz_BANDS, float)) + intercept)
    obs   = np.array([row[c] for c in GHz_COLS], float)
    resid = obs / model  # >1 bump, <1 dip
    max_idx = int(np.nanargmax(np.abs(resid - 1.0)))
    which = GHz_BANDS[max_idx]
    return pd.Series({
        "alpha": slope, "r2": r**2,
        "resid_2695": resid[0], "resid_4995": resid[1],
        "resid_8800": resid[2], "resid_15400": resid[3],
        "max_abs_resid": float(np.nanmax(np.abs(resid - 1.0))),
        "which_band": which
    })

def bfield_from_freq_mhz(f_mhz: float, s: int) -> float:
    # f_ce ≈ 2.8 MHz/G * B  ⇒  B[G] ≈ f_MHz / (2.8 * s)
    return float(f_mhz) / (2.8 * s)

def main(infile: str, outfile: str):
    df = load_year(infile)
    fits = df.apply(fit_day, axis=1)
    out = pd.concat([df[["date"] + ALL_COLS], fits], axis=1)
    # rank
    ranked = out.dropna(subset=["max_abs_resid"]).sort_values("max_abs_resid", ascending=False)
    # quick B at the band with largest residual (s=2 and s=3)
    ranked["B_G_s2"] = ranked["which_band"].apply(lambda f: round(bfield_from_freq_mhz(f, 2), 1))
    ranked["B_G_s3"] = ranked["which_band"].apply(lambda f: round(bfield_from_freq_mhz(f, 3), 1))
    # write CSV
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(outfile, index=False)
    print(f"Wrote {outfile} with {len(ranked)} rows.")

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", help="Yearly RSTN file (e.g., data/rstn/2018.txt)")
    ap.add_argument("-o", "--out", required=True, help="Output CSV path")
    args = ap.parse_args()
    main(args.infile, args.out)
