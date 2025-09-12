#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
P0 analysis for RSTN Palehua yearly data.

Inputs  (wide CSV):
  data/10_2018parsed.csv  with columns:
    date, station, sfu_245MHz, sfu_410MHz, ..., sfu_15400MHz

Outputs:
  results/p0_summary.txt
  results/figures/p0_slope.png
  results/figures/p0_spectrum.png
"""

import sys
import re
import math
from pathlib import Path

import numpy as np
import pandas as pd

# Try scipy for nicer fit; fall back to numpy polyfit
try:
    from scipy.stats import linregress
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

import matplotlib.pyplot as plt


# -------------------------- config ---------------------------------

QUIET_PERCENTILE = 0.10     # 10th percentile as "quiet-Sun"
HI_FREQ_MIN_MHZ  = 4995     # use >= this MHz for slope fit
LOW_TRIPLE        = (245, 410, 610)  # for plasma cutoff check
GYRO_TARGET_MHZ   = 2695    # check bump near 2.695 GHz
GYRO_CONT_POINTS  = (1415, 4995)  # continuum anchor points for 2.695

INPUT_FILE_NAME   = "10_2018parsed.csv"  # change if needed

# -------------------------------------------------------------------


def root_dir_from_here() -> Path:
    """Resolve project root as parent of 'scripts' folder if this file is in scripts/."""
    here = Path(__file__).resolve()
    scripts_dir = here.parent
    # project root is parent of scripts/ (assuming structure .../Projects/HSF001_Sun/scripts)
    return scripts_dir.parent


def load_wide_csv(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        sys.exit(f"[ERROR] Input CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    required_prefix = "sfu_"
    sfu_cols = [c for c in df.columns if c.startswith(required_prefix)]
    if not sfu_cols:
        sys.exit("[ERROR] No 'sfu_XXXMHz' columns found in the CSV.")
    return df


def wide_to_long(df: pd.DataFrame) -> pd.DataFrame:
    value_cols = [c for c in df.columns if c.startswith("sfu_")]
    long = df.melt(id_vars=["date", "station"], value_vars=value_cols,
                   var_name="band", value_name="flux_sfu")
    # Extract frequency number in MHz from 'sfu_XXXXMHz'
    long["f_MHz"] = long["band"].str.extract(r"(\d+)\s*MHz").astype(int)
    # Clean
    long = long.replace([np.inf, -np.inf], np.nan).dropna(subset=["flux_sfu"])
    long = long[long["flux_sfu"] > 0]
    return long


def quiet_percentiles(long: pd.DataFrame, q=QUIET_PERCENTILE) -> pd.DataFrame:
    qs = long.groupby("f_MHz")["flux_sfu"].quantile(q).reset_index()
    return qs


def fit_slope_hi_freq(qs: pd.DataFrame, fmin=HI_FREQ_MIN_MHZ):
    hi = qs[qs["f_MHz"] >= fmin].copy()
    if len(hi) < 2:
        return None  # not enough points
    x = np.log10(hi["f_MHz"].to_numpy())
    y = np.log10(hi["flux_sfu"].to_numpy())
    if HAVE_SCIPY:
        slope, intercept, r, p, se = linregress(x, y)
        return dict(slope=slope, intercept=intercept, R=r, stderr=se, n=len(hi))
    # fallback: numpy polyfit
    slope, intercept = np.polyfit(x, y, 1)
    # simple R estimate
    yhat = intercept + slope * x
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    R = 0.0 if ss_tot == 0 else math.sqrt(max(0.0, 1.0 - ss_res / ss_tot))
    return dict(slope=slope, intercept=intercept, R=R, stderr=float("nan"), n=len(hi))


def plasma_ne_bound(qs: pd.DataFrame, triple=LOW_TRIPLE, ratio_thresh=0.7):
    f_low, f1, f2 = triple
    have = set(qs["f_MHz"])
    if not {f_low, f1, f2}.issubset(have):
        return None
    flux_low = qs.loc[qs["f_MHz"] == f_low, "flux_sfu"].iloc[0]
    interp = np.interp(
        f_low,
        [f1, f2],
        qs.set_index("f_MHz").loc[[f1, f2], "flux_sfu"].to_numpy()
    )
    ratio = float(flux_low) / float(interp)
    if ratio < ratio_thresh:
        fp_hz = f_low * 1e6
        ne_cm3 = (fp_hz / 8980.0) ** 2
        return dict(ratio=ratio, fp_MHz=f_low, ne_cm3=ne_cm3)
    return dict(ratio=ratio, fp_MHz=None, ne_cm3=None)


def gyro_excess(qs: pd.DataFrame, f_target=GYRO_TARGET_MHZ, anchors=GYRO_CONT_POINTS):
    have = set(qs["f_MHz"])
    if f_target not in have or not set(anchors).issubset(have):
        return None
    cont = np.interp(
        f_target,
        list(anchors),
        qs.set_index("f_MHz").loc[list(anchors), "flux_sfu"].to_numpy()
    )
    flux = qs.loc[qs["f_MHz"] == f_target, "flux_sfu"].iloc[0]
    excess = float(flux) - float(cont)
    # B (Gauss) from f = s * 27.992493 GHz/T * B
    B_G_s2 = ((f_target * 1e6) / 2) / (27.992493e9) * 1e4
    B_G_s3 = ((f_target * 1e6) / 3) / (27.992493e9) * 1e4
    return dict(excess_sfu=excess, B_G_s2=B_G_s2, B_G_s3=B_G_s3)


def ensure_dirs(root: Path):
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "results" / "figures").mkdir(parents=True, exist_ok=True)


def save_summary(root: Path, qs: pd.DataFrame, slope_info, plasma_info, gyro_info):
    out = root / "results" / "p0_summary.txt"
    with open(out, "w", encoding="utf-8", newline="") as f:
        f.write("Quiet-Sun 10th-percentile spectrum (sfu):\n")
        for _, r0 in qs.sort_values("f_MHz").iterrows():
            f.write(f"  {int(r0.f_MHz):>5} MHz : {r0.flux_sfu:.2f} sfu\n")

        f.write("\nThin-regime slope (>= %d MHz): " % HI_FREQ_MIN_MHZ)
        if slope_info is None:
            f.write("not enough high-frequency points\n")
        else:
            f.write(
                "alpha = %.2f  (R=%.3f, stderr=%s, n=%d)\n" %
                (
                    slope_info["slope"],
                    slope_info["R"],
                    ("%.2f" % slope_info["stderr"]) if np.isfinite(slope_info["stderr"]) else "NA",
                    slope_info["n"],
                )
            )

        if plasma_info is None:
            f.write("\nPlasma cutoff: insufficient bands to evaluate\n")
        else:
            f.write(
                "245 MHz suppression ratio (obs / interp %d–%d): %.2f\n"
                % (LOW_TRIPLE[1], LOW_TRIPLE[2], plasma_info["ratio"])
            )
            if plasma_info["ne_cm3"] is not None:
                f.write(
                    "-> Interpreted cutoff at %d MHz: fp=%d MHz, ne≈%.2e cm^-3\n"
                    % (LOW_TRIPLE[0], plasma_info["fp_MHz"], plasma_info["ne_cm3"])
                )
            else:
                f.write("-> No strong suppression; cannot claim cutoff at 245 MHz from percentiles.\n")

        if gyro_info is None:
            f.write("\n2.695 GHz excess: insufficient bands to evaluate\n")
        else:
            f.write(
                "\n2.695 GHz excess over continuum: %.2f sfu\n"
                % gyro_info["excess_sfu"]
            )
            f.write(
                "Implied B if gyro-resonant: s=2 -> ~%d G, s=3 -> ~%d G\n"
                % (round(gyro_info["B_G_s2"]), round(gyro_info["B_G_s3"]))
            )
    print(f"[OK] Wrote summary: {out}")


def plot_slope_and_spectrum(root: Path, qs: pd.DataFrame, slope_info):
    figdir = root / "results" / "figures"
    # Full quiet-Sun spectrum (all bands) in log–log
    plt.figure(figsize=(6, 4))
    plt.scatter(np.log10(qs["f_MHz"]), np.log10(qs["flux_sfu"]), s=40)
    plt.xlabel("log10 Frequency [MHz]")
    plt.ylabel("log10 Flux [sfu]")
    plt.title("Quiet-Sun 10th-percentile Spectrum (All Bands)")
    plt.grid(True, which="both", linestyle=":")
    plt.tight_layout()
    out1 = figdir / "p0_spectrum.png"
    plt.savefig(out1, dpi=200)
    plt.close()
    print(f"[OK] Saved {out1}")

    # Slope-fit figure (>= HI_FREQ_MIN_MHZ)
    hi = qs[qs["f_MHz"] >= HI_FREQ_MIN_MHZ].copy()
    if len(hi) >= 2 and slope_info is not None:
        logf = np.log10(hi["f_MHz"].values)
        logs = np.log10(hi["flux_sfu"].values)
        fit_y = slope_info["intercept"] + slope_info["slope"] * logf

        plt.figure(figsize=(6, 4))
        plt.scatter(logf, logs, label="Quiet-Sun (10th pct)", s=40)
        plt.plot(logf, fit_y, label=f"Fit: α = {slope_info['slope']:.2f}", linewidth=2)
        plt.xlabel("log10 Frequency [MHz]")
        plt.ylabel("log10 Flux [sfu]")
        plt.title(f"P0 Thin-Regime Slope (≥ {HI_FREQ_MIN_MHZ/1000:.1f} GHz)")
        plt.legend()
        plt.grid(True, which="both", linestyle=":")
        plt.tight_layout()
        out2 = figdir / "p0_slope.png"
        plt.savefig(out2, dpi=200)
        plt.close()
        print(f"[OK] Saved {out2}")
    else:
        print("[WARN] Not enough high-frequency points to plot slope.")


def main():
    root = root_dir_from_here()
    ensure_dirs(root)

    csv = root / "data" / INPUT_FILE_NAME
    df = load_wide_csv(csv)
    long = wide_to_long(df)
    qs = quiet_percentiles(long, q=QUIET_PERCENTILE)

    slope_info = fit_slope_hi_freq(qs, fmin=HI_FREQ_MIN_MHZ)
    plasma_info = plasma_ne_bound(qs, triple=LOW_TRIPLE, ratio_thresh=0.7)
    gyro_info = gyro_excess(qs, f_target=GYRO_TARGET_MHZ, anchors=GYRO_CONT_POINTS)

    save_summary(root, qs, slope_info, plasma_info, gyro_info)
    plot_slope_and_spectrum(root, qs, slope_info)

    print("\n[Done] P0 analysis complete.")


if __name__ == "__main__":
    main()
