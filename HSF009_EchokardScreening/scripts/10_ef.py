#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, argparse, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def area_length_volume_ml(area_cm2, length_cm):
    """Single-plane area–length volume (mL): V = 8/(3π) * A^2 / L"""
    return (8.0/(3.0*math.pi)) * (area_cm2**2) / np.maximum(length_cm, 1e-9)

def load_series(csv_path):
    df = pd.read_csv(csv_path)
    need = {"time_s","area_cm2","length_cm"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    df = df.sort_values("time_s").reset_index(drop=True)
    df = df.dropna(subset=["time_s","area_cm2","length_cm"])
    return df

def subsample_to_fps(df, fps):
    """Nearest-neighbour decimation to a target FPS (no fitting)."""
    tmin, tmax = df["time_s"].iloc[0], df["time_s"].iloc[-1]
    if tmax <= tmin:
        return df.copy().reset_index(drop=True)
    tt = np.arange(tmin, tmax, 1.0/fps)
    idx = np.searchsorted(df["time_s"].values, tt, side="left")
    idx = np.clip(idx, 0, len(df)-1)
    return df.iloc[idx].copy().reset_index(drop=True)

def compute_ef(df):
    V = area_length_volume_ml(df["area_cm2"].values, df["length_cm"].values)
    EDV, ESV = float(V.max()), float(V.min())
    EF = 100.0*(EDV-ESV)/EDV if EDV > 0 else float("nan")
    out = df.copy()
    out["volume_ml"] = V
    return EF, EDV, ESV, out

def main():
    ap = argparse.ArgumentParser(description="Compute EF from real CSV (area_cm2 & length_cm vs time)")
    ap.add_argument("--in-csv", required=True, help="input CSV: time_s,area_cm2,length_cm[,series_id]")
    ap.add_argument("--out", default="../results", help="output folder")
    ap.add_argument("--fps-list", default="60,30,20,15,10,8,5", help="comma-separated FPS to test")
    args = ap.parse_args()
    if not os.path.isfile(args.in_csv):
        raise FileNotFoundError(f"--in-csv must be a CSV file, got: {args.in_csv}")
    os.makedirs(args.out, exist_ok=True)
    base = os.path.splitext(os.path.basename(args.in_csv))[0]

    df = load_series(args.in_csv)

    # Reference (native sampling)
    EF_ref, EDV_ref, ESV_ref, df_ref = compute_ef(df)

    rows = []
    for tok in args.fps_list.split(","):
        tok = tok.strip()
        if not tok:
            continue
        fps = int(float(tok))
        dff = subsample_to_fps(df, fps)
        EF, EDV, ESV, _ = compute_ef(dff)
        rows.append(dict(fps=fps, EF=EF, EDV=EDV, ESV=ESV, EF_err=EF - EF_ref))

    df_fps = pd.DataFrame(rows).sort_values("fps")
    df_fps.to_csv(os.path.join(args.out, f"{base}_ef_vs_fps.csv"), index=False)

    # Volume time series (native)
    plt.figure(figsize=(8,4))
    plt.plot(df_ref["time_s"], df_ref["volume_ml"], lw=2)
    try:
        t_ED = df_ref["time_s"].iloc[df_ref["volume_ml"].idxmax()]
        t_ES = df_ref["time_s"].iloc[df_ref["volume_ml"].idxmin()]
        plt.scatter([t_ED, t_ES], [df_ref["volume_ml"].max(), df_ref["volume_ml"].min()])
    except Exception:
        pass
    plt.xlabel("Time (s)"); plt.ylabel("LV Volume (mL)")
    plt.title(f"Volume time series — EF≈{EF_ref:.1f}% (EDV {EDV_ref:.0f} mL, ESV {ESV_ref:.0f} mL)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.out, f"{base}_volume_timeseries.png"), dpi=160)
    plt.close()

    # EF vs FPS plot
    if not df_fps.empty:
        plt.figure(figsize=(7,3.6))
        plt.plot(df_fps["fps"], df_fps["EF"], marker="o")
        plt.axhline(EF_ref, ls="--", lw=1)
        plt.xlabel("Frame rate (fps)"); plt.ylabel("EF (%)")
        plt.title("EF vs frame rate (no fit)")
        plt.tight_layout()
        plt.savefig(os.path.join(args.out, f"{base}_ef_vs_fps.png"), dpi=160)
        plt.close()

    # Text summary
    with open(os.path.join(args.out, f"{base}_summary.txt"), "w", encoding="utf-8") as f:
        f.write(f"Source CSV: {args.in_csv}\n")
        f.write(f"Reference EF (native sampling): {EF_ref:.2f}%\n")
        f.write(f"EDV={EDV_ref:.1f} mL, ESV={ESV_ref:.1f} mL\n\n")
        if not df_fps.empty:
            f.write(df_fps.to_string(index=False))
            f.write("\n")

if __name__ == "__main__":
    main()
