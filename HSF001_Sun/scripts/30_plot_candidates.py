#!/usr/bin/env python3
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from pathlib import Path

FREQS = [245,410,610,1415,2695,4995,8800,15400]
COLS  = [f"sfu_{f}MHz" for f in FREQS]
GHz    = np.array([2695,4995,8800,15400], float)

def main(parsed_csv, cand_csv, outdir="results/spectra", topn=8):
    df   = pd.read_csv(parsed_csv, parse_dates=["date"])
    cand = pd.read_csv(cand_csv, parse_dates=["date"])
    need = {"max_abs_resid","which_band","alpha","r2",
            "resid_2695","resid_4995","resid_8800","resid_15400"}
    if not need.issubset(set(cand.columns)):
        raise SystemExit("cand_csv missing required columns; run fit_spectrum.py first.")
    Path(outdir).mkdir(parents=True, exist_ok=True)

    top = cand.sort_values("max_abs_resid", ascending=False).head(topn)
    for _, r in top.iterrows():
        d = r["date"].date()
        row = df[df["date"]==pd.Timestamp(d)]
        if row.empty: 
            continue
        y = row.iloc[0][COLS].astype(float).to_numpy()
        nu = np.array(FREQS, float)

        mask = np.isin(nu, GHz) & np.isfinite(y) & (y>0)
        X, Y = np.log10(nu[mask]), np.log10(y[mask])
        m, b = np.polyfit(X, Y, 1)
        model = 10**(m*np.log10(nu) + b)

        plt.figure(figsize=(7,4))
        plt.loglog(nu, y, "o-", label="noon flux")
        plt.loglog(nu, model, "--", label=f"GHz fit α≈{m:.2f}")
        plt.title(f"Palehua {d} — max residual @ {int(r['which_band'])} MHz")
        plt.xlabel("Frequency (MHz)"); plt.ylabel("Flux (sfu)")
        plt.legend(); plt.tight_layout()
        fn = Path(outdir)/f"spectrum_{d}.png"
        plt.savefig(fn, dpi=140); plt.close()
        print("wrote", fn)

if __name__=="__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("parsed_csv")
    ap.add_argument("cand_csv")
    ap.add_argument("--outdir", default="../results/spectra")
    ap.add_argument("--topn", type=int, default=8)
    args = ap.parse_args()
    main(args.parsed_csv, args.cand_csv, args.outdir, args.topn)
