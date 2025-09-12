#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import math

def fmt_date(ts):
    # ts is Timestamp; return 'YYYY-MM-DD'
    try:
        return ts.date().isoformat()
    except Exception:
        return str(ts)  # fallback

def fmt_num(x, nd=2):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return ""
    return f"{x:.{nd}f}"

def main(cand_csv, out_md="../results/report_assets/README_candidates_2018.md", topn=8):
    cand = pd.read_csv(cand_csv, parse_dates=["date"])
    top  = cand.sort_values("max_abs_resid", ascending=False).head(topn)

    lines = []
    lines += [
      "## HSF P0 — Palehua 2018 (Top spectral-curvature days)\n",
      "| Date | Max-resid band (MHz) | α(GHz) | R² | Resid@2.7 | 5.0 | 8.8 | 15.4 | B[G] s=2 | B[G] s=3 |",
      "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for _, r in top.iterrows():
        d   = fmt_date(r["date"])
        wb  = "" if pd.isna(r["which_band"]) else str(int(r["which_band"]))
        row = (
          f"| {d} | {wb} | {fmt_num(r.get('alpha'))} | {fmt_num(r.get('r2'))} | "
          f"{fmt_num(r.get('resid_2695'))} | {fmt_num(r.get('resid_4995'))} | "
          f"{fmt_num(r.get('resid_8800'))} | {fmt_num(r.get('resid_15400'))} | "
          f"{fmt_num(r.get('B_G_s2'),0)} | {fmt_num(r.get('B_G_s3'),0)} |"
        )
        lines.append(row)

    Path(out_md).parent.mkdir(parents=True, exist_ok=True)
    Path(out_md).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote", out_md)

if __name__=="__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("cand_csv")
    ap.add_argument("--out_md", default="../results/report_assets/README_candidates_2018.md")
    ap.add_argument("--topn", type=int, default=8)
    args = ap.parse_args()
    main(args.cand_csv, args.out_md, args.topn)
