#!/usr/bin/env python3
"""
RSTN (USAF/NOAA) noontime daily flux parser.

Input format (one line per day, variable whitespace):
YYMMDD STATION v245 v410 v610 v1415 v2695 v4995 v8800 v15400

Example:
180101 PALE  13 25 34 46 67 107 225 518

Returns a DataFrame with:
- date (datetime.date)
- station (str)
- sfu_245MHz ... sfu_15400MHz (float, NaN if missing)
"""

from __future__ import annotations
import re
from datetime import datetime
from pathlib import Path
from typing import Iterable, Union

import numpy as np
import pandas as pd

RSTN_FREQS = [245, 410, 610, 1415, 2695, 4995, 8800, 15400]
RSTN_COLS  = [f"sfu_{f}MHz" for f in RSTN_FREQS]

def _parse_line(line: str):
    t = re.split(r"\s+", line.strip())
    if len(t) < 2: 
        return None
    yymmdd, station = t[0], t[1]
    # parse numbers that follow; fewer than 8 allowed → pad with NaN
    vals = []
    for tok in t[2:]:
        try:
            vals.append(float(tok))
        except Exception:
            # tolerate blanks/non-numerics
            continue
    if len(vals) < 8:
        vals += [np.nan] * (8 - len(vals))
    else:
        vals = vals[:8]
    # date conversion
    try:
        dt = datetime.strptime(yymmdd, "%y%m%d").date()
        if dt.year < 1950:  # century roll if needed
            dt = dt.replace(year=dt.year + 100)
    except Exception:
        return None
    row = {"date": dt, "station": station}
    row.update({c: v for c, v in zip(RSTN_COLS, vals)})
    return row

def load_year(path: Union[str, Path]) -> pd.DataFrame:
    p = Path(path)
    rows = []
    with p.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip():
                rec = _parse_line(line)
                if rec:
                    rows.append(rec)
    df = pd.DataFrame(rows).sort_values("date").reset_index(drop=True)
    return df

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("file", help="RSTN yearly text file (e.g., pale_noontime-flux_2018.txt)")
    ap.add_argument("-o", "--out", help="Optional CSV path to write parsed table")
    args = ap.parse_args()
    df = load_year(args.file)
    print(df.head(10).to_string(index=False))
    if args.out:
        df.to_csv(args.out, index=False)
