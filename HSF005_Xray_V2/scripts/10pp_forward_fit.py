#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Forward-cone fit for pp elastic scattering (e.g., TOTEM 7 TeV) from HEPData CSVs.

Run from the folder that contains this script, for example:
  python 10pp_forward_fit.py --in ..\data\raw --out ..\results --tmax 0.12 --rho 0.14
"""

from __future__ import annotations
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

HBARC2_MB_GEV2 = 0.389379  # mb·GeV^2


# ---------- robust CSV reading ----------

def _try_read_csv(csv_path: Path) -> pd.DataFrame:
    tries = [
        dict(engine="python", sep=None, comment="#", encoding="utf-8", on_bad_lines="skip"),
        dict(engine="python", sep=None, comment="#", encoding="utf-8-sig", on_bad_lines="skip"),
        dict(engine="python", sep=None, comment="#", encoding="cp1252", on_bad_lines="skip"),
        dict(engine="python", sep=",",  comment="#", encoding="utf-8", on_bad_lines="skip"),
        dict(engine="python", sep=";",  comment="#", encoding="utf-8", on_bad_lines="skip"),
        dict(engine="python", sep="\t", comment="#", encoding="utf-8", on_bad_lines="skip"),
    ]
    last_exc = None
    for kw in tries:
        try:
            return pd.read_csv(csv_path, **kw)
        except Exception as e:
            last_exc = e

    try:
        raw = pd.read_csv(csv_path, header=None, engine="python", encoding="utf-8", on_bad_lines="skip")
        if raw.shape[1] == 1:
            s = raw.iloc[:, 0].astype(str)
            if s.str.contains(";").mean() > 0.2:
                parts = s.str.split(";", expand=True)
            elif s.str.contains(",").mean() > 0.2:
                parts = s.str.split(",", expand=True)
            elif s.str.contains("\t").mean() > 0.2:
                parts = s.str.split("\t", expand=True)
            else:
                raise last_exc or ValueError("Single-column CSV with unknown delimiter.")
            header_row = parts.iloc[0]
            if header_row.apply(lambda x: not re.fullmatch(r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?", str(x).strip())).any():
                parts.columns = [str(x).strip() for x in header_row]
                parts = parts.iloc[1:].reset_index(drop=True)
            return parts
    except Exception:
        pass
    raise last_exc or ValueError(f"Could not parse {csv_path.name}")


def _pick_col(cols, pattern):
    for c in cols:
        if re.search(pattern, c, flags=re.IGNORECASE):
            return c
    return None


def read_hepdata_table(csv_path: Path) -> pd.DataFrame:
    df = _try_read_csv(csv_path)
    df.columns = [str(c).strip() for c in df.columns]

    cols = list(df.columns)
    t_col = (_pick_col(cols, r"\babs\s*\(?t\)?|\|t\||\bt\b.*ge?v") or
             _pick_col(cols, r"\bq\^2\b|\bq2\b"))
    y_col = (_pick_col(cols, r"d\s*\(?sig\w*\)?\s*/\s*dt|dsig\s*/\s*dt|d\s*sigma\s*/\s*dt") or
             _pick_col(cols, r"\bd.?sd.?t\b|\by\b|value"))
    if t_col is None:
        t_col = _pick_col(cols, r"independent.*value")
    if y_col is None:
        y_col = _pick_col(cols, r"dependent.*value")
    if t_col is None or y_col is None:
        raise ValueError(f"Could not find |t|/dσ/dt in {csv_path.name}. Columns: {cols[:6]}...")

    t = pd.to_numeric(df[t_col], errors="coerce").values
    y = pd.to_numeric(df[y_col], errors="coerce").values

    estat_col = _pick_col(cols, r"stat.*(err|unc)|unc.*stat")
    esys_col  = _pick_col(cols, r"sys.*(err|unc)|unc.*sys")
    estat = pd.to_numeric(df.get(estat_col), errors="coerce").values if estat_col else np.full_like(y, np.nan, float)
    esys  = pd.to_numeric(df.get(esys_col),  errors="coerce").values if esys_col  else np.full_like(y, np.nan, float)
    err = np.sqrt(np.where(np.isfinite(estat), estat, 0.0)**2 + np.where(np.isfinite(esys), esys, 0.0)**2)
    bad = ~np.isfinite(err) | (err <= 0)
    err[bad] = np.abs(y[bad]) * 0.01

    out = pd.DataFrame({"t": t, "dsdt": y, "err": err})
    out = out.replace([np.inf, -np.inf], np.nan).dropna()
    return out


# ---------- models ----------

def model_baselineA(t, A, B):
    """Single exponential (textbook forward cone)."""
    return A * np.exp(-B * np.abs(t))

def model_baselineB(t, sigma_tot, B, rho):
    """Optical-theorem constrained normalization (fits sigma_tot, B; rho fixed)."""
    A_from_sigma = ( (1.0 + rho**2) / (16.0 * np.pi * HBARC2_MB_GEV2) ) * (sigma_tot**2)
    return A_from_sigma * np.exp(-B * np.abs(t))

def model_hsf(t, A, B, C):
    """Example 'HSF' extension: curvature in exponent."""
    tt = np.abs(t)
    return A * np.exp(-B * tt - C * tt*tt)


# ---------- fitting & metrics ----------

def _mask_clean(t, y, e):
    m = np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (y > 0) & (e > 0)
    return t[m], y[m], e[m]

def _chi2(y, f, e):
    r = (y - f) / e
    return float(np.sum(r*r))

def _rms(y, f):
    return float(np.sqrt(np.mean((y - f)**2)))

def _aic_bic(chi2, n, k):
    # For Gaussian errors (up to additive constants that cancel across models)
    aic = chi2 + 2*k
    bic = chi2 + k*np.log(n)
    return float(aic), float(bic)

def fit_baselineA(t, y, e):
    tt, yy, ee = _mask_clean(t, y, e)
    if len(tt) < 6: raise ValueError("Not enough points.")
    if HAVE_SCIPY:
        popt, pcov = curve_fit(model_baselineA, tt, yy, sigma=ee, absolute_sigma=True, p0=(yy.max(), 20.0), maxfev=10000)
        A, B = popt
        return dict(name="BaselineA", params=dict(A=A, B=B), cov=pcov)
    # weighted log-linear fallback
    w = 1.0 / np.maximum((ee / np.maximum(yy, 1e-300))**2, 1e-30)
    X = np.vstack([np.ones_like(tt), -np.abs(tt)]).T
    WX = X * np.sqrt(w[:, None])
    Wy = np.log(np.maximum(yy, 1e-300)) * np.sqrt(w)
    beta, *_ = np.linalg.lstsq(WX, Wy, rcond=None)
    lnA, B = beta[0], beta[1]
    A = float(np.exp(lnA))
    cov = np.linalg.inv(WX.T @ WX)
    covAB = np.zeros((2,2)); covAB[0,0]=(A**2)*cov[0,0]; covAB[1,1]=cov[1,1]
    return dict(name="BaselineA", params=dict(A=A, B=float(B)), cov=covAB)

def fit_baselineB(t, y, e, rho):
    if not HAVE_SCIPY:
        raise RuntimeError("BaselineB requires SciPy (nonlinear fit for sigma_tot).")
    tt, yy, ee = _mask_clean(t, y, e)
    popt, pcov = curve_fit(lambda T, sig, B: model_baselineB(T, sig, B, rho),
                           tt, yy, sigma=ee, absolute_sigma=True,
                           p0=(100.0, 20.0), bounds=(0, np.inf), maxfev=20000)
    sigma_tot, B = popt
    return dict(name="BaselineB", params=dict(sigma_tot=float(sigma_tot), B=float(B)), cov=pcov)

def fit_hsf(t, y, e):
    if not HAVE_SCIPY:
        raise RuntimeError("HSF fit requires SciPy.")
    tt, yy, ee = _mask_clean(t, y, e)
    popt, pcov = curve_fit(model_hsf, tt, yy, sigma=ee, absolute_sigma=True,
                           p0=(yy.max(), 20.0, 1.0), bounds=(0, np.inf), maxfev=20000)
    A, B, C = popt
    return dict(name="HSF", params=dict(A=float(A), B=float(B), C=float(C)), cov=pcov)

def evaluate_model(name, params, t, y, e, rho):
    if name == "BaselineA":
        f = model_baselineA(t, params["A"], params["B"])
        k = 2
    elif name == "BaselineB":
        f = model_baselineB(t, params["sigma_tot"], params["B"], rho)
        k = 2
    elif name == "HSF":
        f = model_hsf(t, params["A"], params["B"], params["C"])
        k = 3
    else:
        raise ValueError(name)
    chi2 = _chi2(y, f, e)
    rms  = _rms(y, f)
    aic, bic = _aic_bic(chi2, len(t), k)
    return dict(chi2=chi2, rms=rms, aic=aic, bic=bic, k=k, yhat=f)


# ---------- original helpers kept for continuity ----------

def optical_sigma_tot(dsdt0, rho):
    return np.sqrt(dsdt0 * 16.0 * np.pi * HBARC2_MB_GEV2 / (1.0 + rho**2))

def sigma_elastic_single_exp(A, B):
    return A / B


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="Forward-cone fit on HEPData CSVs.")
    ap.add_argument("--in", dest="indir", type=str, required=True, help="Folder with CSVs (searched recursively)")
    ap.add_argument("--out", dest="outdir", type=str, required=True, help="Output folder")
    ap.add_argument("--tmax", type=float, default=0.12, help="Use points with |t| <= tmax (GeV^2)")
    ap.add_argument("--rho", type=float, default=0.14, help="rho = Re/Im at t->0")
    ap.add_argument("--title", type=str, default="pp elastic (forward cone)", help="Plot title")
    ap.add_argument("--hold", type=float, default=0.08, help="Hold-out split: fit on |t|<=hold, test on hold<|t|<=tmax")
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    indir = Path(args.indir);  indir = (script_dir / indir).resolve() if not indir.is_absolute() else indir
    outdir = Path(args.outdir); outdir = (script_dir / outdir).resolve() if not outdir.is_absolute() else outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # Find CSVs
    csvs = [p for p in indir.rglob("*.csv") if p.name.lower().startswith("table1")]
    if not csvs:
        raise SystemExit(f"No CSVs found under: {indir}")
    frames = []
    for p in sorted(csvs):
        try:
            frames.append(read_hepdata_table(p))
        except Exception as exc:
            print(f"[WARN] Skipping {p.name}: {exc}")
    if not frames:
        raise SystemExit("No usable tables after parsing.")
    data = pd.concat(frames, ignore_index=True).replace([np.inf, -np.inf], np.nan).dropna()
    use  = data[np.abs(data["t"]) <= args.tmax].copy().sort_values("t")
    if len(use) < 8:
        raise SystemExit(f"Not enough points for |t| <= {args.tmax} (have {len(use)}).")

    # Split fit/test
    fit_mask  = (np.abs(use["t"].values) <= args.hold)
    test_mask = (np.abs(use["t"].values) >  args.hold)
    tf, yf, ef = use["t"].values[fit_mask],  use["dsdt"].values[fit_mask],  use["err"].values[fit_mask]
    tt, yt, et = use["t"].values[test_mask], use["dsdt"].values[test_mask], use["err"].values[test_mask]

    # ----- Fit models -----
    results = []

    # Baseline A (single exp — your original)
    fitA = fit_baselineA(tf, yf, ef)
    evalA_fit  = evaluate_model(fitA["name"], fitA["params"], tf, yf, ef, args.rho)
    evalA_test = evaluate_model(fitA["name"], fitA["params"], tt, yt, et, args.rho) if len(tt)>0 else None
    results.append(("BaselineA", fitA, evalA_fit, evalA_test))

    # Baseline B (optical-theorem constrained)
    try:
        fitB = fit_baselineB(tf, yf, ef, args.rho)
        evalB_fit  = evaluate_model(fitB["name"], fitB["params"], tf, yf, ef, args.rho)
        evalB_test = evaluate_model(fitB["name"], fitB["params"], tt, yt, et, args.rho) if len(tt)>0 else None
        results.append(("BaselineB", fitB, evalB_fit, evalB_test))
    except Exception as e:
        print("[WARN] BaselineB skipped:", e)

    # HSF extension (curved exponent). Replace with your own HSF model if desired.
    try:
        fitH = fit_hsf(tf, yf, ef)
        evalH_fit  = evaluate_model(fitH["name"], fitH["params"], tf, yf, ef, args.rho)
        evalH_test = evaluate_model(fitH["name"], fitH["params"], tt, yt, et, args.rho) if len(tt)>0 else None
        results.append(("HSF", fitH, evalH_fit, evalH_test))
    except Exception as e:
        print("[WARN] HSF fit skipped:", e)

    # ----- Original single-exp summary (keep your legacy outputs) -----
    A = fitA["params"]["A"]; B = fitA["params"]["B"]
    dsdt0 = A
    sigma_tot = optical_sigma_tot(dsdt0, args.rho)
    sigma_el  = sigma_elastic_single_exp(A, B)
    A_err = float(np.sqrt(fitA["cov"][0,0])) if (fitA["cov"] is not None and fitA["cov"].shape==(2,2) and np.isfinite(fitA["cov"][0,0])) else float("nan")
    B_err = float(np.sqrt(fitA["cov"][1,1])) if (fitA["cov"] is not None and fitA["cov"].shape==(2,2) and np.isfinite(fitA["cov"][1,1])) else float("nan")

    # Legacy plot (baseline A only)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.errorbar(use["t"], use["dsdt"], yerr=use["err"], fmt="o", ms=3, lw=1, label="HEPData (forward)")
    tt_lin = np.linspace(0, args.tmax, 300)
    ax.plot(tt_lin, model_baselineA(tt_lin, A, B), label=r"Fit $A\,e^{-B|t|}$")
    ax.set_xlabel(r"$|t|\;[{\rm GeV}^2]$")
    ax.set_ylabel(r"$d\sigma/dt\;[{\rm mb/GeV}^2]$")
    ax.set_yscale("log")
    ax.set_title(args.title)
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend()
    plot_path = outdir / "pp_forward_fit.png"
    fig.tight_layout(); fig.savefig(plot_path, dpi=200); plt.close(fig)

    res = pd.DataFrame([{
        "t_max_GeV2": args.tmax,
        "rho": args.rho,
        "A_dsigdt_t0_mb_per_GeV2": dsdt0,
        "A_err": A_err,
        "B_GeV_minus2": B,
        "B_err": B_err,
        "sigma_tot_mb": sigma_tot,
        "sigma_el_mb": sigma_el,
        "n_points_used": int(len(use))
    }])
    res_path = outdir / "pp_results.csv"
    res.to_csv(res_path, index=False)

    md = outdir / "pp_result.md"
    md.write_text(
f"""# pp elastic: forward-cone fit

**Model**: $d\\sigma/dt = A\\,e^{{-B|t|}}$, fitted for $|t| \\le {args.tmax:.3f}\\,\\mathrm{{GeV}}^2$.

- $A=(d\\sigma/dt)_{{t=0}} = {dsdt0:.3g}\\;\\mathrm{{mb/GeV^2}}$  {"(±{:.2g})".format(A_err) if np.isfinite(A_err) else ""}  
- $B = {B:.3g}\\;\\mathrm{{GeV^{-2}}}$  {"(±{:.2g})".format(B_err) if np.isfinite(B_err) else ""}  
- $\\rho = {args.rho}$ (assumed)

Optical theorem:
\\[
\\sigma_\\mathrm{{tot}} = \\sqrt{{\\frac{{16\\pi\\,0.389379\\;A}}{{1+\\rho^2}}}} = {sigma_tot:.2f}\\;\\mathrm{{mb}}.
\\]

Single-exponential elastic estimate:
\\[
\\sigma_\\mathrm{{el}} \\simeq \\frac{{A}}{{B}} = {sigma_el:.2f}\\;\\mathrm{{mb}}.
\\]

Points used: {len(use)}

Artifacts:
- Plot: `{plot_path.name}`
- Table: `{res_path.name}`
""",
        encoding="utf-8"
    )

    # ----- New: multi-model comparison outputs -----
    # Comparison CSV
    rows = []
    for name, fit, ev_fit, ev_test in results:
        row = dict(model=name,
                   k=ev_fit["k"],
                   t_fit_max=args.hold,
                   t_test_max=args.tmax,
                   chi2_fit=ev_fit["chi2"],
                   rms_fit=ev_fit["rms"],
                   aic_fit=ev_fit["aic"],
                   bic_fit=ev_fit["bic"],
                   n_fit=len(tf))
        if ev_test is not None and len(tt)>0:
            row.update(dict(chi2_test=ev_test["chi2"], rms_test=ev_test["rms"], n_test=len(tt)))
        else:
            row.update(dict(chi2_test=np.nan, rms_test=np.nan, n_test=0))
        # store key parameters for transparency
        row.update({f"param_{k}": v for k, v in fit["params"].items()})
        rows.append(row)
    cmp_df = pd.DataFrame(rows)
    cmp_path = outdir / "pp_model_compare.csv"
    cmp_df.to_csv(cmp_path, index=False)
# --- Residuals plot (fit & hold-out) ---
# Helper to get yhat for any model on given t-array
    def _yhat_for(model_name, fitdict, t_array):
        if model_name == "BaselineA":
            return model_baselineA(t_array, fitdict["params"]["A"], fitdict["params"]["B"])
        elif model_name == "BaselineB":
            return model_baselineB(t_array, fitdict["params"]["sigma_tot"], fitdict["params"]["B"], args.rho)
        elif model_name == "HSF":
            return model_hsf(t_array, fitdict["params"]["A"], fitdict["params"]["B"], fitdict["params"]["C"])
        else:
            raise ValueError(model_name)

    fig, ax = plt.subplots(figsize=(7.6, 4.4))

    # Baseline A residuals (fit)
    yhatA_fit = _yhat_for("BaselineA", fitA, tf)
    ax.errorbar(tf, yf - yhatA_fit, yerr=ef, fmt=".", label="Baseline A residual (fit)")

    # Baseline A residuals (hold-out), if present
    if len(tt) > 0:
        yhatA_test = _yhat_for("BaselineA", fitA, tt)
        ax.errorbar(tt, yt - yhatA_test, yerr=et, fmt=".", label="Baseline A residual (test)")

    # HSF residuals (fit), if we have an HSF fit
    have_hsf = any(n == "HSF" for (n, _, _, _) in results)
    if have_hsf:
        # reuse fitH from earlier try-block
        yhatH_fit = _yhat_for("HSF", fitH, tf)
        ax.plot(tf, yf - yhatH_fit, ".", label="HSF residual (fit)")
        if len(tt) > 0:
            yhatH_test = _yhat_for("HSF", fitH, tt)
            ax.plot(tt, yt - yhatH_test, ".", label="HSF residual (test)")

    # Optional: Baseline B residuals if available
    have_b = any(n == "BaselineB" for (n, _, _, _) in results)
    if have_b:
        yhatB_fit = _yhat_for("BaselineB", fitB, tf)
        ax.plot(tf, yf - yhatB_fit, ".", label="Baseline B residual (fit)")
        if len(tt) > 0:
            yhatB_test = _yhat_for("BaselineB", fitB, tt)
            ax.plot(tt, yt - yhatB_test, ".", label="Baseline B residual (test)")

    ax.axhline(0.0, ls=":", lw=1)
    ax.set_xlabel(r"$|t|\;[{\rm GeV}^2]$")
    ax.set_ylabel(r"Residual  $y - \hat y$  [mb/GeV$^2$]")
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend()
    plt.tight_layout()
    plt.savefig(outdir / "pp_residuals_v2.png", dpi=200)
    plt.close(fig)

    # Multi-model plot
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ax.errorbar(use["t"], use["dsdt"], yerr=use["err"], fmt="o", ms=3, lw=1, label="HEPData (forward)")
    tt_plot = np.linspace(0, args.tmax, 600)
    # Visualize hold-out region (test window) as a light band
    if args.hold < args.tmax:
        ax.axvspan(args.hold, args.tmax, alpha=0.1, zorder=0, label="hold-out")

    # Recompute yhat on full window for plotting
    for name, fit, _, _ in results:
        if name == "BaselineA":
            yhat = model_baselineA(tt_plot, fit["params"]["A"], fit["params"]["B"])
            lab = "Baseline A: A·e^{-B|t|}"
        elif name == "BaselineB":
            yhat = model_baselineB(tt_plot, fit["params"]["sigma_tot"], fit["params"]["B"], args.rho)
            lab = "Baseline B: σ_tot,B (optical)"
        else:
            yhat = model_hsf(tt_plot, fit["params"]["A"], fit["params"]["B"], fit["params"]["C"])
            lab = "HSF (curved exp)"
        ax.plot(tt_plot, yhat, lw=1.6, label=lab)
    ax.set_xlabel(r"$|t|\;[{\rm GeV}^2]$")
    ax.set_ylabel(r"$d\sigma/dt\;[{\rm mb/GeV}^2]$")
    ax.set_yscale("log")
    ax.set_title(f"{args.title}  (fit: |t|≤{args.hold:g},  test: {args.hold:g}–{args.tmax:g})")
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend()
    plot_path2 = outdir / "pp_forward_fit_v2.png"
    fig.tight_layout(); fig.savefig(plot_path2, dpi=200); plt.close(fig)

    print("[OK] Wrote:")
    print(" -", plot_path)
    print(" -", res_path)
    print(" -", md)
    print(" -", plot_path2)
    print(" -", cmp_path)


if __name__ == "__main__":
    main()
