#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HSF008 · Heart EF Screener (no data, no fit)

Berechnet EF aus einem synthetischen A4C-Zyklus mit der Area–Length-Formel
und erzeugt:
  (1) EF vs Frame-Rate (Subsampling, ohne Fit)
  (2) EF-Fehler vs Segmentierungs-Area-Bias (uniform)
  (3) optional: EF-Shift bei ED/ES-ASYMMETRISCHEN Biases (rein analytisch)
  (4) optional: EF-Shift bei frameweisem Area-Noise (ohne Fit)

Alles ist deterministisch steuerbar über CLI-Parameter.
"""

import os, re, math, json, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from textwrap import dedent

# ---------- kleine Helfer ----------
def area_length_volume_ml(area_cm2, length_cm):
    """Single-plane area–length: V = 8/(3π) * A^2 / L  (mL)."""
    return (8.0/(3.0*math.pi)) * (area_cm2**2) / np.maximum(length_cm, 1e-9)

def parse_num_list(s: str):
    """Robust: extrahiert nur Zahlen (float) aus einem String (auch negative).
       Dadurch sind angeklebte Tokens wie '--seed 42' unschädlich."""
    toks = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', str(s))
    return [float(x) for x in toks]

def parse_int_list(s: str):
    return [int(round(float(x))) for x in parse_num_list(s)]

def build_cycle(T, fps_ref, a_ed, b_ed, amp_a, amp_b, phase):
    """Erzeuge einen Beat mit glatter Systole/Diastole (cos-Bump)."""
    t = np.linspace(0.0, T, int(T*fps_ref), endpoint=False)
    a_t = a_ed - amp_a * (0.5*(1 - np.cos(2*math.pi*(t/T + phase))))
    b_t = b_ed - amp_b * (0.5*(1 - np.cos(2*math.pi*(t/T + phase))))
    L_t = 2.0 * a_t
    A_t = math.pi * a_t * b_t
    V_t = area_length_volume_ml(A_t, L_t)
    return t, A_t, L_t, V_t

def subsample_V(V_ref, T, fps_ref, fps):
    """Subsampling ohne Fit: picke Max/Min aus der subsample-Kurve."""
    tt  = np.linspace(0.0, T, int(T*fps), endpoint=False)
    idx = (tt * fps_ref).astype(int) % len(V_ref)
    Vn  = V_ref[idx]
    EDV = float(Vn.max()); ESV = float(Vn.min())
    EF  = 100.0*(EDV-ESV)/EDV
    return EF, EDV, ESV

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser()
    # Geometrie (cm)
    ap.add_argument("--a-ed",   type=float, default=6.0, help="Halbe Langachse a in ED (cm)")
    ap.add_argument("--b-ed",   type=float, default=3.0, help="Halbe Kurzachse b in ED (cm)")
    ap.add_argument("--amp-a",  type=float, default=1.0, help="Systolische Reduktion a (cm)")
    ap.add_argument("--amp-b",  type=float, default=0.8, help="Systolische Reduktion b (cm)")
    ap.add_argument("--phase",  type=float, default=0.05, help="Phasenverschiebung (Zyklen, 0..1)")

    # Timing
    ap.add_argument("--hr",       type=float, default=60.0, help="Herzfrequenz (bpm)")
    ap.add_argument("--fps-ref",  type=int,   default=1000, help="Referenz-FPS (dichte Abtastung)")
    ap.add_argument("--fps-list", default="1000,60,30,20,15,10,8,5",
                    help="Kommagetrennte FPS-Liste (String)")

    # Uniformer Area-Bias (beide Zeitpunkte gleich skaliert) -> EF kürzt sich i. A. weg
    ap.add_argument("--area-bias-list",
                    default="-0.15,-0.10,-0.05,-0.02,0,0.02,0.05,0.10,0.15",
                    help="Kommaliste relativer Flächen-Bias (z.B. -0.1,0,0.1)")

    # Asymmetrische ED/ES-Biases (rein analytische EF-Verschiebung, kein Fit)
    ap.add_argument("--ed-area-bias", type=float, default=0.0, help="rel. Flächen-Bias nur in ED (z.B. 0.05)")
    ap.add_argument("--es-area-bias", type=float, default=0.0, help="rel. Flächen-Bias nur in ES")
    ap.add_argument("--ed-l-bias",    type=float, default=0.0, help="rel. Längen-Bias nur in ED")
    ap.add_argument("--es-l-bias",    type=float, default=0.0, help="rel. Längen-Bias nur in ES")

    # Frameweises Rauschen (optional, no-fit)
    ap.add_argument("--area-noise-sd", type=float, default=0.0,
                    help="Std der relativen Frame-zu-Frame-Flächenstörung (z.B. 0.02)")

    # Sonstiges
    ap.add_argument("--seed", type=int, default=123, help="Seed für deterministische Noise-Tests")
    ap.add_argument("--out",  required=True, help="Output-Ordner")
    ap.add_argument("--config", help="Optional: JSON-Config; CLI überschreibt Einträge")

    args = ap.parse_args()

    # Optional: Config einlesen (CLI hat Vorrang)
    if args.config and os.path.isfile(args.config):
        with open(args.config, "r", encoding="utf-8") as f:
            cfg = json.load(f)
        for k,v in cfg.items():
            if hasattr(args, k):
                # nur übernehmen, wenn Flag noch Default ist
                if str(getattr(args, k)) == str(ap.get_default(k)):
                    setattr(args, k, v)

    os.makedirs(args.out, exist_ok=True)

    # Abgeleitete Größen
    T        = 60.0/args.hr
    fps_list = parse_int_list(args.fps_list)
    area_uni = parse_num_list(args.area_bias_list)

    # 1) Referenz-Beat
    t, A_t, L_t, V_t = build_cycle(T, args.fps_ref, args.a_ed, args.b_ed, args.amp_a, args.amp_b, args.phase)
    ED_idx = int(np.argmax(V_t)); ES_idx = int(np.argmin(V_t))
    EDV    = float(V_t[ED_idx]); ESV    = float(V_t[ES_idx])
    EF_ref = 100.0*(EDV-ESV)/EDV

    # 2) EF vs Frame-Rate (ohne Fit)
    rows_fps = []
    for fps in fps_list:
        EF, ed, es = subsample_V(V_t, T, args.fps_ref, fps)
        rows_fps.append(dict(fps=fps, EF_est=EF, EDV=ed, ESV=es, EF_err=EF - EF_ref))
    df_fps = pd.DataFrame(rows_fps)

    # 3) EF-Fehler vs uniformer Area-Bias (beide Zeitpunkte gleich skaliert)
    rows_bias = []
    for e in area_uni:
        A_b = (1.0 + e) * A_t
        V_b = area_length_volume_ml(A_b, L_t)
        ed  = float(V_b.max()); es = float(V_b.min())
        ef  = 100.0*(ed-es)/ed
        rows_bias.append(dict(area_bias_pct=100*e, EF_est=ef, EDV=ed, ESV=es, EF_err=ef - EF_ref))
    df_bias = pd.DataFrame(rows_bias).sort_values("area_bias_pct")

    # 4) Asymmetrische ED/ES-Biases (rein analytisch, no-fit)
    EDV_biased = EDV * (1 + args.ed_area_bias)**2 / (1 + args.ed_l_bias)
    ESV_biased = ESV * (1 + args.es_area_bias)**2 / (1 + args.es_l_bias)
    EF_biased  = 100.0 * (EDV_biased - ESV_biased) / EDV_biased
    EF_shift   = EF_biased - EF_ref

    # 5) Option: frameweises Area-Noise (no-fit)
    noise_block = ""
    if args.area_noise_sd > 0:
        rng     = np.random.default_rng(args.seed)
        A_noise = A_t * (1.0 + rng.normal(0.0, args.area_noise_sd, size=A_t.shape))
        V_noise = area_length_volume_ml(A_noise, L_t)
        EDV_n, ESV_n = float(V_noise.max()), float(V_noise.min())
        EF_n   = 100.0 * (EDV_n - ESV_n) / EDV_n
        noise_block = f"\n### Frame-Noise (no-fit)\nσ_area={args.area_noise_sd:.1%} → EF={EF_n:.2f}% (Shift {EF_n-EF_ref:+.2f} pp)\n"

    # ---------- Plots ----------
    # Volumenkurve
    plt.figure(figsize=(7.5,4.2))
    plt.plot(t, V_t, lw=2)
    plt.scatter([t[ED_idx]],[V_t[ED_idx]], s=50)
    plt.scatter([t[ES_idx]],[V_t[ES_idx]], s=50)
    plt.xlabel("Zeit (s)"); plt.ylabel("LV-Volumen (mL)"); plt.title("Synthetisches LV-Volumen (1 Beat)")
    plt.tight_layout(); p1 = os.path.join(args.out, "lv_volume_timeseries.png"); plt.savefig(p1, dpi=160); plt.close()

    # EF vs fps
    plt.figure(figsize=(6.4,4.2))
    plt.plot(df_fps["fps"], df_fps["EF_est"], marker="o")
    plt.xlabel("Frame rate (fps)"); plt.ylabel("EF (%)"); plt.title("EF vs Frame rate (no fit)")
    plt.tight_layout(); p2 = os.path.join(args.out, "ef_vs_fps.png"); plt.savefig(p2, dpi=160); plt.close()

    # EF-error vs area bias
    plt.figure(figsize=(6.4,4.2))
    plt.plot(df_bias["area_bias_pct"], df_bias["EF_err"], marker="s")
    plt.xlabel("Segmentation area bias (%)"); plt.ylabel("EF error (pp)"); plt.title("EF error vs area bias (uniform)")
    plt.tight_layout(); p3 = os.path.join(args.out, "ef_error_vs_area_bias.png"); plt.savefig(p3, dpi=160); plt.close()

    # ---------- Outputs ----------
    df_fps.to_csv(os.path.join(args.out, "ef_vs_fps.csv"), index=False)
    df_bias.to_csv(os.path.join(args.out, "ef_error_vs_area_bias.csv"), index=False)

    report = dedent(f"""
    # HSF008 · Heart EF Screener (no data, no fit)

    **Area–Length:** V = 8/(3π)·A²/L,  **EF** = 100·(EDV−ESV)/EDV

    ## Baseline (@ {args.fps_ref} fps, HR={args.hr:.0f} bpm)
    EDV ≈ {EDV:.1f} mL, ESV ≈ {ESV:.1f} mL → **EF ≈ {EF_ref:.1f}%**

    ## EF vs frame rate (Subsampling)
    {df_fps.to_string(index=False)}

    ## EF error vs uniform area bias
    {df_bias.to_string(index=False)}

    ### ED/ES-asymmetrischer Bias (analytisch, no-fit)
    ED-area={args.ed_area_bias:+.2%}, ES-area={args.es_area_bias:+.2%}, ED-L={args.ed_l_bias:+.2%}, ES-L={args.es_l_bias:+.2%}
    → EF_biased = {EF_biased:.2f}%  (Shift {EF_shift:+.2f} pp)
    """).strip()+"\n" + noise_block + \
    "\n**Files:** lv_volume_timeseries.png, ef_vs_fps.png, ef_error_vs_area_bias.png, ef_vs_fps.csv, ef_error_vs_area_bias.csv\n"

    with open(os.path.join(args.out, "summary.md"), "w", encoding="utf-8") as f:
        f.write(report)

    print(f"[OK] Results written to: {args.out}")

if __name__ == "__main__":
    main()
