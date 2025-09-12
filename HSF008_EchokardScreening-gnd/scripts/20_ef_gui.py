#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HSF008 • EF no-data (no-fit) — GUI
Two-column form (Bias & Noise on the RIGHT) + splitter between Plots and Data.
"""

import os, re, math, time, tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# ---------- math helpers (no fit) ----------
def parse_num_list(s: str):
    toks = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', str(s))
    return [float(x) for x in toks]

def parse_int_list(s: str):
    return [int(round(float(x))) for x in parse_num_list(s)]

def area_length_volume_ml(area_cm2, length_cm):
    return (8.0/(3.0*math.pi)) * (area_cm2**2) / np.maximum(length_cm, 1e-9)

def build_cycle(T, fps_ref, a_ed, b_ed, amp_a, amp_b, phase):
    t = np.linspace(0.0, T, int(T*fps_ref), endpoint=False)
    a_t = a_ed - amp_a * (0.5*(1 - np.cos(2*math.pi*(t/T + phase))))
    b_t = b_ed - amp_b * (0.5*(1 - np.cos(2*math.pi*(t/T + phase))))
    L_t = 2.0 * a_t
    A_t = math.pi * a_t * b_t
    V_t = area_length_volume_ml(A_t, L_t)
    return t, A_t, L_t, V_t

def subsample_V(V_ref, T, fps_ref, fps):
    tt  = np.linspace(0.0, T, int(T*fps), endpoint=False)
    idx = (tt * fps_ref).astype(int) % len(V_ref)
    Vn  = V_ref[idx]
    EDV = float(Vn.max()); ESV = float(Vn.min())
    EF  = 100.0*(EDV-ESV)/EDV
    return EF, EDV, ESV

# ---------- GUI ----------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("HSF008 · EF (no-data, no-fit) — GUI")
        self.geometry("1180x900")

        # defaults
        self.vars = {
            "a_ed": tk.StringVar(value="6.0"),
            "b_ed": tk.StringVar(value="3.0"),
            "amp_a": tk.StringVar(value="1.0"),
            "amp_b": tk.StringVar(value="0.8"),
            "phase": tk.StringVar(value="0.05"),
            "hr": tk.StringVar(value="60"),
            "fps_ref": tk.StringVar(value="1000"),
            "fps_list": tk.StringVar(value="1000,60,30,20,15,10,8,5"),
            "seed": tk.StringVar(value="42"),
            "area_bias_list": tk.StringVar(value="-0.15,-0.10,-0.05,-0.02,0,0.02,0.05,0.10,0.15"),
            "ed_area_bias": tk.StringVar(value="0.0"),
            "es_area_bias": tk.StringVar(value="0.0"),
            "ed_l_bias": tk.StringVar(value="0.0"),
            "es_l_bias": tk.StringVar(value="0.0"),
            "area_noise_sd": tk.StringVar(value="0.0"),
            "out_base": tk.StringVar(value=os.path.join("..", "results"))
        }

        self.canvas_widgets = []
        self.fig_refs = []

        self.build_two_col_form()   # top form: left/right columns
        self.build_buttons()        # calculate / end
        self.build_main_splitter()  # plots (top pane) / data (bottom pane)

    # ----- two-column form -----
    def build_two_col_form(self):
        wrap = ttk.Frame(self)
        wrap.pack(side=tk.TOP, fill=tk.X, padx=10, pady=8)

        left  = ttk.Frame(wrap);  right = ttk.Frame(wrap)
        left.grid(row=0, column=0, sticky="nsew", padx=(0,12))
        right.grid(row=0, column=1, sticky="nsew")
        wrap.columnconfigure(0, weight=1)
        wrap.columnconfigure(1, weight=1)

        # LEFT: Geometry + Timing & FPS + Misc
        lf_geom = ttk.LabelFrame(left, text="Geometry (cm)")
        lf_geom.pack(fill=tk.X, padx=0, pady=(0,6))
        self._add_entry(lf_geom, "a_ed",   "a_ed",   hint="half long-axis (ED)")
        self._add_entry(lf_geom, "b_ed",   "b_ed",   hint="half short-axis (ED)")
        self._add_entry(lf_geom, "amp_a",  "amp_a",  hint="systolic reduction of a")
        self._add_entry(lf_geom, "amp_b",  "amp_b",  hint="systolic reduction of b")
        self._add_entry(lf_geom, "phase",  "phase",  hint="phase offset (0..1 cycles)")

        lf_time = ttk.LabelFrame(left, text="Timing & FPS")
        lf_time.pack(fill=tk.X, padx=0, pady=(0,6))
        self._add_entry(lf_time, "hr (bpm)",   "hr")
        self._add_entry(lf_time, "fps_ref",    "fps_ref")
        self._add_entry(lf_time, "fps_list",   "fps_list", width=44, hint='comma list (e.g., "60,30,20,15,10")')

        lf_misc = ttk.LabelFrame(left, text="Misc")
        lf_misc.pack(fill=tk.X, padx=0, pady=(0,0))
        self._add_entry(lf_misc, "seed",     "seed")
        self._add_entry(lf_misc, "out_base", "out_base", width=44, hint="base output folder")

        # RIGHT: Bias & Noise (moved)
        lf_bias = ttk.LabelFrame(right, text="Bias & Noise")
        lf_bias.pack(fill=tk.BOTH, expand=True, padx=0, pady=(0,0))
        self._add_entry(lf_bias, "area_bias_list", "area_bias_list", width=44, hint="uniform area scale list")
        self._add_entry(lf_bias, "ed_area_bias",   "ed_area_bias",   hint="ED only (e.g., 0.05 = +5%)")
        self._add_entry(lf_bias, "es_area_bias",   "es_area_bias",   hint="ES only")
        self._add_entry(lf_bias, "ed_l_bias",      "ed_l_bias",      hint="ED only (length bias)")
        self._add_entry(lf_bias, "es_l_bias",      "es_l_bias",      hint="ES only (length bias)")
        self._add_entry(lf_bias, "area_noise_sd",  "area_noise_sd",  hint="per-frame noise sd (e.g., 0.02)")

    # NOTE: fixed signature — hint first, then width (int).
    def _add_entry(self, parent, label, key, hint=None, width=12):
        row = ttk.Frame(parent); row.pack(fill=tk.X, pady=2)
        ttk.Label(row, text=label, width=14).pack(side=tk.LEFT)
        ent = ttk.Entry(row, textvariable=self.vars[key], width=width)
        ent.pack(side=tk.LEFT)
        if hint:
            ttk.Label(row, text=hint, foreground="#666").pack(side=tk.LEFT, padx=(8,0))

    # ----- buttons -----
    def build_buttons(self):
        bar = ttk.Frame(self)
        bar.pack(side=tk.TOP, fill=tk.X, padx=10, pady=(0,6))
        ttk.Button(bar, text="Calculate", command=self.on_calculate).pack(side=tk.LEFT)
        ttk.Button(bar, text="End", command=self.destroy).pack(side=tk.RIGHT)

    # ----- splitter with plots (top) and data (bottom) -----
    def build_main_splitter(self):
        self.split = ttk.Panedwindow(self, orient=tk.VERTICAL)
        self.split.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=(0,10))

        # --- TOP pane: scrollable plots ---
        plots_pane = ttk.Frame(self.split)
        self.split.add(plots_pane, weight=3)

        self.plot_canvas = tk.Canvas(plots_pane)
        pv = ttk.Scrollbar(plots_pane, orient="vertical", command=self.plot_canvas.yview)
        ph = ttk.Scrollbar(plots_pane, orient="horizontal", command=self.plot_canvas.xview)
        self.plot_canvas.configure(yscrollcommand=pv.set, xscrollcommand=ph.set)
        pv.pack(side=tk.RIGHT, fill=tk.Y)
        ph.pack(side=tk.BOTTOM, fill=tk.X)
        self.plot_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.plot_inner = ttk.Frame(self.plot_canvas)
        self.plot_canvas.create_window((0,0), window=self.plot_inner, anchor="nw")
        self.plot_inner.bind("<Configure>", lambda e: self.plot_canvas.configure(scrollregion=self.plot_canvas.bbox("all")))

        # --- BOTTOM pane: scrollable data/report ---
        data_pane = ttk.Frame(self.split)
        self.split.add(data_pane, weight=2)

        self.data_text = tk.Text(data_pane, wrap=tk.NONE)
        dv = ttk.Scrollbar(data_pane, orient="vertical", command=self.data_text.yview)
        dh = ttk.Scrollbar(data_pane, orient="horizontal", command=self.data_text.xview)
        self.data_text.configure(yscrollcommand=dv.set, xscrollcommand=dh.set)
        dv.pack(side=tk.RIGHT, fill=tk.Y)
        dh.pack(side=tk.BOTTOM, fill=tk.X)
        self.data_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.data_text.insert(tk.END, "Ready. Adjust parameters and click Calculate.\n")

    # ----- helpers -----
    def _clear_plots(self):
        for w in self.plot_inner.winfo_children():
            try: w.destroy()
            except Exception: pass
        self.canvas_widgets.clear()
        self.fig_refs.clear()

    def _embed_fig(self, fig: Figure, title: str, save_path: str):
        fig.savefig(save_path, dpi=160, bbox_inches="tight")
        lbl = ttk.Label(self.plot_inner, text=title, font=("Segoe UI", 10, "bold"))
        lbl.pack(anchor="w", pady=(10,2))
        canvas = FigureCanvasTkAgg(fig, master=self.plot_inner)
        widget = canvas.get_tk_widget()
        widget.pack(fill=tk.BOTH, expand=True)
        canvas.draw()
        self.canvas_widgets += [lbl, widget]
        self.fig_refs.append(fig)

    def _set_data_text(self, text: str):
        self.data_text.delete("1.0", tk.END)
        self.data_text.insert(tk.END, text)

    # ----- compute -----
    def on_calculate(self):
        try:
            # params
            a_ed   = float(self.vars["a_ed"].get())
            b_ed   = float(self.vars["b_ed"].get())
            amp_a  = float(self.vars["amp_a"].get())
            amp_b  = float(self.vars["amp_b"].get())
            phase  = float(self.vars["phase"].get())
            hr     = float(self.vars["hr"].get())
            fps_ref= int(float(self.vars["fps_ref"].get()))
            fps_list  = parse_int_list(self.vars["fps_list"].get())
            seed      = int(float(self.vars["seed"].get()))
            area_list = parse_num_list(self.vars["area_bias_list"].get())
            ed_ab = float(self.vars["ed_area_bias"].get())
            es_ab = float(self.vars["es_area_bias"].get())
            ed_lb = float(self.vars["ed_l_bias"].get())
            es_lb = float(self.vars["es_l_bias"].get())
            area_noise_sd = float(self.vars["area_noise_sd"].get())
            out_base = self.vars["out_base"].get() or "results"

            # output dir
            ts = time.strftime("%Y%m%d_%H%M%S")
            outdir = os.path.join(out_base, f"HSF008_gui_{ts}")
            os.makedirs(outdir, exist_ok=True)

            # physics
            T = 60.0/hr
            t, A_t, L_t, V_t = build_cycle(T, fps_ref, a_ed, b_ed, amp_a, amp_b, phase)
            ED_idx = int(np.argmax(V_t)); ES_idx = int(np.argmin(V_t))
            EDV, ESV = float(V_t[ED_idx]), float(V_t[ES_idx])
            EF_ref = 100.0*(EDV-ESV)/EDV

            # EF vs fps
            rows_fps = []
            for fps in fps_list:
                EF, ed, es = subsample_V(V_t, T, fps_ref, fps)
                rows_fps.append(dict(fps=fps, EF_est=EF, EDV=ed, ESV=es, EF_err=EF - EF_ref))
            df_fps = pd.DataFrame(rows_fps)
            df_fps.to_csv(os.path.join(outdir, "ef_vs_fps.csv"), index=False)

            # EF error vs uniform area bias
            rows_bias = []
            for e in area_list:
                A_b = (1.0 + e) * A_t
                V_b = area_length_volume_ml(A_b, L_t)
                ed  = float(V_b.max()); es = float(V_b.min())
                ef  = 100.0*(ed-es)/ed
                rows_bias.append(dict(area_bias_pct=100*e, EF_est=ef, EDV=ed, ESV=es, EF_err=ef - EF_ref))
            df_bias = pd.DataFrame(rows_bias).sort_values("area_bias_pct")
            df_bias.to_csv(os.path.join(outdir, "ef_error_vs_area_bias.csv"), index=False)

            # Asymmetric ED/ES (analytic)
            EDV_b = EDV * (1 + ed_ab)**2 / (1 + ed_lb)
            ESV_b = ESV * (1 + es_ab)**2 / (1 + es_lb)
            EF_b  = 100.0*(EDV_b-ESV_b)/EDV_b
            EF_shift = EF_b - EF_ref

            # optional noise (report only)
            noise_line = ""
            if area_noise_sd > 0:
                rng = np.random.default_rng(seed)
                A_noise = A_t * (1.0 + rng.normal(0.0, area_noise_sd, size=A_t.shape))
                V_noise = area_length_volume_ml(A_noise, L_t)
                EDV_n, ESV_n = float(V_noise.max()), float(V_noise.min())
                EF_n = 100.0*(EDV_n-ESV_n)/EDV_n
                noise_line = f"\nFrame noise σ_area={area_noise_sd:.1%}: EF={EF_n:.2f}% (shift {EF_n-EF_ref:+.2f} pp)"

            # ---- plots (top pane) ----
            self._clear_plots()

            fig1 = Figure(figsize=(9.6, 3.8))
            ax1 = fig1.add_subplot(111)
            ax1.plot(t, V_t, lw=2)
            ax1.scatter([t[ED_idx]],[V_t[ED_idx]], s=40)
            ax1.scatter([t[ES_idx]],[V_t[ES_idx]], s=40)
            ax1.set_xlabel("Time (s)"); ax1.set_ylabel("LV Volume (mL)")
            ax1.set_title("Synthetic LV volume over one beat")
            self._embed_fig(fig1, "Figure 1 — LV volume time series",
                            os.path.join(outdir, "lv_volume_timeseries.png"))

            fig2 = Figure(figsize=(9.0, 3.6))
            ax2 = fig2.add_subplot(111)
            ax2.plot(df_fps["fps"], df_fps["EF_est"], marker="o")
            ax2.set_xlabel("Frame rate (fps)"); ax2.set_ylabel("EF (%)")
            ax2.set_title("EF vs frame rate (no fit)")
            self._embed_fig(fig2, "Figure 2 — EF vs frame rate",
                            os.path.join(outdir, "ef_vs_fps.png"))

            fig3 = Figure(figsize=(9.0, 3.6))
            ax3 = fig3.add_subplot(111)
            ax3.plot(df_bias["area_bias_pct"], df_bias["EF_err"], marker="s")
            ax3.set_xlabel("Segmentation area bias (%)"); ax3.set_ylabel("EF error (pp)")
            ax3.set_title("EF error vs area bias (uniform)")
            self._embed_fig(fig3, "Figure 3 — EF error vs area bias",
                            os.path.join(outdir, "ef_error_vs_area_bias.png"))

            # ---- data (bottom pane) ----
            lines = []
            lines.append(f"Output folder: {os.path.abspath(outdir)}")
            lines.append(f"Baseline @ {fps_ref} fps, HR={hr:.0f} bpm")
            lines.append(f"  EDV ≈ {EDV:.1f} mL, ESV ≈ {ESV:.1f} mL → EF ≈ {EF_ref:.1f}%")
            lines.append("")
            lines.append("EF vs frame rate"); lines.append(df_fps.to_string(index=False))
            lines.append("")
            lines.append("EF error vs uniform area bias"); lines.append(df_bias.to_string(index=False))
            lines.append("")
            lines.append("Asymmetric ED/ES bias (analytic, no fit)")
            lines.append(f"  ED-area={ed_ab:+.2%}, ES-area={es_ab:+.2%}, ED-L={ed_lb:+.2%}, ES-L={es_lb:+.2%}")
            lines.append(f"  → EF_biased = {EF_b:.2f}%  (shift {EF_shift:+.2f} pp)")
            if noise_line: lines.append(noise_line)
            lines.append("")
            lines.append("Saved files:")
            lines.append(f"  {os.path.join(outdir,'lv_volume_timeseries.png')}")
            lines.append(f"  {os.path.join(outdir,'ef_vs_fps.png')}")
            lines.append(f"  {os.path.join(outdir,'ef_error_vs_area_bias.png')}")
            lines.append(f"  {os.path.join(outdir,'ef_vs_fps.csv')}")
            lines.append(f"  {os.path.join(outdir,'ef_error_vs_area_bias.csv')}")
            self._set_data_text("\n".join(lines))

        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    App().mainloop()
