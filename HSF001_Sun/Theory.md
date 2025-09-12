Main Theory pdf https://doi.org/10.5281/zenodo.16953782
### Data inclusion (freeze)

* **Instruments:** RSTN (245 MHz–15.4 GHz), NoRP (1/2/3.75/9.4/17/35 GHz), EOVSA (≈1–18 GHz).
* **Dates:** any **30** quiet-Sun days within **2017–2025**.
* **Quiet-Sun rule:** GOES 1–8 Å flux **< 1e-6 W/m²** (no C-class+ flare) during the daily window.
* **Required metadata:** frequency centers & bandwidths; cadence ≥ **1 s** (or daily averages for P0); calibrated units.
* **Quality gates:** SNR ≥ **10** at carriers; continuous band ±**3×** linewidth; no obvious RFI at candidate $f_{\rm DM}$.

### PSD & sideband search (freeze)

* Detrend (linear), mean-center.
* **Welch PSD:** Hann window, segment length $N=2^{15}$, **50 %** overlap.
* Target resolution: $\delta f/f \approx 10^{-6}$ (or best achievable).
* $f_{\rm DM}$ grid: log-spaced **1 kHz–10 GHz**, step $\Delta\log_{10} f = 10^{-3}$.
* **Metrics:** $S=(P_{+1}-P_{-1})/(P_{+1}+P_{-1})$; $R=$ (ec sideband power)/(plasma sideband power).
* **Significance:** local-noise test on combined bins within ±1 FWHM; **FDR q=0.05** across the $f_{\rm DM}$ grid.
* **Persistence:** same $f_{\rm DM}$ (within one grid bin) on **≥3 consecutive days**.

### KK (causality) test (freeze)

* **Model A:** absorptive bump only.
* **Model B:** absorptive + **Hilbert-pair** dispersive component (same center/width).
* **Pass:** $\Delta\mathrm{BIC}=\mathrm{BIC}_A-\mathrm{BIC}_B \ge 6$.
* Fit band: ±**5 FWHM** around the feature; baseline linear.

### Numeric pass bands (freeze)

* **P0 slope:** $[-2.2,\,-1.8]$.
* **P1 symmetry:** $|S|\le 0.05$.
* **P1 ec\:plasma ratio:** $1.6 \le R \le 2.4$.
* **P1 width:** $\Delta f/f \in [3\times10^{-7},\,3\times10^{-6}]$.
* **P2 KK:** $\Delta\mathrm{BIC} \ge 6$.
* **P4 GR check:** $|\Delta f/f + \Phi/c^2| \le 5\times10^{-6}$.

### Controls & exclusions (freeze)

* Drop $f_{\rm DM}$ coincident with known instrument clocks/IFs; require cross-instrument consistency when available.
* Same windowing for all; reject within one bin of DC/Nyquist.
* Negative control (off-Sun/night): must be null at the same $f_{\rm DM}$.

### Outputs to publish (freeze)

* `tables/sidebands.csv`: date, carrier, $f_{\rm DM}$, $P_{+1}$, $P_{-1}$, $S$, $R$, q-value, flags.
* `tables/kk_bic.csv`: $\Delta$BIC per feature.
* `figures/`: PSD zooms, symmetry scatter, KK-residual plots.
* `REPORT.md`: Pass/Fail for **P0–P5** and, if null, **$d_e(f_{\rm DM})$** exclusion (state conversion you use).
* **Seed:** 1729 (fixed wherever randomness appears).

