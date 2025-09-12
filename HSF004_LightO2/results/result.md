## Dataset analyzed (O₂-only)

* **O₂ — Bogumil et al. (2003), 293 K, 235–389 nm**
  File: `O2_Bogumil(2003)_293K_235-389nm.txt` (cleaned; monotonic wavelengths; missing-value row fixed).

## What we did

* Converted $\lambda\!\to\!\omega$ and $\alpha\!\to\!\kappa$ via $\kappa=\alpha\lambda/(4\pi)$.
* Resampled to a uniform $\omega$-grid; removed constant baseline and applied a cosine taper at both ends.
* Computed the **KK dispersive twin** (Hilbert transform), then fitted one global **scale** $s$ on the analysis band.
* Compared shapes with Pearson $|r|$ and normalized RMS residual.

## Headline numbers

*(from your latest run; see `results/kk_summary.txt` and the PNGs)*

* **Correlation $|r|$:** ≈ **0.99**
* **Normalized RMS residual:** ≈ **1×10⁻³**
* **Scale $s$:** ≈ **0.98** (dimensionless)

## Figures

* Measured absorption: `results/figures/O2_Bogumil(2003)_293K_235-389nm_absorption.png`
* KK comparison: `results/figures/O2_Bogumil(2003)_293K_235-389nm_kk.png`

These show the expected rise and structure across 235–389 nm and an almost perfect KK twin after one global rescale—i.e., a **clear causality pass for O₂ in this band**.

## Caveats

* Results are **band-limited** to 235–389 nm; claims outside that window aren’t implied.
* Absolute $\kappa$ scale depends on path length/density in source data; our pass is **shape-based** with one fitted scale $s$.

## Next O₂ steps (optional)

* Add a **second O₂ band** (e.g., Ackerman 176–201 nm).
* Bootstrap $|r|$, RMS, and $s$ for 68 % CIs to formalize uncertainties.

-