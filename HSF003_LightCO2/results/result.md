# HSF003_LightCO2 — Kramers–Kronig Consistency Check (Alpha)

**Dataset**  
Lewis & Carver (1983), CO₂ absorption at **370 K**, **121–193 nm** (uniformly resampled to 2048 points).

**Goal**  
Check if the measured absorption κ(λ) is consistent with the dispersion derived via a Kramers–Kronig (KK) transform over the band, and bound any extra dispersive component.

---

## Results (summary)

- Uniform grid points: **2048**  
- Optimal scale (best-fit amplitude factor) **s ≈ 0.981**  
- RMS mismatch (κ vs. scaled KK[κ]): **3.211 × 10⁻³**  
- Linear correlation **|corr(κ, KK)| ≈ 0.989**  
- Practical tip: cropping to a single hump (e.g., **135–185 nm**) reduces edge influence.

**Bound:** From the residuals, any additional dispersive contribution is **≲ 0.3 %** of the band amplitude.

---

## Figures

- `results/figures/kk_absorption.png` — absorption κ(λ) and KK-derived curve (scaled by s)  
- `results/figures/kk_dispersion.png` — residuals and dispersion view

---

## Method (very brief)

1. Convert wavelength to angular frequency grid; resample κ to a uniform grid (N = 2048).  
2. Compute the KK transform to obtain the implied dispersion.  
3. Fit a scalar **s** minimizing RMS between κ and KK[κ] (accounts for overall scale).  
4. Report RMS and Pearson correlation; inspect residuals to bound extra dispersion.

---

## Data source & license

- Lewis, B.R., Carver, J.H. (1983): CO₂ VUV absorption, 370 K, 121–193 nm. *(Add full citation/DOI if available.)*  
- This repository includes derived figures and summary; raw text files are in `data/raw/`.

---

*Reproducibility:* Windows-only. Install deps (`pip install -r requirements.txt`), then run `scripts\10_exe.bat`. Outputs appear in `results\`.
