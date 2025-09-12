# MOL-01 · H₂O Linewidth T/P Sweep (Voigt) — Predictions

**Target line:** ν₀ ≈ **7181.1557 cm⁻¹**  
**Reference broadening:** γ_ref = **0.0997 cm⁻¹** at **T_ref = 296 K**, **P_ref = 1 atm**  
**Temperature exponent:** n_air = **0.760**

## Expected behaviors
1. **Doppler width**  
   \(\mathrm{FWHM}_G(T) = \nu_0\,\sqrt{\frac{8 k_B T \ln 2}{m c^2}}\) → **slope 0.5** in log–log vs. T.

2. **Collisional width**  
   \(\gamma(T,P) = \gamma_\mathrm{ref}\,(T_\mathrm{ref}/T)^n\,(P/P_\mathrm{ref})\) → **slope −n** in log–log vs. T; **linear in P** at fixed T.  
   Expected pressure slope at temperature \(T\):  
   \(\frac{\mathrm d\gamma}{\mathrm d P} = \gamma_\mathrm{ref} (T_\mathrm{ref}/T)^n\).

## Quick numbers
- **At 296 K:** γ ≈ **0.0997 cm⁻¹** (by definition, 1 atm).  
- **At 375 K:** expected \(\mathrm d\gamma/\mathrm dP\) ≈ **0.08329 cm⁻¹/atm**.

## Pass/Fail
- Doppler slope: **0.50 ± 0.05**  
- Collisional slope: **−n_air ± 0.05**  
- γ vs P: linear with **R² > 0.995**, slope within a few %.  
- Fit quality: χ²/ndf ≈ 1, residuals unstructured.

## How to run
```bash
python scripts/23_linewidth_TP_sweep.py \
  --hitran data/raw/H2O_lines_7150-7220.csv \
  --auto-strongest \
  --T 250,300,350,400,500,700 \
  --P 0.2,0.5,1.0 \
  --instr-fwhm 0.10 \
  --win 1.0 \
  --out results/MOL-01_H2O_TP_Voigt
```