# Prediction (pre-fit expectations)

---

## Dataset and window

- Process: \(pp \to pp\) elastic, √s = 7 TeV (TOTEM).
- Observable: single-differential cross-section \(\mathrm{d}\sigma/\mathrm{d}t\).
- Fit window: \(|t| \le t_\text{max}\) with a typical choice \(t_\text{max}=0.12~\mathrm{GeV}^2\).
- Very first few bins can receive Coulomb–nuclear interference (CNI); if needed we drop the first point(s).

## Expected behaviour

- **Exponential forward cone:** \(\mathrm{d}\sigma/\mathrm{d}t \approx A \exp(-B|t|)\).
- **Slope:** \(B \in [19, 21]~\mathrm{GeV}^{-2}\) at √s = 7 TeV for the window above.
- **Goodness of fit:** \(\chi^2/\mathrm{ndf}\sim 1\) if experimental uncertainties are used as reported.
- **Stability:** \(B\) stable when scanning \(t_\text{max}\) within \(0.08\text{–}0.14~\mathrm{GeV}^2\).

## Pass/Fail criteria (for a quick sanity check)

- **P1 (shape):** Residuals (pulls) are structureless; \(\chi^2/\mathrm{ndf} \in [0.6, 1.4]\).
- **P2 (slope):** \(B\) within the expected 7 TeV band (19–21 GeV\(^{-2}\)).
- **P3 (stability):** \(B\) varies by \(\lesssim 5\%\) if \(t_\text{max}\) is changed by ±0.02 GeV\(^2\).

If all three hold, the dataset is forward-cone consistent and suitable for downstream use (e.g., extrapolations or cross-checks).