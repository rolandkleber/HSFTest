# theory.md
Main Theory pdf https://doi.org/10.5281/zenodo.16953782
## Scope

This project checks whether **oxygen (O₂)** optical absorption in the 235–389 nm band obeys **causality** via the **Kramers–Kronig (KK)** relations. We restrict all analysis and claims to O₂ data.

## Objects and relations

Complex refractive index

$$
\tilde n(\omega)=n(\omega)+i\,\kappa(\omega).
$$

Absorption coefficient $\alpha$ and extinction $\kappa$ are linked by

$$
\alpha(\omega)=\frac{4\pi}{\lambda}\,\kappa(\omega)
=\frac{2\omega}{c}\,\kappa(\omega),
\qquad
\kappa(\omega)=\frac{\alpha(\omega)\,\lambda}{4\pi}.
$$

**KK (causality)**

$$
n(\omega)-1=\frac{2}{\pi}\,\mathcal{P}\!\int_{0}^{\infty}
\frac{\omega'\,\kappa(\omega')}{\omega'^2-\omega^2}\,d\omega'.
$$

Because measurements cover a **finite band**, we compare **shapes** (after removing a constant baseline and applying a smooth taper) and allow one global **scale factor** $s$ when matching the KK-predicted dispersive “twin”.

## Pipeline (O₂ only)

1. **Parse** raw two-column text $(\lambda,\,\text{value})$ in nm and the reported absorption unit (auto-detected as $\alpha$ here).
2. **Convert** $\{\lambda\}\!\to\!\{\omega\}$ and $\alpha\!\to\!\kappa$.
3. **Uniformize** on an evenly spaced $\omega$-grid; drop bad rows/NaNs.
4. **Baseline + taper** (cosine, \~7.5 % edges).
5. **Hilbert transform** to compute the **KK dispersive twin** from $\kappa(\omega)$.
6. **Fit one scale $s$** and evaluate:

   * Pearson correlation $|r|$ between measured curve and KK twin (normalized),
   * **RMS residual** on unit-variance signals,
   * the fitted **scale** $s$.

**Pass criterion:** $|r|\gtrsim 0.98$ and small RMS (≈10⁻³–10⁻² after normalization) **within the O₂ band**.

## Why this matters

KK is a **model-independent** consequence of causality for passive media. If O₂ spectra violated KK beyond small numerical/edge effects, that would be a red flag for the measurement or for new physics. A good match sets a clean, quantitative bound on any such anomaly *in this band for O₂*.

---
