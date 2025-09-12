## Title
Main Theory pdf https://doi.org/10.5281/zenodo.16953782
HSF P1 Null Test with Cosmic-Ray Muon Time Series

## Objective

Test the HSF “Track-B” prediction (P1) that a coherent field could imprint a **narrow, coherent modulation** at frequency $f_{\mathrm{DM}}$ on an otherwise stochastic count rate. If no line is observed, report a **95% confidence upper limit** on the fractional modulation amplitude $A(f)$.

## Prediction (P1)

Under HSF, a weak coherent field $\phi(t)=\phi_0\cos(2\pi f_{\mathrm{DM}} t)$ can yield a small periodic modulation of a rate observable $R(t)$:

$$
R(t) = R_0\,[1 + A \cos(2\pi f_{\mathrm{DM}} t + \varphi)] + \text{noise},
$$

with $A \ll 1$. In a power spectrum (periodogram/PSD) of a suitably detrended time series, this appears as a **narrow peak** at $f_{\mathrm{DM}}$ (and possibly small harmonics if the response is mildly nonlinear).

**Null expectation (no HSF):** the spectrum is consistent with Poisson + environmental noise, with **no unresolved narrow lines** above the statistical floor.

## Observable and Data Choice

* **Observable:** counts of **cosmic-ray muons** in fixed time bins $\Delta t$ (here 300 s).
* **Rationale:** Muon counts provide long, simple time series with many events; if a global, weak coherent effect exists, a narrow modulation might be detectable (or constrained) in the PSD.
* **Caveat:** Standard physics already explains the *mean* rate; we are only probing **extra, coherent modulations**.

## Analysis Pipeline (Reproducible)

1. **Parse events** from raw logs (CosmicWatch format variants).
2. **Build rate series:** resample events into fixed bins $\Delta t = 300$ s → $N_k$ counts/bin.
3. **Detrend:** remove slow drifts via a centered rolling mean (window $\sim 30$ min) to suppress low-frequency systematics.
4. **Spectrum:** compute a periodogram/PSD (FFT on detrended series). For gapped data, a Lomb–Scargle variant can be used.
5. **Search band:** frequencies $f \in [f_{\min}, f_{\mathrm{Nyq}}]$ with
   $f_{\min}\approx 1/T$ (total duration $T$) and $f_{\mathrm{Nyq}} \approx 1/(2\Delta t)$.
6. **Decision rule:** no significant unresolved peaks ⇒ compute a **95% C.L. upper limit** $A_{95\%}(f)$ on a hypothetical sinusoidal line, using the measured noise floor (heuristic or Poisson-aware model).

## Reported Quantity

* **Fractional modulation amplitude** $A$: the amplitude of a sinusoidal perturbation relative to the mean rate $R_0$.
* **Upper-limit curve** $A_{95\%}(f)$: the smallest amplitude we would have detected at 95% C.L.; values above this are excluded (given our assumptions and preprocessing).

## Assumptions & Limitations

* The detector response is linear in the small-signal regime ($A\ll 1$).
* Detrending does not suppress a *narrow* line in the search band.
* Environmental effects (pressure/temperature) are not explicitly corrected here; residuals contribute to the noise floor → **limits are conservative**.
* Sidebands around a “carrier” (strict P1 sideband picture) are not applicable to pure count rates; we test **narrow coherent lines** as a proxy.

## Outcome Interpretation

* **Positive detection:** a stable, unresolved peak exceeding the statistical expectation at a fixed $f$. Requires extensive vetting (systematics, instrument artifacts, look-elsewhere).
* **Null:** no such peak; quote $A_{95\%}(f)$. This **constrains** models that predict stronger coherent modulations of muon rates in the searched band.

---

