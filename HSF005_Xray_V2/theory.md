
Main Theory pdf https://doi.org/10.5281/zenodo.16953782
# Theory: Forward pp Elastic Scattering (7 TeV)

This note explains the minimal theory used to analyse the **forward elastic** proton–proton (pp) scattering data from TOTEM (√s = 7 TeV). It is intentionally compact and only contains the pieces we actually use in the fitting code.

---

## 1) Kinematics and variable choice

- Four–momentum transfer: \(t = (p_\text{out} - p_\text{in})^2 \simeq -p^2\theta^2\) for small scattering angle \(\theta\).
- In the **forward cone** (\(|t| \lesssim 0.1\,\mathrm{GeV}^2\)) the elastic differential cross-section is observed to be close to an exponential in \(|t|\).

We analyse the measured differential cross-section \( \mathrm{d}\sigma/\mathrm{d}t \) as a function of \(|t|\).

---

## 2) Forward amplitude and optical theorem

Let \(f(t)\) denote the spin-averaged elastic scattering amplitude. Near \(t=0\) the **optical theorem** relates the imaginary part at \(t=0\) to the **total** cross-section:
\[
\sigma_\text{tot} = \frac{4\pi}{p_\text{cm}}\,\mathrm{Im}\,f(0)\,,
\]
where \(p_\text{cm}\) is the proton momentum in the centre-of-mass frame.

Define the **\(\rho\)-parameter** as the ratio of the real to imaginary parts at \(t=0\):
\[
\rho \equiv \frac{\mathrm{Re}\,f(0)}{\mathrm{Im}\,f(0)}\,.
\]

TOTEM measurements near 7 TeV give \(\rho \sim 0.14\). In our baseline fits we **fix** \(\rho\) to the value passed by the user (default 0.14) to avoid degeneracy with the overall normalisation in the limited \(|t|\) window.

---

## 3) Exponential forward cone

Empirically, in the forward region the modulus of the hadronic amplitude decreases approximately exponentially with \(|t|\). Correspondingly we fit the data with
\[
\frac{\mathrm{d}\sigma}{\mathrm{d}t} \;=\; A\,e^{-B\,|t|}\,,
\]
with two parameters
- \(A\) — normalisation (mb/GeV\(^2\)),
- \(B\) — **slope** (GeV\(^{-2}\)).

The parameter \(B\) governs how rapidly the cross-section falls with \(|t|\). For 7 TeV pp elastic scattering, one expects \(B \approx 19\text{–}21~\mathrm{GeV}^{-2}\) in the range \(|t|\lesssim 0.1~\mathrm{GeV}^2\).

> Why this is sufficient here: within \(|t|\lesssim 0.12~\mathrm{GeV}^2\) the hadronic term dominates and Coulomb–nuclear interference (CNI) corrections only shift the very first bins. A pure exponential captures the physics well enough to test forward consistency and to compare with published slopes.

---

## 4) From the forward fit to consistency checks

- **Slope comparison:** \(B\) extracted from the exponential fit is compared with literature values (forward-cone slopes published by TOTEM for the same dataset).
- **Goodness of fit:** \(\chi^2/\mathrm{ndf}\) close to 1 indicates that the simple exponential is adequate in the chosen window.
- **Stability w.r.t. \(t_\text{max}\):** repeating the fit for several \(t_\text{max}\) values checks for curvature or onset of the diffractive minimum outside the fit window.

---

## 5) Units and conventions

- \(|t|\): GeV\(^2\)
- \(\mathrm{d}\sigma/\mathrm{d}t\): mb/GeV\(^2\)
- \(B\): GeV\(^{-2}\)

All fits are performed in linear space but diagnostics are displayed on a **logarithmic** \(y\)-axis to make the exponential behaviour a straight line.

---

## 6) Optional refinements (not enabled by default)

- Include CNI terms at very small \(|t|\) and fit \(\rho\) simultaneously with \(A, B\).
- Allow a small curvature: \( \ln(\mathrm{d}\sigma/\mathrm{d}t) = \ln A - B|t| + C |t|^2 \).
- Global normalisation nuisance with prior when combining tables from different optics settings.

These are unnecessary for the quick forward-cone validation we perform here.
