# pp elastic: forward-cone fit

**Model**: $d\sigma/dt = A\,e^{-B|t|}$, fitted for $|t| \le 0.120\,\mathrm{GeV}^2$.

- $A=(d\sigma/dt)_{t=0} = 511\;\mathrm{mb/GeV^2}$  (±1.7)  
- $B = 20.1\;\mathrm{GeV^-2}$  (±0.076)  
- $\rho = 0.14$ (assumed)

Optical theorem:
\[
\sigma_\mathrm{tot} = \sqrt{\frac{16\pi\,0.389379\;A}{1+\rho^2}} = 99.01\;\mathrm{mb}.
\]

Single-exponential elastic estimate:
\[
\sigma_\mathrm{el} \simeq \frac{A}{B} = 25.43\;\mathrm{mb}.
\]

Points used: 49

Artifacts:
- Plot: `pp_forward_fit.png`
- Table: `pp_results.csv`
