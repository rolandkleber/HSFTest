# prediction.md

## Prior expectations (before running)

* The O₂ spectrum in **Bogumil et al. (2003), 293 K, 235–389 nm** should yield a **high-quality KK match**:

  * smooth baseline with structured band features,
  * $|r| \approx 0.99$ after one global scale $s=\mathcal{O}(1)$,
  * normalized RMS residual $\lesssim 10^{-3}$–$10^{-2}$.
* Edge effects: the very ends of the band may show mild ringing; tapering should confine residuals to the outer \~10 % if they appear.
* Sensitivity: down-sampling or a few outliers should not change $|r|$ by more than a few $10^{-3}$.

## Falsification

* A persistent low correlation ($|r|<0.95$) **after** cleaning and tapering would contradict the KK expectation in this band and would motivate checks of file integrity, units, or (if robust) deeper investigation.

## Extensions (still O₂-only)

* Repeat on another O₂ dataset (e.g., Ackerman 300 K 176–201 nm or Greenblatt 650–789 nm).
* Report the intersection of passes: consistent high $|r|$ across disjoint O₂ bands gives a stronger statement for O₂ as a species.

---


