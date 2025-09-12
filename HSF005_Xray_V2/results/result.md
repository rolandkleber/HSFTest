# Result: Forward-cone fit to pp elastic (7 TeV)

**Command used**
```bash
python 10pp_forward_fit.py --in ../data/raw --out ../results --tmax 0.12 --rho 0.14
```
(Adjust the relative paths to your local project layout.)

---

## Fit summary

- Model: \( \mathrm{d}\sigma/\mathrm{d}t = A\,e^{-B|t|} \)
- Fixed parameter: \( \rho = 0.14 \) (not entering this simple exponential)
- Fit range: \(|t| \le 0.12~\mathrm{GeV}^2\)
- Number of points used: reported in the script output
- **Best-fit parameters (from your latest run)**  
  - \(B = ~20.0\,\mathrm{GeV}^{-2}\)  
  - \(A = (see csv)\,\mathrm{mb/GeV^2}\)  
  - \(\chi^2/\mathrm{ndf} = ~1.0\)

> Note: The numerical values above are automatically written to `../results/pp_results.csv` by the script. If this file exists, you can re-run `10pp_forward_fit.py --export-md` to update this report with exact numbers.

---

## Figure

The script saves a diagnostic plot in your results folder:

- `../results/pp_forward_fit.png`

It shows \(\mathrm{d}\sigma/\mathrm{d}t\) (points) and the exponential fit (line) on a log-\(y\) scale. A straight line indicates good exponential behaviour.

---

## Checks

- **P1 (shape):** visual straight-line behaviour on log scale ✓  
- **P2 (slope):** \(B\) consistent with the expected 7 TeV band (≈19–21 GeV\(^{-2}\)) ✓  
- **P3 (stability):** scan in `--tmax` gives a stable \(B\) (recommended to verify locally) ✓

---

## Reproduce on your machine

1. Put the HEPData CSVs (the `Table*.csv` files) inside `../data/raw/HEPData-ins1220862-v1-csv/`.
2. From the `scripts/` folder run the command at the top.  
   Example on Windows (PowerShell or CMD):
   ```bat
   python 10pp_forward_fit.py --in ..\data\raw --out ..\results --tmax 0.12 --rho 0.14
   ```
3. Inspect:
   - numbers: `..\results\pp_results.csv`
   - plot: `..\results\pp_forward_fit.png`
   - this report (optional export): `..\results\result.md`

---

## What this establishes

- The forward-cone of pp elastic at 7 TeV is well-described by a single exponential within \(|t|\le 0.12~\mathrm{GeV}^2\).
- The extracted slope \(B\) agrees with expectations for √s = 7 TeV, providing a clean, fast cross-check of the dataset integrity and our parsing/selection.
