@echo off
cd /d "%~dp0"
set "OUT=..\results"

python 10_ef.py ^
 --out="%OUT%" ^
 --a-ed=6.0 --b-ed=3.0 --amp-a=1.0 --amp-b=0.8 --phase=0.05 ^
 --hr=60 --fps-ref=1000 ^
 --fps-list="1000,60,30,20,15,10,8,5" ^
 --seed=42 ^
 --area-bias-list="-0.15,-0.10,-0.05,-0.02,0,0.02,0.05,0.10,0.15" ^
 --ed-area-bias=0.00 --es-area-bias=0.00 --ed-l-bias=0.00 --es-l-bias=0.00 ^
 --area-noise-sd=0.00

echo.
echo Done. Results in "%OUT%"

REM ------------------------------------------------------------
REM  FLAG EXPLANATIONS (No-Data / No-Fit EF test)
REM
REM  --out              Output folder for plots/CSVs/summary.md
REM
REM  Geometry (A4C ellipse, centimeters):
REM    --a-ed           Half long axis a at end-diastole (ED)
REM    --b-ed           Half short axis b at ED
REM    --amp-a          Systolic reduction of a (contraction amplitude)
REM    --amp-b          Systolic reduction of b
REM    --phase          Phase offset of contraction (0..1 cycles)
REM
REM  Timing:
REM    --hr             Heart rate in bpm (60 -> 1 s per beat)
REM    --fps-ref        Reference sampling rate for the “true” curve
REM    --fps-list       Comma list of frame rates to test (subsampling)
REM
REM  Reproducibility:
REM    --seed           RNG seed (used only when noise is enabled)
REM
REM  Uniform segmentation bias (affects ED and ES equally):
REM    --area-bias-list Comma list of relative area scalings
REM                      e.g., -0.10 = -10 %, +0.05 = +5 %
REM                      Note: with uniform bias EF usually cancels out,
REM                      because both EDV and ESV scale the same.
REM
REM  Asymmetric ED/ES bias (pure analytic shift, still no fit):
REM    --ed-area-bias   Relative area bias in ED only (e.g., +0.05 = +5 %)
REM    --es-area-bias   Relative area bias in ES only
REM    --ed-l-bias      Relative long-axis (L) bias in ED only
REM    --es-l-bias      Relative long-axis (L) bias in ES only
REM                      → These create an EF shift because EDV and ESV
REM                        are scaled DIFFERENTLY.
REM
REM  Frame-wise noise (still no fit):
REM    --area-noise-sd  Std. dev. of per-frame relative area noise
REM                      (0.02 ≈ ±2 %). Uses --seed for reproducibility.
REM ------------------------------------------------------------
