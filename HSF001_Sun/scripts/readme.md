# HSF001_Sun (Alpha, Windows-only)

This is an early preview. The repo includes **results/** so you can inspect figures immediately.
The end-to-end pipeline is **partially implemented**; some stages are placeholders.

## What works
- `30_exe.bat`: builds candidate spectra from the included RSTN sample and plots them.
- `35_exe.bat`: generates a small README snippet for top candidates.

## Not yet finalized
- `10_exe.bat` / `S10_rstn_parse.py`: parser still being finished.
- `40_exe.bat` / `40_analysis.py`: final analysis & summary still being finished.

## Quick start (Windows)
1. Install Python 3.12 (64-bit).
2. `pip install -r requirements.txt`
3. Double-click `00_all.bat`  *(skips unfinished steps; won’t error)*

Outputs land in `results\` (figures and snippet are already present).

## License
MIT
