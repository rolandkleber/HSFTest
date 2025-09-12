# HSF009 · Real Echo EF — Results

## Theory predictions (HSF-consistent)
- EF is computed with the **single‑plane area–length identity**: $V=\tfrac{8}{3\pi}\tfrac{A^2}{L}$; **EF=100·(EDV−ESV)/EDV**. No fitting.
- With **ED/ES only**, EF is **independent of frame rate** (subsampling cannot change extrema).
- A **uniform scale** in area or length rescales volumes but **cancels in EF**; absolute volumes may look high with a major‑axis length proxy.
- Expect EF vs FPS to be **flat above ~15–20 fps**; with ED/ES only it is flat at all fps.

## Per‑series results (single‑view 2CH)
| series | EF (%) | EDV (mL) | ESV (mL) | note |
|---|---:|---:|---:|---|
| patient10_2CH | 31.46 | 320.0 | 220.0 | ED/ES only; single view (2CH) |
| patient11_2CH | 17.60 | 408.0 | 337.0 | ED/ES only; single view (2CH) |
| patient12_2CH | 24.93 | 424.0 | 318.0 | ED/ES only; single view (2CH) |
| patient13_2CH | 9.00 | 331.0 | 302.0 | ED/ES only; single view (2CH) |
| patient14_2CH | 23.07 | 233.6 | 179.7 | ED/ES only; single view (2CH) |
| patient15_2CH | 39.61 | 334.6 | 202.1 | ED/ES only; single view (2CH) |

## Remarks
- Volumes are **single‑view estimates** using a **major‑axis** length proxy; EF is robust but volumes can be inflated.
- For more physiologic volumes (and often a different EF), either measure **apex→annulus length** or combine **2CH + 4CH** (biplane area–length).
- Inputs used: `data/patientXX_2CH/lv_trace.csv` (two rows: ED & ES).
