# HSF009_EchokardScreening (Alpha, Windows)

**Was drin ist:** Beispiel-Serien (`data/patientXX_2CH/lv_trace.csv`) und die erzeugten Resultate (`results/`).
**Ziel:** EF aus Fläche & Länge (Single-Plane Area–Length).

## Schnellstart
1) `pip install -r requirements.txt`
2) `cd scripts` → `10_exe.bat`
   - liest `data\*\lv_trace.csv`
   - schreibt je Serie in `results\<serie>\*`

> Hinweis: Die NIfTI-Extraktion (`05_extractdata.py`) ist optional und benötigt `nibabel` & `scikit-image`.
