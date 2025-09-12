#!/usr/bin/env python3
import os, glob, argparse, sys
import numpy as np, pandas as pd
try:
    import nibabel as nib
    from skimage.measure import regionprops
except Exception:
    print("Install deps:  pip install nibabel scikit-image", file=sys.stderr)
    raise

def area_cm2(mask, sx_mm, sy_mm):
    return float(mask.astype(bool).sum() * sx_mm * sy_mm * 0.01)

def long_axis_cm(mask, sx_mm, sy_mm):
    props = regionprops(mask.astype(np.uint8))
    if not props: return 0.0
    s_mm = (sx_mm * sy_mm) ** 0.5
    return float(props[0].major_axis_length * s_mm / 10.0)

def one_row(img_path, msk_path, sid, t):
    img, msk = nib.load(img_path), nib.load(msk_path)
    sx_mm, sy_mm = map(float, img.header.get_zooms()[:2])
    arr = np.asanyarray(msk.dataobj)
    if arr.ndim == 3: arr = arr[..., 0]
    arr = arr > 0
    return dict(time_s=float(t),
                area_cm2=area_cm2(arr, sx_mm, sy_mm),
                length_cm=long_axis_cm(arr, sx_mm, sy_mm),
                series_id=sid)

def find4(pdir, view):
    """Return ED_img, ES_img, ED_gt, ES_gt using glob patterns."""
    def g(pat): 
        hits = sorted(glob.glob(os.path.join(pdir, pat)))
        return hits[0] if hits else None
    ed_img = g(f"*_{view}_ED.nii.gz")
    es_img = g(f"*_{view}_ES.nii.gz")
    ed_gt  = g(f"*_{view}_ED_gt.nii.gz")
    es_gt  = g(f"*_{view}_ES_gt.nii.gz")
    return ed_img, es_img, ed_gt, es_gt

def convert_view(pdir, pat, view, out_root):
    ed_img, es_img, ed_gt, es_gt = find4(pdir, view)
    if not all([ed_img, es_img, ed_gt, es_gt]):
        missing = ["ED_img","ES_img","ED_gt","ES_gt"][ [ed_img,es_img,ed_gt,es_gt].index(None) ]
        return False, f"missing {missing} in {pdir}"
    # Series ID = folder name + view
    series_id = f"{pat}_{view}"
    out_dir = os.path.join(out_root, series_id)
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, "lv_trace.csv")
    rows = [one_row(ed_img, ed_gt, series_id, 0.0),
            one_row(es_img, es_gt, series_id, 1.0)]
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    return True, out_csv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)   # ...\data\raw
    ap.add_argument("--out", required=True)   # ...\data
    a = ap.parse_args()

    patients = sorted([p for p in glob.glob(os.path.join(a.raw, "patient*")) if os.path.isdir(p)])
    print(f"Found {len(patients)} patient folders.")
    done = 0
    for pdir in patients:
        pat = os.path.basename(pdir)
        for view in ("2CH", "4CH"):
            ok, info = convert_view(pdir, pat, view, a.out)
            print(("[OK] " if ok else "[SKIP] ") + f"{pat} {view} -> {info}")
            if ok: done += 1
    print(f"Finished. Wrote {done} series CSV(s).")

if __name__ == "__main__":
    main()
