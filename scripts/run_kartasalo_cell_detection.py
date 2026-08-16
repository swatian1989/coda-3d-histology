#!/usr/bin/env python
"""CODA stage 3, partial: two independent nucleus detectors, compared to each other.

    python scripts/run_kartasalo_cell_detection.py

WHAT THIS CAN AND CANNOT ESTABLISH

The published validation is precision and recall against two human annotators
at a 2 um matching tolerance, with a 90 percent acceptance bar. No manual
annotation exists for this material and none can be manufactured, so precision
and recall against human ground truth are NOT reported here and the 90 percent
bar is not tested. Writing a number against an imagined annotator would be worse
than reporting nothing.

What can be established without human labels is method-to-method agreement. Two
detectors that share no code path are run on the same fields and matched to each
other. Agreement is a necessary condition for correctness and not a sufficient
one: two detectors can agree and both be wrong, most obviously by both missing
the same faint nuclei. It is reported as agreement, never as accuracy.

THE TWO DETECTORS

  A  threshold and watershed. Haematoxylin optical density by colour
     deconvolution, Otsu threshold, distance transform, watershed split of
     touching nuclei. This is the classical route and the one the project's own
     cell detection follows.

  B  Laplacian-of-Gaussian blob detection. Multi-scale second-derivative
     response with non-maximum suppression, which finds nuclei as local scale
     space extrema rather than by segmenting a foreground mask.

They share the haematoxylin channel and nothing else: one segments regions and
splits them, the other detects points at a scale. A false positive from a stain
artefact would have to fool both mechanisms to survive the comparison.

Matching uses the Hungarian algorithm at a 2 um tolerance, which is the
tolerance the published validation used, so the agreement figure is on the same
footing as the published accuracy figure even though it measures a different
thing.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coda_my.deconv import deconvolve  # noqa: E402
from coda_my.loaders.kartasalo import NATIVE_MPP_UM, section_paths  # noqa: E402

DATA = ROOT / "data/raw/kartasalo/extracted/Data_to_IDA"
OUT = ROOT / "results/kartasalo"
MATCH_TOL_UM = 2.0          # [PAPER] the tolerance used for the published validation
logger = logging.getLogger("kartasalo_cells")


def haematoxylin(rgb: np.ndarray) -> np.ndarray:
    ch = deconvolve(np.clip(rgb, 0, 255).astype(np.uint8))
    h = ch["hematoxylin"]
    lo, hi = np.percentile(h, 1), np.percentile(h, 99.5)
    return np.clip((h - lo) / max(hi - lo, 1e-6), 0, 1)


def detect_watershed(h: np.ndarray, mpp: float,
                     min_um: float = 2.5, max_um: float = 15.0) -> np.ndarray:
    """Detector A: Otsu threshold, distance transform, watershed."""
    from skimage.feature import peak_local_max
    from skimage.filters import threshold_otsu
    from skimage.segmentation import watershed

    sm = ndimage.gaussian_filter(h, max(0.5 / mpp, 0.6))
    try:
        thr = threshold_otsu(sm)
    except ValueError:
        return np.empty((0, 2))
    mask = sm > thr
    mask = ndimage.binary_opening(mask, np.ones((3, 3)))
    if mask.sum() < 10:
        return np.empty((0, 2))
    dist = ndimage.distance_transform_edt(mask)
    min_d = max(int(round(min_um / mpp)), 2)
    peaks = peak_local_max(dist, min_distance=min_d, labels=mask)
    if not len(peaks):
        return np.empty((0, 2))
    markers = np.zeros(h.shape, dtype=int)
    for i, (y, x) in enumerate(peaks, start=1):
        markers[y, x] = i
    lab = watershed(-dist, markers, mask=mask)
    min_px = np.pi * (min_um / 2 / mpp) ** 2
    max_px = np.pi * (max_um / 2 / mpp) ** 2
    out = []
    for idx, sl in enumerate(ndimage.find_objects(lab), start=1):
        if sl is None:
            continue
        m = lab[sl] == idx
        a = int(m.sum())
        if min_px <= a <= max_px:
            ys, xs = np.nonzero(m)
            out.append([ys.mean() + sl[0].start, xs.mean() + sl[1].start])
    return np.array(out) if out else np.empty((0, 2))


def detect_log(h: np.ndarray, mpp: float,
               min_um: float = 2.5, max_um: float = 15.0) -> np.ndarray:
    """Detector B: multi-scale Laplacian-of-Gaussian blob detection."""
    from skimage.feature import blob_log

    lo = max(min_um / 2 / mpp, 1.0)
    hi = max(max_um / 2 / mpp, lo + 1.0)
    blobs = blob_log(h, min_sigma=lo, max_sigma=hi, num_sigma=6,
                     threshold=0.06, overlap=0.4)
    return blobs[:, :2] if len(blobs) else np.empty((0, 2))


def match(a: np.ndarray, b: np.ndarray, tol_px: float) -> tuple[int, float]:
    """Hungarian matching within a tolerance. Returns (n_matched, median distance)."""
    if not len(a) or not len(b):
        return 0, float("nan")
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    big = tol_px * 1000.0
    cost = np.where(d <= tol_px, d, big)
    ri, ci = linear_sum_assignment(cost)
    keep = d[ri, ci] <= tol_px
    return int(keep.sum()), float(np.median(d[ri, ci][keep])) if keep.any() else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fields", type=int, default=12)
    ap.add_argument("--field-px", type=int, default=900,
                    help="field side at native resolution")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.FileHandler(ROOT / "logs/kartasalo_cells.log"),
                                  logging.StreamHandler(sys.stdout)])
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None

    mpp = NATIVE_MPP_UM                    # native resolution: nuclei need it
    tol_px = MATCH_TOL_UM / mpp
    paths = section_paths(DATA / "liver")
    picks = np.linspace(0, len(paths) - 1, args.n_fields).astype(int)
    logger.info("stage 3 partial: %d fields of %d px (%.0f um) at %.2f um/px, "
                "match tolerance %.1f um = %.1f px",
                len(picks), args.field_px, args.field_px * mpp, mpp,
                MATCH_TOL_UM, tol_px)

    rows = []
    for k, i in enumerate(picks):
        with Image.open(paths[i]) as im:
            W, H = im.size
            cy, cx = H // 2, W // 2
            s = args.field_px
            crop = np.asarray(im.crop((cx - s // 2, cy - s // 2,
                                       cx + s // 2, cy + s // 2)).convert("RGB"))
        h = haematoxylin(crop)
        if h.std() < 0.02:
            logger.info("  section %d: field is blank, skipped", i + 1)
            continue
        A = detect_watershed(h, mpp)
        B = detect_log(h, mpp)
        n_m, med = match(A, B, tol_px)
        area_mm2 = (args.field_px * mpp / 1000.0) ** 2
        rows.append({"section": int(i + 1), "n_watershed": len(A), "n_log": len(B),
                     "n_matched": n_m,
                     "agreement_vs_watershed": n_m / max(len(A), 1),
                     "agreement_vs_log": n_m / max(len(B), 1),
                     "median_offset_um": med * mpp if np.isfinite(med) else np.nan,
                     "density_watershed_per_mm2": len(A) / area_mm2,
                     "density_log_per_mm2": len(B) / area_mm2})
        logger.info("  section %2d: watershed %4d, LoG %4d, matched %4d "
                    "(%.0f%% / %.0f%%)", i + 1, len(A), len(B), n_m,
                    100 * rows[-1]["agreement_vs_watershed"],
                    100 * rows[-1]["agreement_vs_log"])

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "stage3_cell_detection.csv", index=False)
    if df.empty:
        logger.warning("no usable fields")
        return

    f1 = 2 * df.n_matched.sum() / max(df.n_watershed.sum() + df.n_log.sum(), 1)
    logger.info("")
    logger.info("STAGE 3  method-to-method agreement over %d fields", len(df))
    logger.info("  watershed detected %d nuclei, LoG detected %d, matched %d",
                df.n_watershed.sum(), df.n_log.sum(), df.n_matched.sum())
    logger.info("  agreement F1 (symmetric)     : %.3f", f1)
    logger.info("  matched fraction of watershed: %.3f", df.agreement_vs_watershed.mean())
    logger.info("  matched fraction of LoG      : %.3f", df.agreement_vs_log.mean())
    logger.info("  median centroid offset       : %.2f um (tolerance %.1f)",
                df.median_offset_um.median(), MATCH_TOL_UM)
    logger.info("  density, watershed           : %.0f per mm2 (median)",
                df.density_watershed_per_mm2.median())
    logger.info("  density, LoG                 : %.0f per mm2 (median)",
                df.density_log_per_mm2.median())

    summary = {
        "n_fields": int(len(df)), "mpp_um": mpp, "field_um": args.field_px * mpp,
        "match_tolerance_um": MATCH_TOL_UM,
        "n_watershed": int(df.n_watershed.sum()), "n_log": int(df.n_log.sum()),
        "n_matched": int(df.n_matched.sum()),
        "agreement_f1": float(f1),
        "median_offset_um": float(df.median_offset_um.median()),
        "density_watershed_per_mm2": float(df.density_watershed_per_mm2.median()),
        "density_log_per_mm2": float(df.density_log_per_mm2.median()),
        "IS_ACCURACY": False,
        "note": ("method-to-method agreement between two detectors sharing no code "
                 "path beyond the haematoxylin channel. NOT precision or recall: no "
                 "human annotation exists for this material, so the published 90 "
                 "percent acceptance bar is not tested and is not claimed."),
    }
    (OUT / "stage3_summary.json").write_text(json.dumps(summary, indent=2))
    logger.info("wrote %s", OUT / "stage3_summary.json")


if __name__ == "__main__":
    main()
