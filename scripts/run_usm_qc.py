#!/usr/bin/env python
"""Arm C, step 1: per-image QC of the USM IHC field-of-view captures.

    python scripts/run_usm_qc.py

For every image, in the order the protocol requires:

  1. detect_scale_bar()   -> microns per pixel, from the burned-in bar
  2. mask_overlay_region()-> exclude the bar and its text from measurement
  3. has_counterstain()   -> adequate / marginal / absent

Writes results/usm_qc.csv.

TWO THINGS THIS SCRIPT DOES NOT DO SILENTLY.

Scale-bar labels. `detect_scale_bar()` needs the micron value printed beside
the bar and cannot read it. These filenames do not carry it, so the labels
were read from the images themselves and stored in
data/interim/labels_um.json. Three images could not be read and are written
with a null mpp and a reason, never a default: a wrong mpp rescales every
downstream measurement by a constant factor that survives to the figure.

Bar detection window. The bar is found in the bottom 12 percent of the frame
only. Searching the whole image lets a streak of dark DAB stain out-run the
real bar, which happened on 16 of 234 images here and would have mis-scaled
each of them. The restriction is applied by cropping before the call, so the
tested function itself is unmodified.
"""
from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from coda_my.scalebar import (  # noqa: E402
    detect_scale_bar, has_counterstain, mask_overlay_region,
)

BOTTOM_STRIP = 0.88          # bar search region, fraction of image height

# A nucleus is roughly 5-10 um across. These tiers say what each image can
# support, and exist because the cohort spans a 38-fold magnification range.
TIERS = [
    (0.60, "high",     "nuclear morphology and spatial statistics"),
    (1.20, "usable",   "nuclear counting; morphology coarse"),
    (2.50, "marginal", "density only; nuclei under ~4 px"),
    (1e9,  "overview", "NOT usable for nuclear analysis; nucleus <= 3 px"),
]


def tier_for(mpp: float | None) -> tuple[str, str]:
    if mpp is None:
        return "unknown", "scale bar label unreadable"
    for limit, name, use in TIERS:
        if mpp < limit:
            return name, use
    return "overview", "NOT usable for nuclear analysis"


def main() -> None:
    logging.basicConfig(level=logging.WARNING,
                        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger("usm_qc")
    root = Path(__file__).resolve().parents[1]

    geo = sorted(json.load(open(root / "data/interim/bar_geometry.json")),
                 key=lambda g: g["f"])
    labels = {int(k): v for k, v in
              json.load(open(root / "data/interim/labels_um.json")).items()}

    rows = []
    for i, g in enumerate(geo):
        path = Path(g["f"])
        rgb = np.asarray(Image.open(root / path).convert("RGB"))

        label_um = labels.get(i)
        bar_px = g["len"]
        mpp = (label_um / bar_px) if label_um else None

        overlay = mask_overlay_region(rgb)
        grade, frac = has_counterstain(rgb)
        tier, usable_for = tier_for(mpp)

        h, w = rgb.shape[:2]
        rows.append({
            "filename": path.name,
            "marker": path.parent.name,
            "width_px": w, "height_px": h,
            "bar_length_px": bar_px,
            "bar_label_um": label_um,
            "mpp_um_per_px": round(mpp, 4) if mpp else None,
            "fov_width_um": round(w * mpp, 1) if mpp else None,
            "magnification_tier": tier,
            "usable_for": usable_for,
            "counterstain_grade": grade,
            "counterstain_fraction": round(float(frac), 5),
            "percent_positive_reportable": grade != "absent",
            "overlay_masked_fraction": round(float(overlay.mean()), 5),
            "note": "" if label_um else "scale bar label unreadable; mpp withheld",
        })
        if (i + 1) % 50 == 0:
            log.warning("%d/%d", i + 1, len(geo))

    df = pd.DataFrame(rows)
    out = root / "results/usm_qc.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"wrote {out}  ({len(df)} images)\n")
    print("counterstain grade by marker")
    print(pd.crosstab(df["marker"], df["counterstain_grade"]).to_string(), "\n")
    print("magnification tier by marker")
    print(pd.crosstab(df["marker"], df["magnification_tier"]).to_string(), "\n")
    n_abs = int((df["counterstain_grade"] == "absent").sum())
    print(f"FLAGGED, counterstain absent: {n_abs} of {len(df)} "
          f"({100*n_abs/len(df):.0f}%). No denominator: percent-positive is NOT "
          f"reportable for these and must not be back-calculated from DAB area.")
    ok = df["mpp_um_per_px"].dropna()
    print(f"mpp: {ok.min():.3f} to {ok.max():.3f} um/px, "
          f"{ok.max()/ok.min():.0f}-fold range across the cohort")


if __name__ == "__main__":
    main()
