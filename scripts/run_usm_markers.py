#!/usr/bin/env python
"""Arm C, step 2: marker quantification on the USM IHC captures.

    python scripts/run_usm_markers.py

Reads results/usm_qc.csv so every image is measured at ITS OWN mpp and under
the gates that image passed. Routing follows the protocol exactly:

    ER, PR, Ki67 -> score_marker(), hotspot_vs_average(), to_point_pattern()
    HER2         -> membrane_completeness() only

HER2 is never sent to score_marker(). It is membranous; per-nucleus DAB
scoring produces a confident meaningless number, and the library raises on it
by design.

THREE GATES, APPLIED NOT ASSUMED.

Counterstain. Where has_counterstain() graded "absent" there are no visible
negative nuclei, so there is no denominator. percent_positive is written as
null for those images and `percent_reportable` is False. Positive density per
mm2 and spatial arrangement are still valid and are still computed. Nothing is
back-calculated from DAB area.

Magnification. Images coarser than 2.5 um/px cannot resolve a nucleus and are
skipped for nuclear analysis with the reason recorded.

Scale bar. Images whose bar label could not be read carry no mpp, so no
measurement in microns is possible. They are skipped rather than defaulted.
"""
from __future__ import annotations

import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from coda_my.her2 import HER2Config, membrane_completeness  # noqa: E402
from coda_my.ihc import (  # noqa: E402
    IHCConfig, hotspot_vs_average, score_marker, to_point_pattern,
)

MAX_MPP_FOR_NUCLEI = 2.5      # coarser than this, a nucleus is under ~3 px
HOTSPOT_WINDOW_UM = 500.0     # the reporting window a pathologist would use

# Folder names to the library's marker vocabulary, which is case sensitive:
# ihc.MARKERS is ('ER', 'PR', 'HER2', 'Ki67').
MARKER_MAP = {"ER": "ER", "PR": "PR", "HER2": "HER2", "KI67": "Ki67"}


def main() -> None:
    logging.basicConfig(level=logging.ERROR,
                        format="%(asctime)s %(levelname)s %(message)s")
    warnings.filterwarnings("ignore")
    root = Path(__file__).resolve().parents[1]

    qc = pd.read_csv(root / "results/usm_qc.csv")
    rows, points, skipped = [], [], []

    for i, r in qc.iterrows():
        folder = r["marker"]
        marker = MARKER_MAP[folder]
        path = root / "data/raw/usm" / folder / r["filename"]
        mpp = r["mpp_um_per_px"]

        if pd.isna(mpp):
            skipped.append({**r[["filename", "marker"]].to_dict(),
                            "reason": "scale bar label unreadable, no mpp"})
            continue
        if mpp > MAX_MPP_FOR_NUCLEI and marker != "HER2":
            skipped.append({**r[["filename", "marker"]].to_dict(),
                            "reason": f"mpp {mpp:.2f} too coarse for nuclei"})
            continue

        rgb = np.asarray(Image.open(path).convert("RGB"))
        rec = {"filename": r["filename"], "marker": marker,
               "mpp_um_per_px": mpp,
               "counterstain_grade": r["counterstain_grade"],
               "magnification_tier": r["magnification_tier"]}

        try:
            if marker == "HER2":
                res = membrane_completeness(rgb, HER2Config(mpp=float(mpp)))
                rec.update({
                    "membrane_area_fraction": res["membrane_area_fraction"],
                    "n_enclosed_cells": res["n_enclosed_cells"],
                    "mean_membrane_completeness": res["mean_completeness"],
                    "frac_complete_gt80": res["frac_complete_gt80"],
                    "median_cell_area_um2": res["median_cell_area_um2"],
                    "chicken_wire_index": res["chicken_wire_index"],
                    "percent_reportable": False,
                    "note": "membrane completeness only; never 0/1+/2+/3+",
                })
            else:
                cfg = IHCConfig(mpp=float(mpp))
                nuclei, summary = score_marker(rgb, marker, cfg)
                reportable = r["counterstain_grade"] != "absent"
                area_mm2 = (rgb.shape[0] * rgb.shape[1]) * (mpp ** 2) / 1e6
                rec.update({
                    "n_nuclei_detected": summary.get("n_nuclei"),
                    "n_positive": summary.get("n_positive"),
                    "positive_density_per_mm2": (summary.get("n_positive", 0)
                                                 / area_mm2) if area_mm2 else None,
                    "fov_area_mm2": round(area_mm2, 4),
                    "dab_threshold_used": summary.get("dab_threshold_used"),
                    "percent_positive": (summary.get("percent_positive")
                                         if reportable else None),
                    "percent_reportable": bool(reportable),
                    "note": "" if reportable else
                            "counterstain absent: no denominator, percent withheld",
                })
                if marker == "Ki67" and len(nuclei):
                    hs = hotspot_vs_average(nuclei, window_um=HOTSPOT_WINDOW_UM)
                    rec.update({f"ki67_{k}": v for k, v in hs.items()})
                if reportable and len(nuclei):
                    pp = to_point_pattern(nuclei)
                    pp["filename"] = r["filename"]
                    pp["marker"] = marker
                    points.append(pp)
        except Exception as exc:
            rec["error"] = f"{type(exc).__name__}: {exc}"
        rows.append(rec)

        if (i + 1) % 40 == 0:
            print(f"  {i+1}/{len(qc)}", flush=True)

    out = root / "results"
    df = pd.DataFrame(rows)
    df.to_csv(out / "usm_markers.csv", index=False)
    pd.DataFrame(skipped).to_csv(out / "usm_skipped.csv", index=False)
    if points:
        allp = pd.concat(points, ignore_index=True)
        allp.to_parquet(root / "data/interim/usm_point_patterns.parquet", index=False)

    print(f"\nwrote results/usm_markers.csv  ({len(df)} images analysed)")
    print(f"wrote results/usm_skipped.csv  ({len(skipped)} skipped)")
    if points:
        print(f"wrote data/interim/usm_point_patterns.parquet "
              f"({len(allp):,} nuclei from {allp['filename'].nunique()} images)")

    if "error" in df:
        n_err = int(df["error"].notna().sum())
        if n_err:
            print(f"\n{n_err} images errored:")
            print(df.loc[df["error"].notna(), ["filename", "marker", "error"]]
                    .head(8).to_string(index=False))

    print("\nanalysed per marker:")
    print(df["marker"].value_counts().to_string())
    rep = df[df.get("percent_reportable", False) == True]  # noqa: E712
    print(f"\npercent-positive reportable on {len(rep)} of {len(df)} images")
    if "ki67_hotspot_minus_average" in df:
        k = df[df["marker"] == "Ki67"]["ki67_hotspot_minus_average"].dropna()
        if len(k):
            print(f"\nKi67 hotspot minus average, n={len(k)} images:")
            print(f"  median {k.median():.1f} pp, IQR {k.quantile(.25):.1f} to "
                  f"{k.quantile(.75):.1f}, max {k.max():.1f}")


if __name__ == "__main__":
    main()
