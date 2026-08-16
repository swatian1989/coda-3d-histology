#!/usr/bin/env python
"""Arm C, step 4: measured nuclear diameters and the 2D to 3D count correction.

    python scripts/run_usm_3d_extrapolation.py

WHAT THIS IS, AND WHAT IT IS NOT.

This is NOT a three-dimensional reconstruction. No volume is built, because
building one requires serial sections and these are single fields. The
applicability gate refuses stages 1, 2, 5 and 6 on this dataset for exactly
that reason.

This IS the stereological part of CODA's quantification that a single section
can support. A nucleus is counted whenever any part of it intersects the
section, so the effective sampling depth is the section thickness plus the
nuclear diameter, not the thickness alone. Counting nuclei per unit area of a
section therefore overestimates the number per unit volume of tissue, and
overestimates it most for the largest nuclei, which biases comparisons between
cell populations systematically rather than randomly. The correction is

    C_3D = C_2D * k * T / (T + D)

with T the section thickness and D the nuclear diameter.

TWO PARAMETERS ARE DELIBERATELY NOT TAKEN FROM THE PAPER'S DEFAULTS.

k, the skipped-section factor, is 3 in CODA because every third section was
stained and each stained section stands for three sections of tissue. These
are single fields with no stack, so nothing is being extrapolated across
skipped sections and k is 1 here. Leaving k at 3 would inflate every count
threefold.

D, the nuclear diameter, defaults in the library to CODA's PANCREAS
measurements. The protocol is explicit that these must be measured in the
tissue at hand, because the correction scales counts directly. Diameters are
therefore measured here from the segmented nuclei of this cohort, separately
for marker-positive and marker-negative populations, and reported.

T, the section thickness, is NOT measured from the images and cannot be. The
value used is the paper's 4 um, and every result below scales with it, so it
is reported as an explicit assumption rather than buried.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coda_my.ihc import IHCConfig, score_marker  # noqa: E402
from coda_my.reconstruct import ReconstructionConfig, extrapolate_3d_cell_count  # noqa: E402

SECTION_THICKNESS_UM = 4.0     # CONFIRMED cutting thickness for these blocks
N_DIAMETER_IMAGES = 25         # images sampled to measure diameters
MIN_NUCLEI = 50


def main() -> None:
    logging.basicConfig(level=logging.ERROR)
    warnings.filterwarnings("ignore")

    qc = pd.read_csv(ROOT / "results/usm_qc.csv")
    mk = pd.read_csv(ROOT / "results/usm_markers.csv")

    # ---------------------------------------------------------- measure D
    # Sample the highest-resolution images, where a nuclear boundary is best
    # resolved, and only those with a counterstain so negatives are visible.
    cand = qc[(qc["counterstain_grade"] != "absent")
              & (qc["mpp_um_per_px"] < 0.60)
              & (qc["marker"].isin(["KI67", "ER", "PR"]))]
    cand = cand.sort_values("mpp_um_per_px").head(N_DIAMETER_IMAGES)
    print(f"measuring nuclear diameters on {len(cand)} high-resolution images "
          f"({cand['mpp_um_per_px'].min():.3f} to {cand['mpp_um_per_px'].max():.3f} um/px)")

    pos_d, neg_d, per_image = [], [], []
    marker_map = {"ER": "ER", "PR": "PR", "KI67": "Ki67"}
    for _, r in cand.iterrows():
        path = ROOT / "data/raw/usm" / r["marker"] / r["filename"]
        rgb = np.asarray(Image.open(path).convert("RGB"))
        try:
            nuclei, _ = score_marker(rgb, marker_map[r["marker"]],
                                     IHCConfig(mpp=float(r["mpp_um_per_px"])))
        except Exception:
            continue
        if len(nuclei) < MIN_NUCLEI or "area_um2" not in nuclei:
            continue
        # equivalent circular diameter from the segmented area
        d = 2.0 * np.sqrt(nuclei["area_um2"].to_numpy() / np.pi)
        p = d[nuclei["positive"].to_numpy()]
        n = d[~nuclei["positive"].to_numpy()]
        pos_d.append(p); neg_d.append(n)
        per_image.append({"filename": r["filename"], "marker": r["marker"],
                          "mpp": r["mpp_um_per_px"], "n_nuclei": len(nuclei),
                          "median_d_positive_um": float(np.median(p)) if len(p) else np.nan,
                          "median_d_negative_um": float(np.median(n)) if len(n) else np.nan})

    pos = np.concatenate(pos_d) if pos_d else np.array([])
    neg = np.concatenate(neg_d) if neg_d else np.array([])
    D_pos, D_neg = float(np.median(pos)), float(np.median(neg))
    print(f"\nMEASURED nuclear diameters (equivalent circular, this cohort)")
    print(f"  marker positive : median {D_pos:.2f} um  IQR "
          f"{np.percentile(pos,25):.2f}-{np.percentile(pos,75):.2f}  n={len(pos):,}")
    print(f"  marker negative : median {D_neg:.2f} um  IQR "
          f"{np.percentile(neg,25):.2f}-{np.percentile(neg,75):.2f}  n={len(neg):,}")
    print(f"  library default (CODA pancreas ductal epithelium) was 4.20 um")

    pd.DataFrame(per_image).to_csv(ROOT / "results/usm_nuclear_diameters.csv", index=False)

    # -------------------------------------------------- 2D to 3D correction
    cfg = ReconstructionConfig(
        section_thickness_um=SECTION_THICKNESS_UM,
        sections_skipped=1,                      # single fields; see module docstring
        nuclear_diameters={"positive": D_pos, "negative": D_neg},
    )
    T = cfg.section_thickness_um
    f_pos, f_neg = T / (T + D_pos), T / (T + D_neg)
    print(f"\nAbercrombie correction factors at T={T:.1f} um")
    print(f"  positive nuclei : {f_pos:.3f}  ({100*(1-f_pos):.0f}% of 2D counts are "
          f"section-plane artefact)")
    print(f"  negative nuclei : {f_neg:.3f}")

    rows = []
    nuc = mk[mk["marker"].isin(["ER", "PR", "Ki67"])]
    for _, r in nuc.iterrows():
        if pd.isna(r.get("n_positive")) or pd.isna(r.get("fov_area_mm2")):
            continue
        n_pos = float(r["n_positive"])
        n_tot = float(r.get("n_nuclei_detected") or np.nan)
        n_neg = (n_tot - n_pos) if np.isfinite(n_tot) else np.nan
        area_mm2 = float(r["fov_area_mm2"])
        vol_mm3 = area_mm2 * (T / 1000.0)          # section volume = area x thickness

        c3d = extrapolate_3d_cell_count(
            {"positive": int(n_pos),
             "negative": int(n_neg) if np.isfinite(n_neg) else 0}, cfg)

        rows.append({
            "filename": r["filename"], "marker": r["marker"],
            "counterstain_grade": r["counterstain_grade"],
            "n_positive_2d": n_pos,
            "density_2d_per_mm2": n_pos / area_mm2 if area_mm2 else np.nan,
            "n_positive_3d_corrected": c3d["positive"],
            "density_3d_per_mm3": c3d["positive"] / vol_mm3 if vol_mm3 else np.nan,
            "correction_factor": f_pos,
            "section_thickness_um_assumed": T,
            "n_negative_3d_corrected": (c3d["negative"]
                                        if np.isfinite(n_neg) else np.nan),
        })

    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/usm_3d_extrapolation.csv", index=False)

    print(f"\nwrote results/usm_3d_extrapolation.csv ({len(out)} images)")
    print("\n2D areal density vs 3D volumetric density, by marker")
    for m, g in out.groupby("marker"):
        print(f"  {m:5s} n={len(g):3d}  2D median {g['density_2d_per_mm2'].median():8.0f}/mm2"
              f"   3D median {g['density_3d_per_mm3'].median():10.0f}/mm3")

    meta = {"section_thickness_um_ASSUMED": T,   # key kept for compatibility
            "section_thickness_provenance": "confirmed by the laboratory",
            "sections_skipped_factor": cfg.sections_skipped,
            "measured_diameter_positive_um": round(D_pos, 3),
            "measured_diameter_negative_um": round(D_neg, 3),
            "correction_factor_positive": round(f_pos, 4),
            "correction_factor_negative": round(f_neg, 4),
            "n_images_for_diameter": len(per_image),
            "n_nuclei_measured": int(len(pos) + len(neg)),
            "IS_A_RECONSTRUCTION": False,
            "note": ("stereological correction of 2D counts to volumetric density; "
                     "no volume was built and none can be from single sections")}
    (ROOT / "results/usm_3d_extrapolation_meta.json").write_text(json.dumps(meta, indent=2))
    print("\nwrote results/usm_3d_extrapolation_meta.json")
    print(f"\nSection thickness {T:.1f} um is the confirmed cutting thickness for "
          f"these blocks, not an inherited default. Volumetric densities scale "
          f"linearly with it.")


if __name__ == "__main__":
    main()
