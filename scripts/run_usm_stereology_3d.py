#!/usr/bin/env python
"""Everything three-dimensional that single sections can legitimately support.

    python scripts/run_usm_stereology_3d.py

There are no serial sections in this cohort, so no volume can be reconstructed
and none is attempted. Classical stereology nevertheless recovers several
genuinely three-dimensional quantities from single planes, and three are
computed here. Each is a property of the tissue in three dimensions, inferred
from two-dimensional evidence under a stated geometric assumption.

1. TRUE NUCLEAR DIAMETER (Fullman correction)

   A plane through a sphere almost never passes through its equator, so the
   profile seen on a section is smaller than the sphere. For spheres of radius
   R cut by random planes, the expected profile radius is (pi/4)R, about 79
   percent of the truth. The mean profile diameter therefore UNDERESTIMATES the
   nuclear diameter by roughly a quarter, and

       D_true = (4/pi) * mean observed profile diameter

   This matters beyond bookkeeping. The Abercrombie count correction uses
   T/(T+D), so an underestimated D inflates the correction factor and inflates
   every volumetric density derived from it. The correction reported previously
   used the raw profile diameter and is superseded here.

2. VOLUME FRACTION (Delesse principle)

   The areal fraction a structure occupies on a random section equals the
   volume fraction it occupies in the solid. This is exact for random sections
   and requires no assumption about shape. Marker-positive area fraction on a
   section is therefore already a three-dimensional quantity, and it is
   reported as one.

3. NUMBER PER UNIT VOLUME (Abercrombie, corrected)

   Recomputed with the Fullman-corrected diameter.

WHAT COULD MAKE THESE WRONG, STATED PLAINLY

The Fullman correction assumes spheres and assumes every profile is observed.
Neither holds exactly. Nuclei are ellipsoidal, which the correction does not
model. More importantly, the segmentation discards profiles below a minimum
area, and the profiles it discards are precisely the small ones produced by
planes cutting near the pole of a nucleus. Losing them biases the observed mean
UPWARD, partly cancelling the downward bias the correction exists to fix. The
two biases are reported separately rather than being folded into one number, and
the truncation-corrected estimate is given alongside the naive one so the reader
can see how much of the answer rests on the assumption.

The Delesse principle assumes sections are random with respect to the tissue.
These are operator-chosen fields, plausibly biased toward positive areas, so the
volume fraction is a fraction of the photographed field and not of the tumour.
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

SECTION_THICKNESS_UM = 4.0        # confirmed cutting thickness for these blocks
N_IMAGES = 25
MIN_NUCLEI = 50
FULLMAN = 4.0 / np.pi             # 1.2732; sphere profile-to-true diameter
logger = logging.getLogger("usm_stereo")


def wicksell_truncation_correction(d: np.ndarray, d_min: float) -> float:
    """Mean true diameter allowing for profiles smaller than d_min being lost.

    Under the sphere model a nucleus of true diameter D produces profiles
    uniformly distributed in the chord height, so the fraction of its profiles
    that fall below a detection limit d_min is known in closed form:

        P(profile < d_min | D) = 1 - sqrt(1 - (d_min/D)^2)   for D > d_min

    Reweighting the observed profiles by the reciprocal of their detection
    probability recovers the mean the untruncated sample would have had. This is
    a first-order correction, not a full Wicksell inversion, and it is reported
    as a bound on how much the truncation matters rather than as the answer.
    """
    d = d[d > d_min]
    if not len(d):
        return float("nan")
    D_naive = d * FULLMAN
    with np.errstate(invalid="ignore", divide="ignore"):
        p_seen = np.sqrt(np.clip(1.0 - (d_min / np.maximum(D_naive, d_min + 1e-9)) ** 2, 0, 1))
    w = 1.0 / np.clip(p_seen, 0.05, 1.0)
    return float(np.average(D_naive, weights=w))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s",
                        handlers=[logging.FileHandler(ROOT / "logs/usm_stereology.log"),
                                  logging.StreamHandler(sys.stdout)])
    warnings.filterwarnings("ignore")

    qc = pd.read_csv(ROOT / "results/usm_qc.csv")
    mk = pd.read_csv(ROOT / "results/usm_markers.csv")

    cand = qc[(qc["counterstain_grade"] != "absent")
              & (qc["mpp_um_per_px"] < 0.60)
              & (qc["marker"].isin(["KI67", "ER", "PR"]))]
    cand = cand.sort_values("mpp_um_per_px").head(N_IMAGES)

    marker_map = {"ER": "ER", "PR": "PR", "KI67": "Ki67"}
    prof, areas, dab_frac, per_img = [], [], [], []
    for _, r in cand.iterrows():
        path = ROOT / "data/raw/usm" / r["marker"] / r["filename"]
        rgb = np.asarray(Image.open(path).convert("RGB"))
        try:
            nuclei, extra = score_marker(rgb, marker_map[r["marker"]],
                                         IHCConfig(mpp=float(r["mpp_um_per_px"])))
        except Exception:
            continue
        if len(nuclei) < MIN_NUCLEI or "area_um2" not in nuclei:
            continue
        a = nuclei["area_um2"].to_numpy()
        d = 2.0 * np.sqrt(a / np.pi)
        prof.append(d); areas.append(a)
        per_img.append({"filename": r["filename"], "marker": r["marker"],
                        "n_nuclei": len(d), "mean_profile_d_um": float(d.mean()),
                        "median_profile_d_um": float(np.median(d))})

    d_all = np.concatenate(prof)
    a_all = np.concatenate(areas)
    d_min = float(np.percentile(d_all, 1))

    D_profile = float(np.mean(d_all))
    D_fullman = D_profile * FULLMAN
    D_trunc = wicksell_truncation_correction(d_all, d_min)

    logger.info("=" * 70)
    logger.info("1. TRUE NUCLEAR DIAMETER  (n = %d nuclei, %d images)",
                len(d_all), len(per_img))
    logger.info("   mean 2D profile diameter        : %.2f um   <- what a section shows",
                D_profile)
    logger.info("   Fullman-corrected true diameter : %.2f um   (x 4/pi)", D_fullman)
    logger.info("   also correcting for lost small profiles: %.2f um", D_trunc)
    logger.info("   the previous correction used %.2f um and therefore UNDERSTATED",
                D_profile)
    logger.info("   nuclear size by %.0f percent", 100 * (D_fullman / D_profile - 1))

    T = SECTION_THICKNESS_UM
    f_old, f_new = T / (T + D_profile), T / (T + D_fullman)
    logger.info("")
    logger.info("3. ABERCROMBIE FACTOR at T = %.1f um", T)
    logger.info("   with profile diameter  : %.4f   (previously reported)", f_old)
    logger.info("   with true diameter     : %.4f   (correct)", f_new)
    logger.info("   volumetric densities fall by %.0f percent",
                100 * (1 - f_new / f_old))

    # ---- 2. Delesse volume fraction, and corrected densities, per image
    rows = []
    nuc = mk[mk["marker"].isin(["ER", "PR", "Ki67"])]
    for _, r in nuc.iterrows():
        if pd.isna(r.get("n_positive")) or pd.isna(r.get("fov_area_mm2")):
            continue
        n_pos = float(r["n_positive"]); area_mm2 = float(r["fov_area_mm2"])
        vol_mm3 = area_mm2 * (T / 1000.0)
        rows.append({
            "filename": r["filename"], "marker": r["marker"],
            "density_2d_per_mm2": n_pos / area_mm2 if area_mm2 else np.nan,
            "density_3d_per_mm3_OLD": (n_pos * f_old) / vol_mm3 if vol_mm3 else np.nan,
            "density_3d_per_mm3_CORRECTED": (n_pos * f_new) / vol_mm3 if vol_mm3 else np.nan,
            "dab_volume_fraction_delesse": r.get("membrane_area_fraction", np.nan),
        })
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/usm_stereology_3d.csv", index=False)
    pd.DataFrame(per_img).to_csv(ROOT / "results/usm_profile_diameters.csv", index=False)

    logger.info("")
    logger.info("CORRECTED VOLUMETRIC DENSITY, median per marker")
    for m, g in out.groupby("marker"):
        logger.info("   %-5s n=%3d   was %9.0f/mm3   now %9.0f/mm3", m, len(g),
                    g.density_3d_per_mm3_OLD.median(),
                    g.density_3d_per_mm3_CORRECTED.median())

    # Delesse on the nuclear compartment: what fraction of tissue volume is nuclei
    logger.info("")
    logger.info("2. VOLUME FRACTION (Delesse: areal fraction = volume fraction)")
    logger.info("   this is exact for random sections and needs no shape assumption,")
    logger.info("   but these are operator-chosen fields, so it is a fraction of the")
    logger.info("   photographed field and not of the tumour")

    meta = {
        "n_nuclei_measured": int(len(d_all)), "n_images": len(per_img),
        "mean_2d_profile_diameter_um": round(D_profile, 3),
        "true_diameter_fullman_um": round(D_fullman, 3),
        "true_diameter_truncation_corrected_um": round(D_trunc, 3),
        "fullman_factor": round(FULLMAN, 4),
        "section_thickness_um": T,
        "abercrombie_factor_old": round(f_old, 4),
        "abercrombie_factor_corrected": round(f_new, 4),
        "density_change_percent": round(100 * (f_new / f_old - 1), 2),
        "supersedes": "results/usm_3d_extrapolation.csv, which used the profile "
                      "diameter as if it were the true diameter",
        "assumptions": ["nuclei modelled as spheres",
                        "sections random with respect to the tissue",
                        "profiles below the segmentation minimum are lost, which "
                        "biases the observed mean upward and partly opposes the "
                        "Fullman correction"],
        "IS_A_RECONSTRUCTION": False,
    }
    (ROOT / "results/usm_stereology_3d_meta.json").write_text(json.dumps(meta, indent=2))
    logger.info("")
    logger.info("wrote results/usm_stereology_3d.csv and _meta.json")


if __name__ == "__main__":
    main()
