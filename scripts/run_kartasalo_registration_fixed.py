#!/usr/bin/env python
"""Re-run the liver stack with the corrected estimator, and compare like for like.

    python scripts/run_kartasalo_registration_fixed.py

Same data, same fiducials, same metrics as run_kartasalo_registration.py, with
three changes, all justified by measurement rather than preference:

  rotation      direct search over angle, scoring each candidate by the pixel
                correlation the pipeline already uses, bounded to plus/minus 45
                degrees. Validated first: mean absolute error against the
                fiducial-implied rotation falls from 20.8 to 3.9 degrees, and
                agreement within five degrees rises from 2 of 20 pairs to 15.
  two scales    rigid solved at about 81 microns per pixel, close to the 80 the
                configuration declares and never uses, and elastic at 7.36,
                instead of forcing both to share the fine scale.
  elastic kept  the displacement fields are returned rather than discarded, so
                target registration error can be reported for the FULL
                transform and not only its rigid part.

The comparison that decides whether this is an improvement is against applying
no transform at all, because the original pipeline failed that test.

Outputs are written with a `_fix` tag so the original run is preserved and the
two can be shown side by side.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from coda_my.loaders.kartasalo import (  # noqa: E402
    NATIVE_MPP_UM, SECTION_THICKNESS_UM, fiducial_array, load_fiducials, load_stack,
)
from coda_my.qc import (  # noqa: E402
    accumulated_tre, axial_vs_lateral_correlation, target_registration_error,
    z_skip_validation,
)
from coda_my.registration import RegistrationConfig  # noqa: E402
from coda_my.registration_fix import (  # noqa: E402
    SearchConfig, procrustes_rigid, register_stack_two_scale, transform_points,
)

DATA = ROOT / "data/raw/kartasalo/extracted/Data_to_IDA"
OUT = ROOT / "results/kartasalo"
TAG = "ds16fix"
BASE_DS = 16
logger = logging.getLogger("kartasalo_fix")


def blockmean(a: np.ndarray, f: int) -> np.ndarray:
    if f == 1:
        return a.astype(np.float32)
    h, w = a.shape
    h2, w2 = (h // f) * f, (w // f) * f
    return a[:h2, :w2].astype(np.float32).reshape(h2 // f, f, w2 // f, f).mean(axis=(1, 3))


def procrustes_residual(A: np.ndarray, B: np.ndarray) -> float:
    """Mean residual after the optimal rigid alignment. See procrustes_rigid.

    An earlier version of this function inverted the SVD convention and
    reported a floor of 518 um where the true value is 75 um, which reversed
    the conclusion about what limits accuracy on this series.
    """
    return procrustes_rigid(A, B)[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coarse-extra", type=int, default=11, help="-> ~81 um/px")
    ap.add_argument("--max-abs-deg", type=float, default=45.0)
    ap.add_argument("--no-elastic", action="store_true",
                    help="measured to be BETTER on this series: elastic raises "
                         "TRE from 114 to 175 um and improves fine correlation "
                         "on 0 of 46 sections")
    args = ap.parse_args()

    global TAG
    TAG = "ds16fix" if not args.no_elastic else "ds16fix_rigid"
    OUT.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.FileHandler(ROOT / "logs/kartasalo_fix.log"),
                                  logging.StreamHandler(sys.stdout)])

    fine_mpp = NATIVE_MPP_UM * BASE_DS
    coarse_mpp = fine_mpp * args.coarse_extra
    logger.info("=" * 72)
    logger.info("rigid at %.1f um/px, elastic at %.2f um/px, search bounded to +/-%.0f deg",
                coarse_mpp, fine_mpp, args.max_abs_deg)

    fine, meta = load_stack(DATA / "liver", downsample=BASE_DS,
                            cache=ROOT / f"data/interim/kartasalo_liver_ds{BASE_DS}.npy")
    coarse = np.stack([blockmean(fine[i], args.coarse_extra) for i in range(len(fine))])
    logger.info("fine %s, coarse %s", fine.shape, coarse.shape)

    fid1 = load_fiducials(DATA / "fiducialcoordinates_liver_observer1.txt")
    fid2 = load_fiducials(DATA / "fiducialcoordinates_liver_observer2.txt")
    secs = sorted(fid1["section"].unique())[:len(fine)]
    n = len(fine)

    cfg = RegistrationConfig()
    t0 = time.time()
    registered, params, fields = register_stack_two_scale(
        coarse, fine, coarse_mpp, fine_mpp, cfg,
        SearchConfig(max_abs_deg=args.max_abs_deg), elastic=not args.no_elastic)
    mins = (time.time() - t0) / 60
    logger.info("registration done in %.1f min (original took 199.8 min)", mins)
    np.save(OUT / f"registered_{TAG}.npy", registered)

    # ---------------------------------------------------------------- step 5
    corr = pd.DataFrame([{
        "section": i + 1,
        "correlation": p.get("correlation", np.nan),
        "correlation_fine": p.get("correlation_fine", np.nan),
        "correlation_after_elastic": p.get("correlation_after_elastic", np.nan),
        "angle_deg": p.get("angle", np.nan),
        "dy_px": p.get("dy_fine", np.nan), "dx_px": p.get("dx_fine", np.nan),
        "reference_offset": p.get("reference_offset", np.nan),
    } for i, p in enumerate(params)])
    corr["flagged_below_min"] = corr["correlation"] < cfg.min_correlation
    corr.to_csv(OUT / f"step5_correlation_{TAG}.csv", index=False)
    c = corr["correlation"].dropna()
    n_flag = int(corr["flagged_below_min"].sum())
    logger.info("STEP 5  coarse correlation median %.4f  (was 0.2698); flagged %d/%d "
                "(was 29/47)", c.median(), n_flag, n)

    # ------------------------------------------------ step 6, FULL transform
    shape = fine.shape[1:]
    moved = []
    for i, s in enumerate(secs):
        p = params[i]
        pts = fiducial_array(fid1, s, BASE_DS)
        moved.append(transform_points(pts, shape, p.get("angle", 0.0),
                                      p.get("dy_fine", 0.0), p.get("dx_fine", 0.0),
                                      fields[i]))
    tre = pd.DataFrame([
        {**target_registration_error(moved[i], moved[i + 1], mpp=fine_mpp),
         "pair": f"{secs[i]}-{secs[i+1]}"} for i in range(n - 1)])
    tre.to_csv(OUT / f"step6_tre_pairwise_{TAG}.csv", index=False)

    centre = n // 2
    atre = accumulated_tre(moved, reference_index=centre, mpp=fine_mpp)
    atre.to_csv(OUT / f"step6_atre_{TAG}.csv", index=False)

    identity = np.array([np.linalg.norm(fiducial_array(fid1, secs[i]) -
                                        fiducial_array(fid1, secs[i + 1]),
                                        axis=1).mean() * NATIVE_MPP_UM
                         for i in range(n - 1)])
    floor = np.array([procrustes_residual(fiducial_array(fid1, secs[i]),
                                          fiducial_array(fid1, secs[i + 1]))
                      * NATIVE_MPP_UM for i in range(n - 1)])
    inter = np.concatenate([
        np.linalg.norm(fiducial_array(fid1, s) - fiducial_array(fid2, s), axis=1)
        * NATIVE_MPP_UM for s in secs])

    new_tre = float(tre["tre_mean_um"].mean())
    logger.info("STEP 6  TRE (FULL transform incl. elastic): mean %.0f  median %.0f um",
                new_tre, tre["tre_median_um"].median())
    logger.info("        identity (no transform)      : %.0f um", identity.mean())
    logger.info("        ORIGINAL pipeline            : 2544 um")
    logger.info("        rigid floor (Procrustes)     : %.0f um", floor.mean())
    logger.info("        annotation floor             : %.1f um", np.median(inter))
    verdict = ("BEATS the do-nothing baseline" if new_tre < identity.mean()
               else "still WORSE than doing nothing")
    logger.info("        VERDICT: %s", verdict)
    logger.info("STEP 6  ATRE mean %.0f um (was 2473)", atre["atre_mean_um"].mean())

    # ---------------------------------------------------------------- step 7
    from run_kartasalo_registration import detect_holes, straightness, track_holes
    per_raw = [detect_holes(fine[i], fine_mpp) for i in range(n)]
    per = []
    for i, pts in enumerate(per_raw):
        if not len(pts):
            per.append(pts); continue
        p = params[i]
        per.append(transform_points(pts, shape, p.get("angle", 0.0),
                                    p.get("dy_fine", 0.0), p.get("dx_fine", 0.0),
                                    fields[i]))
    tracks = track_holes(per, max_jump_px=2500.0 / fine_mpp)
    straight = straightness(tracks, fine_mpp, SECTION_THICKNESS_UM)
    straight.to_csv(OUT / f"step7_hole_straightness_{TAG}.csv", index=False)
    if len(straight):
        logger.info("STEP 7  hole straightness mean %.0f um (was 1331), tracked on "
                    "%d-%d sections", straight["residual_mean_um"].mean(),
                    straight["n_sections_found"].min(),
                    straight["n_sections_found"].max())

    # ---------------------------------------------------------------- step 8
    ax = axial_vs_lateral_correlation(registered, mpp=fine_mpp,
                                      section_um=SECTION_THICKNESS_UM)
    ax.to_csv(OUT / f"step8_axial_lateral_{TAG}.csv", index=False)
    zs = z_skip_validation(registered, section_um=SECTION_THICKNESS_UM)
    zs.to_csv(OUT / f"step8_zskip_{TAG}.csv", index=False)
    logger.info("STEP 8  axial vs lateral:\n%s", ax.to_string(index=False))
    logger.info("STEP 8  z-skip:\n%s", zs.to_string(index=False))

    summary = {
        "dataset": "Kartasalo liver, Etsin c76335fa-cdcf-4ddc-ab1c-1882bad82861, CC BY 4.0",
        "n_sections": n, "downsample": BASE_DS,
        "mpp_um": fine_mpp, "coarse_mpp_um": coarse_mpp,
        "mpp_provenance": "paper text; TIFF carries no calibration",
        "section_thickness_um": SECTION_THICKNESS_UM,
        "elastic": not args.no_elastic, "max_abs_deg": args.max_abs_deg,
        "runtime_min": round(mins, 1),
        "correlation_median": float(c.median()), "n_flagged": n_flag,
        "tre_full_mean_um": new_tre,
        "tre_full_median_um": float(tre["tre_median_um"].median()),
        "atre_mean_um": float(atre["atre_mean_um"].mean()),
        "identity_mean_um": float(identity.mean()),
        "original_pipeline_mean_um": 2543.7,
        "rigid_floor_mean_um": float(floor.mean()),
        "interobserver_median_um": float(np.median(inter)),
        "hole_deviation_mean_um": (float(straight["residual_mean_um"].mean())
                                   if len(straight) else None),
        "beats_identity": bool(new_tre < identity.mean()),
        "tre_note": "FULL transform including elastic; fields are returned by the driver",
    }
    (OUT / f"summary_{TAG}.json").write_text(json.dumps(summary, indent=2))
    logger.info("wrote %s", OUT / f"summary_{TAG}.json")


if __name__ == "__main__":
    main()
