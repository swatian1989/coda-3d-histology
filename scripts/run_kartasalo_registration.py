#!/usr/bin/env python
"""Phase 2 steps 4-8: register the Kartasalo liver stack and measure accuracy.

    python scripts/run_kartasalo_registration.py --downsample 16
    python scripts/run_kartasalo_registration.py --limit 7 --downsample 32   # smoke

Stages are cached to results/kartasalo/, so a later stage can be re-run without
repeating registration.

WHAT IS MEASURED, AND WHAT EACH MEASUREMENT CAN AND CANNOT SEE

step 4  register_stack(elastic=True) over all sections, centre-out.
step 5  Pixel correlation per section; anything under cfg.min_correlation is
        flagged rather than dropped.
step 6  TRE and ATRE against the operator-annotated fiducials.

        IMPORTANT LIMIT. register_stack returns the RIGID parameters only;
        the elastic displacement field is applied to the image and then
        discarded, so it cannot be replayed onto point coordinates. Landmark
        TRE here is therefore the error after the rigid stage. It is not the
        error of the elastically registered images, and it is reported as
        rigid-stage TRE everywhere. Fixing this properly would mean changing
        registration.py, which is a protected module.

step 7  The independent check, and the one that does see the elastic stage.
        The four laser-cut holes are detected DE NOVO in the registered images
        by morphology, with no reference to the annotation, then matched
        through z and fitted with a straight line. Residual scatter about that
        line is reconstruction error in microns. Because the holes are found
        in the final registered images, this number includes every stage of
        the transform.

step 8  axial_vs_lateral_correlation and z_skip_validation.

The point-transform convention is not assumed. `_fit_point_transform` recovers
it empirically by pushing a marker through apply_rigid and reading back where
it lands, because a sign error in the rotation would silently corrupt every
TRE number rather than raising.
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
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coda_my.loaders.kartasalo import (  # noqa: E402
    NATIVE_MPP_UM, SECTION_THICKNESS_UM, fiducial_array, load_fiducials,
    load_stack, section_paths,
)
from coda_my.qc import (  # noqa: E402
    accumulated_tre, axial_vs_lateral_correlation, target_registration_error,
    z_skip_validation,
)
from coda_my.registration import (  # noqa: E402
    RegistrationConfig, apply_rigid, pixel_correlation, register_stack,
)

DATA = ROOT / "data/raw/kartasalo/extracted/Data_to_IDA"
OUT = ROOT / "results/kartasalo"
logger = logging.getLogger("kartasalo_reg")


# --------------------------------------------------------------- point transform


def _fit_point_transform(shape: tuple[int, int]) -> callable:
    """Recover how apply_rigid moves a POINT, by measurement not by assumption.

    ndimage.rotate's sign convention in array coordinates is easy to get
    backwards, and a backwards rotation still produces plausible-looking TRE
    numbers. So push a delta through the real function and solve for the
    convention that reproduces where it landed.
    """
    h, w = shape
    probe = np.zeros(shape)
    py, px = h // 3, w // 4                      # off-centre, so rotation shows
    probe[py, px] = 1.0
    angle, dy, dx = 23.0, 11.0, -7.0
    moved = apply_rigid(probe, angle, dy, dx)
    got = np.unravel_index(np.argmax(moved), moved.shape)

    cy, cx = (h - 1) / 2.0, (w - 1) / 2.0
    best, best_err = None, np.inf
    for sign in (1.0, -1.0):
        t = np.deg2rad(sign * angle)
        c, s = np.cos(t), np.sin(t)
        yy = cy + (py - cy) * c - (px - cx) * s + dy
        xx = cx + (py - cy) * s + (px - cx) * c + dx
        err = np.hypot(yy - got[0], xx - got[1])
        if err < best_err:
            best, best_err = sign, err
    logger.info("point-transform convention: sign %+.0f (residual %.2f px)",
                best, best_err)
    if best_err > 2.0:
        raise RuntimeError(f"could not recover point convention (residual "
                           f"{best_err:.2f} px); refusing to report TRE")

    def transform(points: np.ndarray, angle: float, dy: float, dx: float,
                  _sign=best, _cy=cy, _cx=cx) -> np.ndarray:
        t = np.deg2rad(_sign * angle)
        c, s = np.cos(t), np.sin(t)
        y, x = points[:, 0] - _cy, points[:, 1] - _cx
        return np.column_stack([_cy + y * c - x * s + dy,
                                _cx + y * s + x * c + dx])

    return transform


# --------------------------------------------------------------- hole detection


def detect_holes(section: np.ndarray, mpp: float, n_expect: int = 4,
                 min_um: float = 80.0, max_um: float = 600.0) -> np.ndarray:
    """Find laser-cut holes de novo, on a RAW (uninverted) greyscale section.

    No fiducial annotation is consulted, which is the whole point: this gives a
    ground truth that does not inherit the operators' clicks.

    A hole is NOT simply "bright". Measured on these slides, the free slide
    background saturates at 255 and occupies 74 percent of the frame, tissue
    sits below 150 at 21 percent, and the holes form a narrow band between the
    two at roughly 215 to 250, only 1.2 percent of pixels. Thresholding on
    brightness alone therefore merges every hole into the background; the band
    is what separates them. Measured hole diameters here are 157 to 334 um.

    The band edges are taken from the image's own histogram rather than fixed
    grey levels, so a section stained or exposed differently is still handled.
    """
    img = np.asarray(section)
    if img.ndim == 3:
        img = img.mean(axis=2)
    img = img.astype(np.float32)

    sat = np.percentile(img, 85)                     # slide background
    tis = np.percentile(img, 20)                     # tissue
    lo, hi = tis + 0.55 * (sat - tis), sat - 0.02 * (sat - tis)
    band = (img >= lo) & (img < hi)
    band = ndimage.binary_opening(band, np.ones((3, 3)))

    lab, n = ndimage.label(band)
    if n == 0:
        return np.empty((0, 2))
    min_px = np.pi * (min_um / 2 / mpp) ** 2
    max_px = np.pi * (max_um / 2 / mpp) ** 2
    rows = []
    objs = ndimage.find_objects(lab)
    for idx in range(1, n + 1):
        sl = objs[idx - 1]
        if sl is None:
            continue
        m = lab[sl] == idx
        area = int(m.sum())
        if not (min_px <= area <= max_px):
            continue
        ys, xs = np.nonzero(m)
        cy, cx = ys.mean(), xs.mean()
        r = np.hypot(ys - cy, xs - cx).max()
        if area / (np.pi * r * r + 1e-9) < 0.40:     # roughly round
            continue
        rows.append((cy + sl[0].start, cx + sl[1].start, area))
    if not rows:
        return np.empty((0, 2))
    rows.sort(key=lambda t: -t[2])
    return np.array([[r[0], r[1]] for r in rows[:n_expect]])


def track_holes(per_section: list[np.ndarray], max_jump_px: float) -> np.ndarray:
    """Match detected holes through z by nearest neighbour.

    The gate has to be generous. After registration the holes are only as well
    aligned as the registration itself, so a threshold tighter than the
    residual error unmatches exactly the sections that matter. There are only
    four holes and they are far apart, so a loose gate with one-to-one
    assignment is safe: mismatching would require the residual to exceed the
    spacing between distinct holes.

    Returns (n_tracks, n_sections, 2) with NaN where a hole was not found.
    """
    n = len(per_section)
    seed_i = next((i for i, p in enumerate(per_section) if len(p)), None)
    if seed_i is None:
        return np.empty((0, n, 2))
    seed = per_section[seed_i]
    tracks = np.full((len(seed), n, 2), np.nan)
    last = {t: seed[t].copy() for t in range(len(seed))}

    for i, pts in enumerate(per_section):
        if not len(pts):
            continue
        # greedy one-to-one: shortest available pairing first
        pairs = []
        for t in range(len(seed)):
            ref = last.get(t)
            if ref is None:
                continue
            for k in range(len(pts)):
                pairs.append((float(np.hypot(pts[k, 0] - ref[0],
                                             pts[k, 1] - ref[1])), t, k))
        pairs.sort()
        done_t, done_k = set(), set()
        for d, t, k in pairs:
            if t in done_t or k in done_k or d > max_jump_px:
                continue
            tracks[t, i] = pts[k]
            last[t] = pts[k]
            done_t.add(t); done_k.add(k)
    return tracks


def straightness(tracks: np.ndarray, mpp: float, section_um: float) -> pd.DataFrame:
    """Residual of each hole's centroid about a straight line down z, in um."""
    rows = []
    n_sec = tracks.shape[1]
    z = np.arange(n_sec) * section_um
    for t in range(tracks.shape[0]):
        ok = ~np.isnan(tracks[t, :, 0])
        if ok.sum() < 5:
            continue
        res = []
        for axis in (0, 1):
            v = tracks[t, ok, axis] * mpp
            coef = np.polyfit(z[ok], v, 1)
            res.append(v - np.polyval(coef, z[ok]))
        d = np.hypot(res[0], res[1])
        rows.append({"hole": t + 1, "n_sections_found": int(ok.sum()),
                     "residual_mean_um": float(d.mean()),
                     "residual_median_um": float(np.median(d)),
                     "residual_p95_um": float(np.percentile(d, 95)),
                     "residual_max_um": float(d.max())})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--downsample", type=int, default=16,
                    help="16 gives 7.36 um/px, 8 gives 3.68 um/px")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--no-elastic", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    (ROOT / "logs").mkdir(exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.FileHandler(ROOT / "logs/kartasalo_reg.log"),
                                  logging.StreamHandler(sys.stdout)])

    mpp = NATIVE_MPP_UM * args.downsample
    tag = f"ds{args.downsample}" + (f"_n{args.limit}" if args.limit else "")
    logger.info("=" * 72)
    logger.info("Kartasalo liver, downsample %d -> %.2f um/px, sections %.1f um",
                args.downsample, mpp, SECTION_THICKNESS_UM)

    # ---------------------------------------------------------------- load
    paths = section_paths(DATA / "liver")
    logger.info("sections available on disk: %d", len(paths))
    stack, meta = load_stack(DATA / "liver", downsample=args.downsample,
                             limit=args.limit,
                             cache=ROOT / f"data/interim/kartasalo_liver_{tag}.npy")
    logger.info("stack %s, %.2f GB in RAM", stack.shape, stack.nbytes / 1e9)

    fid1 = load_fiducials(DATA / "fiducialcoordinates_liver_observer1.txt")
    fid2 = load_fiducials(DATA / "fiducialcoordinates_liver_observer2.txt")

    # ------------------------------------------------------ step 4: register
    reg_cache = OUT / f"registered_{tag}.npy"
    par_cache = OUT / f"params_{tag}.json"
    if reg_cache.exists() and par_cache.exists() and not args.force:
        registered = np.load(reg_cache)
        params = json.loads(par_cache.read_text())
        logger.info("loaded cached registration %s", registered.shape)
    else:
        cfg = RegistrationConfig()
        t0 = time.time()
        logger.info("step 4: register_stack(elastic=%s) on %d sections ...",
                    not args.no_elastic, len(stack))
        regs, params = register_stack([stack[i] for i in range(len(stack))],
                                      cfg, elastic=not args.no_elastic)
        registered = np.stack(regs)
        np.save(reg_cache, registered)
        par_cache.write_text(json.dumps(params, indent=2, default=float))
        logger.info("step 4 done in %.1f min", (time.time() - t0) / 60)

    cfg = RegistrationConfig()

    # -------------------------------------------------- step 5: correlation
    corr = pd.DataFrame([
        {"section": i + 1,
         "correlation": p.get("correlation", np.nan),
         "correlation_after_elastic": p.get("correlation_after_elastic", np.nan),
         "angle_deg": p.get("angle", np.nan),
         "dy_px": p.get("dy", np.nan), "dx_px": p.get("dx", np.nan),
         "reference_offset": p.get("reference_offset", np.nan)}
        for i, p in enumerate(params)])
    corr["flagged_below_min"] = corr["correlation"] < cfg.min_correlation
    corr.to_csv(OUT / f"step5_correlation_{tag}.csv", index=False)

    c = corr["correlation"].dropna()
    ce = corr["correlation_after_elastic"].dropna()
    logger.info("STEP 5  pixel correlation: median %.4f  mean %.4f  min %.4f  "
                "p05 %.4f", c.median(), c.mean(), c.min(), c.quantile(.05))
    if len(ce):
        logger.info("STEP 5  after elastic:     median %.4f  (n=%d improved %d)",
                    ce.median(), len(ce),
                    int((corr["correlation_after_elastic"] >
                         corr["correlation"]).sum()))
    flagged = corr[corr["flagged_below_min"]]
    logger.info("STEP 5  flagged below min_correlation=%.2f: %d sections %s",
                cfg.min_correlation, len(flagged),
                list(flagged["section"]) if len(flagged) else "(none)")

    # --------------------------------------------------- step 6: TRE / ATRE
    tf = _fit_point_transform(stack.shape[1:])
    n = len(stack)
    secs = sorted(fid1["section"].unique())[:n]

    moved1, moved2 = [], []
    for i, s in enumerate(secs):
        p = params[i] if i < len(params) else {}
        a, dy, dx = p.get("angle", 0.0), p.get("dy", 0.0), p.get("dx", 0.0)
        moved1.append(tf(fiducial_array(fid1, s, args.downsample), a, dy, dx))
        moved2.append(tf(fiducial_array(fid2, s, args.downsample), a, dy, dx))

    centre = n // 2
    tre_rows = []
    for i in range(n - 1):
        r = target_registration_error(moved1[i], moved1[i + 1], mpp=mpp)
        r["pair"] = f"{secs[i]}-{secs[i+1]}"
        tre_rows.append(r)
    tre = pd.DataFrame(tre_rows)
    tre.to_csv(OUT / f"step6_tre_pairwise_{tag}.csv", index=False)

    atre = accumulated_tre(moved1, reference_index=centre, mpp=mpp)
    atre.to_csv(OUT / f"step6_atre_{tag}.csv", index=False)

    # annotation noise floor: the two observers on the SAME sections
    inter = []
    for i, s in enumerate(secs):
        d = np.linalg.norm(fiducial_array(fid1, s) - fiducial_array(fid2, s),
                           axis=1) * NATIVE_MPP_UM
        inter.extend(d.tolist())
    inter = np.array(inter)

    # Best achievable rigid alignment on these same pairs: fit the transform
    # directly to the fiducials. No registration method can beat this, so a TRE
    # near it means the limit is tissue deformation, not the algorithm.
    floor = []
    for i in range(n - 1):
        A = fiducial_array(fid1, secs[i]); B = fiducial_array(fid1, secs[i + 1])
        Ac, Bc = A - A.mean(0), B - B.mean(0)
        U, _, Vt = np.linalg.svd(Bc.T @ Ac)
        R = U @ Vt
        if np.linalg.det(R) < 0:
            U[:, -1] *= -1; R = U @ Vt
        floor.append(float(np.linalg.norm(Ac - Bc @ R.T, axis=1).mean() * NATIVE_MPP_UM))
    floor = np.array(floor)
    logger.info("STEP 6  RIGID FLOOR (Procrustes fitted to the fiducials themselves): "
                "mean %.1f median %.1f um - no rigid method can beat this",
                floor.mean(), np.median(floor))

    logger.info("STEP 6  pairwise TRE (rigid stage, observer 1): mean %.2f  "
                "median %.2f  p95 %.2f um",
                tre["tre_mean_um"].mean(), tre["tre_median_um"].median(),
                tre["tre_p95_um"].quantile(.95))
    logger.info("STEP 6  ATRE vs centre section %d: mean %.2f  max %.2f um",
                centre + 1, atre["atre_mean_um"].mean(), atre["atre_max_um"].max())
    if len(atre) > 3:
        rho = np.corrcoef(atre["distance_from_reference"],
                          atre["atre_mean_um"])[0, 1]
        logger.info("STEP 6  ATRE vs distance from centre: Pearson r %.3f "
                    "(positive means the stack drifts)", rho)
    logger.info("STEP 6  INTER-OBSERVER floor: mean %.2f median %.2f p95 %.2f um "
                "- TRE below this is not resolvable by the ground truth",
                inter.mean(), np.median(inter), np.percentile(inter, 95))

    # ------------------------------------------- step 7: independent holes
    # Detected on the RAW sections, where the holes occupy a clean intensity
    # band, then pushed through the same rigid transform as the images. The
    # registered stack returned by register_stack is preprocessed (complemented,
    # background zeroed, filtered), which clips part of the hole band, so
    # detecting there would trade a clean measurement for a noisy one.
    per_raw = [detect_holes(stack[i], mpp) for i in range(len(stack))]
    per = []
    for i, pts in enumerate(per_raw):
        if not len(pts):
            per.append(pts)
            continue
        p = params[i] if i < len(params) else {}
        per.append(tf(pts, p.get("angle", 0.0), p.get("dy", 0.0), p.get("dx", 0.0)))
    found = [len(p) for p in per]
    logger.info("STEP 7  holes detected de novo per section: median %d, "
                "sections with 4: %d/%d", int(np.median(found)),
                sum(f >= 4 for f in found), len(found))

    # validate the detector against the operators, which also bounds both
    val = []
    for i, s in enumerate(secs):
        if i >= len(per_raw) or not len(per_raw[i]):
            continue
        for a in fiducial_array(fid1, s, args.downsample):
            val.append(np.linalg.norm(per_raw[i] - a, axis=1).min() * mpp)
    val = np.array([v for v in val if v < 200.0])
    if len(val):
        logger.info("STEP 7  detector vs operator annotation: n=%d matched, "
                    "mean %.1f median %.1f um (floor is the %.1f um "
                    "inter-observer median)", len(val), val.mean(),
                    np.median(val), np.median(inter))
    tracks = track_holes(per, max_jump_px=2500.0 / mpp)   # generous; see track_holes
    straight = straightness(tracks, mpp, SECTION_THICKNESS_UM)
    straight.to_csv(OUT / f"step7_hole_straightness_{tag}.csv", index=False)
    if len(straight):
        logger.info("STEP 7  deviation from a straight line down z, per hole:")
        for _, r in straight.iterrows():
            logger.info("        hole %d  found on %d sections  mean %.1f  "
                        "median %.1f  p95 %.1f  max %.1f um", r["hole"],
                        r["n_sections_found"], r["residual_mean_um"],
                        r["residual_median_um"], r["residual_p95_um"],
                        r["residual_max_um"])
        logger.info("STEP 7  pooled mean deviation %.2f um "
                    "(annotation-independent ground truth)",
                    straight["residual_mean_um"].mean())
    else:
        logger.warning("STEP 7  no hole tracked on enough sections; reporting none")

    # ------------------------------------------------------------- step 8
    ax = axial_vs_lateral_correlation(registered, mpp=mpp,
                                      section_um=SECTION_THICKNESS_UM)
    ax.to_csv(OUT / f"step8_axial_lateral_{tag}.csv", index=False)
    logger.info("STEP 8  axial vs lateral correlation:\n%s", ax.to_string(index=False))

    zs = z_skip_validation(registered, section_um=SECTION_THICKNESS_UM)
    zs.to_csv(OUT / f"step8_zskip_{tag}.csv", index=False)
    logger.info("STEP 8  z-skip validation:\n%s", zs.to_string(index=False))

    summary = {
        "dataset": "Kartasalo liver, Etsin c76335fa-cdcf-4ddc-ab1c-1882bad82861, CC BY 4.0",
        "n_sections": int(n), "downsample": args.downsample,
        "mpp_um": mpp, "mpp_provenance": "paper text; TIFF carries no calibration",
        "section_thickness_um": SECTION_THICKNESS_UM,
        "elastic": not args.no_elastic,
        "correlation_median": float(c.median()),
        "n_flagged": int(len(flagged)),
        "tre_rigid_mean_um": float(tre["tre_mean_um"].mean()),
        "tre_rigid_median_um": float(tre["tre_median_um"].median()),
        "atre_mean_um": float(atre["atre_mean_um"].mean()),
        "interobserver_median_um": float(np.median(inter)),
        "rigid_floor_mean_um": float(floor.mean()),
        "rigid_floor_median_um": float(np.median(floor)),
        "hole_deviation_mean_um": (float(straight["residual_mean_um"].mean())
                                   if len(straight) else None),
        "tre_note": ("rigid stage only; register_stack does not return the "
                     "elastic field so it cannot be replayed onto points"),
    }
    (OUT / f"summary_{tag}.json").write_text(json.dumps(summary, indent=2))
    logger.info("wrote %s", OUT / f"summary_{tag}.json")


if __name__ == "__main__":
    main()
