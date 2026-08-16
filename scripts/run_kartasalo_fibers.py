#!/usr/bin/env python
"""CODA stage 7 on the liver volume: fibre alignment, and the sectioning-angle effect.

    python scripts/run_kartasalo_fibers.py

Stage 7 needs an eosin channel, which DAB immunohistochemistry does not have.
The liver series is H&E, so it does.

WHY THIS NEEDS THE VOLUME

Fibre alignment measured on a single section depends on the angle the structure
was cut at. Kiemen et al. could correct for that because they had a volume and
could choose the plane deliberately, and they measured 2.2 to 2.5-fold
differences between longitudinal and axial cuts through the same structures. On
a single section the angle is unknown and uncorrectable, which is why the Arm C
protocol restricts itself to reporting a distribution.

Having reconstructed a volume, the same comparison becomes available here: the
cutting plane (xy) is what the microtome produced, and the two orthogonal planes
(xz and yz) are views no individual section contains. Anisotropy is measured in
all three and compared. This is the one place in this study where sectioning
angle can be treated as a variable rather than as noise.

WHAT IS AND IS NOT COMPARABLE TO THE PAPER

The paper compared longitudinal against axial cuts of identified anatomical
structures, selected by hand from the volume. Here the three orthogonal planes
of the reconstruction are compared instead, without identifying structures,
because that identification needs the trained multi-class segmentation that this
work does not have. The comparison therefore tests whether cutting angle changes
the measurement, which is the methodological claim, and not whether a named
structure is more aligned along its axis, which is the anatomical one.

Voxels are anisotropic in a way that matters here and is not hidden: 7.36 um in
plane against 5 um between sections. The xz and yz planes are resampled to the
in-plane pixel size before measurement so that a window is the same physical
size in every orientation, because the anisotropy index is computed over a fixed
window in microns.
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
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coda_my.deconv import deconvolve  # noqa: E402
from coda_my.fibers import FiberConfig, tiled_anisotropy  # noqa: E402
from coda_my.loaders.kartasalo import (  # noqa: E402
    NATIVE_MPP_UM, SECTION_THICKNESS_UM, load_stack, section_paths,
)
from coda_my.registration import RegistrationConfig, apply_rigid  # noqa: E402
from coda_my.registration_fix import (  # noqa: E402
    SearchConfig, register_stack_two_scale,
)

DATA = ROOT / "data/raw/kartasalo/extracted/Data_to_IDA"
OUT = ROOT / "results/kartasalo"
BASE_DS = 16
logger = logging.getLogger("kartasalo_fibers")


def blockmean(a: np.ndarray, f: int) -> np.ndarray:
    if f == 1:
        return a.astype(np.float32)
    h, w = a.shape
    h2, w2 = (h // f) * f, (w // f) * f
    return a[:h2, :w2].astype(np.float32).reshape(h2 // f, f, w2 // f, f).mean(axis=(1, 3))


def eosin_of(rgb: np.ndarray) -> np.ndarray:
    """Eosin optical density, normalised to [0, 1]."""
    ch = deconvolve(np.clip(rgb, 0, 255).astype(np.uint8))
    e = ch["eosin"]
    lo, hi = np.percentile(e, 1), np.percentile(e, 99)
    return np.clip((e - lo) / max(hi - lo, 1e-6), 0, 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coarse-extra", type=int, default=11)
    ap.add_argument("--n-planes", type=int, default=40,
                    help="orthogonal planes sampled per orientation")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.FileHandler(ROOT / "logs/kartasalo_fibers.log"),
                                  logging.StreamHandler(sys.stdout)])
    fine_mpp = NATIVE_MPP_UM * BASE_DS
    coarse_mpp = fine_mpp * args.coarse_extra

    # grayscale drives the registration; colour is what fibre analysis needs
    grey, _ = load_stack(DATA / "liver", downsample=BASE_DS,
                         cache=ROOT / f"data/interim/kartasalo_liver_ds{BASE_DS}.npy")
    coarse = np.stack([blockmean(grey[i], args.coarse_extra) for i in range(len(grey))])
    n = len(grey)

    rgb_cache = ROOT / f"data/interim/kartasalo_liver_ds{BASE_DS}_rgb.npy"
    if rgb_cache.exists():
        rgb = np.load(rgb_cache)
    else:
        from PIL import Image
        Image.MAX_IMAGE_PIXELS = None
        paths = section_paths(DATA / "liver")
        first = np.asarray(Image.open(paths[0]).reduce(BASE_DS).convert("RGB"))
        rgb = np.zeros((n, *first.shape), dtype=np.uint8)
        rgb[0] = first
        for i, p in enumerate(paths[1:], start=1):
            with Image.open(p) as im:
                a = np.asarray(im.reduce(BASE_DS).convert("RGB"))
            h = min(a.shape[0], first.shape[0]); w = min(a.shape[1], first.shape[1])
            rgb[i, :h, :w] = a[:h, :w]
            if (i + 1) % 10 == 0:
                logger.info("loaded RGB %d/%d", i + 1, n)
        np.save(rgb_cache, rgb)
    logger.info("RGB stack %s (%.2f GB)", rgb.shape, rgb.nbytes / 1e9)

    t0 = time.time()
    _, params, _ = register_stack_two_scale(
        coarse, grey, coarse_mpp, fine_mpp, RegistrationConfig(),
        SearchConfig(max_abs_deg=45.0), elastic=False)
    vol = np.zeros_like(rgb)
    for i in range(n):
        p = params[i]
        for ch in range(3):
            vol[i, :, :, ch] = np.clip(
                apply_rigid(rgb[i, :, :, ch].astype(np.float32), p.get("angle", 0.0),
                            p.get("dy_fine", 0.0), p.get("dx_fine", 0.0)), 0, 255)
    logger.info("colour volume built in %.1f min", (time.time() - t0) / 60)

    # The window has to fit in the THINNEST dimension of the volume, and that is
    # z: 47 sections at 5 um is a block only 235 um deep against 9.8 mm in
    # plane. A window larger than the block depth yields no orthogonal windows
    # at all. The same window is then used in every orientation, because a
    # comparison between planes measured with different window sizes would be
    # measuring the window. CODA used 50 um at 0.5 um/px, which is 100 px; the
    # same pixel count is impossible here and the constraint is the stack depth,
    # not the implementation.
    depth_um = n * SECTION_THICKNESS_UM
    window_um = float(np.floor(depth_um * 0.75 / fine_mpp) * fine_mpp)
    cfg = FiberConfig(mpp=fine_mpp, window_um=window_um)
    logger.info("block depth %.0f um; fibre window %.0f um = %d px at %.2f um/px "
                "(same window in all orientations)",
                depth_um, cfg.window_um, int(cfg.window_um / cfg.mpp), cfg.mpp)

    rows = []
    # --- xy: the plane the microtome actually produced
    idx = np.linspace(0, n - 1, min(args.n_planes, n)).astype(int)
    for i in idx:
        a = tiled_anisotropy(eosin_of(vol[i]), cfg)
        a = a[np.isfinite(a)]
        rows += [{"plane": "xy (cutting plane)", "index": int(i), "anisotropy": float(v)}
                 for v in a]

    # --- xz and yz: views no single section contains.
    # Resample the z axis to the in-plane pixel size so a window is square in
    # microns; without this the window is 7.36 um wide and 5 um tall and the
    # index would measure the resampling, not the tissue.
    zoom_z = SECTION_THICKNESS_UM / fine_mpp
    for name, axis in (("xz (orthogonal)", 1), ("yz (orthogonal)", 2)):
        size = vol.shape[axis]
        picks = np.linspace(size * 0.25, size * 0.75, args.n_planes).astype(int)
        for j in picks:
            plane = vol[:, j, :, :] if axis == 1 else vol[:, :, j, :]
            if plane.shape[0] < 8:
                continue
            pl = ndimage.zoom(plane.astype(np.float32), (zoom_z, 1, 1), order=1)
            a = tiled_anisotropy(eosin_of(pl), cfg)
            a = a[np.isfinite(a)]
            rows += [{"plane": name, "index": int(j), "anisotropy": float(v)}
                     for v in a]

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "stage7_fiber_anisotropy.csv", index=False)
    logger.info("%d windows measured across %d planes", len(df), df["index"].nunique())

    summary = {}
    logger.info("STAGE 7  fibre anisotropy index by plane orientation")
    for name, g in df.groupby("plane"):
        v = g["anisotropy"]
        summary[name] = {"n_windows": int(len(v)), "mean": float(v.mean()),
                         "median": float(v.median()),
                         "q1": float(v.quantile(.25)), "q3": float(v.quantile(.75))}
        logger.info("  %-20s n=%6d  mean %.4f  median %.4f  IQR %.4f-%.4f",
                    name, len(v), v.mean(), v.median(),
                    v.quantile(.25), v.quantile(.75))

    xy = df.loc[df["plane"].str.startswith("xy"), "anisotropy"]
    tests = {}
    for name in ("xz (orthogonal)", "yz (orthogonal)"):
        o = df.loc[df["plane"] == name, "anisotropy"]
        if len(o) < 10:
            continue
        u = stats.mannwhitneyu(xy, o, alternative="two-sided")
        fold = float(o.mean() / xy.mean()) if xy.mean() else float("nan")
        # Cliff's delta on a subsample, since the full pairing is large
        rng = np.random.default_rng(0)
        a = rng.choice(xy.to_numpy(), size=min(2000, len(xy)), replace=False)
        b = rng.choice(o.to_numpy(), size=min(2000, len(o)), replace=False)
        delta = float((np.sign(a[:, None] - b[None, :]).mean()))
        tests[name] = {"U": float(u.statistic), "p": float(u.pvalue),
                       "fold_vs_xy": fold, "cliffs_delta": delta,
                       "n_xy": int(len(xy)), "n_other": int(len(o))}
        logger.info("  %s vs xy: %.2f-fold, Mann-Whitney U=%.0f p=%.3g, "
                    "Cliff's delta %.3f", name, fold, u.statistic, u.pvalue, delta)

    out = {"config": {"mpp_um": fine_mpp, "window_um": cfg.window_um,
                      "section_thickness_um": SECTION_THICKNESS_UM,
                      "n_sections": n},
           "by_plane": summary, "tests": tests,
           "coda_reference_fold": [2.2, 2.5],
           "note": ("orthogonal planes are views no single section contains; the "
                    "comparison tests whether cutting angle changes the measurement, "
                    "not whether a named structure is aligned along its axis")}
    (OUT / "stage7_summary.json").write_text(json.dumps(out, indent=2))
    logger.info("wrote %s", OUT / "stage7_summary.json")


if __name__ == "__main__":
    main()
