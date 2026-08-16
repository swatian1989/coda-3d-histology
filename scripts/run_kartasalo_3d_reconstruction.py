#!/usr/bin/env python
"""CODA stages 5 and 6 on the liver stack: a real volume, and 2D versus 3D counting.

    python scripts/run_kartasalo_3d_reconstruction.py

This is the part of CODA that a single section cannot support at all, and it is
only reachable because the corrected registration brought target error to 175
microns from 2544. It is a genuine three-dimensional reconstruction: 47
registered sections stacked into a volume with real z spacing, not a
stereological correction applied to one plane.

THE HEADLINE MEASUREMENT

CODA's central claim is that counting objects on single sections overestimates
how many objects exist, because one structure crossing several sections is
counted once per section. In pancreas they measured a 12.3-fold mean
overestimate. The measurement is a ratio:

    overcounting = (sum over sections of objects seen in that section)
                   / (connected components in the volume)

A structure spanning k sections contributes k to the numerator and 1 to the
denominator. The ratio therefore reports how much single-section counting
inflates object number, and it cannot be computed without a volume.

WHAT IS BEING COUNTED, AND WHY IT IS NOT THE PAPER'S TEN CLASSES

CODA counts objects from a DeepLab semantic segmentation trained on annotated
histology into ten tissue classes. No such annotation or trained model exists
here, and inventing one would be worse than not having it. Instead the objects
counted are vascular lumina: in liver these are enclosed, background-bright
spaces inside dark parenchyma, separable by the same intensity-band logic that
located the laser-cut fiducial holes, whose accuracy was independently checked
against two human annotators. That is a narrower target than the paper's, and
it is labelled as such everywhere. The overcounting geometry it measures is the
same.

Z SPACING IS ANISOTROPIC AND IS NOT PRETENDED OTHERWISE

Sections are 5 um apart and pixels are 7.36 um, so the voxel is roughly
isotropic in this particular case, which is unusually convenient. Connectivity
is computed in 3D with that spacing recorded, and the volume is never resampled
to fake isotropy.
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
    NATIVE_MPP_UM, SECTION_THICKNESS_UM, load_stack,
)
from coda_my.registration import RegistrationConfig, apply_rigid  # noqa: E402
from coda_my.registration_fix import (  # noqa: E402
    SearchConfig, register_stack_two_scale,
)

DATA = ROOT / "data/raw/kartasalo/extracted/Data_to_IDA"
OUT = ROOT / "results/kartasalo"
BASE_DS = 16
CODA_PANCREAS_OVERCOUNT = 12.3      # [PAPER] Kiemen et al., for reference only
logger = logging.getLogger("kartasalo_3d")


def blockmean(a: np.ndarray, f: int) -> np.ndarray:
    if f == 1:
        return a.astype(np.float32)
    h, w = a.shape
    h2, w2 = (h // f) * f, (w // f) * f
    return a[:h2, :w2].astype(np.float32).reshape(h2 // f, f, w2 // f, f).mean(axis=(1, 3))


def segment_lumina(volume: np.ndarray, mpp: float,
                   min_um: float = 40.0, max_um: float = 1200.0) -> np.ndarray:
    """Binary mask of enclosed bright spaces inside tissue, per section.

    Same logic that located the fiducial holes: the free slide background is
    saturated and occupies most of the frame, tissue is dark, and lumina form a
    band between the two. Thresholding on brightness alone merges lumina into
    the background, so the band and the enclosure test are both required.
    """
    out = np.zeros(volume.shape, dtype=bool)
    min_px = np.pi * (min_um / 2 / mpp) ** 2
    max_px = np.pi * (max_um / 2 / mpp) ** 2
    for i in range(len(volume)):
        img = volume[i]
        if not np.isfinite(img).any() or img.max() <= 0:
            continue
        sat = np.percentile(img, 85)
        tis = np.percentile(img, 20)
        if sat - tis < 1e-6:
            continue
        band = (img >= tis + 0.55 * (sat - tis)) & (img < sat - 0.02 * (sat - tis))
        band = ndimage.binary_opening(band, np.ones((3, 3)))
        lab, n = ndimage.label(band)
        if n == 0:
            continue
        sizes = ndimage.sum(band, lab, range(1, n + 1))
        keep = {j + 1 for j, s in enumerate(sizes) if min_px <= s <= max_px}
        if keep:
            out[i] = np.isin(lab, list(keep))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coarse-extra", type=int, default=11)
    ap.add_argument("--max-abs-deg", type=float, default=45.0)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.FileHandler(ROOT / "logs/kartasalo_3d.log"),
                                  logging.StreamHandler(sys.stdout)])

    fine_mpp = NATIVE_MPP_UM * BASE_DS
    coarse_mpp = fine_mpp * args.coarse_extra
    fine, _ = load_stack(DATA / "liver", downsample=BASE_DS,
                         cache=ROOT / f"data/interim/kartasalo_liver_ds{BASE_DS}.npy")
    coarse = np.stack([blockmean(fine[i], args.coarse_extra) for i in range(len(fine))])
    n = len(fine)

    logger.info("=" * 72)
    logger.info("STAGE 5: building a volume from %d sections, %.2f um/px, %.1f um apart",
                n, fine_mpp, SECTION_THICKNESS_UM)

    cfg = RegistrationConfig()
    t0 = time.time()
    _, params, fields = register_stack_two_scale(
        coarse, fine, coarse_mpp, fine_mpp, cfg,
        SearchConfig(max_abs_deg=args.max_abs_deg), elastic=True)

    # Rebuild at NATURAL intensity: the driver returns preprocessed images
    # (complemented, background zeroed), which destroys the bright/dark
    # relationship the luminal segmentation depends on.
    from coda_my.registration import apply_elastic
    vol = np.zeros((n, *fine.shape[1:]), dtype=np.float32)
    for i in range(n):
        p = params[i]
        img = apply_rigid(fine[i].astype(np.float32), p.get("angle", 0.0),
                          p.get("dy_fine", 0.0), p.get("dx_fine", 0.0))
        if fields[i] is not None:
            img = apply_elastic(img, *fields[i])
        vol[i] = img
    logger.info("volume %s built in %.1f min", vol.shape, (time.time() - t0) / 60)
    np.save(OUT / "volume_natural_ds16fix.npy", vol.astype(np.uint8))

    # ------------------------------------------------------------- stage 6
    logger.info("STAGE 6: segmenting vascular lumina and counting in 2D vs 3D")
    mask = segment_lumina(vol, fine_mpp)
    frac = mask.mean()
    logger.info("lumen voxels: %.3f%% of the volume", 100 * frac)

    # 2D: objects counted independently on each section, as single-section work does
    per_section = []
    for i in range(n):
        lab2, n2 = ndimage.label(mask[i])
        per_section.append(n2)
    total_2d = int(np.sum(per_section))

    # 3D: connected components through the stack
    lab3, n3 = ndimage.label(mask)
    logger.info("2D object count summed over sections: %d", total_2d)
    logger.info("3D connected components in the volume : %d", n3)
    ratio = total_2d / max(n3, 1)
    logger.info("OVERCOUNTING RATIO: %.2f-fold  (CODA reported %.1f-fold in pancreas)",
                ratio, CODA_PANCREAS_OVERCOUNT)

    # how many sections each 3D object spans, which is what drives the ratio
    spans, vols = [], []
    if n3:
        objs = ndimage.find_objects(lab3)
        for idx, sl in enumerate(objs, start=1):
            if sl is None:
                continue
            spans.append(sl[0].stop - sl[0].start)
            vols.append(int((lab3[sl] == idx).sum()))
    spans = np.array(spans); vols = np.array(vols)
    voxel_um3 = fine_mpp * fine_mpp * SECTION_THICKNESS_UM

    pd.DataFrame({"section": np.arange(1, n + 1), "objects_2d": per_section}).to_csv(
        OUT / "stage6_2d_counts.csv", index=False)
    pd.DataFrame({"object": np.arange(1, len(spans) + 1),
                  "sections_spanned": spans,
                  "volume_voxels": vols,
                  "volume_um3": vols * voxel_um3}).to_csv(
        OUT / "stage6_3d_objects.csv", index=False)

    logger.info("sections spanned per 3D object: mean %.2f  median %.0f  max %d",
                spans.mean() if len(spans) else 0,
                np.median(spans) if len(spans) else 0,
                spans.max() if len(spans) else 0)
    single = int((spans == 1).sum())
    logger.info("objects confined to ONE section: %d of %d (%.0f%%) - these are the "
                "ones single-section counting gets right", single, len(spans),
                100 * single / max(len(spans), 1))

    # 2D projections for the figure, so it never has to load the whole volume.
    # Projecting the LUMEN MASK, not the intensity: after rigid warping the area
    # outside the tissue is filled with zeros, so a minimum-intensity projection
    # through the stack returns zero almost everywhere and shows nothing.
    ys, xs = np.nonzero(mask.any(axis=0))
    if len(ys):
        y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    else:
        y0, y1, x0, x1 = 0, mask.shape[1], 0, mask.shape[2]
    np.savez_compressed(
        OUT / "stage6_projections.npz",
        section_mid=vol[n // 2].astype(np.uint8),
        lumen_xz=mask[:, y0:y1, x0:x1].sum(axis=1).astype(np.uint16),
        lumen_yz=mask[:, y0:y1, x0:x1].sum(axis=2).astype(np.uint16),
        lumen_xy=mask.sum(axis=0).astype(np.uint16),
    )
    logger.info("wrote projections for plotting")

    summary = {
        "n_sections": n, "mpp_um": fine_mpp,
        "section_thickness_um": SECTION_THICKNESS_UM,
        "voxel_um3": voxel_um3,
        "is_true_3d_reconstruction": True,
        "objects_counted": "vascular lumina by intensity band, NOT the paper's "
                           "10-class DeepLab segmentation",
        "lumen_volume_fraction": float(frac),
        "total_2d_object_count": total_2d,
        "n_3d_objects": int(n3),
        "overcounting_ratio": float(ratio),
        "coda_pancreas_reference": CODA_PANCREAS_OVERCOUNT,
        "mean_sections_spanned": float(spans.mean()) if len(spans) else None,
        "median_sections_spanned": float(np.median(spans)) if len(spans) else None,
        "max_sections_spanned": int(spans.max()) if len(spans) else None,
        "objects_in_one_section_only": single,
        "median_object_volume_um3": float(np.median(vols) * voxel_um3) if len(vols) else None,
        "registration_tre_um": 174.6,
        "caveat": ("mouse liver, not breast; lumina segmented by intensity rather than "
                   "by a trained multi-class model; the ratio depends on what is counted"),
    }
    (OUT / "stage6_summary.json").write_text(json.dumps(summary, indent=2))
    logger.info("wrote %s", OUT / "stage6_summary.json")


if __name__ == "__main__":
    main()
