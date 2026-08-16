#!/usr/bin/env python
"""Why the liver registration failed: the rigid stage ran at the wrong scale.

    python scripts/diagnose_registration_scale.py

RegistrationConfig declares

    global_mpp: float = 80.0        # [PAPER] register at 80 um/px

and nothing in registration.py ever reads it. The global rigid stage therefore
runs at whatever resolution the caller supplies. The production run supplied
7.36 um/px, because that is the scale the elastic stage needs, which is 10.9
times finer than the value the paper specifies for the rigid stage.

That is not a cosmetic mismatch. The rigid stage estimates rotation from the
Radon transform, and the paper deliberately does it coarse: at 80 um/px a liver
section is a compact blob whose principal axes dominate the transform, while at
7.36 um/px the transform is driven by fine parenchymal texture that differs
between sections. Fine detail is noise for this purpose.

This sweeps the rigid stage across scales, rigid only, and reports TRE against
the floor that a rigid transform fitted directly to the fiducials achieves. It
reuses the cached stack, so it re-decodes nothing.

The module is not modified. The scale is controlled by what is passed in.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coda_my.loaders.kartasalo import (  # noqa: E402
    NATIVE_MPP_UM, fiducial_array, load_fiducials,
)
from coda_my.registration import RegistrationConfig, register_stack  # noqa: E402

DATA = ROOT / "data/raw/kartasalo/extracted/Data_to_IDA"
CACHE = ROOT / "data/interim/kartasalo_liver_ds16.npy"
BASE_DS = 16
logger = logging.getLogger("scale_diag")


def reduce_stack(stack: np.ndarray, factor: int) -> np.ndarray:
    """Block-mean downsample, cropping the remainder."""
    if factor == 1:
        return stack
    n, h, w = stack.shape
    h2, w2 = (h // factor) * factor, (w // factor) * factor
    s = stack[:, :h2, :w2].astype(np.float32)
    return s.reshape(n, h2 // factor, factor, w2 // factor, factor).mean(axis=(2, 4))


def point_transform(shape, angle, dy, dx, pts):
    h, w = shape
    cy, cx = (h - 1) / 2.0, (w - 1) / 2.0
    t = np.deg2rad(angle)
    c, s = np.cos(t), np.sin(t)
    y, x = pts[:, 0] - cy, pts[:, 1] - cx
    return np.column_stack([cy + y * c - x * s + dy, cx + y * s + x * c + dx])


def main() -> None:
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    if not CACHE.exists():
        raise SystemExit(f"need the cached stack at {CACHE}; run the registration first")

    stack = np.load(CACHE)
    fid = load_fiducials(DATA / "fiducialcoordinates_liver_observer1.txt")
    secs = sorted(fid["section"].unique())[:len(stack)]
    n = len(stack)

    # floor: rigid transform fitted straight to the landmarks
    floor = []
    for i in range(n - 1):
        A = fiducial_array(fid, secs[i]); B = fiducial_array(fid, secs[i + 1])
        Ac, Bc = A - A.mean(0), B - B.mean(0)
        U, _, Vt = np.linalg.svd(Bc.T @ Ac)
        R = U @ Vt
        if np.linalg.det(R) < 0:
            U[:, -1] *= -1; R = U @ Vt
        floor.append(np.linalg.norm(Ac - Bc @ R.T, axis=1).mean() * NATIVE_MPP_UM)
    floor = float(np.mean(floor))

    print(f"stack {stack.shape} at {NATIVE_MPP_UM * BASE_DS:.2f} um/px")
    print(f"rigid floor (Procrustes on the fiducials): {floor:.0f} um\n")
    print(f"{'extra':>6} {'um/px':>7} {'shape':>12} {'corr med':>9} "
          f"{'TRE mean':>9} {'TRE med':>8} {'vs floor':>9}")
    print("-" * 66)

    rows = []
    # extra=1 is the production run already measured at 2544 um; re-running it
    # here would cost another 3.3 hours to learn nothing new, so start coarse.
    for extra in (4, 6, 8, 11, 16, 22):
        mpp = NATIVE_MPP_UM * BASE_DS * extra
        s = reduce_stack(stack, extra)
        if min(s.shape[1:]) < 24:
            continue
        cfg = RegistrationConfig()
        regs, params = register_stack([s[i] for i in range(n)], cfg, elastic=False)

        corr = np.array([p.get("correlation", np.nan) for p in params], dtype=float)
        moved = []
        for i, sec in enumerate(secs):
            p = params[i]
            moved.append(point_transform(s.shape[1:], p.get("angle", 0.0),
                                         p.get("dy", 0.0), p.get("dx", 0.0),
                                         fiducial_array(fid, sec, BASE_DS * extra)))
        d = [np.linalg.norm(moved[i] - moved[i + 1], axis=1).mean() * mpp
             for i in range(n - 1)]
        d = np.array(d)
        rows.append({"extra_downsample": extra, "mpp_um": round(mpp, 2),
                     "height": s.shape[1], "width": s.shape[2],
                     "corr_median": round(float(np.nanmedian(corr)), 4),
                     "tre_mean_um": round(float(d.mean()), 1),
                     "tre_median_um": round(float(np.median(d)), 1),
                     "ratio_to_floor": round(float(d.mean() / floor), 2)})
        r = rows[-1]
        print(f"{extra:>6} {mpp:>7.2f} {str(s.shape[1:]):>12} "
              f"{r['corr_median']:>9.4f} {r['tre_mean_um']:>9.1f} "
              f"{r['tre_median_um']:>8.1f} {r['ratio_to_floor']:>8.2f}x")

    df = pd.DataFrame(rows)
    out = ROOT / "results/kartasalo/scale_diagnostic.csv"
    df.to_csv(out, index=False)
    best = df.loc[df["tre_mean_um"].idxmin()]
    print(f"\nbest: {best['mpp_um']:.1f} um/px -> TRE {best['tre_mean_um']:.0f} um "
          f"({best['ratio_to_floor']:.2f}x the floor)")
    print(f"paper specifies 80 um/px for the rigid stage; production run used "
          f"{NATIVE_MPP_UM * BASE_DS:.2f}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
