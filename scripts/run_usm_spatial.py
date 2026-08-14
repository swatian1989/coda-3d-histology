#!/usr/bin/env python
"""Arm C, step 3: spatial statistics of Ki67-positive nuclei.

    python scripts/run_usm_spatial.py

The question this answers is the one a percentage cannot: a clustered 18
percent and a dispersed 18 percent receive the same Ki67 score and the same
treatment decision. These statistics separate them.

Runs the border-corrected estimators from the canvas-brca spatial feature
library, unmodified, on the positive-nucleus point pattern of each image:

    Ripley K and L        border (reduced sample) corrected
    Clark-Evans           Donnelly perimeter corrected
    quadrat dispersion    variance to mean ratio
    KDE hotspot CV        coefficient of variation of kernel density

Border correction matters more here than usual. These are field-of-view
captures, so a large fraction of the field is within one analysis radius of an
edge, and an uncorrected estimator reads that missing area as reduced
clustering.

RADIUS LIMIT. Ripley's K beyond about a quarter of the field width is
unreliable even corrected, so the radii are capped per image at
min(field width, field height) / 4 and the value used is recorded alongside
every statistic. Statistics are computed only where the positive count is
sufficient; below that they are written null with a reason rather than
returned as a number that cannot be supported.
"""
from __future__ import annotations

import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
# The spatial estimators live in the canvas-brca package, which the protocol
# specifies as the consumer of to_point_pattern(). Imported, never copied.
CANVAS_SRC = Path(r"D:\UCAS project\final completed pipline canvas-brca (2)\canvas-brca\src")
sys.path.insert(0, str(CANVAS_SRC))

from canvas_brca.stage5_features.spatial_features import (  # noqa: E402
    FeatureConfig, _clark_evans_donnelly, _kde_cv, _quadrat_vmr, _ripley_k_border,
)

MIN_POSITIVE = 30        # below this, a point pattern statistic is not supportable


def analyse(pts: np.ndarray, cfg: FeatureConfig) -> dict:
    """Border-corrected spatial statistics for one positive-nucleus pattern."""
    from scipy.spatial import cKDTree

    span = pts.max(axis=0) - pts.min(axis=0)
    r_max = float(min(span) / 4.0)                    # the FOV honesty limit
    radii = [r for r in cfg.ripley_radii_um if r <= r_max] or [r_max / 2]

    tree = cKDTree(pts)
    k_vals, l_vals = [], []
    for r in radii:
        k = _ripley_k_border(pts, tree, pts, r)
        k_vals.append(k)
        l_vals.append(np.sqrt(k / np.pi) - r if np.isfinite(k) and k >= 0 else np.nan)

    nnd, _ = tree.query(pts, k=2)
    nnd = nnd[:, 1]
    area = float(max(span[0], 1e-9) * max(span[1], 1e-9))
    lam = len(pts) / area

    return {
        "n_positive": len(pts),
        "r_max_um": round(r_max, 1),
        "radii_used_um": ";".join(str(int(r)) for r in radii),
        "ripley_k_mean": float(np.nanmean(k_vals)),
        "ripley_l_mean": float(np.nanmean(l_vals)),
        "clark_evans_donnelly": _clark_evans_donnelly(nnd, lam, pts, len(pts)),
        "quadrat_vmr": _quadrat_vmr(pts, pts, cfg),
        "kde_hotspot_cv": _kde_cv(pts),
        "mean_nn_distance_um": float(np.mean(nnd)),
        "positive_density_per_mm2": float(len(pts) / (area / 1e6)) if area else np.nan,
    }


def main() -> None:
    logging.basicConfig(level=logging.ERROR)
    warnings.filterwarnings("ignore")

    pp = pd.read_parquet(ROOT / "data/interim/usm_point_patterns.parquet")
    markers = pd.read_csv(ROOT / "results/usm_markers.csv")
    cfg = FeatureConfig()

    rows = []
    for (fn, marker), g in pp.groupby(["filename", "marker"]):
        pos = g[g["habitat"] == 1][["x_um", "y_um"]].to_numpy(float)
        rec = {"filename": fn, "marker": marker, "n_nuclei_total": len(g)}
        if len(pos) < MIN_POSITIVE:
            rec.update({"n_positive": len(pos),
                        "note": f"fewer than {MIN_POSITIVE} positive nuclei; "
                                f"spatial statistics not supportable"})
        else:
            rec.update(analyse(pos, cfg))
            rec["note"] = ""
        rows.append(rec)

    df = pd.DataFrame(rows)
    m = markers[["filename", "ki67_average_percent", "ki67_hotspot_percent",
                 "ki67_hotspot_minus_average"]] \
        if "ki67_average_percent" in markers else markers[["filename"]]
    df = df.merge(m, on="filename", how="left")
    out = ROOT / "results/usm_spatial.csv"
    df.to_csv(out, index=False)

    print(f"wrote {out}  ({len(df)} images)")
    ok = df[df["note"] == ""]
    print(f"statistics computed on {len(ok)}; {len(df)-len(ok)} had too few positives\n")
    if not len(ok):
        return

    k = ok[ok["marker"] == "Ki67"]
    if len(k):
        print(f"Ki67, n={len(k)} images")
        print("  Clark-Evans (Donnelly corrected): <1 clustered, 1 random, >1 regular")
        ce = k["clark_evans_donnelly"].dropna()
        print(f"    median {ce.median():.3f}  IQR {ce.quantile(.25):.3f} to "
              f"{ce.quantile(.75):.3f}  range {ce.min():.3f} to {ce.max():.3f}")
        print(f"    clustered (<1): {int((ce<1).sum())} of {len(ce)} images")
        for col, lbl in [("quadrat_vmr", "quadrat VMR (1 = Poisson)"),
                         ("kde_hotspot_cv", "KDE hotspot CV"),
                         ("ripley_l_mean", "Ripley L (0 = CSR)")]:
            v = k[col].dropna()
            if len(v):
                print(f"  {lbl}: median {v.median():.3f}, IQR {v.quantile(.25):.3f} "
                      f"to {v.quantile(.75):.3f}")


if __name__ == "__main__":
    main()
