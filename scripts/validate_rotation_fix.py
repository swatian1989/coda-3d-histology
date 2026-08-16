#!/usr/bin/env python
"""Head-to-head: Radon rotation estimation versus direct search, against truth.

    python scripts/validate_rotation_fix.py [--pairs 20] [--extra 11]

The ground truth is the rotation implied by the four operator-annotated
fiducials, recovered by Procrustes. Both estimators see identical preprocessed
inputs at identical scale, so the comparison isolates the estimator.

Sign conventions between an image-rotation routine and a Procrustes rotation
are ambiguous, and guessing wrong would make a working estimator look broken.
So each estimator is scored under BOTH conventions and the better one is
reported, with the convention named. That is applied equally to both, so it
cannot flatter one of them.

This is a gate, not a formality. If the search estimator does not beat Radon
here there is no reason to re-run anything.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coda_my.loaders.kartasalo import (  # noqa: E402
    NATIVE_MPP_UM, fiducial_array, load_fiducials,
)
from coda_my.registration import (  # noqa: E402
    RegistrationConfig, estimate_rotation, preprocess,
)
from coda_my.registration_fix import SearchConfig, estimate_rigid_search  # noqa: E402

DATA = ROOT / "data/raw/kartasalo/extracted/Data_to_IDA"
CACHE = ROOT / "data/interim/kartasalo_liver_ds16.npy"
BASE_DS = 16


def blockmean(a: np.ndarray, f: int) -> np.ndarray:
    if f == 1:
        return a.astype(np.float32)
    h, w = a.shape
    h2, w2 = (h // f) * f, (w // f) * f
    return a[:h2, :w2].astype(np.float32).reshape(h2 // f, f, w2 // f, f).mean(axis=(1, 3))


def procrustes_angle(A: np.ndarray, B: np.ndarray) -> float:
    Ac, Bc = A - A.mean(0), B - B.mean(0)
    U, _, Vt = np.linalg.svd(Bc.T @ Ac)
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    return float(np.degrees(np.arctan2(R[1, 0], R[0, 0])))


def wrap(d: np.ndarray) -> np.ndarray:
    return np.abs((d + 180.0) % 360.0 - 180.0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=20)
    ap.add_argument("--extra", type=int, default=11, help="extra downsample of the ds16 cache")
    ap.add_argument("--max-abs-deg", type=float, default=180.0,
                    help="restrict the angular search; 180 = unrestricted")
    args = ap.parse_args()

    stack = np.load(CACHE)
    fid = load_fiducials(DATA / "fiducialcoordinates_liver_observer1.txt")
    secs = sorted(fid["section"].unique())
    cfg = RegistrationConfig()
    search = SearchConfig(max_abs_deg=args.max_abs_deg)
    mpp = NATIVE_MPP_UM * BASE_DS * args.extra
    n = min(args.pairs, len(stack) - 1)

    rows = []
    t_radon = t_search = 0.0
    for i in range(n):
        a = preprocess(blockmean(stack[i], args.extra), cfg)
        b = preprocess(blockmean(stack[i + 1], args.extra), cfg)
        truth = procrustes_angle(fiducial_array(fid, secs[i]),
                                 fiducial_array(fid, secs[i + 1]))

        t0 = time.time()
        rad = estimate_rotation(a, b, cfg)
        t_radon += time.time() - t0

        t0 = time.time()
        res = estimate_rigid_search(a, b, cfg, search)
        t_search += time.time() - t0

        rows.append({"pair": f"{secs[i]}-{secs[i+1]}", "truth_deg": truth,
                     "radon_deg": rad, "search_deg": res["angle"],
                     "search_corr": res["correlation"]})

    df = pd.DataFrame(rows)
    print(f"Kartasalo liver, {n} consecutive pairs at {mpp:.1f} um/px "
          f"(shape {blockmean(stack[0], args.extra).shape})\n")

    summary = []
    for name, col in (("Radon (registration.py)", "radon_deg"),
                      ("direct search (new)", "search_deg")):
        errs = {}
        for sign, lab in ((1.0, "+"), (-1.0, "-")):
            errs[lab] = wrap(sign * df[col].to_numpy() - df["truth_deg"].to_numpy())
        lab = min(errs, key=lambda k: errs[k].mean())
        e = errs[lab]
        summary.append({"estimator": name, "convention": lab,
                        "mean_abs_err_deg": round(float(e.mean()), 2),
                        "median_abs_err_deg": round(float(np.median(e)), 2),
                        "within_5deg": int((e < 5).sum()),
                        "within_10deg": int((e < 10).sum()), "n": len(e)})
        df[f"{col}_err"] = e

    s = pd.DataFrame(summary)
    print(s.to_string(index=False))
    print(f"\ntiming: Radon {t_radon/n:.2f} s/pair, search {t_search/n:.2f} s/pair "
          f"({t_search/max(t_radon,1e-9):.1f}x)")
    print("chance level for a uniform guess is about 90 deg\n")

    out = ROOT / "results/kartasalo/rotation_estimator_comparison.csv"
    df.to_csv(out, index=False)
    s.to_csv(ROOT / "results/kartasalo/rotation_estimator_summary.csv", index=False)

    radon_m = s.loc[0, "mean_abs_err_deg"]
    search_m = s.loc[1, "mean_abs_err_deg"]
    print(f"{'PASS' if search_m < radon_m else 'FAIL'}: search {search_m:.1f} deg vs "
          f"Radon {radon_m:.1f} deg mean absolute error")
    if search_m >= radon_m:
        print("The replacement does not beat the original. Do not re-run the stack.")
    elif search_m > 10:
        print("Better, but still above 10 deg. A re-run may not reach the 727 um "
              "do-nothing baseline; treat any improvement as partial.")
    else:
        print("Good enough to justify re-running the full stack.")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
