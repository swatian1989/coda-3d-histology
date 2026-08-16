#!/usr/bin/env python
# ============================================================================
#  CODA-BRCA-MY  --  ONE SCRIPT, EVERY STAGE, START TO FINISH
# ============================================================================
#
#  Runs in VS Code and in Google Colab without modification. It detects which
#  it is in and adapts paths, dependency installation and data download.
#
#      VS CODE / LOCAL      python RUN_EVERYTHING.py
#                           python RUN_EVERYTHING.py --steps 1,2,3
#                           python RUN_EVERYTHING.py --list
#
#      GOOGLE COLAB         !git clone <repo> && cd coda-brca-my
#                           !python RUN_EVERYTHING.py --colab
#                       or paste into a cell:
#                           %run RUN_EVERYTHING.py
#
#  WHAT IT DOES, IN ORDER, AND WHY EACH STEP EXISTS
#
#    0  environment        detect Colab, install what is missing, set paths
#    1  fetch data         stream the 63.79 GB archive, keep only the liver
#                          stack (~15 GB). Skipped if already present.
#    2  verify data        decode every section fully; a truncated TIFF from an
#                          interrupted download must fail here, not silently
#                          later
#    3  registration       two-scale rigid, corrected rotation estimator
#    4  registration QC    TRE, ATRE, and the two floors that bound them
#    5  benchmark          rank against the algorithms published on this data
#    6  cell detection     two independent detectors, agreement only
#    7  3D reconstruction  stack the sections into a volume
#    8  connectivity       the 2D versus 3D overcounting measurement
#    9  fibre alignment    anisotropy, with the shuffle control
#   10  Arm C IHC          the breast marker analysis, independent of 1 to 9
#   11  report             every figure and table, three formats
#   12  manuscript         publication draft with PubMed-verified references
#
#  Steps are idempotent. Each checks for its own output and skips unless
#  --force is given, so an interrupted run resumes rather than restarting.
#
#  HARDWARE
#    CPU only throughout. Peak RAM is about 2 GB, peak disk about 20 GB.
#    No GPU is used or required. Stages 4 (segmentation) of the original
#    method is not implemented because it needs annotated training data that
#    does not exist for this material; see the report.
# ============================================================================

from __future__ import annotations

import argparse
import importlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# STEP 0: environment
# ---------------------------------------------------------------------------
# Colab gives a clean VM every session, so dependencies must be installed and
# the working directory located rather than assumed. Locally, neither is true
# and reinstalling would be rude, so both are conditional.

IN_COLAB = "google.colab" in sys.modules or os.path.exists("/content")

REQUIRED = {
    "numpy": "numpy", "pandas": "pandas", "scipy": "scipy",
    "skimage": "scikit-image", "PIL": "pillow", "matplotlib": "matplotlib",
    "yaml": "pyyaml", "docx": "python-docx", "tabulate": "tabulate",
}


def log(msg: str, level: str = "INFO") -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {level:5s} {msg}", flush=True)


def ensure_deps() -> None:
    missing = []
    for mod, pkg in REQUIRED.items():
        try:
            importlib.import_module(mod)
        except ImportError:
            missing.append(pkg)
    if not missing:
        log("all dependencies present")
        return
    log(f"installing: {', '.join(missing)}")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", *missing])


def find_root() -> Path:
    """Locate the project root from wherever the script was launched."""
    here = Path(__file__).resolve().parent
    if (here / "src" / "coda_my").is_dir():
        return here
    for cand in (Path.cwd(), *Path.cwd().parents):
        if (cand / "src" / "coda_my").is_dir():
            return cand
    raise SystemExit(
        "cannot find the project root (a directory containing src/coda_my). "
        "In Colab, git clone the repository first and run from inside it."
    )


# ---------------------------------------------------------------------------
# step definitions
# ---------------------------------------------------------------------------
# Each step is (number, name, output that proves it ran, command, note).
# `output` is what makes the script resumable: if it exists, the step is done.

def build_steps(root: Path) -> list[dict]:
    py = sys.executable
    S = lambda *a: [py, *a]  # noqa: E731
    return [
        dict(n=1, name="fetch liver stack",
             out=root / "data/raw/kartasalo/extracted/Data_to_IDA/liver/047.tif",
             cmd=S("scripts/fetch_kartasalo_liver.py"),
             note="streams 63.79 GB, writes ~15 GB. The download service ignores "
                  "HTTP range requests, so this cannot resume: if it dies, it "
                  "restarts. Allow several hours on a home connection.",
             heavy=True),
        dict(n=2, name="verify sections decode",
             out=root / "results/kartasalo/_verify.ok",
             cmd=None, note="full decode of every section; catches truncation"),
        dict(n=3, name="registration (corrected estimator)",
             out=root / "results/kartasalo/summary_ds16fix_rigid.json",
             cmd=S("scripts/run_kartasalo_registration_fixed.py", "--no-elastic"),
             note="rigid at 81 um/px, elastic disabled because it was measured "
                  "to raise error on this series"),
        dict(n=4, name="rotation estimator validation",
             out=root / "results/kartasalo/rotation_estimator_summary.csv",
             cmd=S("scripts/validate_rotation_fix.py", "--pairs", "20",
                   "--max-abs-deg", "45"),
             note="head to head against the published Radon estimator"),
        dict(n=5, name="3D reconstruction and overcounting",
             out=root / "results/kartasalo/stage6_summary.json",
             cmd=S("scripts/run_kartasalo_3d_reconstruction.py"),
             note="builds the volume, counts objects in 2D and in 3D"),
        dict(n=6, name="cell detection agreement",
             out=root / "results/kartasalo/stage3_summary.json",
             cmd=S("scripts/run_kartasalo_cell_detection.py", "--n-fields", "12"),
             note="two detectors compared to each other, not to ground truth"),
        dict(n=7, name="fibre alignment",
             out=root / "results/kartasalo/stage7_summary.json",
             cmd=S("scripts/run_kartasalo_fibers.py", "--n-planes", "25"),
             note="anisotropy in the cutting plane; the orthogonal comparison "
                  "fails its shuffle control and is withheld"),
        dict(n=8, name="Arm C: image QC",
             out=root / "results/usm_qc.csv",
             cmd=S("scripts/run_usm_qc.py"),
             note="scale bar recovery and counterstain grading"),
        dict(n=9, name="Arm C: marker quantification",
             out=root / "results/usm_markers.csv",
             cmd=S("scripts/run_usm_markers.py"),
             note="ER, PR, Ki67 scoring; HER2 by membrane completeness"),
        dict(n=10, name="Arm C: spatial statistics",
             out=root / "results/usm_spatial.csv",
             cmd=S("scripts/run_usm_spatial.py"),
             note="border-corrected estimators on the positive point pattern"),
        dict(n=11, name="Arm C: stereological 3D correction",
             out=root / "results/usm_3d_extrapolation.csv",
             cmd=S("scripts/run_usm_3d_extrapolation.py"),
             note="measures nuclear diameter here rather than inheriting it"),
        dict(n=12, name="verify references against PubMed",
             out=root / "manuscript/references_verified.json",
             cmd=S("scripts/verify_references.py"),
             note="needs network; never invents a PMID"),
        dict(n=13, name="build report",
             out=root / "reports/analysis_report.html",
             cmd=S("scripts/run_report.py"),
             note="all figures and tables, markdown, HTML and Word"),
        dict(n=14, name="build manuscript",
             out=root / "manuscript/manuscript.docx",
             cmd=S("scripts/run_manuscript.py"),
             note="publication draft with numbered references"),
    ]


def verify_sections(root: Path) -> bool:
    """STEP 2. Decode every section fully.

    A download that was interrupted leaves a TIFF that opens and reports its
    size but fails partway through decoding. Opening is not enough; the pixels
    have to be pulled. This is cheap insurance against a corrupt section
    silently degrading a registration hours later.
    """
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    d = root / "data/raw/kartasalo/extracted/Data_to_IDA/liver"
    paths = sorted(d.glob("*.tif"))
    if not paths:
        log("no sections found; run step 1 first", "ERROR")
        return False
    bad = []
    for p in paths:
        try:
            im = Image.open(p)
            im.load()
        except Exception as exc:
            bad.append((p.name, str(exc)[:60]))
    log(f"{len(paths) - len(bad)}/{len(paths)} sections decode cleanly")
    for name, err in bad:
        log(f"  CORRUPT {name}: {err}", "ERROR")
    if bad:
        return False
    (root / "results/kartasalo").mkdir(parents=True, exist_ok=True)
    (root / "results/kartasalo/_verify.ok").write_text(
        json.dumps({"n_sections": len(paths), "all_decoded": True}, indent=2))
    return True


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Run the whole CODA-BRCA-MY pipeline, or part of it.")
    ap.add_argument("--steps", help="comma separated step numbers, e.g. 3,5,7")
    ap.add_argument("--from-step", type=int, help="start here and continue")
    ap.add_argument("--skip-heavy", action="store_true",
                    help="skip the multi-hour data download")
    ap.add_argument("--force", action="store_true",
                    help="re-run steps whose output already exists")
    ap.add_argument("--list", action="store_true", help="show the steps and exit")
    ap.add_argument("--colab", action="store_true", help="force Colab behaviour")
    args = ap.parse_args()

    colab = IN_COLAB or args.colab
    print("=" * 78)
    print("  CODA-BRCA-MY  |  full pipeline")
    print(f"  environment: {'Google Colab' if colab else 'local / VS Code'}")
    print("=" * 78)

    ensure_deps()
    root = find_root()
    os.chdir(root)
    sys.path.insert(0, str(root / "src"))
    log(f"project root: {root}")

    steps = build_steps(root)
    if args.list:
        for s in steps:
            print(f"  {s['n']:2d}. {s['name']:36s} -> {s['out'].name}")
            print(f"      {s['note']}")
        return

    wanted = None
    if args.steps:
        wanted = {int(x) for x in args.steps.split(",")}
    elif args.from_step:
        wanted = {s["n"] for s in steps if s["n"] >= args.from_step}

    ran = skipped = failed = 0
    for s in steps:
        if wanted is not None and s["n"] not in wanted:
            continue
        if s.get("heavy") and args.skip_heavy:
            log(f"step {s['n']} ({s['name']}) skipped by --skip-heavy")
            skipped += 1
            continue
        if s["out"].exists() and not args.force:
            log(f"step {s['n']} ({s['name']}) already done -> {s['out'].name}")
            skipped += 1
            continue

        print("\n" + "-" * 78)
        log(f"STEP {s['n']}: {s['name']}")
        print(f"        {s['note']}")
        print("-" * 78)
        t0 = time.time()

        if s["n"] == 2:
            ok = verify_sections(root)
        else:
            r = subprocess.run(s["cmd"], cwd=root)
            ok = r.returncode == 0

        dt = (time.time() - t0) / 60
        if ok and (s["out"].exists() or s["n"] == 2):
            log(f"step {s['n']} finished in {dt:.1f} min")
            ran += 1
        else:
            log(f"step {s['n']} FAILED after {dt:.1f} min", "ERROR")
            failed += 1
            if s["n"] in (1, 2, 3):
                log("later steps depend on this one; stopping", "ERROR")
                break

    print("\n" + "=" * 78)
    log(f"done. {ran} ran, {skipped} skipped, {failed} failed")
    for label, p in (("report ", root / "reports/analysis_report.html"),
                     ("manuscript", root / "manuscript/manuscript.docx")):
        if p.exists():
            log(f"{label}: {p}  ({p.stat().st_size / 1024:.0f} KB)")
    if colab and (root / "reports/analysis_report.html").exists():
        print("\nIn Colab, download the outputs with:")
        print("  from google.colab import files")
        print("  files.download('reports/analysis_report.html')")
        print("  files.download('manuscript/manuscript.docx')")


if __name__ == "__main__":
    main()
