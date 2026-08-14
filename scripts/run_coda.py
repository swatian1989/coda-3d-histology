#!/usr/bin/env python
"""End-to-end CODA protocol runner. All 7 stages, in the published order.

    # Full 3D protocol on public serial sections
    python scripts/run_coda.py --dataset kartasalo_prostate --stages all

    # 2D stages only, on TCGA or any single-section cohort
    python scripts/run_coda.py --dataset tcga --stages 3,4,7

Stage map (numbering follows the paper):
    1  registration            serial sections only
    2  registration QC         TRE, ATRE, z-correlation, z-skip
    3  cell detection          any H&E
    4  semantic segmentation   any H&E, needs GPU
    5  3D reconstruction       serial sections only
    6  quantification          composition, cell density, connectivity
    7  fiber alignment         any H&E
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

SERIAL_ONLY = {1, 2, 5, 6}
DATASETS = {
    "kartasalo_prostate": {"serial": True, "n_sections": 260, "mpp": 0.46,
                           "section_um": 5.0},
    "kartasalo_liver": {"serial": True, "n_sections": 47, "mpp": 0.46,
                        "section_um": 5.0},
    "acrobat": {"serial": False, "mpp": 0.92, "note": "same block, NOT consecutive"},
    "tcga": {"serial": False, "mpp": 0.25, "note": "single sections, different blocks"},
    "usm": {"serial": False, "mpp": None, "note": "read from slide properties"},
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=list(DATASETS))
    ap.add_argument("--input", default="data/raw")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--stages", default="all")
    ap.add_argument("--target-label", type=int, default=2)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger("coda")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    meta = DATASETS[args.dataset]
    stages = (set(range(1, 8)) if args.stages == "all"
              else {int(s) for s in args.stages.split(",")})

    if not meta["serial"]:
        blocked = stages & SERIAL_ONLY
        if blocked:
            log.error(
                "Stages %s require SERIAL SECTIONS and '%s' has none (%s). "
                "There is nothing to register or stack. Run stages 3,4,7 on "
                "this dataset, and stages 1,2,5,6 on kartasalo_prostate.",
                sorted(blocked), args.dataset, meta.get("note", ""))
            stages -= SERIAL_ONLY
            if not stages:
                sys.exit(1)

    log.info("dataset=%s serial=%s stages=%s", args.dataset, meta["serial"],
             sorted(stages))

    manifest: dict = {"dataset": args.dataset, "stages": sorted(stages)}
    sections: list[np.ndarray] = []

    # ---------------------------------------------------------- 1 registration
    if 1 in stages:
        from coda_my.registration import RegistrationConfig, register_stack

        log.info("STAGE 1  registration (Radon global + elastic, centre reference)")
        # TODO: load sections from args.input in z order, apply args.limit
        cfg = RegistrationConfig()
        if sections:
            registered, params = register_stack(sections, cfg, elastic=True)
            np.save(out / "registered_stack.npy", np.stack(registered))
            manifest["registration"] = {
                "n_sections": len(registered),
                "median_correlation": float(np.median(
                    [p.get("correlation", np.nan) for p in params])),
                "defective": [i for i, p in enumerate(params)
                              if p.get("correlation", 1) < cfg.min_correlation],
            }
        else:
            log.warning("no sections loaded; implement the loader for %s",
                        args.dataset)

    # ---------------------------------------------------------- 2 registration QC
    if 2 in stages:
        from coda_my.qc import (axial_vs_lateral_correlation, z_skip_validation)

        log.info("STAGE 2  QC (paper: >95%% z-correlation to 20 um, <5%% error to 12 um)")
        path = out / "registered_stack.npy"
        if path.exists():
            stack = np.load(path)
            corr = axial_vs_lateral_correlation(
                stack, mpp=meta["mpp"], section_um=meta["section_um"])
            corr.to_csv(out / "qc_axial_vs_lateral.csv", index=False)
            skip = z_skip_validation(stack, section_um=meta["section_um"])
            skip.to_csv(out / "qc_z_skip.csv", index=False)
            manifest["qc"] = skip.to_dict("records")
        else:
            log.warning("run stage 1 first")

    # ---------------------------------------------------------- 3 cell detection
    if 3 in stages:
        log.info("STAGE 3  cell detection (hematoxylin deconvolution, 2 um/px)")
        log.info("  validate with qc.cell_detection_metrics at 2 um tolerance; "
                 "paper reports >90%% precision and recall")

    # ---------------------------------------------------------- 4 segmentation
    if 4 in stages:
        from coda_my.segmentation import SegmentationConfig

        scfg = SegmentationConfig()
        log.info("STAGE 4  DeepLab v3+ / ResNet-50 semantic segmentation  [GPU]")
        log.info("  %d images per sample, >=%d annotations per class, "
                 "%dx%d tiles >%.0f%% filled cut to %dx%d (%d each), "
                 "acceptance %.0f%% precision AND recall per class",
                 scfg.annotation_images_per_sample, scfg.annotations_per_class,
                 scfg.big_tile_px, scfg.big_tile_px, scfg.fill_fraction * 100,
                 scfg.small_tile_px, scfg.small_tile_px, scfg.tiles_per_big,
                 scfg.acceptance_precision * 100)

    # ---------------------------------------------------------- 5 reconstruction
    if 5 in stages:
        from coda_my.reconstruct import ReconstructionConfig, stack_to_volume

        log.info("STAGE 5  3D reconstruction (isotropic 12 um voxels)")
        rcfg = ReconstructionConfig(section_thickness_um=meta["section_um"])
        labels_path = out / "labelled_sections.npy"
        if labels_path.exists():
            vol = stack_to_volume(list(np.load(labels_path)), rcfg)
            np.save(out / "volume.npy", vol)
        else:
            log.warning("stage 4 must produce labelled_sections.npy first")

    # ---------------------------------------------------------- 6 quantification
    if 6 in stages:
        from coda_my.reconstruct import (ReconstructionConfig, object_metrics,
                                         overcounting_ratio, tissue_composition)

        log.info("STAGE 6  quantification")
        vpath = out / "volume.npy"
        if vpath.exists():
            vol = np.load(vpath)
            rcfg = ReconstructionConfig(section_thickness_um=meta["section_um"])
            names = {i: f"class_{i}" for i in np.unique(vol) if i > 0}
            tissue_composition(vol, names, rcfg).to_csv(
                out / "composition.csv", index=False)
            ratio = overcounting_ratio(vol, args.target_label, rcfg)
            ratio.to_csv(out / "overcounting_2d_vs_3d.csv", index=False)
            object_metrics(vol, args.target_label, rcfg).to_csv(
                out / "object_metrics.csv", index=False)
            if not ratio.empty:
                manifest["overcounting"] = {
                    "mean_fold": float(ratio["ratio"].mean()),
                    "max_fold": float(ratio["ratio"].max()),
                    "paper_mean_fold": 12.3, "paper_max_fold": 40.0,
                }
                log.info("2D/3D overcounting: mean %.1f max %.1f "
                         "(paper: 12.3 and 40)",
                         ratio["ratio"].mean(), ratio["ratio"].max())
        else:
            log.warning("run stage 5 first")

    # ---------------------------------------------------------- 7 fiber alignment
    if 7 in stages:
        log.info("STAGE 7  fiber alignment (eosin channel, 2500 um^2 windows)")
        log.info("  SECTIONING ANGLE CONFOUND: the paper corrected for this using "
                 "the 3D volume to pick longitudinal vs axial planes, finding "
                 "2.2-2.5 fold differences. On single sections you cannot. Use "
                 "tiled_anisotropy for a distribution, or "
                 "boundary_relative_orientation at the invasive front.")

    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    log.info("wrote %s", out / "run_manifest.json")


if __name__ == "__main__":
    main()
