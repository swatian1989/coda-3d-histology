"""One function per report table (T1-T14).

Each returns {"id","title","source","csv_path","df","caption"} and writes a
CSV to results/tables/. A table that cannot be built returns a one-row frame
naming the missing input, so nothing is silently dropped.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "results"


def _save(df: pd.DataFrame, name: str, tdir: str) -> str:
    out = Path(tdir); out.mkdir(parents=True, exist_ok=True)
    p = out / f"{name}.csv"; df.to_csv(p, index=False)
    return str(p)


def _ph(tid, title, needs, unblocks, tdir, name) -> dict:
    df = pd.DataFrame([{"status": "MISSING DATA", "needs": needs, "unblocks": unblocks}])
    return {"id": tid, "title": title, "source": "MISSING DATA",
            "csv_path": _save(df, name, tdir), "df": df,
            "caption": f"Placeholder. Requires {needs}."}


def _load(name):
    p = RESULTS / name
    return pd.read_csv(p) if p.exists() else None


def _cfg():
    return yaml.safe_load((ROOT / "config/coda_params.yaml").read_text())


# ==================================================================== T1-T3


def t1_dataset_inventory(tdir: str) -> dict:
    rows = [
        dict(dataset="Kartasalo mouse prostate", n="260 serial sections", species="mouse",
             tissue="prostate", serial_depth="260 consecutive, 5 um",
             stages_supported="1-7 (only source of 3D)",
             accession="RegBenchmark (software only); images require author request",
             status="NOT ACQUIRED"),
        dict(dataset="Kartasalo mouse liver", n="47 serial sections", species="mouse",
             tissue="liver", serial_depth="47 consecutive",
             stages_supported="1-7, laser-cut holes as independent check",
             accession="as above", status="NOT ACQUIRED"),
        dict(dataset="ACROBAT", n="4,212 WSI / 1,153 patients", species="human",
             tissue="breast", serial_depth="same block, NOT consecutive",
             stages_supported="1,2 (H&E to IHC), 3, 4, 7",
             accession="researchdata.se 2022-190-1 (data use agreement)",
             status="NOT ACQUIRED"),
        dict(dataset="TCGA-BRCA", n="~1,100 WSI", species="human", tissue="breast",
             serial_depth="single sections, different blocks",
             stages_supported="3, 4, 7", accession="GDC portal", status="NOT ACQUIRED"),
        dict(dataset="USM breast IHC", n="234 field-of-view captures", species="human",
             tissue="breast", serial_depth="single fields, no serial depth",
             stages_supported="marker quantification, HER2 membrane, spatial statistics",
             accession="institutional (USM Kota Bharu)", status="ACQUIRED, analysed"),
    ]
    df = pd.DataFrame(rows)
    return {"id": "T1", "title": "Dataset inventory", "source": "STATUS",
            "csv_path": _save(df, "T1_dataset_inventory", tdir), "df": df,
            "caption": "Every dataset the design calls for, the stages it can support, "
                       "and its acquisition status. Only the USM captures are in hand, "
                       "which is why Arms A and B appear as placeholders throughout."}


def t2_locked_parameters(tdir: str) -> dict:
    cfg = _cfg()

    def flat(d, pre=""):
        out = {}
        for k, v in (d or {}).items():
            key = f"{pre}{k}"
            if isinstance(v, dict):
                out.update(flat(v, key + "."))
            else:
                out[key] = v
        return out

    lk = flat(cfg["locked"])
    df = pd.DataFrame([{"section": k.split(".")[0], "parameter": k, "value": str(v)}
                       for k, v in lk.items()]).sort_values(["section", "parameter"])
    return {"id": "T2", "title": f"All {len(df)} locked parameters",
            "source": "CONFIG (SHA-256 verified)",
            "csv_path": _save(df, "T2_locked_parameters", tdir), "df": df,
            "caption": f"The {len(df)} parameters transcribed from the CODA Online "
                       f"Methods, grouped by section. guard.verify() hashes this block "
                       f"and fails the run naming any key that drifts."}


def t3_deviations(tdir: str) -> dict:
    devs = _cfg().get("deviations", [])
    df = pd.DataFrame(devs)
    return {"id": "T3", "title": "Declared deviations from the published method",
            "source": "CONFIG", "csv_path": _save(df, "T3_deviations", tdir), "df": df,
            "caption": "Every deliberate departure from the paper, with reason and "
                       "expected impact. A declared deviation is a methods sentence; "
                       "a silent one is an irreproducible result."}


# ==================================================================== T4-T10 Arms A/B


def t4_registration_accuracy(tdir):
    return _ph("T4", "Registration accuracy vs the published benchmark",
               "Kartasalo stack and fiducials", "Arm A stage 2", tdir, "T4_registration")


def t5_z_skip(tdir):
    return _ph("T5", "z-skip validation", "Kartasalo serial stack", "Arm A stage 2",
               tdir, "T5_z_skip")


def t6_cell_detection(tdir):
    return _ph("T6", "Cell detection precision, recall, F1",
               "manual annotations from two annotators", "stage 3 validation",
               tdir, "T6_cell_detection")


def t7_segmentation(tdir):
    return _ph("T7", "Segmentation per-class precision and recall vs the 90% gate",
               "annotated H&E and a GPU", "stage 4", tdir, "T7_segmentation")


def t8_composition_density(tdir):
    return _ph("T8", "Tissue composition and 3D cell density per class",
               "reconstructed volume plus measured nuclear diameters", "Arm A stage 6",
               tdir, "T8_composition")


def t9_overcounting(tdir):
    return _ph("T9", "Overcounting, 2D vs 3D per section",
               "reconstructed volume with 3D connectivity", "Arm A stage 6",
               tdir, "T9_overcounting")


def t10_object_morphology(tdir):
    return _ph("T10", "Per-object 3D morphology summary", "reconstructed volume",
               "Arm A stage 6", tdir, "T10_morphology")


# ==================================================================== T11-T13 Arm C


def t11_usm_qc(tdir: str) -> dict:
    qc = _load("usm_qc.csv")
    if qc is None:
        return _ph("T11", "USM IHC quality control", "results/usm_qc.csv", "Arm C",
                   tdir, "T11_usm_qc")
    df = qc[["filename", "marker", "mpp_um_per_px", "magnification_tier",
             "counterstain_grade", "counterstain_fraction",
             "percent_positive_reportable", "note"]].copy()
    df = df.rename(columns={"percent_positive_reportable": "reportable"})
    n_abs = int((qc["counterstain_grade"] == "absent").sum())
    return {"id": "T11", "title": "Arm C quality control, per image",
            "source": f"REAL (USM IHC, n={len(df)})",
            "csv_path": _save(df, "T11_usm_qc", tdir), "df": df,
            "caption": f"Filename, recovered microns per pixel, counterstain grade and "
                       f"whether a percentage is reportable, for all {len(df)} images. "
                       f"{n_abs} are graded counterstain absent: no negative nuclei are "
                       f"visible, so there is no denominator and percent positive is "
                       f"withheld rather than back-calculated from DAB area."}


def t12_marker_results(tdir: str) -> dict:
    mk = _load("usm_markers.csv")
    if mk is None:
        return _ph("T12", "Marker results per image", "results/usm_markers.csv",
                   "Arm C", tdir, "T12_markers")
    keep = [c for c in ["filename", "marker", "mpp_um_per_px", "counterstain_grade",
                        "n_nuclei_detected", "n_positive", "positive_density_per_mm2",
                        "percent_positive", "percent_reportable",
                        "n_enclosed_cells", "mean_membrane_completeness",
                        "median_cell_area_um2", "note"] if c in mk.columns]
    df = mk[keep]
    return {"id": "T12", "title": "Marker results per image",
            "source": f"REAL (USM IHC, n={len(df)})",
            "csv_path": _save(df, "T12_markers", tdir), "df": df,
            "caption": "Per-image marker quantification. ER, PR and Ki67 carry "
                       "per-nucleus DAB results; HER2 carries membrane completeness only "
                       "and no percentage, because per-nucleus scoring of a membranous "
                       "marker is the wrong operation."}


def t13_ki67_hotspot(tdir: str) -> dict:
    mk = _load("usm_markers.csv")
    if mk is None or "ki67_hotspot_minus_average" not in (mk.columns if mk is not None else []):
        return _ph("T13", "Ki67 hotspot vs average and the gap",
                   "results/usm_markers.csv with Ki67 rows", "Arm C Ki67",
                   tdir, "T13_ki67")
    k = mk[(mk["marker"] == "Ki67") & mk["ki67_hotspot_minus_average"].notna()
           & (mk["percent_reportable"] == True)].copy()  # noqa: E712
    k["crosses_20pc_cutoff"] = (k["ki67_average_percent"] < 20) & \
                               (k["ki67_hotspot_percent"] >= 20)
    df = k[["filename", "mpp_um_per_px", "n_nuclei_detected",
            "ki67_average_percent", "ki67_hotspot_percent",
            "ki67_hotspot_minus_average", "ki67_n_windows",
            "crosses_20pc_cutoff"]] if "ki67_n_windows" in k.columns else \
        k[["filename", "mpp_um_per_px", "n_nuclei_detected", "ki67_average_percent",
           "ki67_hotspot_percent", "ki67_hotspot_minus_average", "crosses_20pc_cutoff"]]
    gap = k["ki67_hotspot_minus_average"]
    flip = int(k["crosses_20pc_cutoff"].sum())
    return {"id": "T13", "title": "Ki67 hotspot versus average, per image",
            "source": f"REAL (USM IHC, Ki67, n={len(df)})",
            "csv_path": _save(df, "T13_ki67", tdir), "df": df,
            "caption": f"Average and hotspot Ki67 index for each of {len(df)} images with "
                       f"a valid denominator, and the gap between them. Median gap "
                       f"{gap.median():.1f} percentage points, maximum {gap.max():.1f}. "
                       f"On {flip} images ({100*flip/len(df):.0f} percent) the average is "
                       f"below the 20 percent cutoff while the hotspot is at or above it, "
                       f"so the choice of scoring method alone changes the treatment "
                       f"decision."}


# ==================================================================== T14


def t14_stage_applicability(tdir: str) -> dict:
    cfg = _cfg()
    rows = []
    stage_names = {1: "nonlinear registration", 2: "registration QC",
                   3: "cell detection", 4: "semantic segmentation",
                   5: "3D reconstruction", 6: "quantification and connectivity",
                   7: "fiber alignment"}
    plan = {
        "A Kartasalo": {s: ("blocked", "dataset not acquired; images require author "
                                       "request") for s in range(1, 8)},
        "B ACROBAT": {**{s: ("blocked", "dataset not acquired; data use agreement")
                         for s in (1, 2, 3, 4, 7)},
                      5: ("blocked", "sections are not consecutive; no volume"),
                      6: ("blocked", "depends on stage 5")},
        "C USM IHC": {1: ("blocked", "single fields, nothing to register"),
                      2: ("blocked", "no registration to assess"),
                      3: ("ran", "positives scored; negatives only where counterstain "
                                 "is present"),
                      4: ("blocked", "a field shows one tissue type, not ten"),
                      5: ("blocked", "no serial sections"),
                      6: ("blocked", "no volume"),
                      7: ("blocked", "needs an eosin channel; DAB-IHC has none")},
    }
    for arm, stages in plan.items():
        for s in sorted(stages):
            status, reason = stages[s]
            rows.append({"arm": arm, "stage": s, "stage_name": stage_names[s],
                         "status": status, "reason": reason})
    extra = [
        {"arm": "C USM IHC", "stage": "+", "stage_name": "DAB marker quantification",
         "status": "ran", "reason": "225 of 234 images"},
        {"arm": "C USM IHC", "stage": "+", "stage_name": "HER2 membrane completeness",
         "status": "ran", "reason": "53 images"},
        {"arm": "C USM IHC", "stage": "+", "stage_name": "spatial statistics of positives",
         "status": "ran", "reason": "64 point patterns, border corrected"},
    ]
    df = pd.DataFrame(rows + extra)
    n_ran = int((df["status"] == "ran").sum())
    return {"id": "T14", "title": "Stage by arm: ran, blocked, and why",
            "source": "RUN STATUS",
            "csv_path": _save(df, "T14_applicability", tdir), "df": df,
            "caption": f"Every CODA stage against every arm, with the reason for each "
                       f"blocked cell. {n_ran} of {len(df)} ran. Blocking is decided by "
                       f"the data: registering non-serial sections produces a transform "
                       f"and a correlation number, and neither means anything."}


ALL_TABLES = [
    t1_dataset_inventory, t2_locked_parameters, t3_deviations,
    t4_registration_accuracy, t5_z_skip, t6_cell_detection, t7_segmentation,
    t8_composition_density, t9_overcounting, t10_object_morphology,
    t11_usm_qc, t12_marker_results, t13_ki67_hotspot, t14_stage_applicability,
]
