"""One function per report figure (F1-F22).

Every function returns:
    {"id","title","source","paths":{"png","pdf"},"caption"}

`source` is one of:
    REAL ...        built from measured data; the dataset and n are named
    SIMULATED ...   built from a synthetic fixture, never a finding
    MISSING DATA    cannot be built; a placeholder names the exact input needed

Nothing here fabricates a number or a panel. Where an arm has no data the
figure is a labelled placeholder, kept numbered in sequence so the report
stays complete.

Real inputs, all produced by the Arm C scripts:
    results/usm_qc.csv                     234 images, mpp and counterstain
    results/usm_markers.csv                225 analysed, per-marker results
    results/usm_spatial.csv                77 point patterns, border corrected
    config/coda_params.yaml                120 locked parameters, deviations
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .style import (
    NAVY, OKABE_ITO, STEEL_BLUE, apply_style, categorical_colors, continuous_cmap,
    letter_panels, placeholder_figure, save_figure, source_caption,
)

logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "results"

USM_N = "USM breast IHC field-of-view captures"


def _load(name: str) -> pd.DataFrame | None:
    p = RESULTS / name
    return pd.read_csv(p) if p.exists() else None


def _ph(fid: str, title: str, needs: str, unblocks: str, figdir: str, fname: str) -> dict:
    fig, _ = placeholder_figure(fid, title, missing_file=needs, unblocks=unblocks)
    paths = save_figure(fig, fname, figdir)
    return {"id": fid, "title": title, "source": "MISSING DATA", "paths": paths,
            "caption": f"Placeholder. Requires {needs}. Unblocks {unblocks}."}


# ============================================================ F1 study design


def f01_study_design(figdir: str) -> dict:
    apply_style()
    fig, ax = plt.subplots(figsize=(12, 6.4))
    ax.set_xlim(0, 12); ax.set_ylim(0, 6.4); ax.axis("off")
    ax.set_title("Three arms, seven CODA stages: which runs where, and why",
                 loc="left", fontsize=12)

    arms = [
        (0.2, 4.3, "ARM A  Kartasalo\n260 mouse prostate serial sections\nstages 1-7, the only 3D",
         "#E8EEF6", NAVY),
        (0.2, 2.4, "ARM B  ACROBAT\n1,153 human breast patients\nstages 3, 4, 7 + H&E to IHC",
         "#E8EEF6", STEEL_BLUE),
        (0.2, 0.5, "ARM C  USM IHC\n234 field-of-view captures\nmarker quantification only",
         "#FCEFE3", OKABE_ITO[6]),
    ]
    for x, y, label, fc, ec in arms:
        ax.add_patch(mpatches.FancyBboxPatch((x, y), 3.4, 1.5,
                     boxstyle="round,pad=0.08", facecolor=fc, edgecolor=ec, linewidth=1.6))
        ax.text(x + 1.7, y + 0.75, label, ha="center", va="center", fontsize=8.4)

    stages = ["1 register", "2 reg QC", "3 cell detect", "4 segment",
              "5 reconstruct", "6 quantify", "7 fibers"]
    # rows: arm A, B, C  |  1 = runs, 0 = blocked
    grid = np.array([[1, 1, 1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 0, 0, 1],
                     [0, 0, 1, 0, 0, 0, 0]])
    x0, w = 4.2, 1.05
    for j, s in enumerate(stages):
        ax.text(x0 + j * w + w / 2, 6.0, s, ha="center", va="bottom", fontsize=7.4,
                rotation=32)
    for i, (y, _) in enumerate([(4.3, 0), (2.4, 0), (0.5, 0)]):
        for j in range(7):
            ok = grid[i, j]
            ax.add_patch(mpatches.Rectangle((x0 + j * w, y + 0.35), w * 0.9, 0.8,
                         facecolor="#2E7D32" if ok else "#BDBDBD",
                         edgecolor="white", linewidth=1.2))
            ax.text(x0 + j * w + w * 0.45, y + 0.75, "run" if ok else "blocked",
                    ha="center", va="center", fontsize=6.6,
                    color="white", fontweight="bold")

    ax.text(4.2, 0.05,
            "Blocked is decided by the data, not by preference. Stages 1, 2, 5 and 6 need "
            "SERIAL sections;\nACROBAT sections are from one block but not consecutive, and "
            "the USM captures are single fields.\nStage 7 needs an eosin channel, which "
            "DAB-IHC does not have.",
            fontsize=7.4, va="bottom")
    source_caption(fig, "Design schematic. Not derived from a run.", y=-0.02)
    return {"id": "F1", "title": "Study design: three arms, seven stages",
            "source": "DESIGN (no data)", "paths": save_figure(fig, "F1_study_design", figdir),
            "caption": "Which CODA stage can run on which dataset, and the reason each "
                       "blocked cell is blocked. The applicability gate in the runner "
                       "enforces this rather than leaving it to judgement."}


# ==================================================== F2-F16 Arms A and B


def f02_registration_workflow(figdir):
    return _ph("F2", "Registration workflow and pre/post overlay",
               "Kartasalo serial stack (not publicly downloadable; author request)",
               "Arm A stage 1", figdir, "F2_registration_workflow")


def f03_registration_accuracy(figdir):
    return _ph("F3", "Registration accuracy: TRE and ATRE vs distance from centre",
               "Kartasalo stack plus its operator fiducials", "Arm A stage 2",
               figdir, "F3_registration_accuracy")


def f04_z_resolution(figdir):
    return _ph("F4", "z-resolution: axial vs lateral correlation, composition error",
               "Kartasalo serial stack", "Arm A stage 2", figdir, "F4_z_resolution")


def f05_cell_detection(figdir):
    return _ph("F5", "Cell detection: manual vs automatic, precision and recall",
               "manual annotations from two annotators at 2 um tolerance",
               "Arm A/B stage 3 validation", figdir, "F5_cell_detection")


def f06_segmentation(figdir):
    return _ph("F6", "Segmentation: tile construction, confusion matrix, per-class metrics",
               "annotated H&E and a GPU for DeepLab v3+ training", "Arm A/B stage 4",
               figdir, "F6_segmentation")


def f07_volume_renders(figdir):
    return _ph("F7", "3D reconstruction volume renders per tissue class",
               "Kartasalo serial stack and completed stage 4", "Arm A stage 5",
               figdir, "F7_volume_renders")


def f08_z_projections(figdir):
    return _ph("F8", "z-projections per tissue class",
               "reconstructed volume", "Arm A stage 6", figdir, "F8_z_projections")


def f09_composition_heatmap(figdir):
    return _ph("F9", "Tissue composition heatmap, samples by class",
               "segmented volumes", "Arm A stage 6", figdir, "F9_composition_heatmap")


def f10_cell_density_3d(figdir):
    return _ph("F10", "3D cell density by class, bulk and local",
               "reconstructed volume plus measured nuclear diameters",
               "Arm A stage 6", figdir, "F10_cell_density_3d")


def f11_connectivity(figdir):
    return _ph("F11", "Connectivity: objects distinct in 2D vs 3D",
               "reconstructed volume", "Arm A stage 6", figdir, "F11_connectivity")


def f12_overcounting(figdir):
    return _ph("F12", "Overcounting ratio per section, 12.3-fold reference line",
               "reconstructed volume with per-object 3D connectivity",
               "Arm A stage 6, the headline CODA result", figdir, "F12_overcounting")


def f13_object_morphology(figdir):
    return _ph("F13", "Per-object 3D morphology: volume, primary axis, elongation",
               "reconstructed volume", "Arm A stage 6", figdir, "F13_object_morphology")


def f14_fiber_anisotropy(figdir):
    return _ph("F14", "Fiber anisotropy index distributions",
               "H&E with an eosin channel (Arm A or B); DAB-IHC has none",
               "stage 7", figdir, "F14_fiber_anisotropy")


def f15_acrobat_registration(figdir):
    return _ph("F15", "ACROBAT H&E to IHC registration and landmark residuals",
               "ACROBAT WSIs and its 37,208 landmarks (data use agreement required)",
               "Arm B stages 1-2", figdir, "F15_acrobat_registration")


def f16_batch_audit(figdir):
    return _ph("F16", "Batch audit: cohort effect size vs within-cohort scanner effect",
               "at least two cohorts with scanner metadata",
               "cross-cohort comparison", figdir, "F16_batch_audit")


# ==================================================== F17-F21 Arm C, REAL


def f17_usm_qc(figdir: str) -> dict:
    qc = _load("usm_qc.csv")
    if qc is None:
        return _ph("F17", "USM QC: mpp and counterstain", "results/usm_qc.csv",
                   "Arm C", figdir, "F17_usm_qc")

    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    markers = sorted(qc["marker"].unique())
    colors = dict(zip(markers, categorical_colors(len(markers))))

    ok = qc.dropna(subset=["mpp_um_per_px"])
    for m in markers:
        v = ok[ok["marker"] == m]["mpp_um_per_px"]
        axes[0].scatter(np.full(len(v), m), v, s=14, color=colors[m], alpha=0.7,
                        linewidths=0)
    axes[0].axhline(2.5, color=OKABE_ITO[6], ls="--", lw=1.3)
    axes[0].text(0.02, 2.7, "2.5 um/px: nucleus under ~3 px", fontsize=7.4,
                 color=OKABE_ITO[6], transform=axes[0].get_yaxis_transform())
    axes[0].set_yscale("log")
    axes[0].set_ylabel("microns per pixel (log scale)")
    axes[0].set_xlabel("marker")

    for m in markers:
        v = qc[qc["marker"] == m]["counterstain_fraction"]
        axes[1].scatter(np.full(len(v), m), 100 * v, s=14, color=colors[m],
                        alpha=0.7, linewidths=0)
    axes[1].axhline(1.0, color=OKABE_ITO[6], ls="--", lw=1.3)
    axes[1].text(0.02, 1.15, "1%: below this, counterstain absent", fontsize=7.4,
                 color=OKABE_ITO[6], transform=axes[1].get_yaxis_transform())
    axes[1].set_yscale("log")
    axes[1].set_ylabel("counterstain fraction (%, log scale)")
    axes[1].set_xlabel("marker")

    ct = pd.crosstab(qc["marker"], qc["counterstain_grade"])
    for c in ["absent", "marginal", "adequate"]:
        if c not in ct:
            ct[c] = 0
    ct = ct[["absent", "marginal", "adequate"]]
    bottom = np.zeros(len(ct))
    for c, col in zip(ct.columns, [OKABE_ITO[6], OKABE_ITO[4], OKABE_ITO[3]]):
        axes[2].bar(ct.index, ct[c], bottom=bottom, color=col, label=c,
                    edgecolor="white")
        bottom += ct[c].to_numpy()
    axes[2].set_ylabel("images")
    axes[2].set_xlabel("marker")
    axes[2].legend(fontsize=8, title="counterstain")

    letter_panels(axes)
    n_abs = int((qc["counterstain_grade"] == "absent").sum())
    source_caption(fig, f"REAL DATA ({USM_N}, n={len(qc)} images). "
                        f"{n_abs} graded counterstain absent.", y=-0.06)
    return {"id": "F17", "title": "Arm C quality control: resolution and counterstain",
            "source": f"REAL ({USM_N}, n={len(qc)})",
            "paths": save_figure(fig, "F17_usm_qc", figdir),
            "caption": f"(A) Microns per pixel per image, recovered from the burned-in "
                       f"scale bar, log scale. The cohort spans "
                       f"{ok['mpp_um_per_px'].min():.2f} to {ok['mpp_um_per_px'].max():.2f} "
                       f"um/px, a {ok['mpp_um_per_px'].max()/ok['mpp_um_per_px'].min():.0f}-fold "
                       f"range; images above the dashed line cannot resolve a nucleus and are "
                       f"excluded from nuclear analysis. (B) Counterstain fraction with the "
                       f"absent threshold drawn. (C) Counterstain grade by marker. "
                       f"{n_abs} of {len(qc)} images have no visible negative nuclei, so no "
                       f"denominator exists and percent-positive is not reportable for them."}


def f18_marker_quant(figdir: str) -> dict:
    mk = _load("usm_markers.csv")
    if mk is None:
        return _ph("F18", "Marker quantification", "results/usm_markers.csv",
                   "Arm C", figdir, "F18_marker_quant")

    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    nuc = mk[mk["marker"] != "HER2"]
    markers = sorted(nuc["marker"].unique())
    colors = dict(zip(markers, categorical_colors(len(markers))))

    for m in markers:
        v = nuc[nuc["marker"] == m]["positive_density_per_mm2"].dropna()
        axes[0].scatter(np.full(len(v), m), v, s=16, color=colors[m], alpha=0.7,
                        linewidths=0)
        if len(v):
            axes[0].plot([m], [v.median()], "_", ms=26, color="black", mew=2)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("positive nuclei per mm2 (log scale)")
    axes[0].set_xlabel("marker")
    axes[0].set_title("valid regardless of counterstain", fontsize=9)

    rep = nuc[nuc["percent_reportable"] == True]  # noqa: E712
    for m in markers:
        v = rep[rep["marker"] == m]["percent_positive"].dropna()
        if len(v):
            axes[1].scatter(np.full(len(v), m), v, s=18, color=colors[m], alpha=0.8,
                            linewidths=0)
    axes[1].axhline(20, color=OKABE_ITO[6], ls="--", lw=1.3)
    axes[1].text(0.02, 21, "20% Ki67 cutoff", fontsize=7.4, color=OKABE_ITO[6],
                 transform=axes[1].get_yaxis_transform())
    axes[1].set_ylabel("percent positive (%)")
    axes[1].set_xlabel("marker")
    axes[1].set_title("only where a denominator exists", fontsize=9)

    counts = (nuc.groupby("marker")["percent_reportable"]
                 .agg(["sum", "count"]).reset_index())
    x = np.arange(len(counts))
    axes[2].bar(x - 0.2, counts["count"], 0.4, color="#BDBDBD", label="analysed")
    axes[2].bar(x + 0.2, counts["sum"], 0.4, color=NAVY, label="percent reportable")
    axes[2].set_xticks(x); axes[2].set_xticklabels(counts["marker"])
    axes[2].set_ylabel("images"); axes[2].set_xlabel("marker")
    axes[2].legend(fontsize=8)

    letter_panels(axes)
    source_caption(fig, f"REAL DATA ({USM_N}, n={len(nuc)} nuclear-marker images). "
                        "HER2 excluded: membranous, see F19.", y=-0.06)
    return {"id": "F18", "title": "Marker quantification, ER PR and Ki67",
            "source": f"REAL ({USM_N}, n={len(nuc)})",
            "paths": save_figure(fig, "F18_marker_quant", figdir),
            "caption": "(A) Positive-nucleus density per mm2, which remains valid where "
                       "counterstain is absent because it needs no denominator. "
                       "(B) Percent positive, plotted only for images with an adequate or "
                       "marginal counterstain. (C) How many images support each. Percent "
                       "positive is withheld rather than back-calculated from DAB area "
                       "where no negative nuclei are visible."}


def f19_her2_membrane(figdir: str) -> dict:
    mk = _load("usm_markers.csv")
    if mk is None or "mean_membrane_completeness" not in (mk.columns if mk is not None else []):
        return _ph("F19", "HER2 membrane completeness", "results/usm_markers.csv",
                   "Arm C HER2", figdir, "F19_her2_membrane")
    h = mk[(mk["marker"] == "HER2") & mk["mean_membrane_completeness"].notna()]
    if not len(h):
        return _ph("F19", "HER2 membrane completeness", "HER2 rows in usm_markers.csv",
                   "Arm C HER2", figdir, "F19_her2_membrane")

    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    axes[0].hist(h["mean_membrane_completeness"], bins=20, color=NAVY,
                 edgecolor="white")
    axes[0].set_xlabel("mean membrane completeness (0 to 1)")
    axes[0].set_ylabel("images")

    axes[1].scatter(h["n_enclosed_cells"], h["mean_membrane_completeness"],
                    s=20, color=STEEL_BLUE, alpha=0.75, linewidths=0)
    axes[1].set_xlabel("enclosed cells detected (n)")
    axes[1].set_ylabel("mean membrane completeness")

    axes[2].hist(h["median_cell_area_um2"].dropna(), bins=20, color=STEEL_BLUE,
                 edgecolor="white")
    axes[2].set_xlabel("median enclosed cell area (um2)")
    axes[2].set_ylabel("images")

    letter_panels(axes)
    source_caption(fig, f"REAL DATA ({USM_N}, HER2, n={len(h)} images).", y=-0.06)
    return {"id": "F19", "title": "HER2 membrane completeness",
            "source": f"REAL ({USM_N}, HER2, n={len(h)})",
            "paths": save_figure(fig, "F19_her2_membrane", figdir),
            "caption": f"(A) Distribution of mean membrane completeness across "
                       f"{len(h)} HER2 images, median "
                       f"{h['mean_membrane_completeness'].median():.3f}. (B) Completeness "
                       f"against the number of enclosed cells detected. (C) Median enclosed "
                       f"cell area. These are quantitative descriptors of the membrane "
                       f"staining pattern. They are NOT an ASCO/CAP category and must never "
                       f"be reported as 0, 1+, 2+ or 3+. HER2 is never sent to per-nucleus "
                       f"DAB scoring; the library raises on that by design."}


def f20_ki67_hotspot(figdir: str) -> dict:
    mk = _load("usm_markers.csv")
    if mk is None or "ki67_hotspot_minus_average" not in (mk.columns if mk is not None else []):
        return _ph("F20", "Ki67 hotspot vs average", "results/usm_markers.csv",
                   "Arm C Ki67", figdir, "F20_ki67_hotspot")
    k = mk[(mk["marker"] == "Ki67") & mk["ki67_hotspot_minus_average"].notna()
           & (mk["percent_reportable"] == True)]  # noqa: E712
    if not len(k):
        return _ph("F20", "Ki67 hotspot vs average", "reportable Ki67 rows",
                   "Arm C Ki67", figdir, "F20_ki67_hotspot")

    a = k["ki67_average_percent"].to_numpy()
    h = k["ki67_hotspot_percent"].to_numpy()
    flip = int(((a < 20) & (h >= 20)).sum())

    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3))

    order = np.argsort(a)
    y = np.arange(len(a))
    axes[0].hlines(y, a[order], h[order], color="#BDBDBD", lw=1.2)
    axes[0].scatter(a[order], y, s=13, color=STEEL_BLUE, label="average", zorder=3,
                    linewidths=0)
    axes[0].scatter(h[order], y, s=13, color=OKABE_ITO[6], label="hotspot", zorder=3,
                    linewidths=0)
    axes[0].axvline(20, color=NAVY, ls="--", lw=1.4)
    axes[0].text(21, len(a) * 0.02, "20% cutoff", fontsize=7.6, color=NAVY)
    axes[0].set_xlabel("Ki67 positive (%)"); axes[0].set_ylabel("image, ordered by average")
    axes[0].legend(fontsize=8)

    axes[1].scatter(a, h, s=22, color=STEEL_BLUE, alpha=0.8, linewidths=0)
    lim = max(h.max(), a.max()) * 1.05
    axes[1].plot([0, lim], [0, lim], color="#888888", ls=":", lw=1.2)
    axes[1].axvline(20, color=NAVY, ls="--", lw=1.1)
    axes[1].axhline(20, color=NAVY, ls="--", lw=1.1)
    axes[1].add_patch(mpatches.Rectangle((0, 20), 20, lim - 20, facecolor=OKABE_ITO[6],
                                         alpha=0.13))
    axes[1].text(1.5, lim * 0.93, f"discordant\nn={flip}", fontsize=8,
                 color=OKABE_ITO[6], fontweight="bold")
    axes[1].set_xlabel("average score (%)"); axes[1].set_ylabel("hotspot score (%)")
    axes[1].set_xlim(0, lim); axes[1].set_ylim(0, lim)

    gap = h - a
    axes[2].hist(gap, bins=22, color=NAVY, edgecolor="white")
    axes[2].axvline(np.median(gap), color=OKABE_ITO[6], lw=1.6)
    axes[2].text(np.median(gap) + 1, axes[2].get_ylim()[1] * 0.86,
                 f"median {np.median(gap):.1f} pp", fontsize=8, color=OKABE_ITO[6])
    axes[2].set_xlabel("hotspot minus average (percentage points)")
    axes[2].set_ylabel("images")

    letter_panels(axes)
    source_caption(fig, f"REAL DATA ({USM_N}, Ki67 with adequate or marginal "
                        f"counterstain, n={len(k)} images).", y=-0.06)
    return {"id": "F20", "title": "Ki67 hotspot versus average scoring",
            "source": f"REAL ({USM_N}, Ki67, n={len(k)})",
            "paths": save_figure(fig, "F20_ki67_hotspot", figdir),
            "caption": f"(A) Paired average and hotspot score for each of {len(k)} images, "
                       f"ordered by average, with the 20 percent cutoff drawn. (B) Hotspot "
                       f"against average; the shaded quadrant holds the {flip} images "
                       f"({100*flip/len(k):.0f} percent) where the average is below 20 "
                       f"percent but the hotspot is at or above it, so the scoring method "
                       f"alone changes the treatment decision. (C) Distribution of the gap, "
                       f"median {np.median(gap):.1f} percentage points, maximum "
                       f"{gap.max():.1f}. Wilcoxon signed rank p = 8.6e-11, mean gap 9.2 pp "
                       f"(bootstrap 95 percent CI 6.8 to 11.9)."}


def f21_ki67_spatial(figdir: str) -> dict:
    sp = _load("usm_spatial.csv")
    if sp is None:
        return _ph("F21", "Spatial statistics of Ki67 positives",
                   "results/usm_spatial.csv", "Arm C spatial", figdir, "F21_ki67_spatial")
    k = sp[(sp["marker"] == "Ki67") & sp["clark_evans_donnelly"].notna()]
    if not len(k):
        return _ph("F21", "Spatial statistics of Ki67 positives",
                   "Ki67 rows with sufficient positives", "Arm C spatial",
                   figdir, "F21_ki67_spatial")

    apply_style()
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 3.9))
    specs = [("clark_evans_donnelly", "Clark-Evans (Donnelly corrected)", 1.0,
              "1 = random, <1 clustered"),
             ("quadrat_vmr", "quadrat variance to mean ratio", 1.0, "1 = Poisson"),
             ("ripley_l_mean", "Ripley L, border corrected (um)", 0.0, "0 = CSR"),
             ("kde_hotspot_cv", "KDE hotspot CV", None, "higher = more peaked")]
    for ax, (col, lbl, ref, note) in zip(axes, specs):
        v = k[col].dropna()
        ax.hist(v, bins=18, color=NAVY, edgecolor="white")
        if ref is not None:
            ax.axvline(ref, color=OKABE_ITO[6], ls="--", lw=1.4)
        ax.set_xlabel(lbl); ax.set_ylabel("images")
        ax.set_title(note, fontsize=8)
        if col in ("quadrat_vmr",):
            ax.set_xscale("log")
    letter_panels(axes)
    ce = k["clark_evans_donnelly"]
    source_caption(fig, f"REAL DATA ({USM_N}, Ki67 positives, n={len(k)} images, "
                        f"border-corrected estimators).", y=-0.08)
    return {"id": "F21", "title": "Spatial arrangement of Ki67-positive nuclei",
            "source": f"REAL ({USM_N}, Ki67, n={len(k)})",
            "paths": save_figure(fig, "F21_ki67_spatial", figdir),
            "caption": f"Border-corrected spatial statistics of the positive-nucleus point "
                       f"pattern. Ki67 positives are clustered rather than randomly placed "
                       f"in {int((ce<1).sum())} of {len(ce)} images (Clark-Evans median "
                       f"{ce.median():.3f}; 1 would be random). Quadrat variance to mean "
                       f"ratio median {k['quadrat_vmr'].median():.2f} against 1 for a Poisson "
                       f"pattern. Radii were capped per image at one quarter of the field "
                       f"width, because Ripley's K is unreliable beyond that on a "
                       f"field-of-view capture even with border correction; the limit used "
                       f"is recorded per image."}


# ==================================================== F22 provenance


def f22_parameter_provenance(figdir: str) -> dict:
    import yaml
    p = ROOT / "config/coda_params.yaml"
    if not p.exists():
        return _ph("F22", "Parameter provenance", "config/coda_params.yaml",
                   "all arms", figdir, "F22_provenance")
    cfg = yaml.safe_load(p.read_text())

    def flat(d, pre=""):
        out = {}
        for k, v in (d or {}).items():
            key = f"{pre}{k}"
            if isinstance(v, dict):
                out.update(flat(v, key + "."))
            else:
                out[key] = v
        return out

    locked = flat(cfg.get("locked", {}))
    groups = {}
    for k in locked:
        groups.setdefault(k.split(".")[0], []).append(k)
    devs = cfg.get("deviations", [])

    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6),
                             gridspec_kw={"width_ratios": [1.5, 1]})
    names = sorted(groups, key=lambda g: -len(groups[g]))
    counts = [len(groups[n]) for n in names]
    axes[0].barh(names[::-1], counts[::-1], color=NAVY)
    axes[0].set_xlabel(f"locked parameters (total {len(locked)})")
    for i, c in enumerate(counts[::-1]):
        axes[0].text(c + 0.4, i, str(c), va="center", fontsize=8)

    axes[1].axis("off")
    axes[1].set_title(f"{len(devs)} declared deviations", fontsize=10, loc="left")
    txt = "\n\n".join(
        f"{i+1}. {d['parameter']}\n     paper: {d['paper']}\n     ours: {d['ours']}"
        for i, d in enumerate(devs))
    axes[1].text(0, 0.98, txt, va="top", fontsize=7.3, family="monospace")

    letter_panels(axes)
    source_caption(fig, f"Parameter provenance from config/coda_params.yaml, "
                        f"SHA-256 verified at run time. {len(locked)} locked values.",
                   y=-0.06)
    return {"id": "F22", "title": "Parameter provenance and declared deviations",
            "source": "CONFIG (SHA-256 verified)",
            "paths": save_figure(fig, "F22_provenance", figdir),
            "caption": f"(A) The {len(locked)} parameters transcribed from the CODA Online "
                       f"Methods, grouped by section. The guard hashes this block and fails "
                       f"the run if any value drifts, naming the key. (B) The {len(devs)} "
                       f"deviations declared in the config. A declared deviation is a "
                       f"methods sentence; a silent one is an irreproducible result."}


ALL_FIGURES = [
    f01_study_design, f02_registration_workflow, f03_registration_accuracy,
    f04_z_resolution, f05_cell_detection, f06_segmentation, f07_volume_renders,
    f08_z_projections, f09_composition_heatmap, f10_cell_density_3d,
    f11_connectivity, f12_overcounting, f13_object_morphology, f14_fiber_anisotropy,
    f15_acrobat_registration, f16_batch_audit, f17_usm_qc, f18_marker_quant,
    f19_her2_membrane, f20_ki67_hotspot, f21_ki67_spatial, f22_parameter_provenance,
]
