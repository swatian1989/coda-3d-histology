"""Report prose. Every paragraph of the CODA report lives here.

Kept separate from the renderers so that editing what the report SAYS never
risks breaking how it is rendered. report.py imports build_sections from here.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Section:
    heading: str
    level: int
    paragraphs: list[str] = field(default_factory=list)
    figure_ids: list[str] = field(default_factory=list)
    table_ids: list[str] = field(default_factory=list)


def _n(tab, tid, default=0):
    try:
        return len(tab[tid]["df"])
    except Exception:
        return default


def build_sections(fig: dict, tab: dict, stats: dict) -> list[Section]:
    n_missing_fig = sum(1 for f in fig.values() if f["source"] == "MISSING DATA")
    n_missing_tab = sum(1 for t in tab.values() if t["source"] == "MISSING DATA")
    n_real_fig = sum(1 for f in fig.values() if f["source"].startswith("REAL"))

    return [
        Section("Summary", 1, [
            "This report reproduces components of CODA (Kiemen et al., Nature Methods "
            "2022) on a Malaysian breast immunohistochemistry series, and states plainly "
            "which parts of the method could not be run and why.",
            "The design has three arms. Arm A is the only place the full seven-stage "
            "pipeline including three-dimensional reconstruction can run, because it is "
            "the only dataset with true serial sections. Arm B supplies human breast "
            "tissue with the same marker panel and 37,208 registration landmarks. Arm C "
            "is the institutional image series, which supports marker quantification and "
            "spatial analysis and nothing else.",
            f"**Only Arm C has data.** The Arm A images are not publicly downloadable: "
            f"the benchmark repository ships the evaluation software, not the sections, "
            f"and the dataset requires a request to the authors. Arm B requires a data "
            f"use agreement. Consequently {n_real_fig} of 22 figures and "
            f"{22 - n_missing_fig - n_real_fig} others are built from real or "
            f"configuration data, while {n_missing_fig} figures and {n_missing_tab} "
            f"tables are labelled placeholders naming the exact input each needs. "
            f"Nothing is fabricated to fill a gap.",
            "**The headline result is from Arm C.** Across 62 Ki67 images with a valid "
            "denominator, scoring the same field by hotspot rather than by average "
            "raises the index by a median of 5.8 percentage points (mean 9.2, bootstrap "
            "95 percent confidence interval 6.8 to 11.9; Wilcoxon signed rank "
            "p = 8.6e-11, maximum gap 52.7 points). On 14 images, 23 percent of the "
            "series, the average lies below the 20 percent cutoff while the hotspot lies "
            "at or above it, so the choice of scoring method alone changes the treatment "
            "decision.",
            "**Spatial arrangement explains part of that discordance, and the scale "
            "matters.** Ki67-positive nuclei are spatially clustered rather than randomly "
            "placed in 51 of 53 images. Coarse-scale clustering, measured as the quadrat "
            "variance to mean ratio over windows comparable to the reporting field, "
            "correlates strongly with the hotspot-minus-average gap (Spearman rho 0.66, "
            "95 percent CI 0.47 to 0.79, q = 4.1e-07). Nearest-neighbour clustering "
            "measured by the Clark-Evans index does not (rho -0.03, q = 0.84). A "
            "statistic computed at the wrong spatial scale carries no information about "
            "the reproducibility problem.",
            "**Three caveats govern everything below.** These are field-of-view captures "
            "rather than whole slides, so no whole-slide inference follows and spatial "
            "statistics are bounded by the field. Counterstain is absent on 150 of 234 "
            "images, which removes the denominator and makes percent-positive "
            "unreportable for those; density and spatial pattern remain valid. No part "
            "of this work performs three-dimensional reconstruction on breast tissue, "
            "and none of it can, because no serial breast sections exist here.",
        ], figure_ids=["F1"], table_ids=["T1", "T14"]),

        Section("Methods", 1, [
            "Parameters are transcribed from the CODA Online Methods into a "
            "configuration file whose locked block is hashed with SHA-256 and verified "
            "at the start of every run; a drifted value fails the run and is named. "
            f"T2 lists all {_n(tab,'T2',120)} locked parameters grouped by Online "
            f"Methods section. T3 lists every declared deviation with its reason and "
            f"expected impact, and each is referenced inline where it applies.",
            "The runner also enforces applicability. Stages 1, 2, 5 and 6 require serial "
            "sections and are refused on any dataset that lacks them, rather than "
            "producing a transform and a correlation coefficient that would look like "
            "results. Registering sections that are not consecutive yields numbers "
            "without meaning, and the gate exists to prevent exactly that.",
        ], figure_ids=["F22"], table_ids=["T2", "T3"]),

        Section("Scale calibration and quality control (Arm C)", 2, [
            "The images are microscope field-of-view captures with a burned-in red scale "
            "bar, and the bar is both the calibration and a contaminant. It gives "
            "microns per pixel, without which every distance is in pixels; it is also a "
            "saturated high-contrast object that a nucleus detector segments and a "
            "spatial statistic reads as a dense corner cluster. The overlay region is "
            "masked before any measurement.",
            "The micron value printed beside the bar is not in the filenames and cannot "
            "be read reliably by the detector, so it was recovered by isolating the "
            "saturated red text and reading it, for 231 of 234 images. Three could not "
            "be read and carry a null value with a stated reason rather than a default. "
            "The recovered calibration reproduces all four independently verified "
            "reference values exactly: 0.222, 0.424, 0.690 and 0.708 microns per pixel.",
            "**One correction was necessary and is worth recording.** The bar detector "
            "takes the longest run of red pixels, and on heavily stained images a streak "
            "of brown diaminobenzidine outran the bar itself, placing the detected bar "
            "in the middle of the tissue. Restricting the search to the bottom twelve "
            "percent of the frame corrected 16 of 234 images. Left uncorrected, each of "
            "those would have been rescaled by a constant factor that survives every "
            "subsequent step and reaches the figures intact.",
            "**The cohort spans a 38-fold range of magnification**, 0.197 to 7.50 microns "
            "per pixel, not the 3-fold range that four sample images had suggested. At "
            "the coarse end a nucleus occupies about one pixel. Images were therefore "
            "tiered, and the 7 coarser than 2.5 microns per pixel are excluded from "
            "nuclear analysis with the reason recorded. Scale-dependent texture "
            "measurements are not comparable across this range at all.",
            "Counterstain was graded by deconvolution rather than an RGB heuristic. It "
            "is absent on 150 of 234 images, and the distribution is uneven in a way "
            "that decides what each marker can support: Ki67 retains a usable "
            "counterstain on 70 of 76 images, whereas ER, PR and HER2 do not on the "
            "large majority. Where counterstain is absent there are no visible negative "
            "nuclei, so no denominator exists. Percent positive is withheld for those "
            "images and is never back-calculated from stained area.",
        ], figure_ids=["F17"], table_ids=["T11"]),

        Section("Marker quantification (Arm C)", 2, [
            "225 of 234 images were analysed; 9 were skipped, 3 for an unreadable scale "
            "bar and 6 for insufficient resolution, each with the reason recorded.",
            "ER, PR and Ki67 were scored per nucleus for diaminobenzidine positivity. "
            "Positive-cell density per square millimetre is reported for every image "
            "because it requires no denominator. Percent positive is reported only where "
            "the counterstain gate permits, which is 9 of 62 ER images and 2 of 39 PR "
            "images.",
            "HER2 was never sent to per-nucleus scoring. It is a membranous marker, and "
            "per-nucleus diaminobenzidine scoring of it produces a confident meaningless "
            "number, which is worse than an error; the library raises on the attempt by "
            "design. Membrane completeness was measured instead across 53 images, median "
            "0.998. These are quantitative descriptors of the staining pattern and are "
            "not an ASCO/CAP category. They must never be reported as 0, 1+, 2+ or 3+.",
        ], figure_ids=["F18", "F19"], table_ids=["T12"]),

        Section("Ki67 hotspot versus average, and the spatial arrangement (Arm C)", 2, [
            "Ki67 scoring is irreproducible in practice because observers disagree on "
            "whether to score a hotspot or an average, and the 20 percent cutoff that "
            "drives chemotherapy decisions sits where that disagreement is worst. Both "
            "scores were computed from the same nuclei on the same image, so the "
            "difference is attributable to the scoring convention alone and to nothing "
            "else.",
            "The hotspot score exceeds the average by a median of 5.8 percentage points "
            "(interquartile range 0.9 to 15.2, maximum 52.7). The mean difference is 9.2 "
            "points with a bootstrap 95 percent confidence interval of 6.8 to 11.9, and "
            "the paired Wilcoxon signed rank test gives p = 8.6e-11. Of 62 images, 31 "
            "fall below the cutoff on both conventions, 17 fall at or above it on both, "
            "and **14 (23 percent) are discordant**, with the average below 20 percent "
            "and the hotspot at or above it.",
            "Positive nuclei are not randomly arranged. Clark-Evans, corrected for edge "
            "effects by Donnelly's perimeter term, has a median of 0.686 and is below 1, "
            "indicating clustering, in 51 of 53 images. The quadrat variance to mean "
            "ratio has a median of 6.63 against 1 for a Poisson pattern. Border "
            "correction matters more than usual here: on a field-of-view capture a large "
            "fraction of the field lies within one analysis radius of an edge, and an "
            "uncorrected estimator reads the missing area as reduced clustering. Radii "
            "were capped per image at one quarter of the field width and the limit is "
            "recorded alongside every value.",
            "**The scale at which clustering is measured determines whether it carries "
            "information.** Quadrat variance to mean ratio, computed over windows "
            "comparable to the reporting field, correlates strongly with the "
            "hotspot-minus-average gap (Spearman rho 0.66, 95 percent CI 0.47 to 0.79, "
            "Benjamini-Hochberg q = 4.1e-07). The kernel density hotspot coefficient of "
            "variation and the border-corrected Ripley L correlate weakly (rho 0.31 and "
            "0.30, q = 0.044 for both). Clark-Evans, which measures nearest-neighbour "
            "spacing at single-cell distances, does not correlate at all (rho -0.03, "
            "q = 0.84). The discordance is produced by large-scale patchiness, not by "
            "whether positive nuclei touch one another, and a statistic computed at the "
            "wrong scale is silent about it.",
            "This is the question a percentage cannot answer. A clustered 18 percent and "
            "a dispersed 18 percent receive the same score and the same treatment "
            "decision, and these statistics separate them.",
        ], figure_ids=["F20", "F21"], table_ids=["T13"]),

        Section("Arms A and B: not run", 2, [
            "Arm A would provide the registration accuracy benchmark and the "
            "two-dimensional versus three-dimensional overcounting result, which is the "
            "headline finding of the original paper. The images are not publicly "
            "downloadable. The benchmark repository contains the evaluation framework "
            "only, and obtaining the sections requires a request to the authors. F2, F3, "
            "F4, F7 through F14 and T4, T5, T8, T9, T10 are therefore placeholders.",
            "Arm B would validate haematoxylin and eosin to immunohistochemistry "
            "registration against 37,208 landmarks on breast tissue with this exact "
            "marker panel, and is the only place fibre alignment could run, because "
            "diaminobenzidine immunohistochemistry has no eosin channel. It requires a "
            "data use agreement. F15, F16, T6 and T7 are placeholders.",
            "Neither is a methodological obstacle. Both are access steps that a human "
            "must complete, and the pipeline that consumes them is written and tested.",
        ], figure_ids=["F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11",
                       "F12", "F13", "F14", "F15", "F16"],
           table_ids=["T4", "T5", "T6", "T7", "T8", "T9", "T10"]),

        Section("Limitations", 1, [
            "**Field of view, not whole slide.** Arm C images are microscope captures "
            "roughly 350 to 1150 microns across. No whole-slide inference follows from "
            "them, and Ripley's K beyond about a quarter of the field width is "
            "unreliable even with border correction.",
            "**Possible non-random field selection.** Photographs taken to document "
            "positive staining are biased toward positive areas. If the fields were "
            "chosen by eye then the sample does not represent the slide. How the fields "
            "were selected is not recorded in the image metadata and should be stated "
            "explicitly by whoever captured them.",
            "**A 38-fold magnification range.** Every measurement is converted to microns "
            "before pooling, but scale-dependent texture features are not comparable "
            "across this range, and 7 images are too coarse for nuclear analysis at all.",
            "**No denominator on 150 of 234 images.** Percent positive is not reportable "
            "for those and has not been estimated by any indirect route.",
            "**No three-dimensional analysis anywhere in this work, and none on breast.** "
            "Arm A is the only source of serial sections and it is mouse prostate and "
            "liver. Even complete, it would validate the pipeline rather than establish "
            "a breast finding.",
            "**Sectioning angle is uncorrected.** Fibre alignment on a single section "
            "depends on the angle the structure was cut at, and the original work could "
            "correct for this only because it had the volume. Stage 7 did not run here, "
            "so the issue does not affect the present results, but it constrains any "
            "future single-section fibre measurement.",
            "**One institution, one scanner, no comparison cohort.** No cross-cohort "
            "comparison is attempted, so no batch audit is reported. Any future "
            "comparison against public cohorts must run the batch sensitivity audit "
            "first, because scanner and protocol differences will otherwise masquerade "
            "as population differences.",
        ]),

        Section("Reproducibility", 1, [
            f"Configuration `{stats['config_path']}`, SHA-256 of the merged locked block "
            f"`{stats['config_hash_sha256']}`, seed {stats['project_seed']}. "
            f"Git SHA **{stats['git_sha']}**. Python {stats['python_version']} on "
            f"{stats['platform']}.",
            "**Software environment.** Every third-party library used, with the version "
            "installed and the role it plays, so the environment can be rebuilt from "
            "this report without reading the source.",
            "SOFTWARE_TABLE_PLACEHOLDER",
            "Runtime on the machine used, CPU only: quality control over 234 images "
            "about 4 minutes; marker quantification about 25 minutes; spatial statistics "
            "about 8 minutes; report generation about 3 minutes. Stage 4 segmentation "
            "would require a GPU and did not run.",
            "Exact commands to regenerate every artefact in this report:",
            "```\n"
            "python -m pytest tests/ -q\n"
            "python scripts/run_usm_qc.py\n"
            "python scripts/run_usm_markers.py\n"
            "python scripts/run_usm_spatial.py\n"
            "python scripts/run_report.py\n"
            "```",
        ]),
    ]
