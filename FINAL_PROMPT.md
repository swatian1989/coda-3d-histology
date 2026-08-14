# FINAL MASTER PROMPT — paste into Claude Code, once

Everything below the line goes into Claude Code. Before pasting:

```bash
unzip coda-brca-my.zip && cd coda-brca-my
python -m venv .venv && source .venv/bin/activate
pip install numpy pandas scipy scikit-learn scikit-image statsmodels pyyaml \
            matplotlib seaborn pytest pillow openslide-python pyarrow tqdm
python -m pytest tests/ -q          # must print 60 passed
git init && git add -A && git commit -m "baseline: 60 tests, CODA params locked"
claude
```

---

You are reproducing CODA (Kiemen A et al., *Nat Methods* 2022;19:1490-1499,
doi:10.1038/s41592-022-01650-9) on public data and on my own breast IHC images.

Read these three files in full before writing any code:
`config/coda_params.yaml`, `docs/FINAL_METHODS.md`, `docs/CODA_PUBLIC_PROTOCOL.md`.

## The design

Three arms. All 7 CODA stages get run; which arm runs which stage is determined
by the data, not by preference.

| Arm | Data | Stages | Gives |
|---|---|---|---|
| A | Kartasalo, 260 mouse prostate serial sections | **1-7, full 2D and full 3D** | Registration accuracy vs published benchmark; 3D reconstruction; 2D-vs-3D overcounting |
| B | ACROBAT, 1,153 human breast patients | 3, 4, 7 + H&E↔IHC registration | 2D on breast with my exact marker panel, validated on 37,208 landmarks |
| C | My USM breast IHC field-of-view captures | marker quantification only | HER2 membrane, Ki67 spatial heterogeneity — results no public data can produce |

3D exists only in Arm A. TCGA, ACROBAT and my images are all single sections or
non-consecutive, so stages 1, 2, 5, 6 cannot run on them. The code enforces this.

## Rules, all phases

1. **Never edit the `locked:` block of coda_params.yaml.** It holds 120
   parameters transcribed from the Online Methods and is SHA-256 hashed.
   `guard.verify()` fails the run and names the drifted key. A deliberate change
   goes in `deviations:` with parameter, paper, ours and reason, then delete
   `config/.coda_params.sha256` to re-register.

2. **Call `guard.verify()` and `guard.check_applicability()` at the start of
   every run.** Never force a stage the gate blocks. Registering non-serial
   sections produces a transform and a correlation number and neither means
   anything.

3. **Do not modify these tested modules:** `registration.py`, `reconstruct.py`,
   `qc.py`, `segmentation.py`, `scalebar.py`, `her2.py`, `deconv.py`,
   `fibers.py`, `ihc.py`, `cohort_compare.py`, `guard.py`. If you think one
   needs changing, say so and wait.

4. **Run `python -m pytest tests/ -q` after every phase.** Baseline is 60. If it
   drops, revert before continuing.

5. **CPU only, 16 GB RAM, no GPU** except where noted. Stream tiles, never load
   a whole WSI, cache anything expensive to disk.

6. **Ask about data units, file layout and column names. Do not guess.** Wrong
   units silently destroy every downstream result while still producing
   clean-looking output.

7. **Label every figure and table REAL (name the dataset and n) or SIMULATED.**
   Never present a synthetic fixture result as a finding.

---

## PHASE 0 — baseline

Write no code. Run:

```
python -m pytest tests/ -q
python -c "import sys; sys.path.insert(0,'src'); from coda_my.guard import verify, snapshot, flatten; p=verify(); snapshot(p); print(len(flatten(p['locked'])), 'locked parameters')"
python scripts/run_coda.py --dataset tcga --stages all
```

Expect 60 tests passed, 120 locked parameters, and the TCGA run dropping stages
1, 2, 5, 6 with an explanation. Report the actual output. Then check `data/raw/`
and tell me exactly what is there.

## PHASE 1 — my IHC images (Arm C)

My images are in `data/raw/usm/`. They are microscope field-of-view captures
with a burned-in red scale bar, at mixed magnifications.

For every image, in this order:

1. `detect_scale_bar()` — pass the label microns. Verified values on my four
   sample images were 0.222, 0.424, 0.690 and 0.708 µm/px. If the label is not
   in the filename, ask me for it. Do not assume a default; a wrong mpp
   mis-scales every measurement by a constant factor that survives to the figure.
2. `mask_overlay_region()` — exclude the bar and its text before any
   measurement. It is a high-contrast object that nucleus detectors segment and
   spatial statistics read as a dense corner cluster.
3. `has_counterstain()` — grade adequate / marginal / absent.

Write `results/usm_qc.csv` with filename, mpp, counterstain grade and fraction.
**Flag every image graded "absent".** Those have no visible negative nuclei, so
there is no denominator: percent-positive is not reportable and must not be
back-calculated from DAB area. Positive density and spatial pattern remain valid.

Then:
- ER, PR, Ki67 → `score_marker()`, `hotspot_vs_average()`, `to_point_pattern()`
  into the border-corrected spatial statistics (Ripley K and L with border
  correction, Donnelly-corrected Clark-Evans, quadrat dispersion, KDE hotspot CV)
- HER2 → `membrane_completeness()` only. Never nuclear scoring; the code raises
  on it by design.

Report the Ki67 hotspot-minus-average gap per image. That gap is the
reproducibility problem this arm exists to quantify.

## PHASE 2 — full 2D and full 3D on Kartasalo (Arm A)

Download the 260-section mouse prostate serial stack and the 47-section liver
stack from https://github.com/BioimageInformaticsTampere/RegBenchmark

Run all seven stages:

1. `register_stack()` — Radon global registration at 80 µm/px, 3 candidate
   references, elastic field from 1.5 mm tiles at 8 µm/px smoothed σ=2 px,
   everything referenced to the **centre** section.
2. `target_registration_error()` and `accumulated_tre()` against the supplied
   fiducials; `axial_vs_lateral_correlation()` and `z_skip_validation()`.
   Compare to the published benchmark. ATRE is the metric that matters — pairwise
   error can look excellent while the stack banana-bends. The liver stack has
   laser-cut holes that must form a straight line after correct alignment; use
   them as an independent check.
3. Cell detection, validated at the 2 µm tolerance.
4. Segmentation — **GPU needed**. Use the paper's data design exactly:
   7 annotation images, ≥50 per class, 9000×9000 tiles >65% filled by overlaying
   the least-represented class, cut to 324 tiles of 500×500, augment rotation +
   scale 0.8-1.2 + hue 0.8-1.2, patience 5, and the >90% precision-AND-recall
   gate that forces retraining rather than reporting whatever comes out.
   If no GPU is available, stop and tell me; do not substitute a smaller model
   without saying so.
5. `stack_to_volume()` — isotropic 12 µm voxels.
6. `tissue_composition()`, `extrapolate_3d_cell_count()`, `overcounting_ratio()`,
   `object_metrics()`.
7. `tiled_anisotropy()` on the eosin channel.

**Measure the nuclear diameters yourself** before running the cell-count
correction. The defaults in `reconstruct.py` are CODA's pancreas values; the
correction scales counts directly, so borrowed diameters bias every density
systematically between cell types rather than randomly.

On the overcounting result: the paper reports mean 12.3-fold and max 40-fold in
human pancreatic precursors. **Do not expect that number.** Mouse prostate has
different architecture. What should reproduce is the direction — 2D overcounts
3D wherever structures branch. Report your number and say plainly that it is not
comparable to the paper's, and why.

## PHASE 3 — 2D on human breast (Arm B)

ACROBAT from https://researchdata.se/en/catalogue/dataset/2022-190-1
4,212 WSIs, 1,153 patients, H&E + ER/PGR/HER2/KI67, 37,208 landmarks.

Register H&E to IHC with `global_register()` and `elastic_field()`, validate
against the landmarks with `target_registration_error()`. Then stages 3, 4, 7.

Sections come from the same block but are **not necessarily consecutive**. Do
not stack them into a volume.

Before comparing anything across cohorts, run `audit_batch_sensitivity()`. Three
scanners are involved in ACROBAT alone. Report the audit table; do not silently
drop the features that fail it.

## PHASE 4 — report

One document, `reports/analysis_report.md` plus a self-contained HTML with
figures embedded. Structure: summary, methods, results by arm, limitations,
reproducibility.

Must include:
- The full deviations table from `coda_params.yaml`
- Which stages ran on which arm, and which were blocked and why
- Effect sizes with confidence intervals, never bare p values
- Wilcoxon rank sum as the paper specifies, Benjamini-Hochberg q alongside p
- The counterstain QC table with every "absent" image flagged
- A blunt limitations section: field-of-view rather than whole slide for Arm C,
  possible non-random field selection, 3-fold magnification range, no denominator
  where counterstain is absent, mouse rather than human for Arm A, and no 3D on
  breast anywhere

## Start

Phase 0 now. Then Phase 1, since it needs no downloads. Report and continue.
