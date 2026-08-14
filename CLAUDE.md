# CLAUDE.md — coda-brca-my

Project instructions for Claude Code. Read fully before writing code.

## What this is

Two arms, one project.

**2D arm (buildable now).** CODA-derived histomorphometry applied to a
Malaysian breast cohort (USM, Kota Bharu) and compared against TCGA-BRCA, plus
quantified spatial analysis of ER/PR/HER2/Ki67 IHC within the Malaysian cohort.

**3D arm (blocked on wet lab).** Full CODA serial-section reconstruction. See
`docs/3D_PROTOCOL_AND_COSTING.md`. It cannot start until blocks are cut. Do not
write 3D reconstruction code until sections exist.

Reference: Kiemen A et al. CODA. *Nat Methods* 2022;19:1490-1499.
Code: github.com/ashleylk/CODA (MATLAB).

## Hard facts that constrain the design

**CODA does not use IHC.** It deliberately avoids it, labelling ten tissue types
from serial H&E alone and describing IHC as an expensive limitation of competing
methods. Do not describe this project as "CODA on IHC".

**TCGA cannot support 3D.** TCGA slides are single sections from different
blocks, never consecutive. No serial data exists in TCGA for any cancer type.
The 3D arm is Malaysian-cohort-only and has no comparison cohort.

**Existing USM scans cannot support 3D.** 200+ scans across four markers and
many patients is breadth, not depth. CODA needs 100+ consecutive sections
through ONE block to build ONE volume.

**Comparisons must be H&E to H&E.** Never compare USM IHC against TCGA H&E. The
difference you measure would be the stain. Every USM case with IHC has a
diagnostic H&E; use that for the cohort comparison, and use the IHC as a
molecular annotation layer within the Malaysian cohort only.

## The threat that decides whether this study is publishable

**Batch effect will masquerade as population difference.** Different scanner,
stain protocol, lab, fixation time. Models can predict TCGA tissue source site
from H&E alone; the site signal is that strong. If you scan in Kota Bharu and
compare to TCGA you WILL find differences and almost none will be biological.

Three defences, all implemented in `cohort_compare.py`, all mandatory:

1. **A third cohort.** CPTAC-BRCA or BCNB. If Malaysia differs from TCGA AND
   from the third cohort in the same direction, that is a signal. If Malaysia is
   simply the odd one out on every feature, that is the scanner.
2. **Clinical matching.** Malaysian breast cancer presents younger, later-stage,
   with more TNBC. Without propensity matching on age, stage, grade and subtype
   you are measuring stage, not population.
3. **Batch sensitivity audit.** `audit_batch_sensitivity()` flags features that
   separate scanners within a cohort more strongly than they separate cohorts.
   Report the audit table. Do not silently drop failures.

If more than 30% of features fail the audit, stain normalisation is not
sufficient. Consider rescanning a subset of TCGA slides on the USM scanner.

## Ethnicity

Malaysia is Malay, Chinese and Indian. **"Malaysian" is not an ancestry group.**
Record ethnicity per case and stratify. Pooling will be flagged by reviewers who
work in this area and can both hide and fabricate effects.

## Methodological constraints from CODA

- Per-image stain vector estimation by k-means over optical densities, not fixed
  textbook vectors. Across cohorts this removes part of the lab-to-lab staining
  difference before any comparison is made.
- Fiber anisotropy index in ~2500 µm² windows (50 × 50 µm at 0.5 µm/px). 1 means
  fully aligned, 0 isotropic.
- Cell detection from the hematoxylin channel, validated at >90% precision and
  recall against two manual annotators.
- Segmentation acceptance bar: >90% precision and recall per class.

## The sectioning-angle caveat, which you must not paper over

Fiber alignment on a 2D section depends on the angle the structure was cut at.
CODA could correct for this because it had the volume and chose the angle
deliberately, finding 2.2 to 2.5-fold differences between longitudinal and axial
sections of the same structures. **On single sections you cannot correct for it.**

Two acceptable responses:
- Restrict measurement to the invasive front, where orientation is defined
  relative to the tumour boundary rather than to a tube (`boundary_relative_orientation`).
- Report the distribution across many windows (`tiled_anisotropy`) and treat
  sectioning angle as noise.

Never report a single per-slide alignment number as a property of the patient.

## HER2

HER2 is membranous. Per-nucleus DAB scoring is the wrong operation and produces
a confident meaningless number, which is worse than an error. `score_marker()`
raises on HER2 by design. Do not remove that guard; implement membrane
completeness separately if HER2 is needed.

## Ki67, the scientific opening

Ki67 scoring is irreproducible because observers disagree on hotspot versus
average assessment, and the 20% cutoff driving chemotherapy decisions sits where
reproducibility is worst. Nobody routinely quantifies HOW positive nuclei are
ARRANGED. A clustered 18% and a dispersed 18% get the same score and the same
treatment decision.

`hotspot_vs_average()` quantifies the discrepancy directly. `to_point_pattern()`
converts scored nuclei into the two-class point pattern consumed by the
canvas-brca spatial statistics (border-corrected Ripley K and L,
Donnelly-corrected Clark-Evans, quadrat dispersion, KDE hotspot CV).

## Code standards

- Python 3.10, type hints, numpy-style docstrings.
- Config-driven, no magic numbers in function bodies.
- Log with `logging`, not `print`.
- Effect sizes with confidence intervals, never bare p values. With hundreds of
  slides trivial differences reach significance; Cliff's delta tells you whether
  it matters.
- Benjamini-Hochberg FDR, report q alongside p.
- Tests use synthetic fixtures. No downloads, no real patient slides in tests.
- Baseline is 21 passing tests. If the count drops, revert.

## Do not modify

`fibers.py`, `deconv.py`, `cohort_compare.py` are tested and correct. If you
think they need improving, say so rather than editing.

## Order of work

1. Verify baseline: `python -m pytest tests/ -q`
2. Batch audit on a small pilot: 10 USM slides, 10 TCGA slides. **Do this before
   processing the full cohort.** If the audit fails badly there is no point
   processing 200 slides.
3. Full feature extraction on USM H&E.
4. TCGA-BRCA download and matched feature extraction.
5. Third cohort (CPTAC or BCNB).
6. Matching, comparison, stratified by ethnicity and subtype.
7. IHC quantification and Ki67 spatial heterogeneity, USM only.
8. Link: does H&E morphology predict measured Ki67 spatial pattern?

## Ask before assuming

- Whether a USM scan is 20× or 40×. Read the slide properties; never assume.
- Which USM cases have H&E, IHC, and complete clinical data. The intersection
  is the real n, and it is usually smaller than any single count.
- Whether ethics approval covers cross-cohort comparison with public data.
