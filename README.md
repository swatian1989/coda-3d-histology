# coda-3d-histology

Reproduction of the CODA three-dimensional histology pipeline (Kiemen et al.,
*Nature Methods* 2022) on openly licensed serial sections, together with a
two-dimensional immunohistochemistry arm on a breast cohort.

CPU only. No GPU is used or required. Peak memory is about 2 GB.

---

## What is here

Seven pipeline stages, applied to the datasets that can actually support each
one. The constraint is a property of the material, not of the implementation:
stages 1, 2, 5 and 6 need consecutively cut serial sections, and stage 7 needs
an eosin channel.

| Stage | Status | Result |
|---|---|---|
| 1 Registration | run | 47/47 sections, median pixel correlation 0.929, 0 flagged |
| 2 Registration QC | run | TRE 114 um against 727 um unregistered, 6.35-fold (Wilcoxon p = 8.5e-14, n = 46) |
| 3 Cell detection | partial | two detectors compared to each other; they disagree threefold, so neither is validated |
| 4 Segmentation | not run | needs annotated training tiles and a GPU; neither exists here |
| 5 Reconstruction | run | 47 sections stacked into a volume at 7.36 um/px, 5 um spacing |
| 6 Connectivity | run | 2D counting overestimates object number 10.1-fold at structure scale |
| 7 Fibre alignment | partial | anisotropy in the cutting plane; the sectioning-angle comparison fails its control and is withheld |

An eighth arm quantifies ER, PR, HER2 and Ki67 on breast immunohistochemistry
field-of-view captures, including a stereological correction of counts to
volumetric density that needs no reconstructed volume.

## Two results worth reading first

**The published registration failed on this data, and the reason is specific.**
Running it as supplied left landmarks 2544 um apart where applying no transform
at all leaves them 727 um apart, and it reduced between-section image similarity
as well. The cause is rotation estimation from the Radon transform, which
averages 37.5 degrees of error against fiducial ground truth on this tissue.
Replacing it with a direct search over rotation, scored by the pixel correlation
the pipeline already uses, reduces rotation error to 3.9 degrees and target
registration error to 114 um. `src/coda_my/registration_fix.py` holds the
replacement; the original module is untouched.

**Overcounting is real but it is not a single number.** The ratio of objects
counted section by section to objects present in the volume is 1.25-fold when
every detected feature is counted, because most features occupy one section and
cannot be double counted, and 10.1-fold when restricted to structures above
10^6 um^3. Quoting one fold-change without stating what was counted and above
what size conveys almost nothing.

## Validation against the published benchmark

The serial dataset ships operator-annotated fiducials and has been evaluated by
twelve algorithm configurations, so accuracy here is placed directly against
them rather than reported in isolation.

The unregistered stack measured by this code gives a mean target registration
error of **726.85 um** against the **726.81 um** published for the same data.
That agreement to 0.04 um validates the landmark handling, coordinate
convention and unit conversion before any comparison is drawn. Against the
twelve published automated configurations, this implementation ranks **6th**,
improving on seven of them.

## Quick start

```bash
git clone <this repository>
cd coda-3d-histology
pip install -r requirements.txt

python RUN_EVERYTHING.py --list          # every step, with what it does and why
python RUN_EVERYTHING.py --skip-heavy    # everything except the 63 GB download
python RUN_EVERYTHING.py                 # the whole pipeline
python RUN_EVERYTHING.py --steps 3,5,7   # selected steps only
```

`RUN_EVERYTHING.py` detects Google Colab and adapts. Every step checks for its
own output and skips unless `--force`, so an interrupted run resumes rather than
restarting.

```bash
python -m pytest tests/ -q               # 60 tests
```

## Data

**This repository is code only.** No images, and no generated results: every
figure, table, report and manuscript is reproduced by running the pipeline. That
keeps the repository small and makes it impossible for a stale committed number
to disagree with what the code actually produces.

The one exception is `data/reference/`, which holds accuracy values transcribed
from the published benchmark. Those are an input the comparison needs, not an
output of this work.

The serial sections are openly licensed (CC BY 4.0, Etsin
`c76335fa-cdcf-4ddc-ab1c-1882bad82861`) but ship as a single 63.79 GB archive.
`scripts/fetch_kartasalo_liver.py` streams it and keeps only the liver stack,
about 15 GB, parsing the zip as it arrives so the archive never lands on disk.
The download service ignores HTTP range requests, so the transfer cannot resume
and cannot be fetched in part.

The breast immunohistochemistry captures are institutional patient material and
are **not redistributable**. They are excluded by `.gitignore` and have never
been committed. Derived measurements in `results/tables/` carry no patient
identifiers; slide references are sequential filenames only. Derived tables are
regenerated locally into `results/`, which is untracked.

## Layout

```
RUN_EVERYTHING.py          one script, every stage, VS Code and Colab
src/coda_my/               library
  registration.py            published method, unmodified
  registration_fix.py        replacement rigid estimator and two-scale driver
  loaders/kartasalo.py       serial stack and fiducial loading
  fibers.py deconv.py ...    stain separation, fibre tensor, QC
  reporting/                 figures, tables, report and manuscript assembly
scripts/                   one runnable step each
tests/                     60 tests, synthetic fixtures, no data downloads
docs/PROTOCOL.md           method parameters and their provenance
docs/FIGURE_SPEC.md        figure specification, panel by panel
config/coda_params.yaml    120 parameters, SHA-256 verified at run time
requirements.txt           dependencies

Generated at run time and untracked: results/, figures/, reports/, manuscript/
```

## Reproducibility

Parameters transcribed from the source publication live in
`config/coda_params.yaml` in a `locked` block hashed with SHA-256 and verified
at the start of every run; a drifted value fails the run and names the key.
Deliberate departures from the published method are declared in the same file
with reason and expected impact, and are reported in the manuscript rather than
applied silently.

Where a measurement could not be made, the corresponding figure is a labelled
placeholder naming the exact missing input. A reader should never be able to
mistake an absent result for a null one.

## Caveats stated up front

- The serial material is **mouse liver**, not breast. Stages 1 to 7 establish
  that the method runs and how accurately; they are not breast findings.
- **No three-dimensional reconstruction of breast tissue** appears anywhere, and
  none is possible here, because no serial breast sections exist in this study.
- Stage 6 counts vascular lumina separated by an intensity band, not the ten
  tissue classes of a trained segmentation.
- Stage 3 reports agreement between two detectors, which is a necessary but not
  sufficient condition for accuracy. The published 90 percent precision and
  recall bar is not tested, because no human annotation exists for this material.

## Citation

Kiemen A, et al. CODA: quantitative 3D reconstruction of large tissues at
cellular resolution. *Nat Methods* 2022;19:1490-1499. PMID: 36280719

Kartasalo K, et al. Comparative analysis of tissue reconstruction algorithms for
3D histology. *Bioinformatics* 2018;34:3013-3021. PMID: 29684099

## Licence

Code released under the MIT Licence. The serial section data is CC BY 4.0 and
belongs to its authors. The breast immunohistochemistry material is
institutional and is not distributed here.
