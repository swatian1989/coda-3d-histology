# coda-brca-my

CODA-derived 2D histomorphometry for a Malaysian breast cohort, compared against
TCGA-BRCA, with quantified ER/PR/HER2/Ki67 spatial analysis.

Reference: Kiemen A et al. CODA. *Nat Methods* 2022;19:1490-1499.

## Two arms

**2D — runs on slides you already have.** No new sectioning, no ethics
amendment, no GPU.

**3D — needs blocks cut.** See `docs/3D_PROTOCOL_AND_COSTING.md`. 7-12 months,
requires ethics approval for destructive block use, sectioning and scanning
budget, and GPU access.

Do the 2D first. It produces the preliminary data that makes the 3D funding
request credible.

## What is built

| Module | Does |
|---|---|
| `deconv.py` | Ruifrok-Johnston deconvolution, per-image stain vector estimation by k-means (CODA method), DAB separation for IHC |
| `fibers.py` | Structure-tensor fiber anisotropy index, orientation maps, tiled distributions, boundary-relative orientation for the invasive front |
| `ihc.py` | Nuclei detection from hematoxylin, per-nucleus DAB scoring, point-pattern export, hotspot-versus-average quantification |
| `cohort_compare.py` | Batch sensitivity audit, propensity matching, Mann-Whitney with Cliff's delta and BH-FDR |

21 tests passing. Verified behaviour: anisotropy 1.00 for aligned fibers, 0.04
for isotropic, rotation-invariant; batch audit correctly flags a
scanner-driven feature as untrustworthy despite it being significant.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install numpy pandas scipy scikit-learn scikit-image statsmodels \
            matplotlib seaborn pytest openslide-python pyarrow
python -m pytest tests/ -q
```

## The first thing to run

Not the full cohort. A batch audit on 10 USM and 10 TCGA slides:

```python
from coda_my.cohort_compare import audit_batch_sensitivity
audit = audit_batch_sensitivity(features, metadata)
print(audit)
```

If most features come back scanner-confounded, processing 200 slides will not
fix it. Find that out on 20 slides, not 200.

## What the IHC is for

Not for comparing to TCGA. TCGA has categorical receptor status from pathology
reports; comparing your DAB images to their H&E measures the stain, not the
population.

The IHC is your within-cohort advantage: continuous, spatially resolved
expression that public data cannot provide. The question worth asking is whether
Ki67 *arrangement* carries information the percentage score discards.

## Known limitations, state them explicitly

- Sectioning angle confounds any single-section fiber alignment measurement.
  Use distributions or boundary-relative orientation.
- HER2 is membranous; nuclear DAB scoring is invalid and the code refuses it.
- Cohort differences may be scanner, not population. The audit is not optional.
- "Malaysian" is not an ancestry group. Stratify by Malay, Chinese, Indian.
