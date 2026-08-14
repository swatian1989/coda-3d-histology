# CODA on Public Data — Protocol

Reproducing Kiemen et al., *Nat Methods* 2022;19:1490-1499 on public datasets,
with no tissue cutting and no ethics application.

Two datasets, two arms. Both start today.

---

## ARM A — Full 3D CODA on the Kartasalo serial stacks

**This is the only public data where the complete CODA protocol runs**, because
it is the only public data with true serial sections.

| | Value |
|---|---|
| Dataset | Kartasalo et al., *Bioinformatics* 2018;34:3013 |
| Benchmark code | github.com/BioimageInformaticsTampere/RegBenchmark |
| Mouse prostate | **260 serial sections**, 20×, 0.46 µm/px, 5 µm thick |
| Mouse liver | 47 serial sections, same resolution |
| Ground truth | Fiducial points marked by human operators on structures visible in both sections, preferentially nuclei split by the sectioning blade |
| Liver extra | Laser-cut holes through the whole tissue, which should form a straight line after correct alignment |

**Why this dataset specifically.** CODA benchmarked its own registration on it
against seven competing methods, and VALIS did the same. So your TRE, ATRE,
RMSE, Jaccard and dA numbers are directly comparable to published values. You
get an objective accuracy result without annotating anything yourself.

### Steps

1. **Register** — `register_stack()`. Radon-transform global registration at
   80 µm/px, three candidate references per section, elastic field from 1.5 mm
   tiles at 8 µm/px smoothed with σ = 2 px, everything referenced to the
   **centre** section.
2. **Validate registration** — `target_registration_error()` and
   `accumulated_tre()` against the supplied fiducials. Compare to the published
   table. ATRE is the metric that matters: pairwise error can look excellent
   while the stack banana-bends.
3. **Detect cells** — hematoxylin deconvolution, 2 µm/px, minima detection.
   Validate with `cell_detection_metrics()` at the 2 µm tolerance CODA used.
4. **Segment tissue** — DeepLab v3+ on ResNet-50, ~50 annotations per class on
   7 evenly spaced sections, 9000×9000 training tiles filled to >65% with
   class-balanced bounding boxes, cut into 500×500 tiles, trained to >90%
   precision and recall per class. **Needs a GPU.**
5. **Reconstruct** — `stack_to_volume()`. Resample 2 × 2 × 12 µm voxels to
   isotropic 12 × 12 × 12 µm. Every volume statistic assumes isotropy.
6. **Quantify** — composition, 3D cell counts with the nuclear-diameter
   correction, connectivity, per-object morphology, z-projections.
7. **Reproduce the headline result** — `overcounting_ratio()`. CODA found 2D
   section counting overcounted true 3D lesion number by a mean of 12.3-fold,
   up to 40-fold, p < 1e-5.
8. **Validate z-resolution** — `z_skip_validation()` and
   `axial_vs_lateral_correlation()`. CODA reported >95% correlation retained
   skipping up to four sections, and <5% error in cell count and composition
   skipping up to two.

### Limitation to state plainly

Mouse, prostate and liver. It validates the pipeline, not a breast finding.

---

## ARM B — CODA 2D components on ACROBAT breast

**The public twin of the USM cohort.**

| | Value |
|---|---|
| Dataset | ACROBAT, Weitz et al., *Med Image Anal* 2024 |
| Access | researchdata.se/en/catalogue/dataset/2022-190-1 |
| Size | 4,212 WSIs, 1,153 breast patients |
| Stains | H&E, **ER, PGR, HER2, KI67** — exactly the USM panel |
| Landmarks | 37,208 manually annotated correspondences |
| Resolution | released at 0.92 µm/px (10×), native 0.23 µm/px |
| Scanners | NanoZoomer S360 and two NanoZoomer XR |
| Source | CHIME study, Södersjukhuset Stockholm, routine diagnostics |

**Critical limitation.** All WSIs of a case come from the same tumour block but
the sections are **not necessarily consecutive**. ACROBAT supports H&E-to-IHC
registration and every 2D measurement. It does **not** support 3D
reconstruction. Do not stack ACROBAT sections into a volume.

### Steps

1. **Register H&E to IHC** — the same `global_register()` and `elastic_field()`
   machinery. Validate against the 37,208 landmarks. This is an objective
   accuracy number on breast tissue with your exact markers.
2. **Quantify IHC** — `score_marker()` for ER, PGR, KI67.
   **HER2 is membranous**; nuclear DAB scoring is invalid and the code refuses
   it by design.
3. **Fiber alignment** — `tiled_anisotropy()` and
   `boundary_relative_orientation()` on the eosin channel of the H&E.
4. **Ki67 spatial heterogeneity** — `hotspot_vs_average()` and
   `to_point_pattern()` into the border-corrected spatial statistics.
5. **Cohort comparison** — ACROBAT (Sweden) and TCGA (USA) become **two**
   independent comparison cohorts for the Malaysian arm. Run
   `audit_batch_sensitivity()` first: three scanners are involved in ACROBAT
   alone.

---

## Sequence

| Step | Data | Ethics | GPU | Start |
|---|---|---|---|---|
| 1. Full 3D CODA | Kartasalo mouse | No | Yes, step 4 only | **Today** |
| 2. 2D + H&E↔IHC | ACROBAT breast | No | No | **Today** |
| 3. Apply to USM | Malaysian slides | **Yes** | No | After approval |
| 4. Three-cohort comparison | USM + ACROBAT + TCGA | Yes | No | After step 3 |

Only step 3 needs ethics. Steps 1, 2 and 4 are public.

Doing steps 1 and 2 first is not a delay. It gives you validated registration
accuracy on landmarked data before you touch patient slides, which is exactly
what a reviewer will ask for, and it is the preliminary data that makes a
serial-sectioning grant application credible.

---

## Parameters that must not drift

All from the CODA Online Methods.

| Parameter | Value |
|---|---|
| Global registration resolution | 80 µm/px |
| Elastic field resolution | 8 µm/px |
| Elastic tile interval | 1.5 mm |
| Field smoothing | Gaussian σ = 2 px |
| Candidate references per section | 3 |
| Registration reference | **Centre section**, never the neighbour |
| Cell detection resolution | 2 µm/px |
| Segmentation resolution | 2 µm/px |
| Cell-match tolerance | 2 µm |
| Segmentation acceptance | >90% precision and recall per class |
| Training tile | 9000×9000, >65% filled, cut to 500×500 |
| Annotations | ~50 per class on 7 sections |
| Isotropic voxel | 12 × 12 × 12 µm |
| 3D cell count | C3D = Σ C_image × 3T / (T + D_subtype) |
| Statistics | Wilcoxon rank sum throughout |

**Measure your own nuclear diameters.** The values in `reconstruct.py` are
CODA's pancreas measurements. The correction scales cell count directly, so
borrowed diameters bias every density you report, and bias it most for the
largest nuclei — a systematic error between cell types, not a random one.
