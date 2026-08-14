# FINAL METHODS — every CODA stage, mapped to every dataset

Measured facts about the four uploaded images come first, because they change
what is possible.

---

## 1. What your images actually are

Read directly from the files, not assumed:

| File | Scale bar | Recovered mpp | Counterstain | Marker (from pattern) |
|---|---|---|---|---|
| viewcapture1 | 135 px = 30 µm | **0.222 µm/px** | marginal, 1.32% | nuclear + cytoplasmic |
| viewcapture32 | 118 px = 50 µm | **0.424 µm/px** | **absent, 0.48%** | nuclear |
| viewcapture49 | 113 px = 80 µm | **0.708 µm/px** | **absent, 0.02%** | **membranous — HER2** |
| viewcapture55 | 116 px = 80 µm | **0.690 µm/px** | marginal, 1.33% | nuclear, blue negatives visible |

Three findings that constrain the analysis:

**(a) These are field-of-view captures, not whole slide images.** ~1611 px wide
with a burned-in scale bar. CODA operates on whole slides at 20×. FOV images
support cell-level and local-neighbourhood measurement; they do not support
whole-slide tissue architecture, and the spatial statistics have a hard extent
limit set by the field, roughly 350 µm to 1150 µm across depending on objective.

**(b) The scale bar is both calibration and contamination.** It gives you mpp,
without which every distance is in pixels. It is also a saturated high-contrast
object that any nucleus detector will segment and any spatial statistic will
read as a dense corner cluster. `mask_overlay_region()` must run before every
measurement.

**(c) Counterstain is marginal to absent.** This is the one that matters most.
Two of four images have essentially no hematoxylin, so **negative nuclei are
invisible**. There is no denominator. A Ki67 index or Allred score cannot be
computed from those images and must not be back-calculated from DAB area.
Positive-cell density and spatial arrangement remain fully valid.

If you have matched images WITH counterstain, or the original slides can be
re-captured with counterstain visible, that single change unlocks
percent-positive analysis across the whole cohort.

**(d) HER2 works, on its own axis.** Measured on viewcapture49: 124 enclosed
cells, mean membrane completeness **0.999**, median cell area 51.9 µm². That is
the complete circumferential pattern. It is a quantitative descriptor, not an
ASCO/CAP category, and must never be reported as 0/1+/2+/3+.

---

## 2. Every CODA stage, against every dataset

| # | CODA stage | Kartasalo serial | ACROBAT breast | TCGA | **Your IHC FOVs** |
|---|---|---|---|---|---|
| 1 | Nonlinear registration | **Yes** | Yes, H&E↔IHC | No | No |
| 2 | Registration QC (TRE, ATRE, z-corr) | **Yes**, has landmarks | Yes, 37,208 landmarks | No | No |
| 3 | Cell detection (hematoxylin deconv) | **Yes** | **Yes** | **Yes** | **Partly** — positives only where counterstain absent |
| 4 | DeepLab v3+ segmentation | **Yes** [GPU] | **Yes** [GPU] | **Yes** [GPU] | Limited — FOV shows one tissue type, not 10 |
| 5 | 3D reconstruction | **Yes** | No | No | No |
| 6 | Quantification, connectivity, 2D-vs-3D | **Yes** | No | No | No |
| 7 | Fiber alignment (eosin channel) | **Yes** | **Yes** | **Yes** | No — DAB-IHC has no eosin |
| + | DAB marker quantification | No | **Yes** | No | **Yes** |
| + | HER2 membrane completeness | No | **Yes** | No | **Yes** |
| + | Spatial statistics of positives | No | **Yes** | No | **Yes**, within FOV extent |

Reading across: the serial dataset is the only place stages 1, 2, 5 and 6 can
run at all. Your images are the only place the marker analyses can run. They are
complementary, not competing.

---

## 3. The three-arm design

### Arm A — full 7-stage CODA, Kartasalo serial sections

260 mouse prostate sections, 0.46 µm/px, 5 µm thick, with operator-marked
fiducials. CODA benchmarked its own registration on this data against seven
competitors, so your TRE and ATRE land on a published scale.

Deliverable: registration accuracy versus the published benchmark, and the
2D-versus-3D overcounting result (paper: mean 12.3-fold, max 40-fold).

Limitation to state: mouse, prostate. Validates the pipeline, not a breast
finding.

### Arm B — stages 1, 2, 3, 4, 7 on ACROBAT breast

4,212 WSIs, 1,153 patients, H&E + ER/PGR/HER2/KI67, 37,208 landmarks. Your
exact marker panel, public, whole slides.

This is where H&E↔IHC registration gets validated objectively on breast tissue,
and where fiber alignment can run because there is an eosin channel.

Limitation: sections are from the same block but **not necessarily consecutive**.
No 3D.

### Arm C — marker analysis on your images

What runs today, with code already written and tested against your files:

1. `detect_scale_bar()` → mpp per image. Verified on all four.
2. `mask_overlay_region()` → exclude the bar and label.
3. `has_counterstain()` → grade adequate / marginal / absent. This gate decides
   whether a percentage is reportable for that image.
4. ER, PR, Ki67 → `score_marker()`, per-nucleus DAB, positive density,
   `hotspot_vs_average()`.
5. HER2 → `membrane_completeness()`, never nuclear scoring.
6. Spatial pattern → `to_point_pattern()` into the border-corrected Ripley K
   and L, Donnelly-corrected Clark-Evans, quadrat dispersion, KDE hotspot CV.

The scientific question Arm C can answer that public data cannot: **does the
spatial arrangement of positive nuclei carry information the percentage
discards?** A clustered 18% and a dispersed 18% receive the same score and the
same treatment decision. Nobody quantifies the difference routinely.

---

## 4. Honest limits of Arm C

State these in the methods, do not work around them:

- **FOV, not WSI.** Spatial statistics are bounded by the field. Ripley's K
  beyond ~1/4 of the field width is unreliable even with border correction.
  Report the radius range analysed.
- **Non-random field selection.** Photographs taken to document positive
  staining are biased toward positive areas. If fields were chosen by eye, the
  sample is not representative of the slide and no whole-slide inference
  follows. If they were chosen systematically, say how.
- **Mixed magnifications.** 0.22 to 0.71 µm/px across four images, a 3-fold
  range. Every measurement must be converted to microns before pooling, and
  scale-dependent texture features are not comparable across objectives at all.
- **No denominator where counterstain is absent.** Density per mm², yes.
  Percent-positive, no.
- **Serial sections do not exist here**, so nothing in stages 1, 2, 5, 6
  applies to your data at any sample size.

---

## 5. Sequence

| Order | Arm | Needs | Start |
|---|---|---|---|
| 1 | C — your images | Nothing. Code is written and tested. | **Now** |
| 2 | B — ACROBAT | Download, ~free | This week |
| 3 | A — Kartasalo 3D | Download + GPU for stage 4 | When GPU available |

Arm C first because it uses data you already have and answers the question only
you can ask. Arms A and B are the methodological validation that makes Arm C
defensible.
