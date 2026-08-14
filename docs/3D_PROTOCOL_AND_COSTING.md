# 3D CODA Protocol and Costing — Breast, USM

Document to take to your supervisor. It states exactly what 3D requires, what it
costs, and what it would produce. Nothing here can be done with slides that are
already cut.

---

## 1. Why the existing 200+ scans cannot be used

CODA (Kiemen et al., *Nat Methods* 2022;19:1490-1499) reconstructs a tissue
volume by registering **consecutive sections cut through one block**. The
published study processed 4,114 sections across 13 samples; one sample alone
had 101 consecutive sections at 4 µm spacing.

The current USM archive is four stains (ER, PR, HER2, Ki67) on separate sections
per patient. Four sections is four sections, however many patients contribute
them. There is no z-axis to reconstruct. Patient count does not substitute for
section depth within a block.

**The requirement is depth in one block, not breadth across patients.**

## 2. What CODA requires, exactly

| Parameter | CODA published value | Notes |
|---|---|---|
| Section thickness | 4 µm | Standard microtome |
| Sectioning | Continuous through the block | No gaps, no re-facing |
| H&E interval | Every 3rd section | Gives 12 µm axial resolution |
| Held-back sections | The 2 between each H&E | **This is where IHC goes** |
| Scan magnification | 20× (0.5 µm/px lateral) | Existing USM scanner is adequate |
| Sections per block | 100-150 H&E | For a ~5 mm block |
| Annotation | ~50 examples per tissue class on 7 sections | Per block |
| Segmentation target | >90% precision and recall per class | CODA's acceptance bar |

CODA validated that skipping two sections is safe: registration holds at >95%
pixel correlation across sections up to five planes apart, and 3D cell count and
tissue composition are preserved to within 5% error at 12 µm spacing. So the
every-third-section design loses almost nothing and cuts the workload by two
thirds.

## 3. The proposal worth making

CODA's own discussion proposes the extension and does not perform it:

> future addition of IHC labeling, spatial 'omics', and gene expression imaging
> to the intervening sections will increase the number of labels

**Serial H&E for 3D architecture, intervening sections for ER/PR/HER2/Ki67.**
Nobody has published this in breast. It uses the exact four markers already
established at USM and the existing IHC quantification pipeline.

### The question it answers

**Is Ki67 heterogeneity a 3D phenomenon that 2D scoring systematically
misrepresents?**

This is not a manufactured question. Ki67 scoring is irreproducible because
observers disagree on hotspot versus average assessment, and the 20% cutoff that
drives adjuvant chemotherapy decisions sits where reproducibility is worst.

CODA showed that 2D counting overcounted pancreatic precursor lesions by up to
40-fold, average 12.3-fold, because lesions that appear separate in one plane are
connected in 3D. If Ki67 hotspots behave the same way, then a hotspot seen on one
section may be a slice through a single connected proliferative domain, or it may
be one of several genuinely independent foci. The 2D image cannot distinguish
these, and they arguably mean different things biologically.

That is a real clinical question with a real answer, and 3D is the only way to
get it.

## 4. Scope: start with two blocks

CODA used 13 samples with a large multi-department team over several years. Two
blocks is the correct pilot.

Suggested selection: one Luminal with low Ki67, one TNBC with high Ki67. The
contrast is the point.

## 5. Costing per block

Fill in local USM rates. Rough structure:

| Item | Quantity per block | Unit cost | Subtotal |
|---|---|---|---|
| Serial sectioning at 4 µm | ~400 sections | | |
| H&E staining | ~130 slides | | |
| IHC staining (4 markers, spread across held-back sections) | ~40 slides | | |
| WSI scanning at 20× | ~170 slides | | |
| Storage (1-2 TB per block) | | | |
| **Total per block** | | | |
| **Two blocks** | | | |

Scanning is normally the dominant line. Get a per-slide quote before committing.

## 6. Compute requirements

| Need | Status at USM | Action |
|---|---|---|
| GPU for DeepLab/ResNet50 training | **Not available** (CPU-only Ryzen, 16 GB) | Required. CODA acknowledges an NVIDIA GPU grant. |
| MATLAB licence | Check | CODA repo is MATLAB: github.com/ashleylk/CODA |
| Storage, 2-4 TB | Check | |
| RAM, 32 GB+ recommended | Currently 16 GB | Registration of large stacks is memory-hungry |

The GPU is not optional. Semantic segmentation of 130 whole slide images per
block cannot be trained or run on CPU in reasonable time.

## 7. Ethics and block access — do this first

**Serial sectioning consumes the block.** Diagnostic material usually cannot be
exhausted without specific approval, because the block may be needed for future
clinical testing.

Requirements before any cutting:
- USM JEPeM ethics approval covering destructive use of archived diagnostic blocks
- Pathology department sign-off on block release
- Confirmation that no clinical need for the blocks remains
- If the TCGA comparison arm proceeds, ethics must also cover cross-cohort
  comparison and international public data use

This is the longest lead item. Start it now, in parallel with everything else.

## 8. Realistic timeline

| Phase | Duration |
|---|---|
| Ethics and block access | 2-4 months |
| Sectioning, staining, scanning (2 blocks) | 1-2 months |
| Registration and segmentation model training | 2-3 months |
| Reconstruction, analysis, IHC integration | 2-3 months |
| **Total** | **7-12 months** |

## 9. Recommended sequence

**Do the 2D study first.** It runs on slides already scanned, needs no ethics
amendment, no sectioning budget and no GPU, and it produces a publishable result
in months rather than a year. It also generates the preliminary data that makes
the 3D funding request credible.

The 2D pipeline is built and tested (`coda-brca-my`, 21 passing tests). It covers
collagen fiber alignment, nuclear morphometry, cell detection, quantified IHC
spatial heterogeneity, and the Malaysia versus TCGA comparison with batch-effect
auditing.

**Then use the 2D result to fund the 3D.** A grant application that says "our 2D
data show Ki67 spatial heterogeneity varies independently of the percentage
score, and we propose to test whether this reflects 3D connectivity invisible in
single sections" is far stronger than one that proposes 3D reconstruction with no
preliminary work.

---

## Questions to answer before committing

1. Are the FFPE blocks available and can they be exhausted?
2. Is there ethics approval for destructive use, or how long to obtain it?
3. What is the per-slide sectioning, staining and scanning cost at USM?
4. Is there GPU access, through USM HPC or otherwise?
5. Is there a MATLAB licence, or is a Python reimplementation acceptable?

If items 1, 2 or 4 are no, the 3D arm cannot proceed regardless of the science.
The 2D arm proceeds either way.
