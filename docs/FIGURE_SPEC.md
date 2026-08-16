# FIGURE SPECIFICATION — mirroring the CODA Extended Data layout

Every figure below reproduces the layout of a specific CODA figure, panel for
panel, on data we can actually run. Where a panel cannot be produced, that is
stated with the reason rather than substituted.

Two figure sets:
- **Set A (3D)** — Kartasalo serial stacks. Mirrors CODA Extended Data 1 to 8.
- **Set C (2D)** — my breast IHC. Mirrors the CODA panel *style*, not its
  content, since CODA has no IHC.

Style throughout, matching the paper: white background, panels lettered lowercase
bold (a, b, c), scale bars burned in bottom-left with the length labelled,
colourblind-safe palettes, class colour legend as a vertical swatch block.

---

## SET A — 3D, Kartasalo serial stacks

### Figure A1 — Registration workflow
*Mirrors CODA Extended Data Fig 1a*

| Panel | Content |
|---|---|
| a-i | Stack schematic: image 1 … centre … image n, with curved arrows showing registration inward to the centre section |
| a-ii | Sample fixed image *n* and sample moving image *n+1*, side by side, raw |
| a-iii | Global registration point of reference marked with a cross on the fixed image |
| a-iv | Filtered fixed image, viridis colormap, after greyscale + background removal + complement + Gaussian |
| a-v | **Radon transform [0,360] degrees**, viridis, the characteristic double-lobe pattern |
| a-vi | **2D cross correlation** heatmap, colorbar 0 to 1, peak marked |
| a-vii | Local registration points of reference: cross markers on a grid over the tissue, 1.5 mm spacing |
| a-viii | Grid representation overlay of elastic registration results: red deformed mesh over the tissue |
| a-ix | Interpolated horizontal and vertical displacement fields, two panels, shared colorbar in pixels |
| a-x | Pre-registration / Global registration / Local registration, three overlays. Pre shows magenta-green offset; after registration both collapse to grey |

The magenta-green overlay is the most persuasive panel in the whole figure. Use
it. Fixed image in magenta, moving in green; perfect alignment reads as neutral
grey.

### Figure A2 — Registration accuracy vs the published benchmark
*Mirrors CODA Extended Data Fig 1b*

| Panel | Content |
|---|---|
| b-i | Image *N* with fiducial markers (green circles), zoom inset showing a nucleus split by the blade |
| b-ii | Image *N+1* with the matched fiducials (blue circles), same zoom |
| b-iii | **Normalised performance scatter**: x axis TRE, ATRE, RMSE, J, dA. Black diamonds = unregistered, red squares = **our** registration, grey circles = the seven other algorithms from Kartasalo 2018. Use the paper's normalisation formulas so higher is always better |

Panel b-iii is the single most important figure in the whole study. It places
our result directly against eight published methods on the same data.

Add a panel the paper does not have, because the liver stack allows it:

| b-iv | **Laser-cut hole deviation**: hole centroids plotted down z with the fitted straight line, deviation in microns annotated. Ground truth independent of the fiducials |

### Figure A3 — Cell detection validation
*Mirrors CODA Extended Data Fig 2a*

| Panel | Content |
|---|---|
| a-i | H&E validation image, 1.5 mm² region, scale bar 1 mm |
| a-ii | Hematoxylin channel, greyscale, same field |
| a-iii | **Manual counts**: two annotators overlaid, green diamonds = person 1, magenta squares = person 2 |
| a-iv | **Automatic counts**: red triangles = ours, yellow circles = HoVer-Net, cyan crosses = QuPath |
| a-v | Precision bar chart, three methods, individual region points overlaid as diamonds (vs person 1) and squares (vs person 2), y axis 0 to 100 |
| a-vi | Recall bar chart, same structure |

Run HoVer-Net and QuPath for real. A two-method comparison where we are the only
automatic method is not a comparison.

### Figure A4 — 2D to 3D cell count extrapolation
*Mirrors CODA Extended Data Fig 2b*

| Panel | Content |
|---|---|
| b-i | **Nuclear diameter table**, our measured values per tissue subtype, alternating row shading. MEASURE THESE. Do not copy the pancreas values |
| b-ii | Annotated H&E crop with individual nuclei measured, diameters labelled in yellow |
| b-iii | Effective thickness schematic: blue circle = nucleus diameter D, black bar = section thickness T, red dashed lines marking effective thickness T + D |
| b-iv | The formula C3D = Σ Σ C_image · 3T/(T + D_subtype) rendered in LaTeX |

### Figure A5 — Segmentation workflow and training data design
*Mirrors CODA Extended Data Fig 3*

| Panel | Content |
|---|---|
| a-i | Stack schematic: training & validation slices, testing slice, unannotated remainder |
| a-ii | Whole section with 50 annotations per class marked in class colours, pencil icon optional |
| a-iii | Bounding box of histology plus the extracted label mask, two stacked panels |
| b-i | **5% filled** tile: sparse annotation crops on black |
| b-ii | **65% filled** tile: the same tile after class-balanced overlay |
| b-iii | Zoom insets from both, green boxes linking them |
| b-iv | **Class percentage table**, showing pixels per class are near even. Ours should land near 12 to 13% each, as the paper's did |
| c-i | Large tile cut into small tiles for training and validation, with the label overlay below |
| c-ii | Unedited histology tiles for the testing set |
| c-iii | Whole slide image cut into tiles for semantic segmentation, with the class colour legend and two scale bars |

### Figure A6 — Segmentation accuracy
*Mirrors CODA Extended Data Fig 4a*

| Panel | Content |
|---|---|
| a | **Confusion matrices**, one per sample, blue colormap, counts in cells, per-class precision down the right edge and recall along the bottom, overall accuracy in the corner |

Mark every class that fails the 90% gate in red. The paper's protocol is to add
annotations and retrain until all classes pass, not to report the failures.

### Figure A7 — z-resolution validation
*Mirrors CODA Extended Data Fig 5a*

| Panel | Content |
|---|---|
| a-i | **Log correlation vs distance in mm**: dashed black = xy (the ceiling), solid black = unregistered, coloured lines = registered at 4, 8, 12, 16, 20 µm spacing |
| a-ii | **% change in cell count from 4 µm** bar chart, x axis 4, 8, 12, 24, 48, 90 µm, blue gradient, error bars |
| a-iii | **% change in composition from 4 µm** bar chart, same x axis, green gradient, error bars |

Draw a horizontal line at 5% on a-ii and a-iii, and at 95% correlation on a-i.
The paper's claim is that error stays under 5% to 12 µm and correlation over 95%
to 20 µm. Show whether ours does.

### Figure A8 — Tissue classes labelled
*Mirrors CODA Extended Data Fig 5b*

Grid of H&E crops, one per tissue class, each with the segmentation boundary
drawn as a **green outline**, class name below, scale bar on the last panel.

### Figure A9 — Fiber alignment and nuclear aspect ratio
*Mirrors CODA Extended Data Fig 5c and Fig 6c*

| Panel | Content |
|---|---|
| a | **Violin plots** of nuclear aspect ratio, axial vs longitudinal, grouped by structure, individual points overlaid, significance brackets |
| b | **Violin plots** of fiber anisotropy index, same grouping, y axis 0 to 0.8 |
| c | Inter-observer panel: person 1 vs person 2, marked ns where not significant |

Use the 3D volume to pick the axial and longitudinal planes, exactly as the
paper did. This is the one place having a volume lets us correct for sectioning
angle. Say so in the caption.

### Figure A10 — 3D reconstruction
*Mirrors CODA Extended Data Figs 6, 7, 8*

| Panel | Content |
|---|---|
| a | **Global reconstruction**: photorealistic H&E volume block on the left, single-class 3D render on the right, both with 0.5 cm scale bars |
| b | **Subregion renders**: five boxed volumes, each showing a different class combination against the H&E base plane |
| c | **Z-projections**: one panel per class on black, class name in the class colour bottom-left, plus an "All tissue" greyscale projection, plus the class colour legend block |

This is the figure that makes the work look like 3D pathology. Match the paper's
presentation exactly: black background for z-projections, class name in colour,
scale bar in grey.

---

## SET C — 2D, my breast IHC

CODA has no IHC figures, so these mirror the *style* and rigour, not the content.
Every caption must say REAL, name the marker and give n.

### Figure C1 — Image QC
| Panel | Content |
|---|---|
| a | Representative field per marker with the detected scale bar boxed and the masked overlay region shaded |
| b | mpp per image, scatter, showing the 3-fold magnification range |
| c | **Counterstain fraction per image** with the "absent" threshold drawn as a red line. Images below it are labelled NOT REPORTABLE for percent-positive |

### Figure C2 — Marker quantification
*Mirrors the style of Extended Data Fig 2a*

| Panel | Content |
|---|---|
| a | Original field, hematoxylin channel, DAB channel, three panels |
| b | Detected nuclei overlaid, positives red, negatives blue, only where counterstain permits |
| c | DAB intensity distribution per nucleus with the threshold marked |
| d | Percent positive per image, **only for images with adequate counterstain** |

### Figure C3 — HER2 membrane
| Panel | Content |
|---|---|
| a | Original HER2 field |
| b | Membrane mask |
| c | Enclosed cells labelled, colour-coded by completeness |
| d | Completeness distribution histogram |

Caption must state: quantitative descriptor, **not** an ASCO/CAP score.

### Figure C4 — Ki67 hotspot versus average
| Panel | Content |
|---|---|
| a | Field with the sliding reporting window drawn at the hotspot position |
| b | Heatmap of window-level percent positive across the field |
| c | Per image: average, hotspot, and the gap, with the 20% clinical cutoff drawn |

Panel c is the point of the whole arm. If cases sit on opposite sides of the 20%
line depending on where the window lands, that is the result.

### Figure C5 — Spatial statistics of positives
| Panel | Content |
|---|---|
| a | Positive nuclei as a point pattern, negatives in grey |
| b | Border-corrected Ripley L vs radius, with the CSR envelope |
| c | Clark-Evans index per image, 1.0 marked as the random line |
| d | KDE hotspot coefficient of variation vs percent positive, scatter |

State the radius range analysed and cap it at one quarter of the field width.
Ripley's K beyond that is unreliable even with border correction.

---

## What cannot be produced, and why

| CODA panel | Why not |
|---|---|
| Any 3D panel on my IHC | Field-of-view captures, not serial sections |
| Fiber alignment on my IHC | DAB-IHC has no eosin channel |
| PanIN 3D phenotypes | Pancreas-specific finding, not a method output |
| 87% collagen composition | Property of one human sample |
| Cancer vs normal cell density | No cancer in the Kartasalo tissue |

Emit a labelled placeholder for each of these in the report, naming the missing
input. Do not substitute something else and let the layout imply equivalence.
