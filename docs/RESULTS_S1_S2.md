# Results

## 1. Pipeline overview and datasets

We implemented the seven stages of CODA and applied them to three datasets, and
we report at the outset which stages each dataset can support, because the
constraint is a property of the material rather than a limitation of the
implementation. Stages 1, 2, 5 and 6, namely registration, registration quality
control, three-dimensional reconstruction and volumetric quantification, require
serial sections cut consecutively through a single block. Stage 7, fibre
alignment, requires an eosin channel. The applicability of every stage to every
dataset is given in Supplementary Table S14, and the runner refuses a stage
rather than producing a transform and a correlation coefficient for material
that cannot support it.

We obtained the mouse liver series of Kartasalo et al. (n = 47 consecutive
sections, 5 um nominal thickness, 0.46 um per pixel at 20x, mouse, liver),
published under CC BY 4.0 as a single 63.79 GB archive. Four laser-cut holes
were introduced into the block before embedding, and their positions were
annotated independently by two observers, giving 4 landmarks on each of 47
sections and 188 landmark positions per observer. The companion mouse prostate
series (n = 260 sections) resides in the same archive and was not retrieved. The
breast whole-slide resource with matched immunohistochemistry and registration
landmarks requires a data use agreement and was not obtained. Arm C comprises
234 breast immunohistochemistry field-of-view captures acquired at Universiti
Sains Malaysia, which are single fields rather than serial sections and support
stage 3 and marker quantification only.

Pixel size and section thickness are not recorded in the image files. The TIFF
headers carry a resolution tag of 72 dots per inch with unit inch, which is the
generic placeholder written by the encoder and not a microscope calibration. The
values used throughout, 0.46 um per pixel and 5 um section thickness, are taken
from the source publication and are propagated to every derived quantity with
that provenance recorded. Every distance we report in microns therefore inherits
that assumption, and a reader who disputes it can rescale all of them by a single
factor.

We distinguish reproduction of the CODA method from reproduction of its
biological findings. The former is a property of the algorithms, is expected to
transfer between tissues, and is what we test here. The latter, including the
12.3-fold two-dimensional to three-dimensional overcounting of pancreatic
precursor lesions, the 2.3-fold cell density decrease between healthy and
invasive pancreas, and the three-phenotype classification of PanIN, are findings
about human pancreatic cancer. They cannot reproduce in mouse liver and we do not
attempt to reproduce them. Where a measurement of ours corresponds to one of
theirs, we state whether we matched it, differed from it, or did not test it, and
why.

## 2. Registration of serial sections

We registered all 47 sections of the liver series, working outward from the
centre section rather than chaining each section to its immediate predecessor, so
that error accumulated over the stack is bounded by distance from the centre
rather than by stack length. Registration was solved in two stages at two
resolutions: a global rigid stage at 80.96 um per pixel and a local elastic stage
at 7.36 um per pixel. The published protocol specifies 80 um per pixel for the
global stage, and we note that the parameter declaring this value in the
configuration is never read by the registration code, so a caller who supplies a
single resolution silently performs the global stage at whatever resolution the
elastic stage requires. We supplied the two resolutions separately.

Rotation was estimated by direct search rather than from the Radon transform. For
each candidate angle we rotated the moving section, recovered translation by
phase cross correlation, and scored the result by the same pixel correlation used
to judge registration quality, retaining the angle with the highest score. The
search was bounded to plus or minus 45 degrees. This bound is a prior about
serial sections mounted by hand rather than a fitted parameter, and it excludes
no true value in this series, where the largest fiducial-implied rotation between
consecutive sections is 34.1 degrees (n = 46 pairs). The reason for replacing the
published estimator is given in Section 3.

Quantification revealed a median pixel correlation between each registered
section and its reference of 0.929 (interquartile range 0.872 to 0.959, minimum
0.473, n = 47 sections). No section fell below the acceptance threshold of 0.30
at which the protocol flags a section as torn, folded or defectively stained
(0 of 47). The centre-out procedure permits each section to register against any
of the three already-registered sections nearest the centre; the immediate
neighbour was selected for 25 of 47 sections, the second nearest for 13 and the
third nearest for 9, indicating that the nearest section is not always the best
reference and that offering alternatives is used in practice rather than being
redundant.

The transforms applied were of the magnitude expected for hand-mounted sections.
Median absolute rotation was 10.0 degrees (maximum 30.0, n = 47) and median
translation was 478 um (maximum 1134 um, n = 47). Registration took 1.0 minute
for the full stack, against 199.8 minutes for the single-resolution
implementation, because solving the rigid stage on a coarse copy of the data is
both cheaper and, as Section 3 shows, more accurate than solving it at the
resolution the elastic stage requires.

We measured the effect of registration on landmark agreement directly. Mean
target registration error between consecutive sections was 114 um after
registration against 727 um for the same landmarks with no transform applied, a
6.35-fold reduction (median 74 um against 681 um, Wilcoxon signed rank test,
W = 2, p = 8.5e-14, n = 46 section pairs). The comparison against no transform is
reported because it is the weakest test a registration must pass, and Section 3
shows that it is not passed automatically.
