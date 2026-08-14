## Title

**Spatial arrangement of Ki67-positive nuclei explains hotspot versus average scoring discordance in breast cancer immunohistochemistry**

## Authors

Akbar Ali^1,\*^

1. Department of Chemical Pathology, School of Medical Sciences, Universiti Sains Malaysia, Health Campus, Kota Bharu, Kelantan, Malaysia

\* Corresponding author. ORCID 0009-0003-6543-3122.

## Abstract

**Background.** Ki67 immunohistochemistry informs chemotherapy decisions in breast cancer, yet scoring is poorly reproducible and the 20 percent cutoff sits where disagreement is worst [9,10,11]. A central source of disagreement is whether to score a hotspot or a field average. Arrangement of positive nuclei is not quantified routinely.

**Methods.** Measurement components of the CODA framework [1] were applied to 234 breast immunohistochemistry field-of-view captures (ER, PR, HER2, Ki67). Resolution was recovered per image from the scale bar and counterstain graded by colour deconvolution [8]. Each Ki67 image was scored twice from the same nuclei, as a field average and as the maximum over a sliding 500 micron window; positives were analysed as a point pattern with border-corrected estimators.

**Results.** Counterstain was absent on 150 of 234 images, removing the denominator. Across 62 evaluable Ki67 images the hotspot score exceeded the average by a median of 5.8 percentage points (mean 9.2, bootstrap 95 percent CI 6.8 to 11.9; Wilcoxon p = 8.6e-11). On 14 images (23 percent) the average fell below the 20 percent cutoff while the hotspot reached it. Positives were clustered in 51 of 53 images. Coarse-scale clustering predicted the gap (quadrat variance to mean ratio, rho 0.66, 95 percent CI 0.47 to 0.79, q = 4.1e-07); nearest-neighbour clustering did not (rho -0.03, q = 0.84).

**Conclusions.** Hotspot versus average discordance is a measurable property of spatial organisation, driven by large-scale patchiness rather than local nucleus arrangement. A statistic computed at the scale of the reporting window flags cases at risk of a scoring-dependent decision. Findings are bounded by field-of-view sampling.

## Keywords

Ki67; breast cancer; immunohistochemistry; spatial statistics; reproducibility; digital pathology; histomorphometry

## Introduction

Quantitative histology has moved from describing tissue to measuring it. CODA reconstructs large tissue volumes from serial sections at cellular resolution and showed that counting structures on a single two-dimensional section overcounts their true three-dimensional number by a mean of 12.3-fold in pancreatic precursor lesions [1]. That result depends on having a volume, which in turn depends on serial sections; it cannot be obtained from a single section at any sample size.

Breast cancer practice depends on a different kind of measurement. Ki67 immunohistochemistry estimates proliferative fraction and informs adjuvant chemotherapy decisions, and a cutoff near 20 percent is widely used. The International Ki67 in Breast Cancer Working Group has repeatedly documented poor reproducibility [9,11], and a formal reproducibility study found substantial variation between laboratories scoring identical material [10]. A recognised contributor is the choice between scoring a hotspot and scoring an average.

What is not measured in routine practice is how the positive nuclei are arranged. Two tumours with an identical index, one with positivity concentrated in a focus and one with positivity dispersed evenly, receive the same score and the same decision. Spatial arrangement of stromal collagen has been shown to carry prognostic information in breast carcinoma [15,14], which establishes that spatial organisation in this tissue is informative, but the analogous question for proliferation marker arrangement is rarely asked.

We therefore applied the measurement components of CODA that a single section can support to a breast immunohistochemistry series, and asked whether the spatial arrangement of Ki67-positive nuclei explains the discordance between the two scoring conventions. We also state explicitly which stages of the framework could not be run and why, because applying a serial-section method to non-serial data produces output that looks like a result and is not one.

## Materials and Methods

### Study material

234 digital field-of-view captures of breast immunohistochemistry were analysed, comprising ER, PR, HER2 and Ki67. Images are microscope captures rather than whole-slide images, approximately 1611 pixels wide, each carrying a burned-in red scale bar. No patient identifiers were used at any stage.

### Parameter provenance and applicability control

All method parameters were transcribed from the CODA Online Methods [1] into a configuration file whose locked block is hashed with SHA-256 and verified at the start of every run; a changed value fails the run and is named. The full set of 120 locked parameters is given in Supplementary Table S2 and every deliberate deviation is listed in Supplementary Table S3 with its reason and expected impact.

Stage applicability is enforced in software rather than left to judgement. Nonlinear registration, registration quality control, three-dimensional reconstruction and volumetric quantification require serial sections and are refused on datasets that lack them. Fibre alignment requires an eosin channel and is refused on diaminobenzidine immunohistochemistry, which has none.

### Scale recovery, overlay masking and counterstain grading

Microns per pixel was recovered for each image from the burned-in scale bar. The bar length in pixels was measured within the lower twelve percent of the frame; searching the whole frame allowed a streak of diaminobenzidine to exceed the bar in length and displaced the measurement on 16 of 234 images, which would have rescaled each of those by a constant factor. The micron value printed beside the bar is not present in the filenames and was read from the image for 231 images; three could not be read and were excluded rather than assigned a default. The recovered calibration reproduced four independently verified reference values exactly (0.222, 0.424, 0.690 and 0.708 microns per pixel).

The scale bar is also a contaminant, being a saturated high-contrast object that nucleus detection segments and spatial statistics read as a dense corner cluster, so the overlay bounding box was masked before any measurement.

Counterstain adequacy was graded by colour deconvolution [8], counting a pixel as counterstained nucleus where haematoxylin concentration exceeded 0.15 and exceeded diaminobenzidine. Where counterstain is absent there are no visible negative nuclei and therefore no denominator; percent positive was withheld for those images and was never derived from stained area.

### Marker quantification

ER, PR and Ki67 were scored per nucleus for diaminobenzidine positivity at the measured resolution of each image. Positive-cell density per square millimetre was computed for every image because it requires no denominator. Percent positive was computed only where the counterstain gate permitted.

HER2 is a membranous marker and per-nucleus diaminobenzidine scoring of it is invalid [12,13]; the implementation raises an error on that operation. Membrane completeness was measured instead, as the fraction of each enclosed cell boundary that is stained. These values are quantitative descriptors of the staining pattern and are not an ASCO/CAP category; they are not reported as 0, 1+, 2+ or 3+.

Images coarser than 2.5 microns per pixel cannot resolve a nucleus and were excluded from nuclear analysis with the reason recorded.

### Hotspot versus average scoring, and spatial statistics

Each image was scored twice from the same detected nuclei. The average score is the positive fraction across the whole field. The hotspot score is the maximum positive fraction over a sliding 500 micron window containing at least 100 nuclei, which approximates the field a pathologist would select. The difference between them is therefore attributable to the scoring convention alone.

Positive nuclei were converted to a labelled point pattern and characterised with border-corrected Ripley K and L, the Clark-Evans index with Donnelly's perimeter correction, quadrat dispersion as a variance to mean ratio, and the coefficient of variation of a kernel density estimate. Border correction is essential on field-of-view captures, where a large fraction of the field lies within one analysis radius of an edge and an uncorrected estimator reads the missing area as reduced clustering. Radii were capped per image at one quarter of the field width and the limit used is recorded with every value.

### Statistics

Paired comparisons used the Wilcoxon signed rank test, and the Wilcoxon rank sum test was the specified test for unpaired comparisons [1]. Effect sizes are reported with confidence intervals rather than p values alone; the mean scoring gap carries a bootstrap 95 percent confidence interval from 2000 resamples, and correlation coefficients carry Fisher z-transformed intervals. Multiple comparisons across the four spatial statistics were controlled by the Benjamini-Hochberg procedure and q values are reported alongside p. No sample size was predetermined and no data were excluded other than for the stated technical gates.

## Results

### Image quality determines what each marker can support

Of 234 images, 225 were analysed and 9 were excluded, three for an unreadable scale bar and six for insufficient resolution. The series spans a 38-fold range of magnification, 0.197 to 7.50 microns per pixel (Figure 1A).

Counterstain was absent on 150 of 234 images, and the distribution across markers determines what each can support (Figure 1B, 1C). Ki67 retained an adequate or marginal counterstain on 70 of 76 images, whereas percent positive was reportable on only 9 of 62 ER images and 2 of 39 PR images. For the remainder, positive-cell density and spatial arrangement remain valid and are reported, while percent positive is withheld.

HER2 membrane completeness was measured on 52 images, median 0.998 (Figure 2).

![Figure 1](D:/UCAS project/new analysis/coda-brca-my/coda-brca-my/figures/F17_usm_qc.png)

**Figure 1. Image quality control.**

![Figure 2](D:/UCAS project/new analysis/coda-brca-my/coda-brca-my/figures/F18_marker_quant.png)

**Figure 2. Marker quantification.**

![Figure 2b](D:/UCAS project/new analysis/coda-brca-my/coda-brca-my/figures/F19_her2_membrane.png)

**Figure 2b. HER2 membrane completeness.**

**Supplementary Table S11** (234 rows)

| filename   | marker   |   mpp_um_per_px | magnification_tier   | counterstain_grade   |   counterstain_fraction | reportable   |   note |
|:-----------|:---------|----------------:|:---------------------|:---------------------|------------------------:|:-------------|-------:|
| 001.png    | ER       |          0.2609 | high                 | absent               |                 3e-05   | False        |    nan |
| 0010.png   | ER       |          0.625  | usable               | marginal             |                 0.0138  | True         |    nan |
| 002.png    | ER       |          0.5714 | high                 | absent               |                 2e-05   | False        |    nan |
| 003.png    | ER       |          0.531  | high                 | absent               |                 5e-05   | False        |    nan |
| 004.png    | ER       |          0.5085 | high                 | absent               |                 0.00015 | False        |    nan |
| 005.png    | ER       |          0.4808 | high                 | absent               |                 6e-05   | False        |    nan |
| 006.png    | ER       |          2.2472 | marginal             | absent               |                 0.00687 | False        |    nan |
| 007.png    | ER       |          0.3788 | high                 | marginal             |                 0.02115 | True         |    nan |
| 008.png    | ER       |          0.6897 | usable               | absent               |                 0.00685 | False        |    nan |
| 009.png    | ER       |          0.4274 | high                 | marginal             |                 0.01425 | True         |    nan |
| 011.png    | ER       |          0.3817 | high                 | absent               |                 0.00646 | False        |    nan |
| 012.png    | ER       |          1.0989 | usable               | absent               |                 0.00818 | False        |    nan |

**Supplementary Table S12** (225 rows)

| filename   | marker   |   mpp_um_per_px | counterstain_grade   |   n_nuclei_detected |   n_positive |   positive_density_per_mm2 |   percent_positive | percent_reportable   |   n_enclosed_cells |   mean_membrane_completeness |   median_cell_area_um2 | note                                                  |
|:-----------|:---------|----------------:|:---------------------|--------------------:|-------------:|---------------------------:|-------------------:|:---------------------|-------------------:|-----------------------------:|-----------------------:|:------------------------------------------------------|
| 001.png    | ER       |          0.2609 | absent               |                   2 |            1 |                   10.351   |           nan      | False                |                nan |                          nan |                    nan | counterstain absent: no denominator, percent withheld |
| 0010.png   | ER       |          0.625  | marginal             |                 619 |          336 |                  606.049   |            54.2811 | True                 |                nan |                          nan |                    nan | nan                                                   |
| 002.png    | ER       |          0.5714 | absent               |                   1 |            0 |                    0       |           nan      | False                |                nan |                          nan |                    nan | counterstain absent: no denominator, percent withheld |
| 003.png    | ER       |          0.531  | absent               |                  23 |            2 |                    4.99769 |           nan      | False                |                nan |                          nan |                    nan | counterstain absent: no denominator, percent withheld |
| 004.png    | ER       |          0.5085 | absent               |                  26 |           13 |                   35.4234  |           nan      | False                |                nan |                          nan |                    nan | counterstain absent: no denominator, percent withheld |
| 005.png    | ER       |          0.4808 | absent               |                   2 |            1 |                    3.04789 |           nan      | False                |                nan |                          nan |                    nan | counterstain absent: no denominator, percent withheld |
| 006.png    | ER       |          2.2472 | absent               |                 474 |          181 |                   25.2536  |           nan      | False                |                nan |                          nan |                    nan | counterstain absent: no denominator, percent withheld |
| 007.png    | ER       |          0.3788 | marginal             |                 293 |           68 |                  333.901   |            23.2082 | True                 |                nan |                          nan |                    nan | nan                                                   |
| 008.png    | ER       |          0.6897 | absent               |                 390 |          219 |                  324.379   |           nan      | False                |                nan |                          nan |                    nan | counterstain absent: no denominator, percent withheld |
| 009.png    | ER       |          0.4274 | marginal             |                 317 |          216 |                  833.131   |            68.1388 | True                 |                nan |                          nan |                    nan | nan                                                   |
| 011.png    | ER       |          0.3817 | absent               |                 214 |          179 |                  865.64    |           nan      | False                |                nan |                          nan |                    nan | counterstain absent: no denominator, percent withheld |
| 012.png    | ER       |          1.0989 | absent               |                2387 |         1959 |                 1143       |           nan      | False                |                nan |                          nan |                    nan | counterstain absent: no denominator, percent withheld |

### Hotspot and average scoring disagree, and disagree decisively at the cutoff

Across 62 Ki67 images with a valid denominator, the median average score was 11.9 percent and the median hotspot score 19.9 percent. The hotspot score exceeded the average by a median of 5.8 percentage points (interquartile range 0.9 to 15.2, maximum 52.7). The mean difference was 9.2 percentage points with a bootstrap 95 percent confidence interval of 6.8 to 11.9, and the paired Wilcoxon signed rank test gave p = 8.6e-11 (Figure 3).

At the 20 percent cutoff, 14 images (23 percent) were discordant, with the average below the cutoff and the hotspot at or above it. Because both scores derive from the same nuclei on the same image, the change in category is attributable to the scoring convention and to nothing else.

![Figure 3](D:/UCAS project/new analysis/coda-brca-my/coda-brca-my/figures/F20_ki67_hotspot.png)

**Figure 3. Ki67 hotspot versus average scoring.**

**Supplementary Table S13** (62 rows)

| filename                          |   mpp_um_per_px |   n_nuclei_detected |   ki67_average_percent |   ki67_hotspot_percent |   ki67_hotspot_minus_average |   ki67_n_windows | crosses_20pc_cutoff   |
|:----------------------------------|----------------:|--------------------:|-----------------------:|-----------------------:|-----------------------------:|-----------------:|:----------------------|
| _20250203150520_viewcapture64.png |          1.0101 |                 468 |               11.3248  |               11.1111  |                    -0.213675 |                6 | False                 |
| _20250203150520_viewcapture65.png |          0.381  |                 190 |                7.89474 |                7.55814 |                    -0.336597 |                1 | False                 |
| _20250203150533_viewcapture66.png |          0.9346 |                9615 |                7.5299  |               10.9627  |                     3.43281  |               23 | False                 |
| _20250203150533_viewcapture67.png |          0.5882 |                 837 |               12.0669  |               16.4835  |                     4.41661  |                8 | False                 |
| _20250203150546_viewcapture69.png |          1.0526 |                 832 |                9.13462 |                9.00474 |                    -0.129876 |               12 | False                 |
| _20250203153511_viewcapture70.png |          0.381  |                2057 |                8.65338 |               15.1316  |                     6.4782   |                6 | False                 |
| _20250203153538_viewcapture71.png |          0.4918 |                1307 |                2.0658  |                2.94985 |                     0.884053 |                6 | False                 |
| _20250204172838_viewcapture62.png |          0.566  |                3418 |               12.639   |               14.8173  |                     2.17835  |                8 | False                 |
| _20250204191949_viewcapture61.png |          0.4425 |                 766 |               12.7937  |               24       |                    11.2063   |                5 | True                  |
| _20250205110933_viewcapture55.png |          0.6897 |                 947 |               53.7487  |               74.4939  |                    20.7452   |                8 | False                 |
| _20250205110946_viewcapture59.png |          0.8929 |                6461 |               16.6538  |               23.3333  |                     6.67956  |               23 | True                  |
| _20250205110946_viewcapture60.png |          0.4505 |                2456 |               14.5765  |               17.6086  |                     3.03202  |                6 | False                 |

### Positive nuclei are clustered, and the scale of clustering carries the signal

Ki67-positive nuclei were spatially clustered rather than randomly distributed in 51 of 53 images, with a median Donnelly-corrected Clark-Evans index of 0.686 against 1 for a random pattern, and a median quadrat variance to mean ratio of 6.63 against 1 for a Poisson pattern (Figure 4).

The spatial statistics differ sharply in whether they explain the scoring gap. The quadrat variance to mean ratio, computed over windows comparable in size to the reporting field, correlated strongly with the hotspot minus average difference (Spearman rho 0.66, 95 percent confidence interval 0.47 to 0.79, q = 4.1e-07). The kernel density hotspot coefficient of variation and the border-corrected Ripley L correlated weakly (rho 0.31 and 0.30 respectively, q = 0.044 for both). The Clark-Evans index, which measures nearest-neighbour spacing at single-cell distances, showed no association (rho -0.03, q = 0.84).

The discordance is therefore generated by large-scale patchiness in the distribution of proliferating cells, not by whether positive nuclei lie adjacent to one another. A spatial statistic evaluated at the wrong scale is silent about the problem even when it correctly reports that clustering exists.

![Figure 4](D:/UCAS project/new analysis/coda-brca-my/coda-brca-my/figures/F21_ki67_spatial.png)

**Figure 4. Spatial arrangement of Ki67-positive nuclei.**

### Stages that could not be run

The serial-section stages of the framework were not run, and could not be. The benchmark serial dataset is not publicly downloadable, the public repository containing only the evaluation software, and the breast whole-slide resource with matched markers and registration landmarks [3] requires a data use agreement. The present material consists of single fields and cannot support registration [4], reconstruction or volumetric quantification at any sample size. Supplementary Table S14 records every stage against every arm with the reason for each block.

**Supplementary Table S14** (24 rows)

| arm         |   stage | stage_name                      | status   | reason                                              |
|:------------|--------:|:--------------------------------|:---------|:----------------------------------------------------|
| A Kartasalo |       1 | nonlinear registration          | blocked  | dataset not acquired; images require author request |
| A Kartasalo |       2 | registration QC                 | blocked  | dataset not acquired; images require author request |
| A Kartasalo |       3 | cell detection                  | blocked  | dataset not acquired; images require author request |
| A Kartasalo |       4 | semantic segmentation           | blocked  | dataset not acquired; images require author request |
| A Kartasalo |       5 | 3D reconstruction               | blocked  | dataset not acquired; images require author request |
| A Kartasalo |       6 | quantification and connectivity | blocked  | dataset not acquired; images require author request |
| A Kartasalo |       7 | fiber alignment                 | blocked  | dataset not acquired; images require author request |
| B ACROBAT   |       1 | nonlinear registration          | blocked  | dataset not acquired; data use agreement            |
| B ACROBAT   |       2 | registration QC                 | blocked  | dataset not acquired; data use agreement            |
| B ACROBAT   |       3 | cell detection                  | blocked  | dataset not acquired; data use agreement            |
| B ACROBAT   |       4 | semantic segmentation           | blocked  | dataset not acquired; data use agreement            |
| B ACROBAT   |       5 | 3D reconstruction               | blocked  | sections are not consecutive; no volume             |

## Discussion

Ki67 scoring irreproducibility is usually framed as observer variability, and the remedies proposed are training, standardised protocols and automated counting [11]. Our results indicate that a substantial part of the disagreement is a property of the tissue rather than of the observer. When proliferation is spatially patchy, a hotspot convention and an average convention are measuring different things, and both are correct measurements of different quantities.

The scale-dependence is the practically important finding. Clark-Evans correctly reported clustering in almost every image, yet carried no information about the scoring gap, because nearest-neighbour spacing operates at single-cell distances. The quadrat variance to mean ratio, evaluated over windows comparable to the reporting field, explained the gap well. Any attempt to use spatial statistics to flag unreliable Ki67 cases must therefore match the statistic to the scale at which the score is formed, and reporting that positive cells are clustered is not by itself useful.

This has a direct clinical reading. Cases with high coarse-scale patchiness are those where the treatment decision is most likely to depend on where the pathologist looks. Such cases could be flagged for a defined scoring protocol or for a second reader, which is a more targeted intervention than applying the same standardisation everywhere.

What did not reproduce is the three-dimensional component. CODA's overcounting result [1] is the strongest argument that single-section counting misrepresents tissue, and it requires serial sections that do not exist for this material and are not publicly available for the benchmark tissue. We report that as a gap rather than substituting a weaker analysis, because the substitutes available would produce numbers without the property that makes the original result meaningful.

The HER2 handling deserves a note. Per-nucleus scoring of a membranous marker produces a confident and meaningless number, which is more dangerous than an obvious error, and the implementation refuses the operation [13]. Membrane completeness is reported instead as a quantitative descriptor and deliberately not mapped onto the ASCO/CAP categories, which are defined by a scoring procedure this measurement does not reproduce.

## Limitations

**Field of view rather than whole slide.** The images are microscope captures of roughly 350 to 1150 microns across. No whole-slide inference follows, and Ripley's K beyond about a quarter of the field width is unreliable even with border correction, which is why radii were capped and the limit recorded.

**Possible non-random field selection.** Fields photographed to document staining are plausibly biased toward positive areas. If selection was by eye the sample does not represent the slide, and the selection procedure is not recorded in the image metadata.

**A wide magnification range.** The series spans 38-fold in microns per pixel. Measurements are converted to microns before pooling, but scale-dependent texture features are not comparable across this range, and 6 images were too coarse for nuclear analysis.

**No denominator on 150 of 234 images.** Percent positive is not reportable for those and was not estimated indirectly. This constrains ER and PR far more than Ki67.

**No three-dimensional analysis anywhere in this work, and none on breast.** The only serial material contemplated is mouse prostate and liver, which even if obtained would validate the pipeline rather than establish a breast finding.

**Single institution, single scanner, no comparison cohort.** No cross-cohort comparison was attempted and no batch sensitivity audit is reported. Any future comparison against public cohorts must run that audit first, because scanner and protocol differences otherwise masquerade as population differences.

**Ki67 detection was not validated against manual counts** in this material. The framework specifies validation at a 2 micron matching tolerance against two annotators [1], which requires annotation effort not yet performed. The paired comparison is internally controlled, since both scores derive from the same detections, but absolute index values should be read with that in mind.

## Conclusion

The disagreement between hotspot and average Ki67 scoring is measurable, large enough to change the treatment category in roughly a quarter of evaluable images, and explained by coarse-scale spatial patchiness of proliferating cells. Spatial statistics computed at the scale of the reporting window identify the cases at risk; the same statistics computed at single-cell scale do not. Quantifying arrangement, not only proportion, is a tractable addition to Ki67 reporting.

## Data and code availability

Analysis code, the parameter configuration with its SHA-256 locked block, and the test suite are available in the project repository. Per-image quality control, marker and spatial results are provided as supplementary tables. The image material is institutional and is not publicly redistributable. Public datasets referenced but not obtained are the serial benchmark stacks [2] and the ACROBAT breast cohort [3]. Software used includes OpenSlide [7] and QuPath [6] for related handling, and colour deconvolution follows [8].

## References

1. Kiemen AL, Braxton AM, Grahn MP, Han KS, Babu JM, Reichel R, et al. CODA: quantitative 3D reconstruction of large tissues at cellular resolution. Nat Methods. 2022;19:1490-1499. PMID: 36280719. doi:10.1038/s41592-022-01650-9

2. Kartasalo K, Latonen L, Vihinen J, Visakorpi T, Nykter M, Ruusuvuori P. Comparative analysis of tissue reconstruction algorithms for 3D histology. Bioinformatics. 2018;34:3013-3021. PMID: 29684099. doi:10.1093/bioinformatics/bty210

3. Weitz P, Valkonen M, Solorzano L, Carr C, Kartasalo K, Boissin C, et al. The ACROBAT 2022 challenge: Automatic registration of breast cancer tissue. Med Image Anal. 2024;97:103257. PMID: 38981282. doi:10.1016/j.media.2024.103257

4. Borovec J, Kybic J, Arganda-Carreras I, Sorokin DV, Bueno G, Khvostikov AV, et al. ANHIR: Automatic Non-Rigid Histological Image Registration Challenge. IEEE Trans Med Imaging. 2020;39:3042-3052. PMID: 32275587. doi:10.1109/TMI.2020.2986331

5. Graham S, Vu QD, Raza SEA, Azam A, Tsang YW, Kwak JT, et al. Hover-Net: Simultaneous segmentation and classification of nuclei in multi-tissue histology images. Med Image Anal. 2019;58:101563. PMID: 31561183. doi:10.1016/j.media.2019.101563

6. Bankhead P, Loughrey MB, Fernández JA, Dombrowski Y, McArt DG, Dunne PD, et al. QuPath: Open source software for digital pathology image analysis. Sci Rep. 2017;7:16878. PMID: 29203879. doi:10.1038/s41598-017-17204-5

7. Goode A, Gilbert B, Harkes J, Jukic D, Satyanarayanan M. OpenSlide: A vendor-neutral software foundation for digital pathology. J Pathol Inform. 2013;4:27. PMID: 24244884. doi:10.4103/2153-3539.119005

8. Ruifrok AC, Johnston DA. Quantification of histochemical staining by color deconvolution. Anal Quant Cytol Histol. 2001;23:291-9. PMID: 11531144.

9. Dowsett M, Nielsen TO, A'Hern R, Bartlett J, Coombes RC, Cuzick J, et al. Assessment of Ki67 in breast cancer: recommendations from the International Ki67 in Breast Cancer working group. J Natl Cancer Inst. 2011;103:1656-64. PMID: 21960707. doi:10.1093/jnci/djr393

10. Polley MY, Leung SC, Gao D, Mastropasqua MG, Zabaglo LA, Bartlett JM, et al. An international study to increase concordance in Ki67 scoring. Mod Pathol. 2015;28:778-86. PMID: 25698062. doi:10.1038/modpathol.2015.38

11. Nielsen TO, Leung SCY, Rimm DL, Dodson A, Acs B, Badve S, et al. Assessment of Ki67 in Breast Cancer: Updated Recommendations From the International Ki67 in Breast Cancer Working Group. J Natl Cancer Inst. 2021;113:808-819. PMID: 33369635. doi:10.1093/jnci/djaa201

12. Wolff AC, Hammond MEH, Allison KH, Harvey BE, Mangu PB, Bartlett JMS, et al. Human Epidermal Growth Factor Receptor 2 Testing in Breast Cancer: American Society of Clinical Oncology/College of American Pathologists Clinical Practice Guideline Focused Update. J Clin Oncol. 2018;36:2105-2122. PMID: 29846122. doi:10.1200/JCO.2018.77.8738

13. Wolff AC, Somerfield MR, Dowsett M, Hammond MEH, Hayes DF, McShane LM, et al. Human Epidermal Growth Factor Receptor 2 Testing in Breast Cancer: ASCO-College of American Pathologists Guideline Update. J Clin Oncol. 2023;41:3867-3872. PMID: 37284804. doi:10.1200/JCO.22.02864

14. Conklin MW, Eickhoff JC, Riching KM, Pehlke CA, Eliceiri KW, Provenzano PP, et al. Aligned collagen is a prognostic signature for survival in human breast carcinoma. Am J Pathol. 2011;178:1221-32. PMID: 21356373. doi:10.1016/j.ajpath.2010.11.076

15. Provenzano PP, Eliceiri KW, Campbell JM, Inman DR, White JG, Keely PJ. Collagen reorganization at the tumor-stromal interface facilitates local invasion. BMC Med. 2006;4:38. PMID: 17190588. doi:10.1186/1741-7015-4-38

16. Encoder-decoder with atrous separable convolution for semantic image segmentation. ECCV 2018, arXiv:1802.02611. **PMID: not found** (not indexed in PubMed).

17. Deep residual learning for image recognition. CVPR 2016, arXiv:1512.03385. **PMID: not found** (not indexed in PubMed).
