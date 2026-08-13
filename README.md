# Transmission acceleration outperforms outbreak threshold for dengue outbreak detection

The repository reproduces every figure, supplementary figure, table, and CSV output of the manuscript from a single weekly dengue case dataset and five sequential R scripts.

## Background

Routine dengue surveillance in more than 130 countries triggers an outbreak alert when weekly case counts exceed the historical mean plus two standard deviations. This endemic-channel rule was originally devised for malaria control under stable transmission and has not been independently validated as a signal for anticipatory action in hyperendemic dengue. We evaluated 14 outbreak-detection methods drawn a priori from four paradigms: surveillance-guideline percentile thresholds, retrospective thresholds, contemporary surveillance algorithms, and transmission-acceleration measures. Performance was assessed against an anticipatory anchor framework that classifies every detector trigger by its temporal relation to the annual epidemic peak. The framework was derived in Quezon City, Philippines (2013–2024), and applied without parameter modification to all 17 Philippine administrative regions and to eight dengue-endemic countries (2016–2024).

The Constant Transmission Acceleration, a hysteretic short-term-average to long-term-average ratio adapted from seismological signal detection, generalised across all three scales as the operational leader on the five primary metrics: True-Alarm Magnitude, the year-mean number of true alarms, sensitivity, mean lead time, and warning persistence. The WHO-mandated retrospective Outbreak Threshold (the rolling week-specific donor-year mean plus two standard deviations) did not.

## Authors

Keanu John Pelitro<sup>1</sup>, Julia Fye Manzano<sup>1</sup>, Troy Owen Matavia<sup>1</sup>, Kylone Soriano<sup>1</sup>, Klara Bilbao<sup>1</sup>, Gereka Marie Garcia<sup>1</sup>, Aira Joy Delos Angeles<sup>1</sup>, Alfredo Mahar Lagmay<sup>1,3</sup>, and DJ Darwin Bandoy<sup>1,2,*</sup>

<sup>1</sup> University of the Philippines Resilience Institute, Quezon City, Philippines.
<sup>2</sup> College of Veterinary Medicine, University of the Philippines, Los Baños, Laguna, Philippines.
<sup>3</sup> National Institute of Geological Sciences, University of the Philippines Diliman, Quezon City, Philippines.

<sup>*</sup> Corresponding author: drbandoy@up.edu.ph

## Repository structure

```
dengue-transmission-acceleration/
├── README.md
├── LICENSE                                          MIT
├── CITATION.cff
├── run_all.R                                         reproduces the full pipeline in order
├── R/                                                shared configuration, theme, and detector code
│   ├── 00_config.R                                   paths, package loader, global seed, eta-threshold handoff
│   ├── 01_publication_theme.R                        shared ggplot2 publication theme
│   ├── 02_detection_framework.R                      shared detector functions (Stages 4–5)
│   └── install_dependencies.R                        one-off CRAN dependency installer
├── scripts/
│   ├── Stage1a_eta_threshold_derivation_QC.R         empirical η_ON, η_OFF derivation
│   ├── Stage1b_anchor_focused_window_selection_QC.R  STA/LTA window-selection robustness check
│   ├── Stage2_outbreak_threshold_drift_QC.R          Figure 1
│   ├── Stage3_QC_analysis.R                          Figures 2 and 3, Tables 1–2
│   ├── Stage4_Regional_Analysis.R                    Figures 4 and 5
│   └── Stage5_Country_Analysis.R                     Figures 6 and 7
├── data/
│   ├── Dengue-Rainfall_Dataset.xlsx                  deidentified weekly data (place here; see Data)
│   └── geometry/                                     optional local shapefile/GeoPackage for the Stage 4 map
└── outputs/                                          created automatically on first run
```

The scripts are designed to run in numerical order via `run_all.R`, which locates the project root, loads `R/00_config.R`, and sources each script in `scripts/` in turn. Stage 1a derives the two hysteresis thresholds used by the Constant Transmission Acceleration detector and writes them to `outputs/Stage1_eta_thresholds_QC/eta_thresholds_derived.R`, which `load_eta_thresholds()` (in `R/00_config.R`) reads back for Stages 3, 4, and 5 — so re-running Stage 1a with different data propagates automatically to every downstream stage rather than requiring the thresholds to be hard-coded. If that file is absent, `load_eta_thresholds()` falls back to the published constants (η_ON = 1.33, η_OFF = 0.73) with a warning. Stage 1b is a self-contained script (base R + `readxl` only) that does not depend on `R/00_config.R` and does not feed any downstream stage. Stages 2 to 5 are mutually independent and can be run in any order once Stage 1a has completed.

> `run_all.R` must be sourced with the working directory set to the project root (the folder containing `R/` and `scripts/`). It stops immediately with a clear message if either is missing.

## Data

The complete weekly dataset is deposited as a single Excel workbook on Zenodo:

*Dengue-Rainfall_Dataset.xlsx.* Zenodo. https://doi.org/10.5281/zenodo.19448854

Download the workbook and place it at `data/Dengue-Rainfall_Dataset.xlsx` in the project root before running the pipeline; `R/00_config.R` looks for it at that exact path and stops with an explicit error if it is missing.

The workbook contains five sheets. The two metadata sheets (Dataset Summary and Data Dictionary) describe provenance and column definitions. The three analytical sheets each correspond to one geographic scale of the analysis.

| Sheet           | Records | Period    | Analytical columns used by the pipeline |
|-----------------|---------|-----------|------------------------------------------|
| `QC Data`       | 832     | 2010–2025 | `YR`, `WN`, `DC_QC`                      |
| `Regional Data` | 7 072   | 2016–2025 | `REGION`, `YR`, `WN`, `DC_DOH`           |
| `Country Data`  | 3 216   | 2016–2025 | `COUNTRY`, `YR`, `WN`, `DC_OPENDENGUE`   |

`YR` is the ISO calendar year, `WN` is the ISO epidemiological week number, and the `DC_*` columns are the weekly suspected dengue case counts at the corresponding scale. Each observation is indexed to the Monday of the corresponding ISO week. `DC_QC` is from the Quezon City Epidemiology and Surveillance Division. `DC_DOH` is from PIDSR-reported regional submissions archived at the Humanitarian Data Exchange (https://data.humdata.org). `DC_OPENDENGUE` is from the OpenDengue repository (https://opendengue.org).

The analytical sheets carry additional auxiliary columns retained for compatibility with related work but not consumed by the present pipeline: rainfall fields (`RF_NASA`, `RF_PAGASA`, `RF_HDX`) and data-quality flags (`FLAG_COVID`, `FLAG_PLAUSIBILITY`, `FLAG_SINGLE_CELL_RF`, `FLAG_DEKADAL_APPROX`, `FLAG_TERMINAL_GAP`). All auxiliary columns are documented in full in the `Data Dictionary` sheet of the workbook.

Cases were classified as suspected dengue under the Philippine Integrated Disease Surveillance and Response (PIDSR) Programme, consistent with the WHO 2009 dengue case definition adopted by PIDSR in 2011. Laboratory confirmation was not required, consistent with passive surveillance practice throughout the study period. All surveillance data are aggregated weekly totals containing no individual-level identifiers.

## Reproducing the analysis

### Requirements

R 4.4.1 or newer under Linux, macOS, or Windows. Approximately 8 GB of RAM is sufficient. `R/00_config.R` sets the global seed (12345) at project-configuration time to guarantee reproducibility of every randomised step (year-cluster bootstraps, ROC bootstraps, and parametric bootstraps inside the η-derivation); each script additionally checks its own dependencies via `require_packages()` and stops with a clear message if any are missing, rather than installing them silently.

Run the one-off installer once per machine:

```r
source("R/install_dependencies.R")
```

**Required:** `readxl`, `dplyr`, `tidyr`, `purrr`, `tibble`, `stringr`, `readr`, `zoo`, `ISOweek`, `scales`, `rlang`, `pROC`, `MASS`, `boot`, `ggplot2`, `ggrepel`, `patchwork`, `cowplot`, `viridisLite`, `grid`, `sf`, `ggspatial`, `prettymapr`, `here`.

**Optional:** `geodata`, `rnaturalearth`, `rnaturalearthdata`, `terra` (Stage 4 map geometry — only needed if no local shapefile/GeoPackage is placed in `data/geometry/`); `irr` (an inter-rater-agreement diagnostic in Stage 3 that degrades gracefully if absent); `rstudioapi` (lets Stage 1b auto-detect its own file path when run interactively in RStudio; Stage 1b falls back to other detection methods if it is absent).

`sf` requires the GDAL, GEOS, and PROJ system libraries; see the comments at the foot of `R/install_dependencies.R` for platform-specific installation notes if it fails to build.

### Run order

```r
# From the project root (the folder containing R/, scripts/, and data/),
# with Dengue-Rainfall_Dataset.xlsx placed at data/Dengue-Rainfall_Dataset.xlsx:

source("run_all.R")
```

or, to run stages individually from the project root:

```r
source("scripts/Stage1a_eta_threshold_derivation_QC.R")        # produces eta_thresholds_derived.R
source("scripts/Stage1b_anchor_focused_window_selection_QC.R") # optional: window-selection robustness check
source("scripts/Stage2_outbreak_threshold_drift_QC.R")         # Figure 1
source("scripts/Stage3_QC_analysis.R")                         # Figures 2, 3 and tables
source("scripts/Stage4_Regional_Analysis.R")                    # Figures 4, 5 and tables
source("scripts/Stage5_Country_Analysis.R")                      # Figures 6, 7 and tables
```

Each script writes its outputs to a single Stage-specific subdirectory of `outputs/`. Subdirectories are created on first run if they do not already exist.

### Approximate runtime on a recent laptop

Stage 1a, four to six minutes (six independent derivation methods, the slowest of which performs B = 2 000 ROC bootstrap replicates). Stage 1b (optional), a few minutes (window-selection robustness check; B = 2 000 season-cluster BCa bootstrap replicates per candidate window). Stage 2, under one minute. Stage 3, six to ten minutes (year-cluster bootstrap with B = 1 000, fourteen detectors, six pairwise Wilcoxon contrasts on the year × detector matrix). Stage 4, eight to fifteen minutes (the same machinery on 17 regions). Stage 5, six to twelve minutes (the same machinery on 8 countries). Total wall-clock time for the complete pipeline is typically 25 to 50 minutes.

## Pipeline overview

### Stage 1a: Empirical η_ON and η_OFF derivation

The activation and deactivation thresholds for the Constant Transmission Acceleration detector were derived empirically from the Quezon City series rather than imported from the seismological literature. Six independent derivation methods are applied: (i) baseline-ratio asymmetric quantiles at the 90th and 50th percentiles; (ii) ROC-style operating points at 95% and 99% sensitivity for the separation of pre-peak from baseline weeks (B = 2 000 bootstrap replicates); (iii) calibration to an average run length to false alarm of 12 weeks under a baseline-only fraction of 30%, with secondary targets at 6 and 26 weeks; (iv) leave-one-year-out cross-validation under an asymmetric utility favouring early triggering; (v) negative-binomial theoretical thresholds at ±1σ; and (vi) parametric bootstrap quantiles at the 90th and 25th percentiles. The median across the six methods is adopted as the headline pair: η_ON = 1.33 for activation, η_OFF = 0.73 for deactivation. The persistent-coverage criterion freezes the long-term-average reference at activation until eight consecutive sub-threshold weeks have elapsed. Adopted thresholds are inherited unchanged by the regional and country stages.

Outputs (`Stage1_eta_thresholds_QC/`): a four-panel diagnostic figure, six per-method derivation tables, the final decision record, and a sourceable R file (`eta_thresholds_for_figure2.R`) used by every downstream stage.

### Stage 1b: STA/LTA window-selection robustness check

An independent check of the prespecified STA/LTA window widths against the two anticipatory anchors (Supplementary Table 21). For each of the Continuous TA (2/8, 3/12, 4/16) and Constant TA (2/13, 3/20, 4/26) families, the prespecified window is tested head-to-head against its scale-matched alternatives using exact paired season-level tests (randomization tests for operational metrics, sign-flip tests for AUC outcomes), each with 2 000 season-cluster BCa bootstrap confidence intervals. The comparison is repeated under the primary anchor specification and across the Actionable-Window and Epidemic-Burden sensitivity grids. This stage does not feed any downstream stage; it is a standalone robustness check on the window choice fixed in Stage 1a.

Outputs (`Stage1b_anchor_focused_window_selection_QC/`): 19 CSV tables spanning the primary and sensitivity head-to-head comparisons, plus a run-summary text file.

### Stage 2: Outbreak threshold drift in Quezon City

Weekly notified dengue case counts in Quezon City for 2013–2025 are overlaid with the rolling week-specific Mean + 2 SD outbreak threshold under two baselines: one that includes the 2020–2021 pandemic years and one that excludes them. The figure documents the structural reactivity of the endemic channel to recent epidemic years and motivates the framework that follows.

Outputs (`Stage2_outbreak_threshold_drift_QC/`): `fig1_dengue_thresholds_vertical.pdf`, `fig1_dengue_thresholds_vertical.png`, `fig1_thresholds_with_pandemic.csv`, `fig1_thresholds_without_pandemic.csv`.

### Stage 3: Quezon City analysis

The fourteen detection methods are evaluated against the anticipatory anchor framework on the Quezon City series for 10 evaluable seasons (2013–2019, 2022–2024; 2020 and 2021 excluded for pandemic surveillance disruption, 2025 excluded as out-of-distribution).

The framework defines two compartments. The Actionable Window (AW) is the four to eight weeks immediately preceding the annual peak: the interval in which vector-control mobilisation is operationally plausible. The Epidemic Burden block (EB) is the smallest contiguous block of weeks containing the annual peak whose cumulative case count reaches at least 70% of the year's total. A trigger fired at week *t* is classified as a true alarm if and only if *t* belongs to AW ∪ EB; all other triggers are false alarms.

Each detector is evaluated on eight headline metrics partitioned across three operational categories: alarm accuracy and epidemic-burden coverage (True-Alarm Magnitude, year-mean number of true alarms, positive predictive value, sensitivity); early-warning timeliness (mean lead time, warning persistence, actionable lead-time yield); and false-alarm volume (year-mean count of triggers landing outside AW ∪ EB). Five of these eight metrics drive the within-stratum dominance and consensus analyses: True-Alarm Magnitude, the year-mean number of true alarms, sensitivity, mean lead time, and warning persistence. Year-cluster bootstrap (B = 1 000) yields 95% confidence intervals around every metric, resampling whole evaluable years to preserve within-year correlation.

Outputs (`Stage3_QC_analysis/`): the Figure 2 dominance dashboard (Panel A dominance matrix; Panel B True-Alarm Magnitude versus mean lead time scatter), the Figure 3 multipanel for the Constant Transmission Acceleration, the Supplementary Figure 3 multipanel for the Continuous Transmission Acceleration, six head-to-head Wilcoxon comparison figures (Constant TA versus Outbreak Threshold; Continuous TA versus Outbreak Threshold; Constant TA versus Continuous TA; and Constant TA versus each of the three contemporary algorithms — Farrington, EWARS, EARS), Tables 1, 1b, 2, 2A, and 2B, and sensitivity tables reporting robustness to AW width, EB cumulative-case fraction, and year-inclusion rules. CSV exports of all figure-underlying data are written to `Stage3_QC_analysis/evaluation_framework/`.

### Stage 4: Regional analysis (17 Philippine regions)

The full fourteen-method framework, the AW/EB anchors, and the η values from Stage 1a are applied to each of BARMM, CAR, MIMAROPA, NCR, and Regions I to XIII (Region IV-B is folded into MIMAROPA per Philippine Statistics Authority convention). A region is retained where the time series carries at least 12 annual peak cases and at least five evaluable years; six to seven evaluable years are available per region.

Within each region the operational leader is identified through a bootstrap dominance probability: the proportion of B = 1 000 year-cluster bootstrap replicates in which a given detector achieves the highest within-region composite mean rank across the five primary metrics. The Figure 4/5 dashboard reports three display-level pairwise contrasts (Constant TA vs Outbreak Threshold, Continuous TA vs Outbreak Threshold, Constant TA vs Continuous TA). The primary regional consensus classification instead uses the complete four-comparator analysis — Constant TA versus the Outbreak Threshold, Farrington, EWARS, and EARS — under a strict cross-region Bonferroni correction (k = 68, i.e. 17 regions × 4 contrasts; corrected α ≈ 7.4 × 10⁻⁴). A region achieves strong consensus if Constant TA wins all four comparisons at the corrected threshold with no losses; partial consensus requires at least one significant win and no losses; lead-only indicates a point-estimate lead with no significant wins or losses; contested indicates mixed wins and losses, or no significant pair at the corrected threshold. Within-region-only Bonferroni (k = 4) is reported as a sensitivity check.

Outputs (`Stage4_Regional_Analysis/`): six per-metric multipanel figures (`Figure4_TAM`, `Figure4_N_True_Alarms`, `Figure4_Sensitivity`, `Figure4_Mean_Lead_Time`, `Figure4_Warning_Persistence`, `Figure4_Method_Summary`); the composite Figure 5 (`Figure5_Detector_Map`) carrying the choropleth detector map and the per-detector dot plot of regional dominance probabilities; and supporting CSVs covering framework metrics with confidence intervals, the dominance matrix, dominance probabilities, bootstrap replicates, the four-comparator Wilcoxon table, and the inclusion audit.

### Stage 5: Country analysis (8 dengue-endemic countries)

The same pipeline is applied to Brazil, Colombia, Mexico, Peru, the Philippines, Singapore, Sri Lanka, and Taiwan. As at the regional scale, the primary country-level consensus classification uses the four-comparator analysis (Constant TA versus the Outbreak Threshold, Farrington, EWARS, and EARS) under a strict cross-country Bonferroni correction (k = 32, i.e. 8 countries × 4 contrasts; corrected α ≈ 1.6 × 10⁻³). Within-country-only Bonferroni (k = 4) is reported as a sensitivity check. The auxiliary cross-country paired Wilcoxon on the per-country dominance probability column is power-limited at *n* = 8 and is reported as a sanity check rather than as a primary inferential output; the within-stratum head-to-head consensus framework is the load-bearing inferential machinery at every scale.

Outputs (`Stage5_Country_Analysis/`): six per-metric multipanel figures (`Figure6_TAM`, `Figure6_N_True_Alarms`, `Figure6_Sensitivity`, `Figure6_Mean_Lead_Time`, `Figure6_Warning_Persistence`, `Figure6_Method_Summary`); the composite two-panel `Figure7_Country_Summary` carrying the per-country dominance dot plot; and supporting CSVs in the same format as the regional stage.

## Methodological notes

### The fourteen detectors

The detectors are grouped a priori into four comparator paradigms.

The **surveillance-guideline percentile thresholds** paradigm (2 detectors) comprises the WHO 75th Percentile Threshold and the WHO 90th Percentile Threshold; both fire when weekly cases exceed the rolling week-specific donor-year case-distribution percentile.

The **retrospective thresholds** paradigm (3 detectors) comprises the Outbreak Threshold (rolling week-specific donor-year mean plus two standard deviations; the WHO-mandated comparator), the Alarm Threshold (the same construction at one standard deviation), and a one-sided walk-forward Cumulative Sum Control chart with k = 0.5σ and h = 5σ, reset at calendar-year boundaries.

The **contemporary surveillance algorithms** paradigm (3 detectors) comprises the Farrington algorithm (quasi-Poisson regression with reweighting to reduce the influence of past outbreaks), the Early Aberration Reporting System (EARS, cumulative sums over short baseline windows), and the WHO/PIDSR Early Warning, Alert and Response System (EWARS, logistic regression for dengue-specific alarm classification).

The **transmission-acceleration** paradigm (6 detectors) comprises the Continuous Transmission Acceleration (a 3-week to 12-week moving-average ratio with no hysteresis); the Constant Transmission Acceleration (a 4-week to 26-week mean ratio with a 2-week guard band, hysteretic activation at η_ON = 1.33 and deactivation at η_OFF = 0.73, the long-term reference frozen at activation and unfrozen only after eight consecutive sub-threshold weeks); the Incidence Gradient (the first difference of the 3-week rolling mean against its rolling 80th percentile); the Critical Transition Indicator (an 8-week to 26-week rolling-variance ratio against its rolling 80th percentile); the Hydrological Inflection Measure (the joint condition that weekly cases exceed the rolling 75th percentile and that the first difference exceeds the rolling 75th percentile of the donor-year rise rate); and the Composite Outbreak Signal (the conjunction of an Outbreak Threshold breach with a Constant Transmission Acceleration on-state).

### Donor-pool construction

All percentile-, standard-deviation-, and historical-baseline-based detectors use a walk-forward donor pool. For evaluable target year *Y*, donors comprise the historical evaluable years preceding *Y* subject to the 2020, 2021, and 2025 exclusions, with each week's baseline statistics computed over donor-year observations of the same calendar week. This precludes look-ahead bias, so every reported metric is honestly prospective.

### Statistical analysis

Each metric is reported with its bootstrap confidence interval resampling whole evaluable years to preserve within-year correlation (B = 2 000 BCa replicates at city scale; B = 1 000 percentile-interval replicates at regional and country scale). Four prespecified pairwise Wilcoxon signed-rank contrasts anchor the primary significance classification at every scale: Constant TA versus the Outbreak Threshold, Farrington, EWARS, and EARS. At the local scale (n = 10 paired years), each contrast is assessed at α = 0.05 without across-pair correction since each tests a distinct hypothesis. At the regional scale, the primary correction applies a strict cross-region Bonferroni factor across all 68 region × contrast cells (4 × 17). At the country scale, the factor is k = 32 (4 × 8). Within-stratum-only Bonferroni (k = 4) is reported at every scale as a sensitivity check. Stage 3 additionally reports two supplementary contrasts (Continuous TA versus the Outbreak Threshold, and Constant TA versus Continuous TA) for the acceleration-family comparison shown in Figure 3.

### Reproducibility, ethics, and reporting

All framework parameters were specified a priori before the head-to-head analyses were run. Sensitivity tests around each (the AW width, the EB cumulative-case fraction, the year-inclusion rules, and the η thresholds) are reported in the supplementary tables. The study is reported following the Standards for Reporting Diagnostic Accuracy Studies (STARD 2015), adapted for binary alarm-algorithm evaluation on aggregated time-series data; STARD was preferred over TRIPOD because the analysis evaluates binary detection algorithms against an externally specified temporal reference (the anticipatory anchor framework) rather than fitting an individual-level prediction model. Ethical approval was granted by the University of the Philippines Research Ethics Committee (Protocol No. 2024-0004-F-FMDS).

## Citation

If you use this code or the dataset, please cite both the manuscript and the dataset deposition:

Pelitro KJ, Manzano JF, Matavia TO, Soriano K, Bilbao K, Garcia GM, Delos Angeles AJ, Lagmay AM, Bandoy DJD. *Transmission acceleration outperforms outbreak threshold for dengue outbreak detection.* DOI to be assigned upon acceptance.

Pelitro KJ, et al. *Dengue-Rainfall_Dataset* (Version 1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.19448854

A `CITATION.cff` file is provided so that the GitHub repository renders a "Cite this repository" widget once the manuscript DOI is assigned.

## Licence

The code in this repository is released under the MIT Licence (see `LICENSE`).

The data are distributed under the Open Data Commons Open Database License (ODC-ODbL) v1.0 (https://opendatacommons.org/licenses/odbl/1.0/), which permits the use, distribution, and adaptation of data, provided that appropriate attribution is given to the UP Resilience Institute–NOAH (UPRI-NOAH) and its contributors and that derivative databases are shared under the same license.

## Funding

Research reported in this publication was supported by the National Institute of Environmental Health Sciences of the National Institutes of Health under Award Number P20ES036118. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

## Acknowledgements

We thank the Quezon City Epidemiology and Surveillance Division for the city-level dengue surveillance data, the Philippine Atmospheric, Geophysical, and Astronomical Services Administration for the meteorological data archived alongside the case series, and the Humanitarian Data Exchange platform for the regional dengue dataset. The country-level case counts were drawn from the OpenDengue repository.

## Contact

For questions about the code or to report a reproducibility issue, please open a GitHub issue. For questions about the underlying surveillance data or the manuscript, contact Keanu John Pelitro at kapelitro@up.edu.ph.
