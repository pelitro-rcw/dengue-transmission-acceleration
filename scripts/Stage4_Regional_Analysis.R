# ==============================================================================
# STAGE 4 — REGIONAL ANALYSIS
# 17 Philippine regions; Figure 4 and 5, plus tables
# ------------------------------------------------------------------------------
# Regional generalisability of dengue outbreak detection methods across the
# 17 administrative regions of the Philippines (2016-2024, excluding the
# 2020-2021 pandemic disruption and the truncated 2025 surveillance year).
#
# OUTPUTS
# -------
# Figure 4 (six multipanel figures, one PDF/PNG per file):
#   Per-metric regional dominance for the five reported operational metrics
#   plus a cross-metric Method Summary figure:
#
#     Figure4_TAM
#     Figure4_N_True_Alarms
#     Figure4_Sensitivity
#     Figure4_Mean_Lead_Time
#     Figure4_Warning_Persistence
#     Figure4_Method_Summary
#
#   Each per-metric figure has three sub-panels:
#     a. Regional dominance matrix  - 11 methods (rows) by 17 regions
#                                    (columns); cell shading by within-region
#                                    min-max-normalised dominance score; right
#                                    column reports the per-method count of
#                                    regions swept at score >= 0.75.
#     b. Per-detector dot plot      - one dot per region per detector;
#                                    bootstrap 95% CI bars; group-mean
#                                    crossbar.
#     c. Wilcoxon detector-paired   - three pairwise contrasts among
#                                    Constant TA, Continuous TA, Outbreak
#                                    Threshold; paired by region (n = 17).
#
# Figure 5 (single composite figure):
#   a. Choropleth detector map of the Philippines (consensus winner per
#      region, four-tier classification scheme).
#   b. Per-region metric table with embedded per-metric significance.
#   c. Per-detector dot plot of regional dominance probabilities.
#
# Tables (CSV):
#   Stage4_Regional_Framework_Metrics.csv
#   Stage4_Regional_Framework_Metrics_with_CIs.csv
#   Stage4_Regional_8Metric_Summary.csv
#   Stage4_Regional_8Metric_Summary_with_CIs.csv
#   Stage4_Regional_Dominance_Matrix.csv
#   Stage4_Regional_Dominance_Probabilities.csv
#   Stage4_Regional_Wilcoxon_PerMetric.csv
#   Stage4_Regional_Wilcoxon_ConstantTA_vs_Comparators.csv  # OT, Farrington, EWARS, EARS
#   Stage4_Method_Summary_Long.csv
#   Stage4_Method_Summary_Aggregate.csv
#   Stage4_Regional_Bootstrap_Replicates.csv
#   Stage4_Detector_Map_RegionTable.csv
#   Stage4_Detector_Map_PerMetricSignificance.csv
#   Stage4_Detector_Map_CrossRegion_Wilcoxon.csv
#   Stage4_Detector_Map_PanelC_DotPlot.csv
#   Stage4_Detector_Map_PanelD_HeatmapData.csv
#   Stage4_Detector_Map_JoinAudit.csv
#   Stage4_Figure4_Legend.txt
#   Stage4_Figure5_Legend.txt
#
# FRAMEWORK
# ---------
# Two-anchor true-alarm rule:
#   T = A1 union A2
#     A1 (Actionable Window):  4-8 weeks pre-peak.
#     A2 (Epidemic Burden):    contiguous weeks centred on the seasonal
#                              peak that cumulate to 70% of seasonal cases.
#
# Two-bucket compartment scheme on True alarms:
#   Actionable: 4 <= lead <= 8 (within A1).
#   Reactive:   any other lead value on a True alarm.
#
# Same-denominator early warning timeliness:
#   Mean_Lead_Time and WP per (region, detector) divide by the count of
#   evaluable years; years that do not produce a qualifying trigger
#   contribute zero rather than NA.
#
# A1-restricted Sensitivity:
#   years_with_A1_true / years_evaluable.
#
# Year-cluster bootstrap (Cameron, Gelbach & Miller 2008):
#   B = 1000 replicates per region; year is the cluster unit because
#   within-year weekly observations are not independent.
#
# Regional inclusion criteria:
#   Annual peak >= 30 cases per (region, year);
#   >= 5 evaluable years per region.
#
# REQUIREMENTS
# ------------
# R packages (auto-installed if missing):
#   readxl, dplyr, tidyr, purrr, ggplot2, zoo, ISOweek, scales, tibble,
#   grid, patchwork, ggrepel, rlang, sf, readr, stringr.
# Optional geometry sources (one of):
#   geodata, rnaturalearth (with rnaturalearthdata).
#
# Input: regional weekly dengue dataset specified by `INPUT_PATH`,
# sheet `Regional Data`, columns REGION, YR, WN, DC_DOH.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------------------------
SCRIPT_TITLE   <- "STAGE 4 - Regional analysis (17 Philippine regions; Figure 4 and 5, plus tables)"
SCRIPT_VERSION <- "1.0 (final)"

cat("\n=============================================================\n")
cat(SCRIPT_TITLE, "\n", sep = "")
cat("Version: ", SCRIPT_VERSION, "\n", sep = "")
cat("=============================================================\n\n")

REQUIRED_PACKAGES <- c(
  "readxl", "dplyr", "tidyr", "purrr", "ggplot2", "zoo", "ISOweek",
  "scales", "tibble", "grid", "patchwork", "ggrepel", "rlang",
  "sf", "readr", "stringr", "ggspatial", "prettymapr", "viridisLite"
)
# Geometry sources for the Figure 5 map. Only needed when no local boundary
# file is present in data/geometry/.
OPTIONAL_GEOM_PACKAGES <- c("geodata", "rnaturalearth", "rnaturalearthdata")

# --- Project bootstrap -------------------------------------------------------
# Portable paths + shared publication theme. Replaces the previous inline
# install.packages() loops and the hard-coded Desktop path.
.bootstrap <- function() {
  root <- Sys.getenv("TA_PROJECT_ROOT", unset = "")
  if (!nzchar(root) || !dir.exists(root)) {
    d <- normalizePath(getwd(), winslash = "/")
    for (i in seq_len(6)) {
      if (dir.exists(file.path(d, "R")) && dir.exists(file.path(d, "scripts"))) {
        root <- d; break
      }
      pp <- dirname(d); if (identical(pp, d)) break; d <- pp
    }
  }
  if (!nzchar(root)) {
    stop("Set the working directory to the project root before sourcing.",
         call. = FALSE)
  }
  source(file.path(root, "R", "00_config.R"))
}
.bootstrap()

require_packages(REQUIRED_PACKAGES, purpose = "Stage 4")
invisible(lapply(REQUIRED_PACKAGES, function(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))))

source(file.path(DIR_R, "01_publication_theme.R"), local = TRUE)

# Shared detection framework: build_ears(), build_farrington(), build_ewars(),
# the anchor/compartment helpers and hh_paired_wilcoxon().
# local = TRUE is REQUIRED -- run_all.R runs each stage with
# sys.source(envir = env); without it these functions land in globalenv and
# cannot see this stage's own constants.
source(file.path(DIR_R, "02_detection_framework.R"), local = TRUE)

# The map labels use ggrepel's text halo (bg.colour / bg.r), added in 0.9.0.
# Checked up front so an outdated install fails immediately with a clear
# message, rather than partway through building Figure 5.
if (utils::packageVersion("ggrepel") < "0.9.0") {
  stop("ggrepel >= 0.9.0 is required (map label halos use bg.colour). ",
       "Installed: ", as.character(utils::packageVersion("ggrepel")),
       "\nRun: source(\"R/install_dependencies.R\")", call. = FALSE)
}

set.seed(GLOBAL_SEED)
options(scipen = 999)

# ------------------------------------------------------------------------------
# 1. USER PARAMETERS
# ------------------------------------------------------------------------------
# Path to the input regional dengue dataset. Edit as appropriate, or use
# `file.choose()` interactively.
INPUT_PATH <- DATA_FILE
SHEET_NAME <- SHEET_REGIONAL

OUTPUT_DIR <- file.path(DIR_OUTPUT, "Stage4_Regional_analysis")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# PH_SHAPEFILE is resolved in R/00_config.R: any .shp/.gpkg/.geojson placed in
# data/geometry/ is picked up automatically, which removes the last network
# dependency and pins the administrative boundary vintage.

base_family_global <- PUB_FAMILY

# Delegates to save_pub(): vector PDF + 600 dpi PNG, dimensions clamped to the
# journal print area.
save_plot_file <- function(stem, plot_obj, width, height, dpi = 600,
                           max_height = NC_H_SUPP) {
  save_pub(stem = sub("\\.(pdf|png)$", "", stem), plot = plot_obj,
           width = width, height = height, dir = OUTPUT_DIR, dpi = dpi,
           max_height = max_height)
}

# ------------------------------------------------------------------------------
# 2. FRAMEWORK CONSTANTS
# ------------------------------------------------------------------------------
# Constant TA (Vaezi-style hysteresis)
# Thresholds are read from Stage 1's derivation rather than hard-coded, so
# re-running Stage 1 propagates here. If Stage 1 output is absent,
# load_eta_thresholds() falls back to the published constants (1.33 / 0.73)
# with a warning, preserving the manuscript's numerical results exactly.
.eta    <- load_eta_thresholds()
ETA_ON  <- .eta$ETA_ON
ETA_OFF <- .eta$ETA_OFF
STA_WIN        <- 4L
LTA_WIN        <- 26L
GUARD          <- 2L

# Continuous TA (3-week / 12-week ratio); same ON threshold as Constant TA
ETA_ON_CLASSIC <- ETA_ON

# Variance-ratio guard (prevents division by tiny variances).
VAR26_GUARD    <- 0.5

# Vaezi STA/LTA OFF-state reset (consecutive-off weeks forcing LTA recompute).
MIN_OFF_RESET  <- 8L

# Anchor parameters (two-anchor framework: A1 and A2; no A3).
A1_LEAD_MIN    <- 4L
A1_LEAD_MAX    <- 8L
# A2_BURDEN_FRAC (0.70) and A2_BURDEN_PCT ("70%") are defined in R/00_config.R.
# Fixed by design: this stage reads no file and depends on no other stage for it.

# Lead-time compartment thresholds.
COMPARTMENT_ACTIONABLE_MIN <- 4L  # lower bound of A1
COMPARTMENT_ACTIONABLE_MAX <- 8L  # upper bound of A1

# Year inclusion
EXCLUDED_YEARS  <- c(2020L, 2021L, 2025L)
TARGET_YEARS    <- 2016L:2025L
EVALUABLE_YEARS <- setdiff(TARGET_YEARS, EXCLUDED_YEARS)  # n = 7

# Regional inclusion criteria
MIN_PEAK_CASES_PER_YEAR        <- 30L
MIN_EVALUABLE_YEARS_PER_REGION <- 5L

# Bootstrap configuration
BOOT_N_CI <- 1000L
# N_BOOTS is an alias used throughout the pairwise-test and dominance-summary
# helpers below; defined here (immediately after BOOT_N_CI) rather than near
# its heaviest use further down, since it is referenced as an eager default
# argument value (B = N_BOOTS) well before that point in the script.
N_BOOTS   <- BOOT_N_CI

# Canonical list of the 17 Philippine administrative regions. Defined here
# (moved up from its original position in the "REGION CANONICALISER" section
# further down) because Supplementary Table 13's construction, earlier in
# this script, references it directly.
CANONICAL_17 <- c("BARMM", "CAR", "MIMAROPA", "NCR",
                  "REGION I", "REGION II", "REGION III", "REGION IV-A",
                  "REGION V", "REGION VI", "REGION VII", "REGION VIII",
                  "REGION IX", "REGION X", "REGION XI", "REGION XII",
                  "REGION XIII")

# Significance and Bonferroni configuration
SIG_LEVEL          <- 0.05
HH_ALPHA           <- 0.05
BONF_LEVEL1_FACTOR <- 5L  # 5 metrics within (region x ordered pair)
BONF_LEVEL2_FACTOR <- 3L  # 3 pairs within metric

# Within-region dominance score threshold
DOMINANCE_THRESHOLD <- 0.75
CONFIDENCE_HI       <- 0.90

# ------------------------------------------------------------------------------
# 3. COLOR SYSTEM AND THEMES
# ------------------------------------------------------------------------------
TYPE_COLORS <- c(
  "National Standard"        = "#F8766D",
  "Retrospective Thresholds" = "#00BA38",
  "Acceleration Measures"    = "#619CFF"
)
type_order <- c(
  "National Standard",
  "Retrospective Thresholds",
  "Acceleration Measures"
)

# Built on the shared publication theme so typography matches every other
# stage; only the dashboard legend placement is overridden.
theme_dashboard <- function(base_size = PUB_BASE, base_family = PUB_FAMILY) {
  theme_pub(base_size, base_family) %+replace%
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "horizontal"
    )
}

# ------------------------------------------------------------------------------
# 4. DATA LOADING
# ------------------------------------------------------------------------------
if (!file.exists(INPUT_PATH)) {
  stop(paste0("Input file not found at:\n", INPUT_PATH))
}

df_raw <- readxl::read_excel(INPUT_PATH, sheet = SHEET_NAME)

required_cols <- c("REGION", "YR", "WN", "DC_DOH")
missing_cols  <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns in sheet '", SHEET_NAME, "': ",
              paste(missing_cols, collapse = ", ")))
}

df_all <- df_raw %>%
  dplyr::mutate(
    REGION = as.character(REGION),
    YR     = suppressWarnings(as.integer(YR)),
    WN     = suppressWarnings(as.integer(WN)),
    DC_DOH = suppressWarnings(as.numeric(DC_DOH))
  ) %>%
  dplyr::filter(!is.na(REGION), !is.na(YR), !is.na(WN),
                WN >= 1, WN <= 53) %>%
  dplyr::mutate(
    ISOweek = sprintf("%d-W%02d", YR, WN),
    Date    = ISOweek::ISOweek2date(paste0(ISOweek, "-1"))
  ) %>%
  dplyr::filter(!is.na(Date)) %>%
  dplyr::arrange(REGION, Date) %>%
  dplyr::filter(YR %in% TARGET_YEARS)

if (nrow(df_all) == 0) {
  stop("No rows remain after filtering to target years.")
}

# Anchor functions reference a column named DC_QC. Add an alias so the
# regional pipeline can use them verbatim.
df_all <- df_all %>% dplyr::mutate(DC_QC = DC_DOH)

# ------------------------------------------------------------------------------
# 5. SURGE GENERATION (PER REGION)
# ------------------------------------------------------------------------------
# Per-region detector pipeline:
#   - Week-specific baselines (mean, sd, p75, p90, alarm = mean+1sd,
#     outbreak = mean+2sd) computed from up to 5 prior donor years.
#   - Rolling P75 / P80 thresholds for dc_diff and ratio_var, donor-year based.
#   - VAR26_GUARD prevents the variance-ratio from blowing up at tiny variances.
#   - Walk-forward CUSUM with year-boundary state reset and rolling baseline
#     mean/sd as the location/scale parameters.
#   - Vaezi STA/LTA hysteresis with year-boundary reset, MIN_OFF_RESET
#     consecutive-off weeks forcing fresh LTA recompute, and >= for ON
#     activation.
#   - Composite = mean+2sd AND (Vaezi OR Critical Transition).
# ------------------------------------------------------------------------------

safe_quantile <- function(x, probs, na.rm = TRUE) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  suppressWarnings(as.numeric(stats::quantile(
    x, probs = probs, na.rm = na.rm, names = FALSE
  )))
}

create_surges <- function(df_region) {
  df_region <- df_region %>% dplyr::arrange(Date)
  n <- nrow(df_region)

  # 5.1 Per-region rolling donor map -----------------------------------------
  region_target_years <- sort(unique(df_region$YR))
  donor_pool_region   <- region_target_years
  rolling_map <- list()
  for (y in region_target_years) {
    donors <- intersect(seq(y - 5L, y - 1L), donor_pool_region)
    if (length(donors) < 3L) {
      donors <- tail(donor_pool_region[donor_pool_region < y], 5L)
    }
    rolling_map[[as.character(y)]] <- sort(unique(donors))
  }

  # 5.2 Rolling weekly baselines (week-specific, donor-year based) -----------
  bl_rows <- list()
  for (y in names(rolling_map)) {
    donors <- rolling_map[[y]]; y_int <- as.integer(y)
    if (length(donors) == 0L) next
    weeks_y <- df_region %>% dplyr::filter(YR == y_int) %>%
      dplyr::pull(WN) %>% unique() %>% sort()
    if (length(weeks_y) == 0L) next
    for (w in weeks_y) {
      vals <- df_region %>% dplyr::filter(YR %in% donors, WN == w) %>%
        dplyr::pull(DC_DOH)
      vals <- vals[is.finite(vals)]
      n_vals <- length(vals)
      if (n_vals == 0L) {
        mean_val <- NA_real_; sd_val <- NA_real_
        p75 <- NA_real_; p90 <- NA_real_
      } else if (n_vals == 1L) {
        mean_val <- vals[1]; sd_val <- NA_real_
        p75 <- vals[1]; p90 <- vals[1]
      } else {
        mean_val <- mean(vals, na.rm = TRUE)
        sd_val   <- stats::sd(vals, na.rm = TRUE)
        p75      <- safe_quantile(vals, 0.75)
        p90      <- safe_quantile(vals, 0.90)
      }
      bl_rows[[length(bl_rows) + 1L]] <- data.frame(
        YR = y_int, WN = as.integer(w),
        bl_mean = mean_val, bl_sd = sd_val,
        bl_p75 = p75, bl_p90 = p90,
        bl_alarm    = if (is.na(sd_val)) mean_val else mean_val + 1 * sd_val,
        bl_outbreak = if (is.na(sd_val)) mean_val else mean_val + 2 * sd_val,
        bl_n = n_vals,
        stringsAsFactors = FALSE
      )
    }
  }
  thresholds_region <- if (length(bl_rows) > 0L) {
    dplyr::bind_rows(bl_rows) %>% dplyr::arrange(YR, WN)
  } else {
    data.frame(
      YR = integer(), WN = integer(),
      bl_mean = numeric(), bl_sd = numeric(),
      bl_p75 = numeric(), bl_p90 = numeric(),
      bl_alarm = numeric(), bl_outbreak = numeric(), bl_n = integer(),
      stringsAsFactors = FALSE
    )
  }
  df_region <- df_region %>%
    dplyr::left_join(thresholds_region, by = c("YR", "WN"))

  # 5.3 Derived series -------------------------------------------------------
  dc_ma3  <- zoo::rollmean(df_region$DC_DOH, 3,  fill = NA, align = "right")
  dc_ma12 <- zoo::rollmean(df_region$DC_DOH, 12, fill = NA, align = "right")
  dc_diff <- c(NA_real_, diff(dc_ma3))

  var8  <- zoo::rollapply(df_region$DC_DOH, 8,
                          function(x) stats::var(x, na.rm = TRUE),
                          fill = NA, align = "right")
  var26 <- zoo::rollapply(df_region$DC_DOH, 26,
                          function(x) stats::var(x, na.rm = TRUE),
                          fill = NA, align = "right")
  ratio_var <- ifelse(
    !is.na(var26) & is.finite(var26) & var26 > VAR26_GUARD,
    var8 / var26, NA_real_
  )

  # 5.4 Rolling derived-series quantiles (donor-year based) ------------------
  rolling_quantile_by_target_local <- function(series, years_vec,
                                               rolling_map, prob) {
    out <- rep(NA_real_, length(series))
    for (y in names(rolling_map)) {
      donors <- rolling_map[[y]]; y_int <- as.integer(y)
      donor_idx <- which(years_vec %in% donors & is.finite(series))
      if (length(donor_idx) < 5L) next
      out[which(years_vec == y_int)] <- safe_quantile(series[donor_idx], prob)
    }
    out
  }
  df_region$rise_thr_rolling  <- rolling_quantile_by_target_local(
    dc_diff,   df_region$YR, rolling_map, 0.75)
  df_region$roc80_thr_rolling <- rolling_quantile_by_target_local(
    dc_diff,   df_region$YR, rolling_map, 0.80)
  df_region$ct80_thr_rolling  <- rolling_quantile_by_target_local(
    ratio_var, df_region$YR, rolling_map, 0.80)

  # 5.5 Walk-forward CUSUM (year-boundary reset + rolling baseline) ----------
  k_mult <- 0.5; h_mult <- 5.0
  cusum_vec <- rep(0, n)
  cusum_h   <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    is_year_boundary <- i > 1L &&
      !is.na(df_region$YR[i]) && !is.na(df_region$YR[i - 1L]) &&
      df_region$YR[i] != df_region$YR[i - 1L]
    mu0    <- df_region$bl_mean[i]
    sigma0 <- df_region$bl_sd[i]
    prev_val <- if (i == 1L || is_year_boundary) 0 else cusum_vec[i - 1L]
    if (!is.na(df_region$DC_DOH[i]) && !is.na(mu0) &&
        !is.na(sigma0) && sigma0 > 0) {
      cusum_vec[i] <- max(0, prev_val +
                            (df_region$DC_DOH[i] - mu0 - k_mult * sigma0))
    } else {
      cusum_vec[i] <- prev_val
    }
    cusum_h[i] <- if (!is.na(sigma0) && sigma0 > 0) h_mult * sigma0 else NA_real_
  }
  df_region$cusum_val <- cusum_vec
  df_region$cusum_h   <- cusum_h

  # 5.6 Primary detector signals (using rolling thresholds) ------------------
  df_region <- df_region %>%
    dplyr::mutate(
      surge_mean_2sd       = as.integer(!is.na(DC_DOH) & !is.na(bl_outbreak) &
                                          DC_DOH > bl_outbreak),
      surge_mean_1sd       = as.integer(!is.na(DC_DOH) & !is.na(bl_alarm) &
                                          DC_DOH > bl_alarm),
      surge_who_75         = as.integer(!is.na(DC_DOH) & !is.na(bl_p75) &
                                          DC_DOH > bl_p75),
      surge_who_90         = as.integer(!is.na(DC_DOH) & !is.na(bl_p90) &
                                          DC_DOH > bl_p90),
      surge_rate_change    = as.integer(!is.na(dc_diff) &
                                          !is.na(roc80_thr_rolling) &
                                          dc_diff > roc80_thr_rolling),
      surge_sta_lta        = as.integer(!is.na(dc_ma3) & !is.na(dc_ma12) &
                                          dc_ma12 > 0 &
                                          (dc_ma3 / dc_ma12) > ETA_ON_CLASSIC),
      surge_cusum          = as.integer(!is.na(cusum_val) & !is.na(cusum_h) &
                                          cusum_val > cusum_h),
      surge_critical_trans = as.integer(!is.na(ratio_var) &
                                          !is.na(ct80_thr_rolling) &
                                          ratio_var > ct80_thr_rolling),
      surge_hydrology      = as.integer(
        !is.na(DC_DOH) & !is.na(dc_diff) & !is.na(bl_p75) &
          !is.na(rise_thr_rolling) &
          DC_DOH > bl_p75 & dc_diff > rise_thr_rolling
      )
    )

  # 5.7 Vaezi-style STA/LTA hysteresis ---------------------------------------
  MIN_T <- STA_WIN + GUARD + LTA_WIN
  dc <- df_region$DC_DOH
  triggered_v <- rep(FALSE, n)
  R_vaezi_v   <- rep(NA_real_, n)
  is_on <- FALSE; frozen_lta <- NA_real_; consec_off <- 0L
  for (t in seq_len(n)) {
    if (t > 1L && !is.na(df_region$YR[t]) &&
        !is.na(df_region$YR[t - 1L]) &&
        df_region$YR[t] != df_region$YR[t - 1L]) {
      is_on <- FALSE; frozen_lta <- NA_real_; consec_off <- 0L
    }
    if (!is_on) consec_off <- consec_off + 1L else consec_off <- 0L
    if (!is_on && consec_off >= MIN_OFF_RESET) frozen_lta <- NA_real_
    if (t < MIN_T) next

    sta_vals <- dc[(t - STA_WIN + 1L):t]
    sta <- if (all(is.na(sta_vals))) NA_real_ else mean(sta_vals, na.rm = TRUE)

    if (!is_on || is.na(frozen_lta)) {
      lta_idx <- (t - STA_WIN - GUARD - LTA_WIN + 1L):(t - STA_WIN - GUARD)
      if (length(lta_idx) == LTA_WIN && min(lta_idx) > 0L) {
        lta_vals <- dc[lta_idx]
        frozen_lta <- if (all(is.na(lta_vals))) NA_real_
                      else mean(lta_vals, na.rm = TRUE)
      } else {
        frozen_lta <- NA_real_
      }
    }
    R_t <- if (!is.na(frozen_lta) && frozen_lta > 0 && !is.na(sta))
      sta / frozen_lta else NA_real_
    R_vaezi_v[t] <- R_t

    if (!is_on && !is.na(R_t) && R_t >= ETA_ON) { is_on <- TRUE; consec_off <- 0L }
    if ( is_on && !is.na(R_t) && R_t < ETA_OFF) { is_on <- FALSE; frozen_lta <- NA_real_ }
    triggered_v[t] <- is_on
  }
  df_region$surge_sta_lta_vaezi <- as.integer(triggered_v)
  df_region$R_vaezi             <- R_vaezi_v

  # 5.8 Composite signal -----------------------------------------------------
  df_region$surge_composite <- as.integer(
    df_region$surge_mean_2sd == 1L &
      (df_region$surge_sta_lta_vaezi == 1L | df_region$surge_critical_trans == 1L)
  )

  # --- Contemporary comparators: Farrington, EARS, EWARS --------------------
  # Built with the shared constructors in R/02_detection_framework.R, so these
  # are the same trigger definitions used by Stage 3 and Stage 6.
  # EWARS here is an EWARS-STYLE alarm-indicator model, not the WHO software.
  .cs <- ifelse(is.na(df_region$DC_DOH), 0, df_region$DC_DOH)
  .rn <- if ("RF_HDX" %in% names(df_region)) df_region$RF_HDX else rep(NA_real_, nrow(df_region))
  .ev_idx <- which(df_region$YR %in% EVALUABLE_YEARS)

  df_region$surge_ears <- build_ears(.cs)
  df_region$surge_farrington <- build_farrington(.cs, df_region$YR, .ev_idx)$alarm
  df_region$surge_ewars <- build_ewars(.cs, .rn, df_region$YR,
                                       eval_years = EVALUABLE_YEARS,
                                       excluded_years = EXCLUDED_YEARS)$alarm

  df_region
}

cat("Generating per-region surge signals...\n")
df_all <- df_all %>%
  dplyr::group_by(REGION) %>%
  dplyr::group_modify(~ create_surges(.x)) %>%
  dplyr::ungroup()

# ------------------------------------------------------------------------------
# 6. METHOD MAPPING
# ------------------------------------------------------------------------------
surge_defs <- c(
  "Continuous Transmission Acceleration" = "surge_sta_lta",
  "Constant Transmission Acceleration"   = "surge_sta_lta_vaezi",
  "Incidence Gradient"                   = "surge_rate_change",
  "Hydrological Inflection Measure"      = "surge_hydrology",
  "Alarm Threshold"                      = "surge_mean_1sd",
  "Composite Outbreak Signal"            = "surge_composite",
  "Outbreak Threshold"                   = "surge_mean_2sd",
  "WHO 90th Percentile Threshold"        = "surge_who_90",
  "Critical Transition Indicator"        = "surge_critical_trans",
  "Cumulative Sum Control"               = "surge_cusum",
  "WHO 75th Percentile Threshold"        = "surge_who_75",
  "Farrington"                           = "surge_farrington",
  "EARS"                                 = "surge_ears",
  "EWARS"                                = "surge_ewars"
)

# NAME-KEYED, not positional: a positional vector here is what produced the
# "Size 14 / Size 11" tibble error when detectors were added in Stage 3.
PARADIGM_LOOKUP <- c(
  "Continuous Transmission Acceleration" = "Acceleration Measures",
  "Constant Transmission Acceleration"   = "Acceleration Measures",
  "Incidence Gradient"                   = "Acceleration Measures",
  "Hydrological Inflection Measure"      = "Acceleration Measures",
  "Critical Transition Indicator"        = "Acceleration Measures",
  "Cumulative Sum Control"               = "Acceleration Measures",
  "Composite Outbreak Signal"            = "Acceleration Measures",
  "Alarm Threshold"                      = "Retrospective Thresholds",
  "Outbreak Threshold"                   = "Retrospective Thresholds",
  "WHO 90th Percentile Threshold"        = "National Standard",
  "WHO 75th Percentile Threshold"        = "National Standard",
  "Farrington"                           = "Contemporary Surveillance Algorithms",
  "EARS"                                 = "Contemporary Surveillance Algorithms",
  "EWARS"                                = "Contemporary Surveillance Algorithms"
)
.untyped <- setdiff(names(surge_defs), names(PARADIGM_LOOKUP))
if (length(.untyped) > 0) {
  stop("Detector(s) with no paradigm: ", paste(.untyped, collapse = ", "),
       call. = FALSE)
}
method_type_map <- tibble::tibble(
  Method   = names(surge_defs),
  Paradigm = unname(PARADIGM_LOOKUP[names(surge_defs)])
)

# Canonical method display order (paradigm-grouped):
#   Acceleration Measures (7), Retrospective Thresholds (2), National Standard (2).
method_order <- c(
  "Continuous Transmission Acceleration",
  "Constant Transmission Acceleration",
  "Incidence Gradient",
  "Hydrological Inflection Measure",
  "Critical Transition Indicator",
  "Cumulative Sum Control",
  "Composite Outbreak Signal",
  "Alarm Threshold",
  "Outbreak Threshold",
  "WHO 90th Percentile Threshold",
  "WHO 75th Percentile Threshold",
  "Farrington",
  "EARS",
  "EWARS"
)
stopifnot(setequal(method_order, names(surge_defs)))
stopifnot(length(method_order) == length(surge_defs))

# Multi-line wrapped detector tick labels.
method_two_line <- c(
  "Farrington"                           = "Farrington",
  "EARS"                                 = "EARS",
  "EWARS"                                = "EWARS",
  "Continuous Transmission Acceleration" = "Continuous\nTransmission\nAcceleration",
  "Constant Transmission Acceleration"   = "Constant\nTransmission\nAcceleration",
  "Incidence Gradient"                   = "Incidence\nGradient",
  "Hydrological Inflection Measure"      = "Hydrological\nInflection\nMeasure",
  "Critical Transition Indicator"        = "Critical\nTransition\nIndicator",
  "Cumulative Sum Control"               = "Cumulative\nSum\nControl",
  "Composite Outbreak Signal"            = "Composite\nOutbreak\nSignal",
  "Alarm Threshold"                      = "Alarm\nThreshold",
  "Outbreak Threshold"                   = "Outbreak\nThreshold",
  "WHO 90th Percentile Threshold"        = "WHO 90th\nPercentile\nThreshold",
  "WHO 75th Percentile Threshold"        = "WHO 75th\nPercentile\nThreshold"
)
.missing_two_line_keys <- setdiff(method_order, names(method_two_line))
if (length(.missing_two_line_keys) > 0L) {
  stop("method_two_line is missing keys for: ",
       paste(.missing_two_line_keys, collapse = ", "),
       ". Update the named vector to include all method_order entries.")
}

# ------------------------------------------------------------------------------
# 7. ANCHORS AND TRIGGER CLASSIFICATION
# ------------------------------------------------------------------------------
peak_index_whichmax <- function(dc_vec) {
  if (length(dc_vec) == 0) return(NA_integer_)
  valid_dc <- ifelse(is.na(dc_vec), -Inf, dc_vec)
  if (all(!is.finite(valid_dc))) return(NA_integer_)
  pk <- which.max(valid_dc)
  if (length(pk) == 0 || !is.finite(valid_dc[pk])) return(NA_integer_)
  as.integer(pk)
}

compute_A1 <- function(df_in, year,
                       lead_min = A1_LEAD_MIN, lead_max = A1_LEAD_MAX) {
  df_y <- df_in %>% dplyr::filter(YR == year)
  if (nrow(df_y) == 0 || all(is.na(df_y$DC_QC)))
    return(list(peak_week = NA_integer_, A1_weeks = integer(0)))
  peak_idx  <- which.max(df_y$DC_QC)
  peak_week <- df_y$WN[peak_idx]
  A1_weeks  <- (peak_week - lead_max):(peak_week - lead_min)
  A1_weeks  <- A1_weeks[A1_weeks >= 1L]
  list(peak_week = as.integer(peak_week),
       A1_weeks  = as.integer(A1_weeks))
}

compute_A2 <- function(df_in, year, burden_frac = A2_BURDEN_FRAC) {
  df_y <- df_in %>% dplyr::filter(YR == year) %>% dplyr::arrange(WN)
  if (nrow(df_y) == 0 || all(is.na(df_y$DC_QC)))
    return(list(start_week = NA_integer_, end_week = NA_integer_,
                A2_weeks = integer(0)))
  cases <- ifelse(is.na(df_y$DC_QC), 0, df_y$DC_QC)
  weeks <- df_y$WN; n_w <- length(weeks); total <- sum(cases)
  if (total <= 0)
    return(list(start_week = NA_integer_, end_week = NA_integer_,
                A2_weeks = integer(0)))
  peak_idx <- which.max(cases); lo <- hi <- peak_idx; S <- cases[peak_idx]
  while (S / total < burden_frac && (lo > 1L || hi < n_w)) {
    L <- if (lo > 1L)  cases[lo - 1L] else -Inf
    R <- if (hi < n_w) cases[hi + 1L] else -Inf
    if (L >= R) { lo <- lo - 1L; S <- S + L }
    else        { hi <- hi + 1L; S <- S + R }
  }
  list(start_week = as.integer(weeks[lo]),
       end_week   = as.integer(weeks[hi]),
       A2_weeks   = as.integer(weeks[lo:hi]))
}

compute_anchors_for_year <- function(df_in, year,
                                     lead_min = A1_LEAD_MIN,
                                     lead_max = A1_LEAD_MAX,
                                     burden_frac = A2_BURDEN_FRAC) {
  A1 <- compute_A1(df_in, year, lead_min, lead_max)
  A2 <- compute_A2(df_in, year, burden_frac)
  list(year = year, peak_week = A1$peak_week,
       A1_weeks = A1$A1_weeks, A2_weeks = A2$A2_weeks)
}

# True alarm iff t in A1 OR t in A2 (T = A1 union A2).
classify_trigger <- function(trigger_week, A1_weeks, A2_weeks) {
  in_A1 <- trigger_week %in% A1_weeks
  in_A2 <- trigger_week %in% A2_weeks
  is_true <- in_A1 || in_A2
  list(is_true = is_true, in_A1 = in_A1, in_A2 = in_A2)
}

# Two-bucket compartment classifier on True alarms.
classify_compartment <- function(lead_time) {
  if (is.na(lead_time)) return(NA_character_)
  if (lead_time >= COMPARTMENT_ACTIONABLE_MIN &&
      lead_time <= COMPARTMENT_ACTIONABLE_MAX) return("Actionable")
  return("Reactive")
}

# ------------------------------------------------------------------------------
# 8. REGIONAL INCLUSION CRITERIA
# ------------------------------------------------------------------------------
peak_per_region_year <- df_all %>%
  dplyr::filter(YR %in% EVALUABLE_YEARS) %>%
  dplyr::group_by(REGION, YR) %>%
  dplyr::summarise(
    annual_peak = suppressWarnings(max(DC_DOH, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    annual_peak = ifelse(is.finite(annual_peak), annual_peak, 0)
  )

evaluable_region_year <- peak_per_region_year %>%
  dplyr::filter(annual_peak >= MIN_PEAK_CASES_PER_YEAR) %>%
  dplyr::select(REGION, YR)

region_year_count <- evaluable_region_year %>%
  dplyr::count(REGION, name = "n_evaluable_years")

regions_in <- region_year_count %>%
  dplyr::filter(n_evaluable_years >= MIN_EVALUABLE_YEARS_PER_REGION) %>%
  dplyr::pull(REGION)

regions_excluded <- setdiff(unique(df_all$REGION), regions_in)

cat("Included regions: ", length(regions_in), "; excluded: ",
    length(regions_excluded), "\n", sep = "")

# ------------------------------------------------------------------------------
# 9. PER-REGION FRAMEWORK METRICS
# ------------------------------------------------------------------------------
build_trigger_detail_for_region <- function(df_region, evaluable_years_for_region) {
  rows <- list()
  for (yr in evaluable_years_for_region) {
    anchors <- compute_anchors_for_year(df_region, yr)
    df_y <- df_region %>% dplyr::filter(YR == yr) %>% dplyr::arrange(WN)
    if (nrow(df_y) == 0) next

    for (method_name in names(surge_defs)) {
      col_name <- surge_defs[[method_name]]
      trig_idx <- which(df_y[[col_name]] == 1L)
      if (length(trig_idx) == 0L) next
      for (k in trig_idx) {
        wk  <- df_y$WN[k]
        cls <- classify_trigger(wk, anchors$A1_weeks, anchors$A2_weeks)
        lead_time <- if (!is.na(anchors$peak_week))
          as.integer(anchors$peak_week - wk) else NA_integer_
        # Compartment is meaningful only for True alarms.
        comp <- if (isTRUE(cls$is_true)) classify_compartment(lead_time)
                else NA_character_
        case_count <- df_y$DC_DOH[k]
        rows[[length(rows) + 1L]] <- data.frame(
          Year = yr, Method = method_name, Week = as.integer(wk),
          Peak_Week = anchors$peak_week, Lead_Time = lead_time,
          Compartment = comp,
          IsTrue = cls$is_true,
          InA1 = cls$in_A1, InA2 = cls$in_A2,
          DC = if (is.na(case_count)) 0 else case_count,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      Year = integer(), Method = character(), Week = integer(),
      Peak_Week = integer(), Lead_Time = integer(),
      Compartment = character(), IsTrue = logical(),
      InA1 = logical(), InA2 = logical(),
      DC = numeric(), stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(rows)
}

compute_yearly_lead_for_region <- function(trig_aug, evaluable_years_for_region) {
  if (nrow(trig_aug) == 0)
    return(data.frame(
      Method = character(), Year = integer(),
      First_A1_True_Week = integer(), Lead_Time_Yr = numeric(),
      stringsAsFactors = FALSE
    ))
  trig_aug %>%
    dplyr::filter(Year %in% evaluable_years_for_region, InA1, IsTrue) %>%
    dplyr::group_by(Method, Year) %>%
    dplyr::slice_min(Week, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      Method,
      Year = as.integer(Year),
      First_A1_True_Week = as.integer(Week),
      Lead_Time_Yr = as.numeric(Lead_Time)
    )
}

compute_method_metrics_region <- function(trig_aug, yearly_lead,
                                          evaluable_years_for_region) {
  rows <- list()
  n_eval <- length(evaluable_years_for_region)
  for (m in names(surge_defs)) {
    trig_sub <- trig_aug %>%
      dplyr::filter(Method == m, Year %in% evaluable_years_for_region)
    lead_sub <- yearly_lead %>%
      dplyr::filter(Method == m, Year %in% evaluable_years_for_region)

    total      <- nrow(trig_sub)
    true_n     <- sum(trig_sub$IsTrue, na.rm = TRUE)
    false_n    <- total - true_n
    reactive_n <- sum(trig_sub$Compartment == "Reactive", na.rm = TRUE)
    truact_n   <- sum(trig_sub$IsTrue & trig_sub$Compartment == "Actionable",
                      na.rm = TRUE)
    truact_lt  <- trig_sub$Lead_Time[
      trig_sub$IsTrue & trig_sub$Compartment == "Actionable"
    ]

    ppv <- if (total > 0)  true_n / total else NA_real_

    # A1-restricted Sensitivity.
    years_with_a1_true <- length(unique(
      trig_sub$Year[trig_sub$IsTrue & trig_sub$Compartment == "Actionable"]
    ))
    sens <- if (n_eval > 0) years_with_a1_true / n_eval else NA_real_

    # Mean Lead Time: same-denominator headline + conditional variant.
    mean_lead_conditional <- if (nrow(lead_sub) > 0)
      mean(lead_sub$Lead_Time_Yr, na.rm = TRUE) else NA_real_
    if (is.nan(mean_lead_conditional)) mean_lead_conditional <- NA_real_
    mean_lead <- if (n_eval > 0)
      sum(lead_sub$Lead_Time_Yr, na.rm = TRUE) / n_eval else NA_real_
    if (is.nan(mean_lead)) mean_lead <- NA_real_

    # Warning Persistence: same-denominator headline + conditional variant.
    wp_conditional <- if (length(truact_lt) > 0)
      mean(truact_lt, na.rm = TRUE) else NA_real_
    if (is.nan(wp_conditional)) wp_conditional <- NA_real_
    if (n_eval > 0L) {
      per_year_wp <- vapply(evaluable_years_for_region, function(y) {
        lt <- trig_sub$Lead_Time[
          trig_sub$Year == y &
            trig_sub$IsTrue &
            trig_sub$Compartment == "Actionable"
        ]
        if (length(lt) > 0L) mean(lt, na.rm = TRUE) else 0
      }, numeric(1))
      wp <- mean(per_year_wp, na.rm = TRUE)
      if (is.nan(wp)) wp <- NA_real_
    } else {
      wp <- NA_real_
    }

    aly       <- if (true_n > 0) truact_n / true_n else NA_real_
    n_true_yr <- if (n_eval > 0) true_n  / n_eval else NA_real_
    n_false_yr <- if (n_eval > 0) false_n / n_eval else NA_real_

    # TAM (True-Alarm Magnitude): per-year T-restricted alarm-on case sum,
    # averaged across evaluable years.
    tam_per_year <- vapply(evaluable_years_for_region, function(y) {
      ix <- which(trig_sub$Year == y & trig_sub$IsTrue)
      if (length(ix) == 0) 0 else sum(trig_sub$DC[ix], na.rm = TRUE)
    }, numeric(1))
    tam <- mean(tam_per_year, na.rm = TRUE)

    rows[[m]] <- data.frame(
      Method = m,
      # Reported metrics
      TAM            = tam,
      N_True_Alarms  = n_true_yr,
      Sensitivity    = sens,
      Mean_Lead_Time = mean_lead,
      WP             = wp,
      # Diagnostic / wide-CSV-only columns
      Mean_Lead_Time_conditional = mean_lead_conditional,
      WP_conditional             = wp_conditional,
      PPV            = ppv,
      ALY            = aly,
      N_False_Alarms = n_false_yr,
      Total_Triggers = total,
      True_Alarms    = true_n,
      False_Alarms   = false_n,
      n_Reactive            = reactive_n,
      n_TrueActionable      = truact_n,
      N_Years_with_A1_True  = years_with_a1_true,
      N_Years_Evaluable     = n_eval,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

# Per-region pipeline driver
regional_metrics_list <- list()
for (reg in regions_in) {
  yrs_reg <- evaluable_region_year %>%
    dplyr::filter(REGION == reg) %>%
    dplyr::pull(YR) %>% sort()
  if (length(yrs_reg) < MIN_EVALUABLE_YEARS_PER_REGION) next
  df_reg     <- df_all %>% dplyr::filter(REGION == reg)
  trig_aug   <- build_trigger_detail_for_region(df_reg, yrs_reg)
  yearly_ld  <- compute_yearly_lead_for_region(trig_aug, yrs_reg)
  metrics    <- compute_method_metrics_region(trig_aug, yearly_ld, yrs_reg)
  metrics$REGION <- reg
  regional_metrics_list[[reg]] <- metrics
}

regional_metrics <- dplyr::bind_rows(regional_metrics_list) %>%
  dplyr::left_join(method_type_map, by = "Method") %>%
  dplyr::mutate(
    Method   = factor(Method,   levels = method_order),
    Paradigm = factor(Paradigm, levels = type_order)
  )

if (nrow(regional_metrics) == 0)
  stop("No regional metrics were computed. Check inclusion criteria and data.")

# ------------------------------------------------------------------------------
# 10. YEAR-CLUSTER BOOTSTRAP (Cameron-Gelbach-Miller, B = 1000)
# ------------------------------------------------------------------------------
# For each region, resample evaluable years with replacement (year is the
# cluster unit because within-year weekly observations are not independent).
# On each replicate, recompute all framework metrics from the cached trigger
# detail subset. Trigger columns and anchors are deterministic given case
# data and are cached outside the bootstrap loop for speed.
# ------------------------------------------------------------------------------
cat("\n=== Year-cluster bootstrap (B =", BOOT_N_CI, "per region) ===\n")
cat("This computation may take a few minutes for ",
    length(regions_in), " regions x ", BOOT_N_CI, " replicates.\n", sep = "")

# Anchor cache (one entry per region/year).
anchor_cache <- list()
for (reg in regions_in) {
  yrs_reg <- evaluable_region_year %>%
    dplyr::filter(REGION == reg) %>%
    dplyr::pull(YR) %>% sort()
  if (length(yrs_reg) < MIN_EVALUABLE_YEARS_PER_REGION) next
  df_reg <- df_all %>% dplyr::filter(REGION == reg)
  anchor_cache[[reg]] <- list()
  for (yr in yrs_reg) {
    anchor_cache[[reg]][[as.character(yr)]] <- compute_anchors_for_year(df_reg, yr)
  }
}

# Bootstrap-specific trigger detail builder using the anchor cache.
build_trigger_detail_cached <- function(df_region, evaluable_years_for_region,
                                        region_name) {
  rows <- list()
  for (yr in evaluable_years_for_region) {
    anchors <- anchor_cache[[region_name]][[as.character(yr)]]
    if (is.null(anchors)) next
    df_y <- df_region %>% dplyr::filter(YR == yr) %>% dplyr::arrange(WN)
    if (nrow(df_y) == 0) next
    for (method_name in names(surge_defs)) {
      col_name <- surge_defs[[method_name]]
      trig_idx <- which(df_y[[col_name]] == 1L)
      if (length(trig_idx) == 0L) next
      for (k in trig_idx) {
        wk  <- df_y$WN[k]
        cls <- classify_trigger(wk, anchors$A1_weeks, anchors$A2_weeks)
        lead_time <- if (!is.na(anchors$peak_week))
          as.integer(anchors$peak_week - wk) else NA_integer_
        comp <- if (isTRUE(cls$is_true)) classify_compartment(lead_time)
                else NA_character_
        case_count <- df_y$DC_DOH[k]
        rows[[length(rows) + 1L]] <- data.frame(
          Year = yr, Method = method_name, Week = as.integer(wk),
          Peak_Week = anchors$peak_week, Lead_Time = lead_time,
          Compartment = comp,
          IsTrue = cls$is_true,
          InA1 = cls$in_A1, InA2 = cls$in_A2,
          DC = if (is.na(case_count)) 0 else case_count,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      Year = integer(), Method = character(), Week = integer(),
      Peak_Week = integer(), Lead_Time = integer(),
      Compartment = character(), IsTrue = logical(),
      InA1 = logical(), InA2 = logical(),
      DC = numeric(), stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(rows)
}

# Pre-compute the FULL (unsampled) trigger detail per region.
full_trig_detail <- list()
for (reg in names(anchor_cache)) {
  df_reg  <- df_all %>% dplyr::filter(REGION == reg)
  yrs_reg <- as.integer(names(anchor_cache[[reg]]))
  full_trig_detail[[reg]] <- build_trigger_detail_cached(df_reg, yrs_reg, reg)
}

# Single-replicate metric computation.
bootstrap_metrics_one_replicate <- function(reg, yrs_resampled) {
  trig_full <- full_trig_detail[[reg]]
  if (nrow(trig_full) == 0L) return(NULL)

  trig_sub <- dplyr::bind_rows(
    lapply(yrs_resampled, function(y) trig_full %>% dplyr::filter(Year == y))
  )

  if (nrow(trig_sub) == 0L) {
    yearly_ld <- data.frame(
      Method = character(), Year = integer(),
      Lead_Time_Yr = numeric(), stringsAsFactors = FALSE
    )
  } else {
    yearly_ld <- trig_sub %>%
      dplyr::filter(InA1, IsTrue) %>%
      dplyr::group_by(Method, Year) %>%
      dplyr::slice_min(Week, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(Method, Year = as.integer(Year),
                       Lead_Time_Yr = as.numeric(Lead_Time))
    yr_mult <- table(yrs_resampled)
    yearly_ld <- yearly_ld %>%
      dplyr::mutate(.mult = as.integer(unname(yr_mult[as.character(Year)]))) %>%
      dplyr::filter(!is.na(.mult), .mult > 0L) %>%
      tidyr::uncount(.mult)
  }

  rows <- list()
  n_eval <- length(yrs_resampled)
  for (m in names(surge_defs)) {
    trig_m <- trig_sub %>% dplyr::filter(Method == m)
    lead_m <- yearly_ld %>% dplyr::filter(Method == m)

    total      <- nrow(trig_m)
    true_n     <- sum(trig_m$IsTrue, na.rm = TRUE)
    false_n    <- total - true_n
    reactive_n <- sum(trig_m$Compartment == "Reactive", na.rm = TRUE)
    truact_n   <- sum(trig_m$IsTrue & trig_m$Compartment == "Actionable",
                      na.rm = TRUE)
    truact_lt  <- trig_m$Lead_Time[
      trig_m$IsTrue & trig_m$Compartment == "Actionable"
    ]

    ppv <- if (total > 0)  true_n / total else NA_real_
    years_with_a1_true <- length(unique(
      trig_m$Year[trig_m$IsTrue & trig_m$Compartment == "Actionable"]
    ))
    sens <- if (n_eval > 0) years_with_a1_true / n_eval else NA_real_

    mean_lead <- if (n_eval > 0)
      sum(lead_m$Lead_Time_Yr, na.rm = TRUE) / n_eval else NA_real_
    if (is.nan(mean_lead)) mean_lead <- NA_real_

    if (n_eval > 0L) {
      per_year_wp <- vapply(yrs_resampled, function(y) {
        lt <- trig_m$Lead_Time[
          trig_m$Year == y &
            trig_m$IsTrue &
            trig_m$Compartment == "Actionable"
        ]
        if (length(lt) > 0L) mean(lt, na.rm = TRUE) else 0
      }, numeric(1))
      wp <- mean(per_year_wp, na.rm = TRUE)
      if (is.nan(wp)) wp <- NA_real_
    } else {
      wp <- NA_real_
    }

    aly       <- if (true_n > 0) truact_n / true_n else NA_real_
    n_true_yr <- if (n_eval > 0) true_n  / n_eval else NA_real_
    n_false_yr <- if (n_eval > 0) false_n / n_eval else NA_real_

    # Cluster-bootstrap TAM: per-year sums computed once on the unduplicated
    # trigger detail, then averaged across the resampled year vector.
    trig_m_unique <- full_trig_detail[[reg]] %>% dplyr::filter(Method == m)
    yr_to_sum <- if (sum(trig_m_unique$IsTrue, na.rm = TRUE) > 0L) {
      tapply(
        trig_m_unique$DC[trig_m_unique$IsTrue],
        trig_m_unique$Year[trig_m_unique$IsTrue],
        sum, na.rm = TRUE
      )
    } else {
      stats::setNames(numeric(0), character(0))
    }
    tam_per_year <- vapply(yrs_resampled, function(y) {
      v <- yr_to_sum[as.character(y)]
      if (is.null(v) || length(v) == 0L || is.na(v)) 0 else as.numeric(v)
    }, numeric(1))
    tam <- mean(tam_per_year, na.rm = TRUE)

    rows[[m]] <- data.frame(
      Method = m,
      TAM = tam, N_True_Alarms = n_true_yr,
      Sensitivity = sens,
      Mean_Lead_Time = mean_lead, WP = wp,
      PPV = ppv, ALY = aly, N_False_Alarms = n_false_yr,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

# Bootstrap loop
boot_results <- list()
set.seed(20260101L)
for (reg in names(anchor_cache)) {
  yrs_reg <- as.integer(names(anchor_cache[[reg]]))
  if (length(yrs_reg) < 2L) {
    cat("  Skipping region '", reg, "' (only ",
        length(yrs_reg), " evaluable year(s)).\n", sep = "")
    next
  }
  reps <- vector("list", BOOT_N_CI)
  for (b in seq_len(BOOT_N_CI)) {
    yrs_resampled <- sample(yrs_reg, size = length(yrs_reg), replace = TRUE)
    reps[[b]] <- bootstrap_metrics_one_replicate(reg, yrs_resampled)
  }
  reps_df <- dplyr::bind_rows(reps, .id = ".rep")

  ci_summary <- reps_df %>%
    dplyr::group_by(Method) %>%
    dplyr::summarise(
      TAM_lo            = safe_quantile(TAM, 0.025),
      TAM_hi            = safe_quantile(TAM, 0.975),
      N_True_Alarms_lo  = safe_quantile(N_True_Alarms, 0.025),
      N_True_Alarms_hi  = safe_quantile(N_True_Alarms, 0.975),
      PPV_lo            = safe_quantile(PPV, 0.025),
      PPV_hi            = safe_quantile(PPV, 0.975),
      Sens_lo           = safe_quantile(Sensitivity, 0.025),
      Sens_hi           = safe_quantile(Sensitivity, 0.975),
      MLT_lo            = safe_quantile(Mean_Lead_Time, 0.025),
      MLT_hi            = safe_quantile(Mean_Lead_Time, 0.975),
      WP_lo             = safe_quantile(WP, 0.025),
      WP_hi             = safe_quantile(WP, 0.975),
      ALY_lo            = safe_quantile(ALY, 0.025),
      ALY_hi            = safe_quantile(ALY, 0.975),
      N_False_Alarms_lo = safe_quantile(N_False_Alarms, 0.025),
      N_False_Alarms_hi = safe_quantile(N_False_Alarms, 0.975),
      .groups = "drop"
    ) %>%
    dplyr::mutate(REGION = reg)

  boot_results[[reg]] <- list(replicates = reps_df, ci = ci_summary)
  cat("  Region '", reg, "' bootstrap complete (",
      BOOT_N_CI, " replicates).\n", sep = "")
}

boot_ci <- dplyr::bind_rows(lapply(boot_results, function(x) x$ci))

# Merge CIs into the metrics table
regional_metrics <- regional_metrics %>%
  dplyr::mutate(REGION = as.character(REGION)) %>%
  dplyr::left_join(
    boot_ci %>% dplyr::mutate(Method = as.character(Method)),
    by = c("REGION", "Method")
  ) %>%
  dplyr::mutate(REGION = factor(REGION, levels = unique(REGION)))

# ------------------------------------------------------------------------------
# 11. DOMINANCE PROBABILITY (composite mean-rank across the 5 reported metrics)
# ------------------------------------------------------------------------------
# For each replicate, rank the three target detectors on each of the five
# reported metrics (rank 1 = best; ties get average rank; all metrics are
# higher-is-better), then compute the mean rank across the five metrics.
# The replicate winner is the detector with the lowest mean rank. The
# Dominance_Probability for a region is the proportion of replicates in
# which the observed point-estimate winner wins on the same composite.
#
# Properties:
#   - Equal weight to each of the five reported metrics (Borda-style).
#   - Insensitive to absolute scale across metrics.
#   - Ties handled by ties.method = "average".
#   - PPV is computed as a diagnostic but is NOT used in this composite.
# ------------------------------------------------------------------------------
COMPOSITE_DOMINANCE_METRICS <- c(
  "TAM", "N_True_Alarms", "Sensitivity", "Mean_Lead_Time", "WP"
)

.composite_winner_per_rep <- function(reps_target,
                                      composite_metrics = COMPOSITE_DOMINANCE_METRICS) {
  ok <- Reduce(`&`, lapply(composite_metrics,
                            function(m) is.finite(reps_target[[m]])))
  reps_target <- reps_target[ok, , drop = FALSE]
  if (nrow(reps_target) == 0L) return(reps_target[0, , drop = FALSE])
  reps_target <- reps_target %>%
    dplyr::group_by(.rep) %>%
    dplyr::mutate(
      r_TAM  = rank(-TAM,            ties.method = "average"),
      r_NTA  = rank(-N_True_Alarms,  ties.method = "average"),
      r_Sens = rank(-Sensitivity,    ties.method = "average"),
      r_MLT  = rank(-Mean_Lead_Time, ties.method = "average"),
      r_WP   = rank(-WP,             ties.method = "average"),
      mean_rank = (r_TAM + r_NTA + r_Sens + r_MLT + r_WP) / 5
    ) %>%
    dplyr::slice_min(mean_rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  reps_target
}

dominance_results <- list()
# Six target detectors: the three original plus the contemporary comparators,
# so the dominance-probability figures cover every method the manuscript
# benchmarks rather than only three of them.
TARGET_DETECTORS <- c("Constant Transmission Acceleration",
                      "Continuous Transmission Acceleration",
                      "Outbreak Threshold",
                      "Farrington", "EARS", "EWARS")

# P_* column name for each target detector, generated rather than hard-coded so
# adding a detector cannot leave a stale set of columns behind.
TARGET_P_COLS <- c(
  "Constant Transmission Acceleration"   = "P_ConstantTA",
  "Continuous Transmission Acceleration" = "P_ContinuousTA",
  "Outbreak Threshold"                   = "P_OutbreakThreshold",
  "Farrington"                           = "P_Farrington",
  "EARS"                                 = "P_EARS",
  "EWARS"                                = "P_EWARS"
)
for (reg in names(boot_results)) {
  reps_df     <- boot_results[[reg]]$replicates
  reps_target <- reps_df %>% dplyr::filter(Method %in% TARGET_DETECTORS)

  rep_winners <- .composite_winner_per_rep(reps_target,
                                           COMPOSITE_DOMINANCE_METRICS)

  obs_pt <- regional_metrics %>%
    dplyr::filter(REGION == reg, Method %in% TARGET_DETECTORS) %>%
    dplyr::select(Method, dplyr::all_of(COMPOSITE_DOMINANCE_METRICS))
  if (nrow(obs_pt) >= 1L &&
      all(vapply(COMPOSITE_DOMINANCE_METRICS,
                 function(m) is.finite(obs_pt[[m]][1]),
                 logical(1)))) {
    obs_pt <- obs_pt %>%
      dplyr::mutate(
        r_TAM  = rank(-TAM,            ties.method = "average"),
        r_NTA  = rank(-N_True_Alarms,  ties.method = "average"),
        r_Sens = rank(-Sensitivity,    ties.method = "average"),
        r_MLT  = rank(-Mean_Lead_Time, ties.method = "average"),
        r_WP   = rank(-WP,             ties.method = "average"),
        mean_rank = (r_TAM + r_NTA + r_Sens + r_MLT + r_WP) / 5
      ) %>%
      dplyr::arrange(mean_rank)
    obs_winner <- as.character(obs_pt$Method[1])
  } else {
    obs_winner <- character(0)
  }

  win_tab    <- table(rep_winners$Method)
  total_reps <- sum(win_tab)

  dom_prob <- if (length(obs_winner) > 0L && obs_winner %in% names(win_tab)) {
    as.numeric(win_tab[obs_winner]) / total_reps
  } else {
    NA_real_
  }

  dominance_results[[reg]] <- data.frame(
    REGION = reg,
    Observed_Winner = if (length(obs_winner) > 0L) obs_winner else NA_character_,
    Dominance_Probability = dom_prob,
    Bootstrap_N    = total_reps,
    Composite_Rule = "mean_rank_5metrics",
    stringsAsFactors = FALSE
  )
  # One P_* column per target detector, added programmatically.
  for (dn in TARGET_DETECTORS) {
    dominance_results[[reg]][[unname(TARGET_P_COLS[dn])]] <-
      if (dn %in% names(win_tab)) as.numeric(win_tab[dn]) / total_reps else 0
  }
}
dominance_df <- dplyr::bind_rows(dominance_results) %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("P_"),
                              ~ ifelse(is.na(.), 0, .)))

cat("\n=== Bootstrap dominance probability (composite mean-rank, 5 metrics) ===\n")
print(as.data.frame(dominance_df), row.names = FALSE)
cat("\n")


# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 13 ADDITIONAL OUTPUT
# Constant TA anchor consensus check against four comparators:
#   1. Outbreak Threshold
#   2. Farrington
#   3. EWARS
#   4. EARS
#
# This block is CSV-only. It does not alter the existing Figure 5 consensus
# analysis, detector map, or the original three-detector inference objects.
# It mirrors the existing regional bootstrap head-to-head logic, but uses four
# Constant-TA-anchored comparisons. Therefore:
#   strict cross-region Bonferroni: k = 4 x 17 = 68
#   within-region Bonferroni:       k = 4
# ------------------------------------------------------------------------------

ST13_COMPARATORS <- c(
  "Outbreak Threshold" = "P_OutbreakThreshold",
  "Farrington"         = "P_Farrington",
  "EWARS"              = "P_EWARS",
  "EARS"               = "P_EARS"
)
ST13_K_WITHIN <- length(ST13_COMPARATORS)
ST13_K_STRICT <- ST13_K_WITHIN * length(unique(as.character(dominance_df$REGION)))

.st13_pair_test <- function(p_const, p_comp, B = N_BOOTS) {
  if (!is.finite(p_const) || !is.finite(p_comp)) {
    return(list(p_raw = NA_real_, winner = NA_character_))
  }
  n_const <- round(p_const * B)
  n_comp  <- round(p_comp  * B)
  n_pair  <- n_const + n_comp
  if (n_pair <= 0L || n_const == n_comp) {
    return(list(p_raw = if (n_pair <= 0L) NA_real_ else 1.0,
                winner = NA_character_))
  }
  if (n_const > n_comp) {
    p_raw <- stats::binom.test(n_const, n_pair, p = 0.5,
                               alternative = "greater")$p.value
    winner <- "Constant Transmission Acceleration"
  } else {
    p_raw <- stats::binom.test(n_comp, n_pair, p = 0.5,
                               alternative = "greater")$p.value
    winner <- names(which.max(c(Comparator = n_comp)))[1]
    winner <- "Comparator"
  }
  list(p_raw = p_raw, winner = winner)
}

st13_rows <- list()
for (i in seq_len(nrow(dominance_df))) {
  rr <- dominance_df[i, , drop = FALSE]
  reg <- as.character(rr$REGION)
  p_const <- as.numeric(rr$P_ConstantTA)

  one_region <- lapply(names(ST13_COMPARATORS), function(comp) {
    p_comp <- as.numeric(rr[[unname(ST13_COMPARATORS[[comp]])]])
    tt <- .st13_pair_test(p_const, p_comp, B = N_BOOTS)
    winner <- if (is.na(tt$winner)) NA_character_ else if (tt$winner == "Comparator") comp else tt$winner
    data.frame(
      REGION = reg,
      Comparator = comp,
      P_ConstantTA = p_const,
      P_Comparator = p_comp,
      Pair_Winner = winner,
      p_raw = tt$p_raw,
      p_bonf_strict = ifelse(is.na(tt$p_raw), NA_real_, min(1, tt$p_raw * ST13_K_STRICT)),
      p_bonf_within = ifelse(is.na(tt$p_raw), NA_real_, min(1, tt$p_raw * ST13_K_WITHIN)),
      stringsAsFactors = FALSE
    )
  })
  st13_rows[[reg]] <- dplyr::bind_rows(one_region)
}

st13_long <- dplyr::bind_rows(st13_rows) %>%
  dplyr::mutate(
    Sig_Strict = !is.na(p_bonf_strict) & p_bonf_strict < HH_ALPHA,
    Sig_Within = !is.na(p_bonf_within) & p_bonf_within < HH_ALPHA,
    Strict_Result = dplyr::case_when(
      is.na(p_bonf_strict) ~ "NA",
      !Sig_Strict ~ "ns",
      Pair_Winner == "Constant Transmission Acceleration" ~ "Constant TA",
      TRUE ~ Comparator
    ),
    Within_Result = dplyr::case_when(
      is.na(p_bonf_within) ~ "NA",
      !Sig_Within ~ "ns",
      Pair_Winner == "Constant Transmission Acceleration" ~ "Constant TA",
      TRUE ~ Comparator
    )
  )

# Wide table designed for direct rendering of Supplementary Table 13.
st13_wide <- st13_long %>%
  dplyr::group_by(REGION) %>%
  dplyr::summarise(
    `Pr_ConstantTA` = dplyr::first(P_ConstantTA),
    `Strict_ConstantTA_vs_OT` = Strict_Result[Comparator == "Outbreak Threshold"][1],
    `Strict_ConstantTA_vs_Farrington` = Strict_Result[Comparator == "Farrington"][1],
    `Strict_ConstantTA_vs_EWARS` = Strict_Result[Comparator == "EWARS"][1],
    `Strict_ConstantTA_vs_EARS` = Strict_Result[Comparator == "EARS"][1],
    `Within_ConstantTA_vs_OT` = Within_Result[Comparator == "Outbreak Threshold"][1],
    `Within_ConstantTA_vs_Farrington` = Within_Result[Comparator == "Farrington"][1],
    `Within_ConstantTA_vs_EWARS` = Within_Result[Comparator == "EWARS"][1],
    `Within_ConstantTA_vs_EARS` = Within_Result[Comparator == "EARS"][1],
    `n_sig_wins_strict` = sum(Sig_Strict & Pair_Winner == "Constant Transmission Acceleration", na.rm = TRUE),
    `n_sig_losses_strict` = sum(Sig_Strict & !is.na(Pair_Winner) & Pair_Winner != "Constant Transmission Acceleration", na.rm = TRUE),
    `n_sig_wins_within` = sum(Sig_Within & Pair_Winner == "Constant Transmission Acceleration", na.rm = TRUE),
    `n_sig_losses_within` = sum(Sig_Within & !is.na(Pair_Winner) & Pair_Winner != "Constant Transmission Acceleration", na.rm = TRUE),
    `Weakest_link_P_strict` = ifelse(
      all(is.na(p_bonf_strict[Pair_Winner == "Constant Transmission Acceleration"])),
      NA_real_,
      max(p_bonf_strict[Pair_Winner == "Constant Transmission Acceleration"], na.rm = TRUE)
    ),
    `Consensus_winner_strict` = dplyr::case_when(
      sum(Sig_Strict & Pair_Winner == "Constant Transmission Acceleration", na.rm = TRUE) == ST13_K_WITHIN ~ "Constant Transmission Acceleration",
      sum(Sig_Strict & Pair_Winner != "Constant Transmission Acceleration", na.rm = TRUE) > 0 ~ "Contested",
      TRUE ~ "No strict consensus"
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    REGION = factor(REGION, levels = CANONICAL_17)
  ) %>%
  dplyr::arrange(REGION) %>%
  dplyr::mutate(REGION = as.character(REGION))

utils::write.csv(
  st13_long,
  file.path(OUTPUT_DIR, "Stage4_Regional_ConstantTA_4Comparator_Consensus_Long.csv"),
  row.names = FALSE
)
utils::write.csv(
  st13_wide,
  file.path(OUTPUT_DIR, "Stage4_Regional_ConstantTA_4Comparator_Consensus.csv"),
  row.names = FALSE
)

cat("Saved Supplementary Table 13 four-comparator consensus outputs:\n")
cat("  Stage4_Regional_ConstantTA_4Comparator_Consensus_Long.csv\n")
cat("  Stage4_Regional_ConstantTA_4Comparator_Consensus.csv\n")
cat("  Strict Bonferroni k = ", ST13_K_STRICT,
    "; within-region Bonferroni k = ", ST13_K_WITHIN, "\n", sep = "")

# Join dominance probability into regional_metrics (region-level property,
# inherited by every method row in that region).
regional_metrics <- regional_metrics %>%
  dplyr::mutate(REGION = as.character(REGION)) %>%
  dplyr::left_join(
    dominance_df %>%
      dplyr::select(REGION, Observed_Winner, Dominance_Probability,
                    dplyr::all_of(unname(TARGET_P_COLS)), Bootstrap_N),
    by = "REGION"
  ) %>%
  dplyr::mutate(REGION = factor(REGION, levels = unique(REGION)))

# ------------------------------------------------------------------------------
# 12. ORDER REGIONS AND METHODS; ESTABLISH METRIC METADATA
# ------------------------------------------------------------------------------
region_ppv_order <- regional_metrics %>%
  dplyr::group_by(REGION) %>%
  dplyr::summarise(mean_ppv = mean(PPV, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(mean_ppv)) %>%
  dplyr::pull(REGION)

regional_metrics <- regional_metrics %>%
  dplyr::mutate(REGION = factor(REGION, levels = rev(region_ppv_order))) %>%
  dplyr::mutate(
    Method   = factor(as.character(Method),   levels = method_order),
    Paradigm = factor(as.character(Paradigm), levels = type_order)
  )

# Five reported head-to-head metrics (Figure 4 figure set).
HH_METRICS_REGIONAL <- c(
  "TAM", "N_True_Alarms", "Sensitivity",
  "Mean_Lead_Time", "WP"
)

hh_metric_full_name <- c(
  TAM            = "True-Alarm Magnitude",
  N_True_Alarms  = "Number of True Alarms",
  Sensitivity    = "Sensitivity",
  Mean_Lead_Time = "Mean Lead Time",
  WP             = "Warning Persistence"
)

hh_metric_y_label <- c(
  TAM            = "True-Alarm Magnitude (cases per year)",
  N_True_Alarms  = "Number of True Alarms (per year)",
  Sensitivity    = "Sensitivity",
  Mean_Lead_Time = "Mean Lead Time (weeks before peak)",
  WP             = "Warning Persistence (weeks before peak)"
)

hh_metric_ci_cols <- list(
  TAM            = c("TAM_lo",            "TAM_hi"),
  N_True_Alarms  = c("N_True_Alarms_lo",  "N_True_Alarms_hi"),
  Sensitivity    = c("Sens_lo",           "Sens_hi"),
  Mean_Lead_Time = c("MLT_lo",            "MLT_hi"),
  WP             = c("WP_lo",             "WP_hi")
)

# All five reported metrics are higher-is-better.
hh_metric_higher_better <- c(
  TAM = TRUE,  N_True_Alarms = TRUE,  Sensitivity = TRUE,
  Mean_Lead_Time = TRUE, WP = TRUE
)

# Operational reference threshold (NA = no fixed threshold).
hh_metric_threshold <- c(
  TAM            = NA_real_,
  N_True_Alarms  = NA_real_,
  Sensitivity    = 1.00,
  Mean_Lead_Time = 4,
  WP             = 4
)

hh_metric_is_proportion <- c(
  TAM = FALSE, N_True_Alarms = FALSE, Sensitivity = TRUE,
  Mean_Lead_Time = FALSE, WP = FALSE
)

# Sub-letter numeric subscript used in per-metric panel headers.
hh_metric_subscript <- c(
  TAM            = "1",
  N_True_Alarms  = "2",
  Sensitivity    = "3",
  Mean_Lead_Time = "4",
  WP             = "5"
)

# Filename slug per metric (Figure 4 outputs).
hh_metric_slug <- c(
  TAM            = "TAM",
  N_True_Alarms  = "N_True_Alarms",
  Sensitivity    = "Sensitivity",
  Mean_Lead_Time = "Mean_Lead_Time",
  WP             = "Warning_Persistence"
)

# Detector-pair definitions for the Wilcoxon side panel.
DETECTOR_PAIRS <- list(
  list(A = "Constant Transmission Acceleration",
       B = "Outbreak Threshold",
       short = "Constant TA  vs  OT"),
  list(A = "Continuous Transmission Acceleration",
       B = "Outbreak Threshold",
       short = "Continuous TA  vs  OT"),
  list(A = "Constant Transmission Acceleration",
       B = "Continuous Transmission Acceleration",
       short = "Constant TA  vs  Continuous TA")
)
HH_ALPHA_BONF <- HH_ALPHA / length(DETECTOR_PAIRS)  # reference value

# Additional cross-region Constant TA head-to-head contrasts requested for
# CSV output only. These do not alter the existing Figure 4 or Figure 5 pair
# definitions. Outbreak Threshold is included here so the dedicated CSV contains
# the complete four-comparator Constant TA analysis.
ADDITIONAL_DETECTOR_PAIRS <- list(
  list(A = "Constant Transmission Acceleration",
       B = "Outbreak Threshold",
       short = "Constant TA  vs  OT"),
  list(A = "Constant Transmission Acceleration",
       B = "Farrington",
       short = "Constant TA  vs  Farrington"),
  list(A = "Constant Transmission Acceleration",
       B = "EWARS",
       short = "Constant TA  vs  EWARS"),
  list(A = "Constant Transmission Acceleration",
       B = "EARS",
       short = "Constant TA  vs  EARS")
)

# ------------------------------------------------------------------------------
# 13. WILCOXON DETECTOR-PAIRED TESTS (PER METRIC)
# ------------------------------------------------------------------------------
# For each metric, three pairwise Wilcoxon signed-rank tests paired by
# REGION (n = 17). Significance is evaluated PER PAIRWISE COMPARISON at
# alpha = 0.05; we do NOT multiply p-values across the three pairs within
# a metric (each pair answers a distinct scientific question).
# ------------------------------------------------------------------------------
compute_detector_paired_wilcoxon <- function(metric_id, regional_metrics_df, detector_pairs = DETECTOR_PAIRS) {
  metric_long <- regional_metrics_df %>%
    dplyr::transmute(
      REGION = as.character(REGION),
      Method = as.character(Method),
      Score  = .data[[metric_id]]
    )
  metric_wide <- metric_long %>%
    tidyr::pivot_wider(names_from = Method, values_from = Score) %>%
    as.data.frame()

  rows <- list()
  for (pp in detector_pairs) {
    a_vals <- if (pp$A %in% names(metric_wide)) metric_wide[[pp$A]]
              else rep(NA_real_, nrow(metric_wide))
    b_vals <- if (pp$B %in% names(metric_wide)) metric_wide[[pp$B]]
              else rep(NA_real_, nrow(metric_wide))
    paired_keep <- !is.na(a_vals) & !is.na(b_vals)
    a_paired <- a_vals[paired_keep]
    b_paired <- b_vals[paired_keep]
    n_pairs  <- length(a_paired)
    diffs    <- a_paired - b_paired
    median_diff <- if (n_pairs > 0L) stats::median(diffs, na.rm = TRUE)
                   else NA_real_
    abs_median_diff <- if (!is.na(median_diff)) abs(median_diff) else NA_real_

    if (n_pairs >= 3L && any(diffs != 0, na.rm = TRUE)) {
      wt <- suppressWarnings(stats::wilcox.test(
        a_paired, b_paired,
        paired = TRUE,
        alternative = "two.sided",
        exact = FALSE
      ))
      v_stat <- as.numeric(wt$statistic)
      p_val  <- as.numeric(wt$p.value)
    } else {
      v_stat <- NA_real_
      p_val  <- NA_real_
    }
    p_pairwise <- p_val
    p_bonf     <- p_val
    sig        <- if (!is.na(p_pairwise)) p_pairwise < HH_ALPHA else NA

    rows[[length(rows) + 1L]] <- data.frame(
      Metric = metric_id,
      Comparison = pp$short,
      Detector_A = pp$A,
      Detector_B = pp$B,
      N_Regions_Paired = n_pairs,
      Median_Diff_AminusB = median_diff,
      Abs_Median_Diff = abs_median_diff,
      V_statistic = v_stat,
      p_value = p_val,
      p_pairwise = p_pairwise,
      p_bonferroni = p_bonf,
      Significant_005 = sig,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

all_wilcoxon_results <- dplyr::bind_rows(
  lapply(HH_METRICS_REGIONAL, function(m) {
    compute_detector_paired_wilcoxon(m, regional_metrics)
  })
)

cat("\n=== Wilcoxon detector-paired (5 metrics x 3 pairs = 15 rows) ===\n")
print(all_wilcoxon_results, row.names = FALSE, digits = 3)
cat("\n")

# Additional CSV-only head-to-head results: Constant TA vs Outbreak Threshold,
# Farrington, EWARS, and EARS. The same five regional metrics, pairing unit,
# two-sided Wilcoxon test, and alpha = 0.05 are retained from the primary
# regional comparison.
additional_wilcoxon_results <- dplyr::bind_rows(
  lapply(HH_METRICS_REGIONAL, function(m) {
    compute_detector_paired_wilcoxon(
      metric_id = m,
      regional_metrics_df = regional_metrics,
      detector_pairs = ADDITIONAL_DETECTOR_PAIRS
    )
  })
)

cat("\n=== Additional Wilcoxon head-to-head: Constant TA vs OT, Farrington, EWARS, EARS (5 metrics x 4 pairs = 20 rows) ===\n")
print(additional_wilcoxon_results, row.names = FALSE, digits = 3)
cat("\n")

# ------------------------------------------------------------------------------
# 14. REGIONAL DOMINANCE MATRIX (LONG)
# ------------------------------------------------------------------------------
# For each metric, compute per-(region, detector) normalized dominance score
# in [0, 1] using min-max normalization within (region, metric) across the
# 11 detectors. All five reported metrics are higher-is-better, so no
# directional flip is needed.
# ------------------------------------------------------------------------------
rescale_01_safe_local <- function(v) {
  vmin <- min(v, na.rm = TRUE)
  vmax <- max(v, na.rm = TRUE)
  if (!is.finite(vmin) || !is.finite(vmax) || vmax == vmin)
    return(rep(NA_real_, length(v)))
  (v - vmin) / (vmax - vmin)
}

compute_regional_dominance_long <- function() {
  rows <- list()
  for (m in HH_METRICS_REGIONAL) {
    higher_better <- hh_metric_higher_better[[m]]
    sub <- regional_metrics %>%
      dplyr::transmute(
        REGION = as.character(REGION),
        Method = as.character(Method),
        Paradigm = as.character(Paradigm),
        Value = .data[[m]]
      )
    sub <- sub %>%
      dplyr::group_by(REGION) %>%
      dplyr::mutate(
        raw_norm = rescale_01_safe_local(Value),
        Dominance_Score = if (higher_better) raw_norm else 1 - raw_norm
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        Metric = m,
        Is_Dominant = !is.na(Dominance_Score) &
          Dominance_Score >= DOMINANCE_THRESHOLD
      )
    rows[[m]] <- sub
  }
  dplyr::bind_rows(rows)
}
regional_dominance_long <- compute_regional_dominance_long()

# Single-color (blue) cell fill scaled by dominance score.
DOMINANCE_BLUE_LIGHT <- "#FFFFFF"
DOMINANCE_BLUE_DARK  <- "#0B2447"
mix_blue_intensity <- function(score) {
  if (is.na(score)) return("#F2F2F2")
  rgb_lo <- grDevices::col2rgb(DOMINANCE_BLUE_LIGHT) / 255
  rgb_hi <- grDevices::col2rgb(DOMINANCE_BLUE_DARK)  / 255
  r <- rgb_lo[1] * (1 - score) + rgb_hi[1] * score
  g <- rgb_lo[2] * (1 - score) + rgb_hi[2] * score
  b <- rgb_lo[3] * (1 - score) + rgb_hi[3] * score
  grDevices::rgb(r, g, b, maxColorValue = 1)
}
regional_dominance_long$cell_fill <- vapply(
  regional_dominance_long$Dominance_Score, mix_blue_intensity, character(1)
)

# ------------------------------------------------------------------------------
# 15. METHOD-ACROSS-METRICS AGGREGATION
# ------------------------------------------------------------------------------
compute_method_summary_long <- function() {
  regional_dominance_long %>%
    dplyr::filter(Metric %in% HH_METRICS_REGIONAL) %>%
    dplyr::group_by(Method, Metric) %>%
    dplyr::summarise(
      Sweep_Count = sum(Is_Dominant, na.rm = TRUE),
      N_Regions   = sum(!is.na(Dominance_Score)),
      .groups     = "drop"
    ) %>%
    dplyr::mutate(
      Method = factor(as.character(Method), levels = method_order),
      Metric = factor(as.character(Metric), levels = HH_METRICS_REGIONAL)
    )
}
method_summary_long <- compute_method_summary_long()

method_summary_agg <- method_summary_long %>%
  dplyr::group_by(Method) %>%
  dplyr::summarise(
    Mean_Sweep = mean(Sweep_Count, na.rm = TRUE),
    SD_Sweep   = stats::sd(Sweep_Count, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    grand_mean = mean(Mean_Sweep, na.rm = TRUE),
    grand_sd   = stats::sd(Mean_Sweep, na.rm = TRUE),
    Z_Sweep    = ifelse(is.finite(grand_sd) & grand_sd > 0,
                        (Mean_Sweep - grand_mean) / grand_sd,
                        NA_real_),
    Sweep_Class = dplyr::case_when(
      is.na(Z_Sweep)  ~ "Insufficient",
      Z_Sweep >=  1.0 ~ "Sweeper",
      Z_Sweep <= -1.0 ~ "Non-sweeper",
      TRUE            ~ "Average"
    )
  ) %>%
  dplyr::arrange(dplyr::desc(Mean_Sweep))

method_summary_order <- as.character(method_summary_agg$Method)

MAX_REGIONS_POSSIBLE <- {
  .mrp <- suppressWarnings(max(method_summary_long$N_Regions, na.rm = TRUE))
  if (!is.finite(.mrp) || .mrp <= 0) 17L else as.integer(.mrp)
}

mix_blue_count <- function(count, max_count = MAX_REGIONS_POSSIBLE) {
  if (is.na(count) || !is.finite(count)) return("#FFFFFF")
  frac <- pmax(0, pmin(1, count / max(1, max_count)))
  mix_blue_intensity(frac)
}
method_summary_long$cell_fill <- vapply(
  method_summary_long$Sweep_Count, mix_blue_count, character(1)
)

# ------------------------------------------------------------------------------
# 16. FIGURE 4 BUILDERS — PER-METRIC PANELS AND METHOD SUMMARY
# ------------------------------------------------------------------------------

# 16a. Regional dominance matrix (per metric) ---------------------------------
build_regional_dominance_matrix <- function(metric_id) {
  metric_label <- hh_metric_full_name[[metric_id]]
  this_metric_df <- regional_dominance_long %>%
    dplyr::filter(Metric == metric_id)

  # Per-method count of regions swept at score >= threshold.
  per_method_count <- this_metric_df %>%
    dplyr::group_by(Method) %>%
    dplyr::summarise(
      N_Regions_Swept = sum(Is_Dominant, na.rm = TRUE),
      N_Regions_Total = sum(!is.na(Dominance_Score)),
      Mean_Score      = mean(Dominance_Score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Sweep_Label = paste0(N_Regions_Swept, "/", N_Regions_Total))

  method_sweep_order <- per_method_count %>%
    dplyr::arrange(dplyr::desc(N_Regions_Swept),
                   dplyr::desc(Mean_Score)) %>%
    dplyr::pull(Method) %>% as.character()

  region_score_order <- this_metric_df %>%
    dplyr::group_by(REGION) %>%
    dplyr::summarise(mean_score = mean(Dominance_Score, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_score)) %>%
    dplyr::pull(REGION) %>% as.character()

  this_metric_df <- this_metric_df %>%
    dplyr::mutate(
      REGION = factor(as.character(REGION), levels = region_score_order),
      Method = factor(as.character(Method), levels = rev(method_sweep_order))
    )
  per_method_count <- per_method_count %>%
    dplyr::mutate(Method = factor(as.character(Method),
                                  levels = rev(method_sweep_order)))

  N_REGIONS     <- length(region_score_order)
  SWEEP_COUNT_X <- N_REGIONS + 1L

  ggplot2::ggplot(this_metric_df) +
    ggplot2::geom_tile(
      ggplot2::aes(x = REGION, y = Method),
      fill = this_metric_df$cell_fill,
      colour = "white", linewidth = 0.5, show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = per_method_count,
      ggplot2::aes(x = SWEEP_COUNT_X, y = Method, label = Sweep_Label),
      size = 2.6, fontface = "bold", family = base_family_global,
      colour = "#0B2447", inherit.aes = FALSE, show.legend = FALSE
    ) +
    ggplot2::geom_vline(xintercept = SWEEP_COUNT_X - 0.5,
                        linetype = "solid",
                        linewidth = 0.50, colour = "grey25") +
    ggplot2::annotate(
      "text", x = SWEEP_COUNT_X, y = length(method_sweep_order) + 1.20,
      label = paste0("Regions Swept\n(score \u2265 ",
                     sprintf("%.2f", DOMINANCE_THRESHOLD), ")"),
      size = 2.2, family = base_family_global, fontface = "bold",
      colour = "#0B2447", lineheight = 0.92, vjust = 0
    ) +
    ggplot2::scale_x_discrete(
      limits = c(region_score_order, "__SWEEP_COUNT__"),
      labels = {
        raw <- as.character(c(region_score_order, " "))
        ifelse(is.na(raw) | raw == "", " ", raw)
      },
      position = "top",
      expand = ggplot2::expansion(add = c(0.04, 0.50))
    ) +
    ggplot2::scale_y_discrete(
      limits = rev(method_sweep_order),
      labels = {
        raw <- as.character(unname(method_two_line[rev(method_sweep_order)]))
        ifelse(is.na(raw) | raw == "", " ", raw)
      },
      expand = ggplot2::expansion(add = c(0.30, 1.45))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = paste0("a", hh_metric_subscript[[metric_id]],
                     "    Regional dominance matrix on ", metric_label,
                     " (rows = methods, columns = regions)"),
      x = NULL, y = NULL
    ) +
    theme_dashboard(base_size = PUB_BASE, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(color = "black", size = 6.8,
                                              hjust = 0.5, vjust = 0,
                                              margin = ggplot2::margin(b = 4)),
      axis.text.y     = ggplot2::element_text(color = "black", size = 6.8,
                                              lineheight = 0.85,
                                              hjust = 1, vjust = 0.5,
                                              margin = ggplot2::margin(r = 4)),
      panel.grid    = ggplot2::element_blank(),
      panel.border  = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin   = ggplot2::margin(40, 14, 8, 8)
    )
}

# 16b. Dominance score legend strip -------------------------------------------
build_dominance_legend <- function() {
  legend_score_seq <- seq(0, 1, length.out = 100L)
  legend_df <- data.frame(
    x = legend_score_seq,
    y = 1L,
    fill_col = vapply(legend_score_seq, mix_blue_intensity, character(1))
  )
  ggplot2::ggplot(legend_df) +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y),
                       fill = legend_df$fill_col,
                       width = 1 / nrow(legend_df),
                       height = 1) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 0.25, 0.50, DOMINANCE_THRESHOLD, 1.00),
      labels = c("0.00\nweak", "0.25", "0.50\nmoderate",
                 sprintf("%.2f\ndominant", DOMINANCE_THRESHOLD),
                 "1.00\nsweep"),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      title = paste0("Dominance score (min-max normalized within region; ",
                     "darker = stronger dominance)"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_void(base_family = base_family_global) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 7.5,
                                          hjust = 0,
                                          margin = ggplot2::margin(b = 4)),
      axis.text.x = ggplot2::element_text(size = PUB_AXIS_TXT, lineheight = 0.95,
                                          margin = ggplot2::margin(t = 2)),
      plot.margin = ggplot2::margin(10, 12, 8, 8)
    )
}

# 16c. Per-detector dot plot --------------------------------------------------
build_regional_dot_panel <- function(metric_id) {
  metric_label  <- hh_metric_full_name[[metric_id]]
  ylab          <- hh_metric_y_label[[metric_id]]
  higher_better <- hh_metric_higher_better[[metric_id]]
  threshold     <- hh_metric_threshold[[metric_id]]
  is_proportion <- hh_metric_is_proportion[[metric_id]]
  ci_cols       <- hh_metric_ci_cols[[metric_id]]

  plot_df <- regional_metrics %>%
    dplyr::transmute(
      Method, REGION,
      Value = .data[[metric_id]]
    ) %>%
    dplyr::mutate(REGION_chr = as.character(REGION),
                  Method_chr = as.character(Method))

  if (all(ci_cols %in% names(boot_ci))) {
    plot_df <- plot_df %>%
      dplyr::left_join(
        boot_ci %>%
          dplyr::transmute(
            REGION_chr = as.character(REGION),
            Method_chr = as.character(Method),
            ci_lo = .data[[ci_cols[1]]],
            ci_hi = .data[[ci_cols[2]]]
          ),
        by = c("REGION_chr", "Method_chr")
      )
  } else {
    plot_df$ci_lo <- NA_real_
    plot_df$ci_hi <- NA_real_
  }

  fav_layer <- NULL
  ref_line  <- NULL
  if (!is.na(threshold)) {
    if (higher_better) {
      fav_layer <- ggplot2::annotate(
        "rect", xmin = -Inf, xmax = Inf,
        ymin = threshold, ymax = Inf,
        fill = "#cfe2f3", alpha = 0.30
      )
    } else {
      fav_layer <- ggplot2::annotate(
        "rect", xmin = -Inf, xmax = Inf,
        ymin = -Inf, ymax = threshold,
        fill = "#cfe2f3", alpha = 0.30
      )
    }
    ref_line <- ggplot2::geom_hline(
      yintercept = threshold, linetype = "dashed",
      linewidth = 0.35, colour = "grey45"
    )
  }

  pj <- ggplot2::position_jitter(width = 0.18, height = 0, seed = 12345)
  p  <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Method, y = Value))
  if (!is.null(fav_layer)) p <- p + fav_layer
  if (!is.null(ref_line))  p <- p + ref_line

  p <- p +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      position = pj,
      colour = "grey55", alpha = 0.55, linewidth = 0.30, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = Method),
      shape = 21, position = pj,
      colour = "grey20", stroke = 0.30, size = 2.2, alpha = 0.92,
      show.legend = FALSE, na.rm = TRUE
    ) +
    ggplot2::stat_summary(
      fun = mean, geom = "crossbar",
      width = 0.50, linewidth = 0.40,
      colour = "black", fatten = 0,
      fill = NA, na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(rep("#FFFFFF", length(method_order)),
                               as.character(method_order)),
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      limits = method_order,
      labels = {
        raw <- as.character(unname(method_two_line[as.character(method_order)]))
        ifelse(is.na(raw) | raw == "", " ", raw)
      }
    ) +
    ggplot2::labs(
      title = paste0("b", hh_metric_subscript[[metric_id]],
                     "    Per-detector regional values on ", metric_label),
      x = NULL, y = ylab
    ) +
    theme_dashboard(base_size = PUB_BASE, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = PUB_AXIS_TXT, lineheight = 0.85,
                                          hjust = 0.5, vjust = 1,
                                          margin = ggplot2::margin(t = 2)),
      axis.text.y = ggplot2::element_text(size = PUB_AXIS_TXT),
      axis.title.y = ggplot2::element_text(size = PUB_AXIS_TIT,
                                           margin = ggplot2::margin(r = 4)),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 12, 22, 8)
    )

  if (is_proportion) {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.05),
      expand = ggplot2::expansion(mult = c(0, 0))
    )
  }
  p
}

# 16d. Wilcoxon detector-paired sig bars --------------------------------------
build_regional_sig_panel <- function(metric_id, all_wilcoxon_results) {
  res <- all_wilcoxon_results %>%
    dplyr::filter(Metric == metric_id) %>%
    dplyr::mutate(Comparison = factor(Comparison,
                                      levels = vapply(DETECTOR_PAIRS,
                                                      function(p) p$short,
                                                      character(1))))
  res <- res %>%
    dplyr::mutate(
      p_label = ifelse(
        is.na(p_bonferroni),
        "n/a",
        ifelse(p_bonferroni < 0.001,
               "p < 0.001",
               paste0("p = ", sprintf("%.3f", p_bonferroni)))
      ),
      sig_star = ifelse(!is.na(Significant_005) & Significant_005, "*", ""),
      bar_fill = ifelse(!is.na(Significant_005) & Significant_005,
                        "#1F2D5C", "#C7CCD9")
    )

  x_max <- max(res$Abs_Median_Diff, na.rm = TRUE)
  if (!is.finite(x_max) || x_max <= 0) x_max <- 1
  x_lim_hi <- x_max * 1.45

  ggplot2::ggplot(res, ggplot2::aes(y = Comparison, x = Abs_Median_Diff)) +
    ggplot2::geom_col(
      fill = res$bar_fill,
      width = 0.55, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(p_label, " ", sig_star),
        x = ifelse(is.na(Abs_Median_Diff),
                   x_lim_hi * 0.05,
                   pmin(Abs_Median_Diff + x_lim_hi * 0.04, x_lim_hi * 0.95))
      ),
      hjust = 0,
      family = base_family_global, size = pub_text_size(PUB_ANNOT),
      fontface = "bold",
      colour = "black", na.rm = TRUE   # p-values in bold black
    ) +
    # CROPPING FIX. The p-value labels are drawn hjust = 0 starting as far right
    # as 0.95 * x_lim_hi, so with zero expansion their text ran past the panel
    # edge and was cut. The right-hand expansion reserves room for them.
    ggplot2::scale_x_continuous(
      limits = c(0, x_lim_hi),
      expand = ggplot2::expansion(mult = c(0, 0.30))
    ) +
    ggplot2::scale_y_discrete(
      limits = rev(vapply(DETECTOR_PAIRS, function(p) p$short, character(1)))
    ) +
    ggplot2::labs(
      title = paste0("c", hh_metric_subscript[[metric_id]],
                     "    Wilcoxon detector-paired (n = ",
                     {
                       n_max <- suppressWarnings(max(res$N_Regions_Paired,
                                                     na.rm = TRUE))
                       if (!is.finite(n_max)) "?" else as.character(n_max)
                     },
                     " regions)"),
      x = "|Median diff|",
      y = NULL
    ) +
    theme_dashboard(base_size = PUB_BASE - 0.5, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = PUB_AXIS_TXT, face = "bold",
                                           lineheight = 0.95),
      axis.text.x  = ggplot2::element_text(size = PUB_AXIS_TXT),
      axis.title.x = ggplot2::element_text(size = PUB_AXIS_TIT),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "none"
    )
}

# 16e. Compose one per-metric multipanel figure -------------------------------
build_regional_multipanel <- function(metric_id) {
  top          <- build_regional_dominance_matrix(metric_id)
  legend_strip <- build_dominance_legend()
  bot_l        <- build_regional_dot_panel(metric_id)
  bot_r        <- build_regional_sig_panel(metric_id, all_wilcoxon_results)
  bottom       <- (bot_l | bot_r) + patchwork::plot_layout(widths = c(7, 3))

  (top / legend_strip / bottom) +
    patchwork::plot_layout(heights = c(1.55, 0.18, 1.00)) +
    patchwork::plot_annotation(
      title = paste0("Regional generalisability: ",
                     hh_metric_full_name[[metric_id]]),
      # Panel letters are INLINE in each sub-panel's own title
      # ("a<sub>metric</sub>    Regional dominance matrix on ...").
      # No patchwork tag is used: tag_levels would also letter the legend
      # strip, and a separate tag would print a second letter above the title.
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = PUB_TITLE + 1,
                                           family = base_family_global,
                                           margin = ggplot2::margin(b = 6)),
        plot.tag = ggplot2::element_text(size = PUB_TAG, face = "bold",
                                         colour = "black",
                                         family = base_family_global,
                                       hjust = 0, vjust = 1),
      # Same left edge as the title: tags anchor to the plot, and
      # plot.title.position = "plot" moves the title off the panel edge to match.
      plot.tag.position     = "topleft",
      plot.title.position   = "plot",
        plot.margin = ggplot2::margin(10, 12, 8, 8)
      )
    ) &
    # Applied to EVERY sub-panel with `&`. plot.title.position = "plot" keeps
    # each inline-lettered title flush with the plot's left edge rather than
    # indented past the y-axis labels.
    ggplot2::theme(
      plot.margin           = ggplot2::margin(10, 12, 8, 8),
      plot.tag              = ggplot2::element_text(size = PUB_TAG,
                                                    face = "bold",
                                                    colour = "black",
                                                    family = base_family_global,
                                                    hjust = 0, vjust = 1),
      plot.tag.position     = "topleft",
      plot.title.position   = "plot"
    )
}

# 16f. Method-across-metrics summary panels -----------------------------------
build_method_summary_matrix <- function() {
  df <- method_summary_long %>%
    dplyr::mutate(
      Method = factor(as.character(Method), levels = rev(method_summary_order)),
      Metric = factor(as.character(Metric), levels = HH_METRICS_REGIONAL)
    )
  agg <- method_summary_agg %>%
    dplyr::mutate(
      Method = factor(as.character(Method), levels = rev(method_summary_order)),
      Summary_Label = sprintf("%.1f (z=%+.2f)", Mean_Sweep, Z_Sweep)
    )

  N_METRICS <- length(HH_METRICS_REGIONAL)
  SUMMARY_X <- N_METRICS + 1L
  N_METHODS <- length(method_summary_order)
  HEADER_Y  <- N_METHODS + 1L

  metric_full_names <- vapply(HH_METRICS_REGIONAL,
                              function(m) hh_metric_full_name[[m]],
                              character(1))
  metric_cell_labels <- c(
    "TAM"            = "True-Alarm\nMagnitude",
    "N_True_Alarms"  = "Number of\nTrue Alarms",
    "Sensitivity"    = "Sensitivity",
    "Mean_Lead_Time" = "Mean\nLead Time",
    "WP"             = "Warning\nPersistence"
  )
  for (m in HH_METRICS_REGIONAL) {
    if (is.na(metric_cell_labels[m]) || is.null(metric_cell_labels[[m]])) {
      metric_cell_labels[m] <- metric_full_names[[m]]
    }
  }

  header_df <- data.frame(
    x_pos = c(seq_len(N_METRICS), SUMMARY_X),
    y_pos = HEADER_Y,
    label = c(unname(metric_cell_labels[HH_METRICS_REGIONAL]),
              "Mean Sweep\n(z-score)"),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(df) +
    ggplot2::geom_tile(
      ggplot2::aes(x = as.integer(Metric), y = as.integer(Method)),
      fill = df$cell_fill, colour = "white", linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = as.integer(Metric), y = as.integer(Method),
                   label = sprintf("%d", Sweep_Count),
                   colour = Sweep_Count >= round(MAX_REGIONS_POSSIBLE * 0.55)),
      size = 2.6, family = base_family_global, fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "grey15", "TRUE" = "white"),
      guide  = "none"
    ) +
    ggplot2::geom_text(
      data = agg,
      ggplot2::aes(x = SUMMARY_X, y = as.integer(Method), label = Summary_Label),
      size = 2.5, fontface = "bold", family = base_family_global,
      colour = "#0B2447", inherit.aes = FALSE
    ) +
    # Header background tiles removed (no background shading).
    ggplot2::geom_text(
      data = header_df,
      ggplot2::aes(x = x_pos, y = y_pos, label = label),
      size = 2.4, family = base_family_global, fontface = "bold",
      colour = "#0B2447", lineheight = 0.85, inherit.aes = FALSE
    ) +
    ggplot2::geom_vline(xintercept = SUMMARY_X - 0.5,
                        linetype = "solid",
                        linewidth = 0.50, colour = "grey25") +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      limits = c(0.5, SUMMARY_X + 0.5),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(seq_len(N_METHODS), HEADER_Y),
      labels = c(unname(method_two_line[rev(method_summary_order)]), ""),
      limits = c(0.5, HEADER_Y + 0.5),
      expand = ggplot2::expansion(add = c(0.05, 0.30))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Regional dominance matrix on 5 operational metrics",
      x = NULL, y = NULL
    ) +
    theme_dashboard(base_size = PUB_BASE, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_blank(),
      axis.text.x     = ggplot2::element_blank(),
      axis.ticks.x    = ggplot2::element_blank(),
      axis.text.y     = ggplot2::element_text(color = "black", size = 6.8,
                                              lineheight = 0.85,
                                              hjust = 1, vjust = 0.5,
                                              margin = ggplot2::margin(r = 4)),
      panel.grid    = ggplot2::element_blank(),
      panel.border  = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin   = ggplot2::margin(20, 14, 8, 8)
    )
}

build_method_summary_dot <- function() {
  df_dots <- method_summary_long %>%
    dplyr::mutate(
      Method = factor(as.character(Method), levels = rev(method_summary_order))
    )
  df_means <- method_summary_agg %>%
    dplyr::mutate(
      Method = factor(as.character(Method), levels = rev(method_summary_order))
    )

  ggplot2::ggplot() +
    ggplot2::geom_jitter(
      data = df_dots,
      ggplot2::aes(x = Sweep_Count, y = Method),
      width = 0, height = 0.18,
      shape = 21, fill = "#3673B6", colour = "white",
      stroke = 0.4, size = 2.4, alpha = 0.85
    ) +
    ggplot2::geom_segment(
      data = df_means,
      ggplot2::aes(x = Mean_Sweep, xend = Mean_Sweep,
                   y    = as.numeric(Method) - 0.30,
                   yend = as.numeric(Method) + 0.30),
      colour = "#0B2447", linewidth = 1.0
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, MAX_REGIONS_POSSIBLE + 0.5),
      breaks = pretty(c(0, MAX_REGIONS_POSSIBLE), n = 5),
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    ) +
    ggplot2::scale_y_discrete(
      limits = rev(method_summary_order),
      labels = unname(method_two_line[rev(method_summary_order)])
    ) +
    ggplot2::labs(
      title = "Per-method sweep counts across the 5 metrics",
      x = "Regions swept (count)", y = NULL
    ) +
    theme_dashboard(base_size = PUB_BASE, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black", size = 6.8,
                                          lineheight = 0.85,
                                          hjust = 1, vjust = 0.5,
                                          margin = ggplot2::margin(r = 4)),
            panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(10, 12, 8, 8)
    )
}

# Composite Wilcoxon panel for the method summary figure (DUAL-STATISTIC).
#
#   p-value      : paired Wilcoxon signed-rank on within-region Dominance_Score
#                  paired by (REGION x METRIC). n_pairs ~ 17 x 5 = 85.
#                  Properly powered; p < 0.05 is reachable.
#   Bar length   : |median(diff in regions-swept count)| across the 5 reported
#                  metrics. Interpretable in operational units of "regions".
#
# The dual specification answers two distinct questions in the same panel:
# (i) is the difference statistically real (n=85 powered test); and (ii) how
# big is the practical difference in operational regions (n=5 sweep-count
# median). Both numbers are also written to the console for auditability.
build_method_summary_sig <- function() {
  TARGETS <- c("Constant Transmission Acceleration",
               "Continuous Transmission Acceleration",
               "Outbreak Threshold")
  PAIRS_LIST <- list(
    list(a = "Constant Transmission Acceleration",
         b = "Continuous Transmission Acceleration",
         label = "Constant TA vs Continuous TA"),
    list(a = "Constant Transmission Acceleration",
         b = "Outbreak Threshold",
         label = "Constant TA vs Outbreak Threshold"),
    list(a = "Continuous Transmission Acceleration",
         b = "Outbreak Threshold",
         label = "Continuous TA vs Outbreak Threshold")
  )

  pivot_score <- regional_dominance_long %>%
    dplyr::filter(as.character(Method) %in% TARGETS,
                  as.character(Metric) %in% HH_METRICS_REGIONAL) %>%
    dplyr::select(REGION, Metric, Method, Dominance_Score) %>%
    dplyr::mutate(Method = as.character(Method)) %>%
    tidyr::pivot_wider(names_from = Method, values_from = Dominance_Score)

  pivot_count <- method_summary_long %>%
    dplyr::filter(as.character(Method) %in% TARGETS,
                  as.character(Metric) %in% HH_METRICS_REGIONAL) %>%
    dplyr::select(Method, Metric, Sweep_Count) %>%
    dplyr::mutate(Method = as.character(Method)) %>%
    tidyr::pivot_wider(names_from = Method, values_from = Sweep_Count)

  rows <- list()
  for (pr in PAIRS_LIST) {
    a_score <- pivot_score[[pr$a]]; b_score <- pivot_score[[pr$b]]
    keep_s  <- !is.na(a_score) & !is.na(b_score)
    a_score <- a_score[keep_s]; b_score <- b_score[keep_s]
    if (length(a_score) >= 2L && any(a_score != b_score)) {
      tt <- suppressWarnings(stats::wilcox.test(a_score, b_score,
                                                paired = TRUE,
                                                exact = FALSE))
      pval           <- tt$p.value
      med_score_diff <- stats::median(a_score - b_score, na.rm = TRUE)
    } else {
      pval           <- NA_real_
      med_score_diff <- NA_real_
    }
    n_score_pairs <- length(a_score)

    a_count <- pivot_count[[pr$a]]; b_count <- pivot_count[[pr$b]]
    keep_c  <- !is.na(a_count) & !is.na(b_count)
    a_count <- a_count[keep_c]; b_count <- b_count[keep_c]
    med_count_diff <- if (length(a_count) >= 1L) {
      stats::median(a_count - b_count, na.rm = TRUE)
    } else {
      NA_real_
    }

    rows[[length(rows) + 1L]] <- data.frame(
      Pair               = pr$label,
      a_method           = pr$a,
      b_method           = pr$b,
      Median_Score_Diff  = med_score_diff,
      Median_Count_Diff  = med_count_diff,
      P_Value            = pval,
      N_Score_Pairs      = n_score_pairs,
      N_Metric_Pairs     = length(a_count),
      Sig                = !is.na(pval) & pval < 0.05,
      stringsAsFactors   = FALSE
    )
  }
  res <- do.call(rbind, rows)
  res$Pair_short <- factor(res$Pair, levels = res$Pair)

  res$Bar_Mag_raw <- ifelse(is.na(res$Median_Count_Diff), 0,
                            abs(res$Median_Count_Diff))
  res$Sig_plot    <- ifelse(is.na(res$Sig), FALSE, res$Sig)
  MIN_VISIBLE_BAR <- 0.5
  res$Bar_Mag <- ifelse(res$Sig_plot & res$Bar_Mag_raw < MIN_VISIBLE_BAR,
                        MIN_VISIBLE_BAR, res$Bar_Mag_raw)

  cat("\n--- Method Summary Panel c: dual-statistic Wilcoxon ---\n")
  cat("    Test: paired Wilcoxon on Dominance_Score (17 REGIONS x 5 METRICS, n ~ 85)\n")
  cat("    Bar:  |median(diff in sweep counts)| across the 5 metrics, units 'regions'\n")
  print(res[, c("Pair", "Median_Score_Diff", "Median_Count_Diff",
                "P_Value", "N_Score_Pairs", "N_Metric_Pairs", "Sig")],
        row.names = FALSE)

  n_typical <- {
    .nm <- suppressWarnings(median(res$N_Score_Pairs, na.rm = TRUE))
    if (!is.finite(.nm)) "?" else as.character(round(.nm))
  }
  res$Sig_star <- ifelse(res$Sig_plot, "*", "")

  ggplot2::ggplot(res) +
    ggplot2::geom_col(
      ggplot2::aes(x = Bar_Mag, y = Pair_short, fill = Sig_plot),
      width = 0.65
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = Bar_Mag, y = Pair_short,
                   label = paste0(
                     ifelse(is.na(P_Value), "n/a",
                            ifelse(P_Value < 0.001,
                                   "p < 0.001",
                                   sprintf("p = %.3g", P_Value))),
                     " ", Sig_star
                   )),
      hjust = -0.10, size = 2.5, family = base_family_global,
      fontface = "bold", colour = "grey15"
    ) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "#BFC6CE", "TRUE" = "#3673B6"),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.30))
    ) +
    ggplot2::labs(
      title = paste0("c    Wilcoxon detector-paired (n \u2248 ",
                     n_typical, ", 17 regions \u00D7 5 metrics)"),
      x = "Median diff in regions swept (across 5 metrics)",
      y = NULL
    ) +
    theme_dashboard(base_size = PUB_BASE, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black", size = 7.0,
                                          hjust = 1, vjust = 0.5,
                                          margin = ggplot2::margin(r = 4)),
            panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 14, 8, 8)
    )
}

build_method_summary_multipanel <- function() {
  top          <- build_method_summary_matrix()
  legend_strip <- build_dominance_legend()
  bot_l        <- build_method_summary_dot()
  bot_r        <- build_method_summary_sig()
  bottom       <- (bot_l | bot_r) + patchwork::plot_layout(widths = c(7, 3))

  (top / legend_strip / bottom) +
    patchwork::plot_layout(heights = c(1.50, 0.18, 1.00)) +
    patchwork::plot_annotation(
      title = paste0("Figure 4 | Method-across-metrics summary ",
                     "(11 methods x 5 reported metrics)"),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = PUB_TITLE + 1,
                                           family = base_family_global,
                                           margin = ggplot2::margin(b = 6)),
        plot.margin = ggplot2::margin(10, 12, 8, 8)
      )
    ) &
    # Applied to EVERY sub-panel with `&`. plot.title.position = "plot" keeps
    # each inline-lettered title flush with the plot's left edge rather than
    # indented past the y-axis labels.
    ggplot2::theme(
      plot.margin           = ggplot2::margin(10, 12, 8, 8),
      plot.tag              = ggplot2::element_text(size = PUB_TAG,
                                                    face = "bold",
                                                    colour = "black",
                                                    family = base_family_global,
                                                    hjust = 0, vjust = 1),
      plot.tag.position     = "topleft",
      plot.title.position   = "plot"
    )
}

# ------------------------------------------------------------------------------
# 17. SAVE FIGURE 4 (PER-METRIC + METHOD SUMMARY)
# ------------------------------------------------------------------------------
# Figure 4 set, authored at final print size. Previously 16 x 13.5 in, which
# production reduced by ~0.44x, printing axis text near 3.5 pt.
fig4_w <- NC_W_DOUBLE
# Composite stacks a dominance matrix, a legend strip and two sub-panels;
# the taller supplementary canvas prevents the rows from compressing.
fig4_h <- NC_H_SUPP

cat("\n=== Building 5 per-metric multipanel figures (Figure 4 set) ===\n")
for (m in HH_METRICS_REGIONAL) {
  fig <- build_regional_multipanel(m)
  slug <- hh_metric_slug[[m]]
  save_plot_file(paste0("Figure4_", slug), fig, fig4_w, fig4_h)
}

cat("\n=== Building Method-Across-Metrics Summary figure (Figure 4 set) ===\n")
{
  fig_summary <- build_method_summary_multipanel()
  save_plot_file("Figure4_Method_Summary", fig_summary, fig4_w, fig4_h)
}

# ------------------------------------------------------------------------------
# 18. WRITE FIGURE 4 CSV TABLES
# ------------------------------------------------------------------------------
# Wide regional metrics table (carries the 5 reported metrics plus diagnostic
# columns). RFR / RFR_RX columns dropped (not used by this pipeline).
utils::write.csv(
  regional_metrics %>% dplyr::mutate(REGION = as.character(REGION),
                                     Method = as.character(Method)),
  file.path(OUTPUT_DIR, "Stage4_Regional_Framework_Metrics.csv"),
  row.names = FALSE
)

# 8-metric harmonized regional CSV.
regional_8metric_csv <- regional_metrics %>%
  dplyr::mutate(REGION = as.character(REGION),
                Method = as.character(Method),
                Paradigm = as.character(Paradigm)) %>%
  dplyr::transmute(
    REGION, Method, Paradigm,
    TAM, N_True_Alarms, PPV, Sensitivity,
    Mean_Lead_Time, WP, ALY,
    Mean_Lead_Time_conditional, WP_conditional,
    N_False_Alarms,
    Total_Triggers, True_Alarms, False_Alarms,
    n_TrueActionable, n_Reactive,
    N_Years_with_A1_True, N_Years_Evaluable
  )
utils::write.csv(
  regional_8metric_csv,
  file.path(OUTPUT_DIR, "Stage4_Regional_8Metric_Summary.csv"),
  row.names = FALSE
)

# 8-metric CSV with bootstrap 95% CIs joined per (REGION, Method).
regional_8metric_with_ci <- regional_8metric_csv %>%
  dplyr::left_join(
    boot_ci %>% dplyr::mutate(REGION = as.character(REGION),
                              Method = as.character(Method)) %>%
      dplyr::select(REGION, Method,
                    TAM_lo, TAM_hi,
                    N_True_Alarms_lo, N_True_Alarms_hi,
                    PPV_lo, PPV_hi,
                    Sens_lo, Sens_hi,
                    MLT_lo, MLT_hi,
                    WP_lo, WP_hi,
                    ALY_lo, ALY_hi,
                    N_False_Alarms_lo, N_False_Alarms_hi),
    by = c("REGION", "Method")
  )
utils::write.csv(
  regional_8metric_with_ci,
  file.path(OUTPUT_DIR, "Stage4_Regional_8Metric_Summary_with_CIs.csv"),
  row.names = FALSE
)

# Wide regional metrics table with bootstrap CIs.
metrics_with_ci <- regional_metrics %>%
  dplyr::mutate(REGION = as.character(REGION),
                Method = as.character(Method)) %>%
  dplyr::left_join(
    boot_ci %>%
      dplyr::mutate(REGION = as.character(REGION),
                    Method = as.character(Method)),
    by = c("REGION", "Method")
  )
utils::write.csv(
  metrics_with_ci,
  file.path(OUTPUT_DIR, "Stage4_Regional_Framework_Metrics_with_CIs.csv"),
  row.names = FALSE
)

# Dominance probability (region-level).
utils::write.csv(
  dominance_df,
  file.path(OUTPUT_DIR, "Stage4_Regional_Dominance_Probabilities.csv"),
  row.names = FALSE
)

# Regional Dominance Matrix (long).
utils::write.csv(
  regional_dominance_long %>%
    dplyr::transmute(
      REGION = as.character(REGION),
      Method = as.character(Method),
      Paradigm = as.character(Paradigm),
      Metric = Metric,
      Value = Value,
      Dominance_Score = Dominance_Score,
      Is_Dominant_at_threshold_075 = Is_Dominant
    ),
  file.path(OUTPUT_DIR, "Stage4_Regional_Dominance_Matrix.csv"),
  row.names = FALSE
)

# Wilcoxon per-metric detector-paired results.
utils::write.csv(
  all_wilcoxon_results,
  file.path(OUTPUT_DIR, "Stage4_Regional_Wilcoxon_PerMetric.csv"),
  row.names = FALSE
)

# Additional cross-region head-to-head comparisons requested for CSV output.
# Kept separate so the existing Figure 4 and Figure 5 inference objects remain
# unchanged.
utils::write.csv(
  additional_wilcoxon_results,
  file.path(OUTPUT_DIR, "Stage4_Regional_Wilcoxon_ConstantTA_vs_Comparators.csv"),
  row.names = FALSE
)

# Method-across-metrics summary CSVs.
utils::write.csv(
  method_summary_long %>%
    dplyr::mutate(
      Method = as.character(Method),
      Metric = as.character(Metric)
    ) %>%
    dplyr::select(Method, Metric, Sweep_Count, N_Regions),
  file.path(OUTPUT_DIR, "Stage4_Method_Summary_Long.csv"),
  row.names = FALSE
)
utils::write.csv(
  method_summary_agg %>%
    dplyr::mutate(Method = as.character(Method)) %>%
    dplyr::select(Method, Mean_Sweep, SD_Sweep, Z_Sweep, Sweep_Class),
  file.path(OUTPUT_DIR, "Stage4_Method_Summary_Aggregate.csv"),
  row.names = FALSE
)

# Bootstrap replicates (long; one row per REGION x .rep x Method).
all_reps_df <- dplyr::bind_rows(
  lapply(names(boot_results), function(reg) {
    reps <- boot_results[[reg]]$replicates
    reps$REGION <- reg
    reps
  })
) %>%
  dplyr::select(REGION, .rep, Method,
                TAM, N_True_Alarms, Sensitivity,
                Mean_Lead_Time, WP, PPV, ALY, N_False_Alarms)
utils::write.csv(
  all_reps_df,
  file.path(OUTPUT_DIR, "Stage4_Regional_Bootstrap_Replicates.csv"),
  row.names = FALSE
)
cat("Saved bootstrap replicates: ",
    file.path(OUTPUT_DIR, "Stage4_Regional_Bootstrap_Replicates.csv"),
    "\n", sep = "")

# ==============================================================================
# FIGURE 5 — DETECTOR MAP, METRIC TABLE, AND DOT PLOT
# ==============================================================================
# Composes a single composite figure summarising regional dominance:
#   a. Choropleth detector map of the 17 administrative regions, encoding the
#      consensus winner under a four-tier classification (strong, partial,
#      lead_only, contested) with Bonferroni-adjusted head-to-head testing.
#   b. Per-region metric table with embedded per-metric significance below
#      each cell value.
#   c. Per-detector dot plot showing the bootstrap dominance probabilities
#      for each region and detector.
#
# Per-region per-metric significance is also computed for six ORDERED
# detector pairs (each ordered direction is its own one-sided test of the
# explicit signed difference X - Y). The full 6 x 5 x 17 table is saved as
# Stage4_Detector_Map_PanelD_HeatmapData.csv for downstream reporting; it
# is NOT rendered as a figure panel.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 19. REGION CANONICALISER
# ------------------------------------------------------------------------------
.normalize_key <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s <- gsub("[().,]", " ", s)
  s <- gsub("\\s+", " ", s)
  trimws(s)
}

region_recoder <- c(
  "NATIONAL CAPITAL REGION" = "NCR", "METROPOLITAN MANILA" = "NCR",
  "METRO MANILA" = "NCR", "MANILA" = "NCR", "NCR" = "NCR",
  "CORDILLERA ADMINISTRATIVE REGION" = "CAR",
  "CORDILLERA ADMINISTRATIVE REGION CAR" = "CAR",
  "CORDILLERA" = "CAR", "CAR" = "CAR",
  "ILOCOS REGION" = "REGION I", "ILOCOS" = "REGION I",
  "REGION 1" = "REGION I", "REGION I" = "REGION I",
  "CAGAYAN VALLEY" = "REGION II", "REGION 2" = "REGION II",
  "REGION II" = "REGION II",
  "CENTRAL LUZON" = "REGION III", "REGION 3" = "REGION III",
  "REGION III" = "REGION III",
  "CALABARZON" = "REGION IV-A", "REGION 4-A" = "REGION IV-A",
  "REGION 4A" = "REGION IV-A", "REGION IV-A" = "REGION IV-A",
  "REGION IVA" = "REGION IV-A", "REGION IV A" = "REGION IV-A",
  "MIMAROPA REGION" = "MIMAROPA", "MIMAROPA" = "MIMAROPA",
  "REGION 4-B" = "MIMAROPA", "REGION 4B" = "MIMAROPA",
  "REGION IV-B" = "MIMAROPA", "REGION IV B" = "MIMAROPA",
  "SOUTHWESTERN TAGALOG REGION" = "MIMAROPA",
  "BICOL REGION" = "REGION V", "BICOL" = "REGION V",
  "REGION 5" = "REGION V", "REGION V" = "REGION V",
  "WESTERN VISAYAS" = "REGION VI", "REGION 6" = "REGION VI",
  "REGION VI" = "REGION VI",
  "CENTRAL VISAYAS" = "REGION VII", "REGION 7" = "REGION VII",
  "REGION VII" = "REGION VII",
  "EASTERN VISAYAS" = "REGION VIII", "REGION 8" = "REGION VIII",
  "REGION VIII" = "REGION VIII",
  "ZAMBOANGA PENINSULA" = "REGION IX", "REGION 9" = "REGION IX",
  "REGION IX" = "REGION IX",
  "NORTHERN MINDANAO" = "REGION X", "REGION 10" = "REGION X",
  "REGION X" = "REGION X",
  "DAVAO REGION" = "REGION XI", "DAVAO" = "REGION XI",
  "REGION 11" = "REGION XI", "REGION XI" = "REGION XI",
  "SOCCSKSARGEN" = "REGION XII", "REGION 12" = "REGION XII",
  "REGION XII" = "REGION XII",
  "CARAGA" = "REGION XIII", "CARAGA REGION" = "REGION XIII",
  "REGION 13" = "REGION XIII", "REGION XIII" = "REGION XIII",
  "AUTONOMOUS REGION IN MUSLIM MINDANAO" = "BARMM",
  "AUTONOMOUS REGION OF MUSLIM MINDANAO" = "BARMM",
  "ARMM" = "BARMM", "BANGSAMORO" = "BARMM",
  "BANGSAMORO AUTONOMOUS REGION IN MUSLIM MINDANAO" = "BARMM",
  "BARMM" = "BARMM"
)

# Province-to-region lookup for GADM/Natural Earth dissolves
province_to_region <- c(
  # CAR
  "ABRA" = "CAR", "APAYAO" = "CAR", "BENGUET" = "CAR", "IFUGAO" = "CAR",
  "KALINGA" = "CAR", "MOUNTAIN PROVINCE" = "CAR", "BAGUIO" = "CAR",
  "BAGUIO CITY" = "CAR", "CITY OF BAGUIO" = "CAR",
  # Region I
  "ILOCOS NORTE" = "REGION I", "ILOCOS SUR" = "REGION I",
  "LA UNION" = "REGION I", "PANGASINAN" = "REGION I",
  # Region II
  "BATANES" = "REGION II", "CAGAYAN" = "REGION II", "ISABELA" = "REGION II",
  "NUEVA VIZCAYA" = "REGION II", "QUIRINO" = "REGION II",
  # Region III
  "AURORA" = "REGION III", "BATAAN" = "REGION III", "BULACAN" = "REGION III",
  "NUEVA ECIJA" = "REGION III", "PAMPANGA" = "REGION III",
  "TARLAC" = "REGION III", "ZAMBALES" = "REGION III",
  # Region IV-A
  "BATANGAS" = "REGION IV-A", "CAVITE" = "REGION IV-A",
  "LAGUNA" = "REGION IV-A", "QUEZON" = "REGION IV-A",
  "RIZAL" = "REGION IV-A",
  # MIMAROPA
  "MARINDUQUE" = "MIMAROPA", "OCCIDENTAL MINDORO" = "MIMAROPA",
  "ORIENTAL MINDORO" = "MIMAROPA", "MINDORO OCCIDENTAL" = "MIMAROPA",
  "MINDORO ORIENTAL" = "MIMAROPA", "PALAWAN" = "MIMAROPA",
  "ROMBLON" = "MIMAROPA",
  # Region V
  "ALBAY" = "REGION V", "CAMARINES NORTE" = "REGION V",
  "CAMARINES SUR" = "REGION V", "CATANDUANES" = "REGION V",
  "MASBATE" = "REGION V", "SORSOGON" = "REGION V",
  # Region VI
  "AKLAN" = "REGION VI", "ANTIQUE" = "REGION VI", "CAPIZ" = "REGION VI",
  "GUIMARAS" = "REGION VI", "ILOILO" = "REGION VI",
  "NEGROS OCCIDENTAL" = "REGION VI",
  # Region VII
  "BOHOL" = "REGION VII", "CEBU" = "REGION VII",
  "NEGROS ORIENTAL" = "REGION VII", "SIQUIJOR" = "REGION VII",
  # Region VIII
  "BILIRAN" = "REGION VIII", "EASTERN SAMAR" = "REGION VIII",
  "LEYTE" = "REGION VIII", "NORTHERN SAMAR" = "REGION VIII",
  "SAMAR" = "REGION VIII", "WESTERN SAMAR" = "REGION VIII",
  "SAMAR WESTERN" = "REGION VIII", "SOUTHERN LEYTE" = "REGION VIII",
  # Region IX
  "ZAMBOANGA DEL NORTE" = "REGION IX", "ZAMBOANGA DEL SUR" = "REGION IX",
  "ZAMBOANGA SIBUGAY" = "REGION IX", "CITY OF ISABELA" = "REGION IX",
  "ISABELA CITY" = "REGION IX",
  # Region X
  "BUKIDNON" = "REGION X", "CAMIGUIN" = "REGION X",
  "LANAO DEL NORTE" = "REGION X", "MISAMIS OCCIDENTAL" = "REGION X",
  "MISAMIS ORIENTAL" = "REGION X",
  # Region XI
  "COMPOSTELA VALLEY" = "REGION XI", "DAVAO DE ORO" = "REGION XI",
  "DAVAO DEL NORTE" = "REGION XI", "DAVAO DEL SUR" = "REGION XI",
  "DAVAO OCCIDENTAL" = "REGION XI", "DAVAO ORIENTAL" = "REGION XI",
  # Region XII
  "COTABATO" = "REGION XII", "NORTH COTABATO" = "REGION XII",
  "COTABATO NORTH" = "REGION XII", "SARANGANI" = "REGION XII",
  "SOUTH COTABATO" = "REGION XII", "COTABATO SOUTH" = "REGION XII",
  "SULTAN KUDARAT" = "REGION XII",
  # Region XIII
  "AGUSAN DEL NORTE" = "REGION XIII", "AGUSAN DEL SUR" = "REGION XIII",
  "DINAGAT ISLANDS" = "REGION XIII", "ISLANDS OF DINAGAT" = "REGION XIII",
  "SURIGAO DEL NORTE" = "REGION XIII", "SURIGAO DEL SUR" = "REGION XIII",
  # BARMM
  "BASILAN" = "BARMM", "LANAO DEL SUR" = "BARMM", "MAGUINDANAO" = "BARMM",
  "MAGUINDANAO DEL NORTE" = "BARMM", "MAGUINDANAO DEL SUR" = "BARMM",
  "SULU" = "BARMM", "TAWI-TAWI" = "BARMM", "TAWI TAWI" = "BARMM",
  # NCR cities
  "QUEZON CITY" = "NCR", "CALOOCAN" = "NCR",
  "LAS PINAS" = "NCR", "LAS PINAS CITY" = "NCR", "MAKATI" = "NCR",
  "MAKATI CITY" = "NCR", "MALABON" = "NCR", "MANDALUYONG" = "NCR",
  "MARIKINA" = "NCR", "MUNTINLUPA" = "NCR", "NAVOTAS" = "NCR",
  "PARANAQUE" = "NCR", "PASAY" = "NCR", "PASIG" = "NCR",
  "SAN JUAN" = "NCR", "TAGUIG" = "NCR", "VALENZUELA" = "NCR",
  "PATEROS" = "NCR"
)

canonical_region <- function(x) {
  s   <- .normalize_key(x)
  out <- unname(region_recoder[s])
  is_unmatched <- is.na(out) | !(out %in% CANONICAL_17)
  if (any(is_unmatched)) {
    prov_lookup <- unname(province_to_region[s[is_unmatched]])
    out[is_unmatched] <- prov_lookup
  }
  out[is.na(out)] <- s[is.na(out)]
  out
}

# Apply canonicaliser to dominance and metrics tables before joining map.
dominance_df_canonical <- dominance_df %>%
  dplyr::mutate(REGION = canonical_region(REGION))
regional_metrics_canonical <- regional_metrics %>%
  dplyr::mutate(REGION = canonical_region(REGION))
boot_ci_canonical <- boot_ci %>%
  dplyr::mutate(REGION = canonical_region(REGION))

# ------------------------------------------------------------------------------
# 20. WINNER-ROW TABLE AND HEAD-TO-HEAD CONSENSUS TIER CLASSIFICATION
# ------------------------------------------------------------------------------
# For each region, the all-pairs head-to-head test is computed on the 2x2
# bootstrap win-counts table (chi-squared, Fisher exact fallback when any
# expected cell < 5), one-sided p-value based on which detector has the
# higher empirical win count. Three pairs per region:
#   ConstantTA  vs OutbreakThr
#   ContinuousTA vs OutbreakThr
#   ConstantTA  vs ContinuousTA
#
# Dual Bonferroni reporting:
#   Primary (strict)        : factor 3 within-region x N_regions cross-region
#                             = total_corr_strict
#   Sensitivity (within-only): factor 3 within-region only
#
# Four-tier consensus classification (using primary p_bonf_strict):
#   strong    : both pairs sig in the Observed_Winner's favour.
#   partial   : exactly ONE pair sig in the Observed_Winner's favour
#               (the other ns), no rival is sig against it.
#   lead_only : Observed_Winner has zero significant wins, zero significant
#               losses, AND P_Observed_Winner >= 0.50.
#   contested : a rival is sig against the Observed_Winner, OR
#               P_Observed_Winner < 0.50.
# ------------------------------------------------------------------------------

winner_rows <- dominance_df_canonical %>%
  dplyr::select(REGION, Observed_Winner, Dominance_Probability,
                P_ConstantTA, P_ContinuousTA, P_OutbreakThreshold) %>%
  dplyr::left_join(
    regional_metrics_canonical %>%
      dplyr::select(REGION, Method, TAM, N_True_Alarms,
                    Sensitivity, Mean_Lead_Time, WP) %>%
      dplyr::mutate(Method = as.character(Method)),
    by = c("REGION" = "REGION", "Observed_Winner" = "Method")
  )

n_missing_metric <- sum(is.na(winner_rows$TAM))
if (n_missing_metric > 0) {
  warning("Could not match metric values for ", n_missing_metric, " region(s).")
}

# Pair test (chi-squared with Fisher fallback).
.pair_test <- function(name_X, p_X, name_Y, p_Y, N) {
  if (any(is.na(c(p_X, p_Y, N))) || N <= 0L)
    return(list(p_one_sided = NA_real_, pair_winner = NA_character_))
  n_X <- as.integer(round(p_X * N))
  n_Y <- as.integer(round(p_Y * N))
  if (n_X == n_Y) return(list(p_one_sided = 1.0, pair_winner = NA_character_))
  if (n_X > n_Y) {
    n_higher <- n_X; n_lower <- n_Y; pair_winner <- name_X
  } else {
    n_higher <- n_Y; n_lower <- n_X; pair_winner <- name_Y
  }
  tab <- matrix(c(n_higher, N - n_higher, n_lower, N - n_lower),
                nrow = 2L, byrow = TRUE)
  expected_cells <- as.numeric(rowSums(tab) %o% colSums(tab) / sum(tab))
  use_fisher <- any(expected_cells < 5)
  if (use_fisher) {
    p_two <- tryCatch(stats::fisher.test(tab)$p.value,
                      error = function(e) NA_real_)
  } else {
    p_two <- tryCatch(
      suppressWarnings(stats::chisq.test(tab, correct = FALSE))$p.value,
      error = function(e) NA_real_)
  }
  if (is.na(p_two)) return(list(p_one_sided = NA_real_,
                                pair_winner = NA_character_))
  p_one <- pmin(pmax(p_two / 2, .Machine$double.eps), 1)
  list(p_one_sided = p_one, pair_winner = pair_winner)
}

N_PAIRS  <- 3L
n_regions <- nrow(winner_rows)
total_corr_strict <- N_PAIRS * n_regions
total_corr_within <- N_PAIRS

p_const_vs_ot_raw   <- numeric(n_regions)
p_cont_vs_ot_raw    <- numeric(n_regions)
p_const_vs_cont_raw <- numeric(n_regions)
w_const_vs_ot       <- character(n_regions)
w_cont_vs_ot        <- character(n_regions)
w_const_vs_cont     <- character(n_regions)
for (i in seq_len(n_regions)) {
  r  <- winner_rows[i, , drop = FALSE]
  t1 <- .pair_test("ConstantTA",   r$P_ConstantTA,
                   "OutbreakThr",  r$P_OutbreakThreshold, N_BOOTS)
  t2 <- .pair_test("ContinuousTA", r$P_ContinuousTA,
                   "OutbreakThr",  r$P_OutbreakThreshold, N_BOOTS)
  t3 <- .pair_test("ConstantTA",   r$P_ConstantTA,
                   "ContinuousTA", r$P_ContinuousTA,      N_BOOTS)
  p_const_vs_ot_raw[i]   <- t1$p_one_sided
  p_cont_vs_ot_raw[i]    <- t2$p_one_sided
  p_const_vs_cont_raw[i] <- t3$p_one_sided
  w_const_vs_ot[i]       <- t1$pair_winner
  w_cont_vs_ot[i]        <- t2$pair_winner
  w_const_vs_cont[i]     <- t3$pair_winner
}

# Helper: long detector names → short codes.
.detector_short <- function(d) {
  dplyr::case_when(
    d == "Constant Transmission Acceleration"   ~ "ConstantTA",
    d == "Continuous Transmission Acceleration" ~ "ContinuousTA",
    d == "Outbreak Threshold"                   ~ "OutbreakThr",
    TRUE                                        ~ NA_character_
  )
}

winner_rows <- winner_rows %>%
  dplyr::mutate(
    p_const_vs_ot_bonf   = pmin(p_const_vs_ot_raw   * total_corr_strict, 1),
    p_cont_vs_ot_bonf    = pmin(p_cont_vs_ot_raw    * total_corr_strict, 1),
    p_const_vs_cont_bonf = pmin(p_const_vs_cont_raw * total_corr_strict, 1),
    p_const_vs_ot_bonf_within   = pmin(p_const_vs_ot_raw   * total_corr_within, 1),
    p_cont_vs_ot_bonf_within    = pmin(p_cont_vs_ot_raw    * total_corr_within, 1),
    p_const_vs_cont_bonf_within = pmin(p_const_vs_cont_raw * total_corr_within, 1),
    pair_winner_const_vs_ot   = w_const_vs_ot,
    pair_winner_cont_vs_ot    = w_cont_vs_ot,
    pair_winner_const_vs_cont = w_const_vs_cont,
    sig_const_vs_ot   = !is.na(p_const_vs_ot_bonf)   & p_const_vs_ot_bonf   < SIG_LEVEL,
    sig_cont_vs_ot    = !is.na(p_cont_vs_ot_bonf)    & p_cont_vs_ot_bonf    < SIG_LEVEL,
    sig_const_vs_cont = !is.na(p_const_vs_cont_bonf) & p_const_vs_cont_bonf < SIG_LEVEL,
    sig_const_vs_ot_within   = !is.na(p_const_vs_ot_bonf_within)   &
      p_const_vs_ot_bonf_within   < SIG_LEVEL,
    sig_cont_vs_ot_within    = !is.na(p_cont_vs_ot_bonf_within)    &
      p_cont_vs_ot_bonf_within    < SIG_LEVEL,
    sig_const_vs_cont_within = !is.na(p_const_vs_cont_bonf_within) &
      p_const_vs_cont_bonf_within < SIG_LEVEL,
    const_dominates =
      sig_const_vs_ot   & pair_winner_const_vs_ot   == "ConstantTA" &
      sig_const_vs_cont & pair_winner_const_vs_cont == "ConstantTA",
    cont_dominates =
      sig_cont_vs_ot    & pair_winner_cont_vs_ot    == "ContinuousTA" &
      sig_const_vs_cont & pair_winner_const_vs_cont == "ContinuousTA",
    ot_dominates =
      sig_const_vs_ot & pair_winner_const_vs_ot == "OutbreakThr" &
      sig_cont_vs_ot  & pair_winner_cont_vs_ot  == "OutbreakThr",
    n_dominators = as.integer(const_dominates) +
                   as.integer(cont_dominates) +
                   as.integer(ot_dominates),
    consensus_dominator = dplyr::case_when(
      n_dominators != 1L ~ NA_character_,
      const_dominates    ~ "Constant Transmission Acceleration",
      cont_dominates     ~ "Continuous Transmission Acceleration",
      ot_dominates       ~ "Outbreak Threshold",
      TRUE               ~ NA_character_
    ),
    p_consensus = dplyr::case_when(
      const_dominates ~ pmax(p_const_vs_ot_bonf, p_const_vs_cont_bonf),
      cont_dominates  ~ pmax(p_cont_vs_ot_bonf,  p_const_vs_cont_bonf),
      ot_dominates    ~ pmax(p_const_vs_ot_bonf, p_cont_vs_ot_bonf),
      TRUE            ~ NA_real_
    ),
    sig_star_dominance = dplyr::case_when(
      is.na(p_consensus)  ~ "ns",
      p_consensus < 0.001 ~ "***",
      p_consensus < 0.01  ~ "**",
      p_consensus < 0.05  ~ "*",
      TRUE                ~ "ns"
    ),
    consensus_R1_pass = !is.na(Observed_Winner) & nzchar(Observed_Winner),
    consensus_R2_pass = !is.na(consensus_dominator),
    consensus_pass    = consensus_R1_pass & consensus_R2_pass &
                        (consensus_dominator == Observed_Winner),
    OW_short = .detector_short(Observed_Winner),
    OW_sig_win_pair_1 = dplyr::case_when(
      OW_short == "ConstantTA"   ~ sig_const_vs_ot   &
        pair_winner_const_vs_ot   == "ConstantTA",
      OW_short == "ContinuousTA" ~ sig_cont_vs_ot    &
        pair_winner_cont_vs_ot    == "ContinuousTA",
      OW_short == "OutbreakThr"  ~ sig_const_vs_ot   &
        pair_winner_const_vs_ot   == "OutbreakThr",
      TRUE                       ~ FALSE
    ),
    OW_sig_win_pair_2 = dplyr::case_when(
      OW_short == "ConstantTA"   ~ sig_const_vs_cont &
        pair_winner_const_vs_cont == "ConstantTA",
      OW_short == "ContinuousTA" ~ sig_const_vs_cont &
        pair_winner_const_vs_cont == "ContinuousTA",
      OW_short == "OutbreakThr"  ~ sig_cont_vs_ot    &
        pair_winner_cont_vs_ot    == "OutbreakThr",
      TRUE                       ~ FALSE
    ),
    OW_sig_loss_pair_1 = dplyr::case_when(
      OW_short == "ConstantTA"   ~ sig_const_vs_ot   &
        pair_winner_const_vs_ot   == "OutbreakThr",
      OW_short == "ContinuousTA" ~ sig_cont_vs_ot    &
        pair_winner_cont_vs_ot    == "OutbreakThr",
      OW_short == "OutbreakThr"  ~ sig_const_vs_ot   &
        pair_winner_const_vs_ot   == "ConstantTA",
      TRUE                       ~ FALSE
    ),
    OW_sig_loss_pair_2 = dplyr::case_when(
      OW_short == "ConstantTA"   ~ sig_const_vs_cont &
        pair_winner_const_vs_cont == "ContinuousTA",
      OW_short == "ContinuousTA" ~ sig_const_vs_cont &
        pair_winner_const_vs_cont == "ConstantTA",
      OW_short == "OutbreakThr"  ~ sig_cont_vs_ot    &
        pair_winner_cont_vs_ot    == "ContinuousTA",
      TRUE                       ~ FALSE
    ),
    OW_n_sig_wins   = as.integer(OW_sig_win_pair_1) + as.integer(OW_sig_win_pair_2),
    OW_n_sig_losses = as.integer(OW_sig_loss_pair_1) + as.integer(OW_sig_loss_pair_2),
    consensus_tier = dplyr::case_when(
      !consensus_R1_pass                                   ~ "contested",
      OW_n_sig_wins == 2L & OW_n_sig_losses == 0L          ~ "strong",
      OW_n_sig_wins == 1L & OW_n_sig_losses == 0L          ~ "partial",
      OW_n_sig_wins == 0L & OW_n_sig_losses == 0L &
        !is.na(Dominance_Probability) &
        Dominance_Probability >= 0.50                      ~ "lead_only",
      TRUE                                                 ~ "contested"
    ),
    Leader_Label_short = dplyr::case_when(
      Observed_Winner == "Constant Transmission Acceleration"   ~ "Constant TA",
      Observed_Winner == "Continuous Transmission Acceleration" ~ "Continuous TA",
      Observed_Winner == "Outbreak Threshold"                   ~ "Outbreak Threshold",
      TRUE                                                       ~ NA_character_
    ),
    consensus_winner = dplyr::case_when(
      consensus_tier == "strong"    ~ Observed_Winner,
      consensus_tier == "partial"   ~ Observed_Winner,
      consensus_tier == "lead_only" ~ "No consensus",
      TRUE                          ~ "No consensus"
    )
  ) %>%
  dplyr::select(-OW_short)

n_strong    <- sum(winner_rows$consensus_tier == "strong",    na.rm = TRUE)
n_partial   <- sum(winner_rows$consensus_tier == "partial",   na.rm = TRUE)
n_lead_only <- sum(winner_rows$consensus_tier == "lead_only", na.rm = TRUE)
n_contested <- sum(winner_rows$consensus_tier == "contested", na.rm = TRUE)

n_strong_within <- sum(
  with(winner_rows, {
    OW_within_w1 <- dplyr::case_when(
      Observed_Winner == "Constant Transmission Acceleration"   ~
        sig_const_vs_ot_within   & pair_winner_const_vs_ot   == "ConstantTA",
      Observed_Winner == "Continuous Transmission Acceleration" ~
        sig_cont_vs_ot_within    & pair_winner_cont_vs_ot    == "ContinuousTA",
      Observed_Winner == "Outbreak Threshold"                   ~
        sig_const_vs_ot_within   & pair_winner_const_vs_ot   == "OutbreakThr",
      TRUE                                                       ~ FALSE
    )
    OW_within_w2 <- dplyr::case_when(
      Observed_Winner == "Constant Transmission Acceleration"   ~
        sig_const_vs_cont_within & pair_winner_const_vs_cont == "ConstantTA",
      Observed_Winner == "Continuous Transmission Acceleration" ~
        sig_const_vs_cont_within & pair_winner_const_vs_cont == "ContinuousTA",
      Observed_Winner == "Outbreak Threshold"                   ~
        sig_cont_vs_ot_within    & pair_winner_cont_vs_ot    == "OutbreakThr",
      TRUE                                                       ~ FALSE
    )
    OW_within_w1 & OW_within_w2
  }),
  na.rm = TRUE
)

cat(sprintf(
  "\nConsensus tiers: strong=%d, partial=%d, lead_only=%d, contested=%d (of %d).\n",
  n_strong, n_partial, n_lead_only, n_contested, n_regions))
cat(sprintf(
  "Within-region-only Bonferroni sensitivity: %d regions meet 'strong' (vs %d under strict).\n",
  n_strong_within, n_strong))

# ------------------------------------------------------------------------------
# 21. PER-REGION PER-METRIC SIGNIFICANCE (six ordered pairs)
# ------------------------------------------------------------------------------
# Each ordered pair X vs Y is its own one-sided test of H1: X > Y on the
# explicit signed difference X - Y. The two ordered directions of an
# unordered comparison give DIFFERENT p-values that sum to ~1.
#
# Test method: paired bootstrap on year-cluster replicates (rigorous).
# Bonferroni: factor 5 within (REGION x ordered_pair). Each ordered pair is
# its own family of 5 metrics; the six ordered pairs are six separate families.
# ------------------------------------------------------------------------------

DETECTOR_PAIRS_FOR_MAP <- list(
  list(A = "Constant Transmission Acceleration",
       B = "Outbreak Threshold",
       short = "Constant TA vs OT",
       code  = "const_vs_ot"),
  list(A = "Continuous Transmission Acceleration",
       B = "Outbreak Threshold",
       short = "Continuous TA vs OT",
       code  = "cont_vs_ot"),
  list(A = "Constant Transmission Acceleration",
       B = "Continuous Transmission Acceleration",
       short = "Constant TA vs Continuous TA",
       code  = "const_vs_cont")
)

ORDERED_PAIRS <- list(
  list(label = "Constant TA vs OT",
       X_full = "Constant Transmission Acceleration",
       Y_full = "Outbreak Threshold",
       X_short = "Constant TA", Y_short = "OT",
       code = "const_vs_ot", sym_code = "ot_vs_const"),
  list(label = "Continuous TA vs OT",
       X_full = "Continuous Transmission Acceleration",
       Y_full = "Outbreak Threshold",
       X_short = "Continuous TA", Y_short = "OT",
       code = "cont_vs_ot", sym_code = "ot_vs_cont"),
  list(label = "Constant TA vs Continuous TA",
       X_full = "Constant Transmission Acceleration",
       Y_full = "Continuous Transmission Acceleration",
       X_short = "Constant TA", Y_short = "Continuous TA",
       code = "const_vs_cont", sym_code = "cont_vs_const"),
  list(label = "Outbreak vs Constant TA",
       X_full = "Outbreak Threshold",
       Y_full = "Constant Transmission Acceleration",
       X_short = "OT", Y_short = "Constant TA",
       code = "ot_vs_const", sym_code = "const_vs_ot"),
  list(label = "Outbreak vs Continuous TA",
       X_full = "Outbreak Threshold",
       Y_full = "Continuous Transmission Acceleration",
       X_short = "OT", Y_short = "Continuous TA",
       code = "ot_vs_cont", sym_code = "cont_vs_ot"),
  list(label = "Continuous TA vs Constant TA",
       X_full = "Continuous Transmission Acceleration",
       Y_full = "Constant Transmission Acceleration",
       X_short = "Continuous TA", Y_short = "Constant TA",
       code = "cont_vs_const", sym_code = "const_vs_cont")
)

CI_LO_HI <- list(
  TAM            = c("TAM_lo",            "TAM_hi"),
  N_True_Alarms  = c("N_True_Alarms_lo",  "N_True_Alarms_hi"),
  Sensitivity    = c("Sens_lo",           "Sens_hi"),
  Mean_Lead_Time = c("MLT_lo",            "MLT_hi"),
  WP             = c("WP_lo",             "WP_hi")
)

# Build long-format replicates frame (with canonical region names) for the
# paired-bootstrap test.
all_reps_canonical <- dplyr::bind_rows(
  lapply(names(boot_results), function(reg) {
    reps <- boot_results[[reg]]$replicates
    reps$REGION <- canonical_region(reg)
    reps
  })
)

compute_level1_paired_bootstrap_ordered <- function(reps_df, regions,
                                                    ordered_pairs, metrics) {
  out <- list()
  for (reg in regions) {
    sub <- reps_df %>% dplyr::filter(REGION == reg)
    if (nrow(sub) == 0L) next
    for (op in ordered_pairs) {
      x_df <- sub %>% dplyr::filter(Method == op$X_full) %>% dplyr::arrange(.rep)
      y_df <- sub %>% dplyr::filter(Method == op$Y_full) %>% dplyr::arrange(.rep)
      if (nrow(x_df) == 0L || nrow(y_df) == 0L) next
      n_common <- min(nrow(x_df), nrow(y_df))
      x_df <- x_df[seq_len(n_common), , drop = FALSE]
      y_df <- y_df[seq_len(n_common), , drop = FALSE]
      for (m in metrics) {
        diffs <- x_df[[m]] - y_df[[m]]
        diffs <- diffs[is.finite(diffs)]
        if (length(diffs) < 50L) {
          out[[length(out) + 1L]] <- data.frame(
            REGION = reg, Ordered_Pair = op$label, Pair_Code = op$code,
            Sym_Pair_Code = op$sym_code,
            Detector_X = op$X_full, Detector_Y = op$Y_full,
            Metric = m, n_used = length(diffs),
            prob_X_better = NA_real_, median_diff_XY = NA_real_,
            p_one_sided = NA_real_, p_two_sided = NA_real_,
            method = "bootstrap_paired", stringsAsFactors = FALSE)
          next
        }
        prX     <- mean(diffs > 0)
        med     <- stats::median(diffs)
        p_floor <- 1 / length(diffs)
        p_one   <- max(1 - prX, p_floor)
        p_two   <- max(2 * min(prX, 1 - prX), p_floor)
        out[[length(out) + 1L]] <- data.frame(
          REGION = reg, Ordered_Pair = op$label, Pair_Code = op$code,
          Sym_Pair_Code = op$sym_code,
          Detector_X = op$X_full, Detector_Y = op$Y_full,
          Metric = m, n_used = length(diffs),
          prob_X_better = prX, median_diff_XY = med,
          p_one_sided = p_one, p_two_sided = p_two,
          method = "bootstrap_paired", stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, out)
}

regions_in_order <- unique(winner_rows$REGION)

cat("\nPer-region per-metric inference: paired bootstrap, 6 ordered pairs.\n")
level1_results_ordered <- compute_level1_paired_bootstrap_ordered(
  all_reps_canonical, regions_in_order, ORDERED_PAIRS, HH_METRICS_REGIONAL
)
level1_method <- "bootstrap_paired"

# Bonferroni adjustment within (REGION x Pair_Code).
if (nrow(level1_results_ordered) > 0L) {
  level1_results_ordered <- level1_results_ordered %>%
    dplyr::group_by(REGION, Pair_Code) %>%
    dplyr::mutate(
      p_bonf = pmin(p_one_sided * BONF_LEVEL1_FACTOR, 1),
      sig_star = dplyr::case_when(
        is.na(p_bonf)  ~ "ns",
        p_bonf < 0.001 ~ "***",
        p_bonf < 0.01  ~ "**",
        p_bonf < 0.05  ~ "*",
        TRUE           ~ "ns"
      )
    ) %>%
    dplyr::ungroup()
}

cat("=== Per-region per-metric significance (6 ordered pairs x 5 metrics x ",
    n_regions, " regions = ", nrow(level1_results_ordered), " tests) ===\n",
    sep = "")
cat("Method:                     ", level1_method, "\n", sep = "")
cat("Test:                       one-sided H1: X > Y per ordered pair\n")
cat("Bonferroni factor:          ", BONF_LEVEL1_FACTOR,
    " (within REGION x ordered pair)\n", sep = "")
cat("Significant at p<0.05:      ",
    sum(level1_results_ordered$sig_star %in% c("*", "**", "***"), na.rm = TRUE),
    " of ", sum(!is.na(level1_results_ordered$p_bonf)), "\n", sep = "")

# ------------------------------------------------------------------------------
# 22. CROSS-REGION PER-METRIC WILCOXON (formatted from in-memory results)
# ------------------------------------------------------------------------------
level2_results <- all_wilcoxon_results %>%
  dplyr::transmute(
    Metric, Comparison,
    Detector_A, Detector_B,
    n_pairs     = N_Regions_Paired,
    median_diff = Median_Diff_AminusB,
    p_value     = p_value
  ) %>%
  dplyr::group_by(Metric) %>%
  dplyr::mutate(
    p_bonf = pmin(p_value * BONF_LEVEL2_FACTOR, 1),
    sig_star = dplyr::case_when(
      is.na(p_bonf)  ~ "ns",
      p_bonf < 0.001 ~ "***",
      p_bonf < 0.01  ~ "**",
      p_bonf < 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  ) %>%
  dplyr::ungroup()

# Per-metric significance for the consensus winner (weakest-link rule among
# the two ordered pairs the winner participates in as X).
build_per_metric_sig_for_region <- function(reg, consensus_winner) {
  sub <- level1_results_ordered %>% dplyr::filter(REGION == reg)
  if (nrow(sub) == 0L) {
    return(data.frame(REGION = reg, Metric = HH_METRICS_REGIONAL,
                      sig_star = "na", p_bonf = NA_real_,
                      reference_pair = NA_character_,
                      stringsAsFactors = FALSE))
  }
  rows <- list()
  for (m in HH_METRICS_REGIONAL) {
    sub_m <- sub %>% dplyr::filter(Metric == m)
    if (nrow(sub_m) == 0L) {
      rows[[length(rows) + 1L]] <- data.frame(
        REGION = reg, Metric = m, sig_star = "na",
        p_bonf = NA_real_, reference_pair = NA_character_,
        stringsAsFactors = FALSE)
      next
    }
    if (consensus_winner == "Constant Transmission Acceleration") {
      ref_codes <- c("const_vs_ot", "const_vs_cont")
    } else if (consensus_winner == "Continuous Transmission Acceleration") {
      ref_codes <- c("cont_vs_ot", "cont_vs_const")
    } else if (consensus_winner == "Outbreak Threshold") {
      ref_codes <- c("ot_vs_const", "ot_vs_cont")
    } else {
      ref_codes <- unique(sub_m$Pair_Code)
    }
    sub_m_ref <- sub_m %>% dplyr::filter(Pair_Code %in% ref_codes)
    if (nrow(sub_m_ref) == 0L) {
      rows[[length(rows) + 1L]] <- data.frame(
        REGION = reg, Metric = m, sig_star = "na",
        p_bonf = NA_real_, reference_pair = NA_character_,
        stringsAsFactors = FALSE)
      next
    }
    if (consensus_winner %in% c("Constant Transmission Acceleration",
                                "Continuous Transmission Acceleration",
                                "Outbreak Threshold")) {
      idx_pick <- which.max(sub_m_ref$p_bonf)[1]
    } else {
      idx_pick <- which.min(sub_m_ref$p_bonf)[1]
    }
    pick <- sub_m_ref[idx_pick, , drop = FALSE]
    rows[[length(rows) + 1L]] <- data.frame(
      REGION = reg, Metric = m,
      sig_star = pick$sig_star, p_bonf = pick$p_bonf,
      reference_pair = pick$Ordered_Pair, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

per_region_per_metric_sig <- do.call(rbind, lapply(seq_len(nrow(winner_rows)),
  function(i) {
    build_per_metric_sig_for_region(
      winner_rows$REGION[i], winner_rows$consensus_winner[i]
    )
  }))

sig_wide <- per_region_per_metric_sig %>%
  dplyr::select(REGION, Metric, sig_star) %>%
  tidyr::pivot_wider(names_from = Metric, values_from = sig_star,
                     names_prefix = "sig_") %>%
  as.data.frame()

winner_rows <- winner_rows %>% dplyr::left_join(sig_wide, by = "REGION")

# ------------------------------------------------------------------------------
# 23. PHILIPPINE GEOMETRY
# ------------------------------------------------------------------------------
get_ph_regions <- function() {
  g <- NULL

  # 1. Local boundary file (preferred: fully offline, pins boundary vintage).
  if (!is.null(PH_SHAPEFILE) && file.exists(PH_SHAPEFILE)) {
    message("[map] Geometry source: local file -> ", basename(PH_SHAPEFILE))
    g <- sf::st_read(PH_SHAPEFILE, quiet = TRUE)

  # 2. GADM via geodata. Cached under data/geometry/ rather than tempdir(), so
  #    a second run does not re-download and the archive becomes reproducible.
  } else if (requireNamespace("geodata", quietly = TRUE)) {
    message("[map] Geometry source: GADM via geodata (cached in data/geometry/)")
    for (lvl in c(1L, 2L)) {
      gg <- try(geodata::gadm(country = "PHL", level = lvl,
                              path = DIR_GEOMETRY), silent = TRUE)
      if (!inherits(gg, "try-error") && !is.null(gg)) {
        g <- sf::st_as_sf(gg)
        break
      }
    }

  # 3. Natural Earth states as a last resort.
  } else if (requireNamespace("rnaturalearth", quietly = TRUE)) {
    message("[map] Geometry source: Natural Earth states")
    gg <- try(rnaturalearth::ne_states(country = "Philippines",
                                       returnclass = "sf"), silent = TRUE)
    if (!inherits(gg, "try-error") && !is.null(gg)) g <- gg
  }

  if (is.null(g)) {
    stop("No Philippines geometry available.\n",
         "  Either place a .shp/.gpkg/.geojson in data/geometry/,\n",
         "  or install 'geodata' (needs internet on first run).",
         call. = FALSE)
  }

  # --- CRS validation ---------------------------------------------------------
  # Sources vary: GADM ships EPSG:4326, some national shapefiles ship PRS92 or
  # a local grid, and a few carry no CRS at all. Normalising to WGS84 here means
  # the later projection to EPSG:3123 is well defined, so the scale bar measures
  # true ground distance rather than degrees.
  if (is.na(sf::st_crs(g))) {
    warning("[map] Geometry has no CRS; assuming EPSG:4326 (WGS84).",
            call. = FALSE)
    sf::st_crs(g) <- CRS_WGS84
  }
  g <- sf::st_transform(g, CRS_WGS84)

  # Repair invalid rings (self-intersections in coastline data are common and
  # otherwise cause st_union/summarise to fail).
  if (!all(sf::st_is_valid(g))) {
    message("[map] Repairing invalid geometries with st_make_valid().")
    g <- sf::st_make_valid(g)
  }
  g
}
ph_geom <- get_ph_regions()

name_candidates <- c("REGION", "REGION_NAM", "Region", "name", "NAME_1",
                     "name_1", "NAME_2", "name_2", "REGNAME",
                     "ADM1_EN", "ADM1_PCODE", "ADM1_NAME",
                     "region", "Region_1")
nm_col <- intersect(name_candidates, names(ph_geom))[1]
if (is.na(nm_col))
  stop("Cannot identify region/province column in geometry.", call. = FALSE)
ph_geom$REGION_RAW <- ph_geom[[nm_col]]
ph_geom$REGION     <- canonical_region(ph_geom$REGION_RAW)
if (length(unique(ph_geom$REGION)) > 17 || anyDuplicated(ph_geom$REGION) > 0) {
  ph_geom <- ph_geom %>%
    dplyr::group_by(REGION) %>%
    dplyr::summarise(.groups = "drop")
}

# Join audit
audit_df <- dplyr::full_join(
  ph_geom %>% sf::st_drop_geometry() %>% dplyr::distinct(REGION) %>%
    dplyr::mutate(in_geometry = TRUE),
  winner_rows %>% dplyr::distinct(REGION) %>% dplyr::mutate(in_csv = TRUE),
  by = "REGION"
) %>%
  dplyr::mutate(
    in_geometry = !is.na(in_geometry),
    in_csv      = !is.na(in_csv),
    status = dplyr::case_when(
      in_geometry &  in_csv  ~ "matched",
      in_geometry & !in_csv  ~ "geometry only",
      !in_geometry &  in_csv ~ "csv only",
      TRUE                   ~ "neither"
    )
  ) %>%
  dplyr::arrange(status, REGION)

map_df <- ph_geom %>%
  dplyr::left_join(winner_rows, by = "REGION") %>%
  dplyr::filter(!is.na(Dominance_Probability),
                REGION %in% CANONICAL_17) %>%
  dplyr::mutate(
    Detector_Label = dplyr::case_when(
      consensus_winner == "Constant Transmission Acceleration"   ~ "Constant TA",
      consensus_winner == "Continuous Transmission Acceleration" ~ "Continuous TA",
      consensus_winner == "Outbreak Threshold"                   ~ "Outbreak Threshold",
      consensus_winner == "No consensus"                         ~ "No consensus",
      TRUE                                                       ~ consensus_winner
    ),
    Confidence_Tier = dplyr::case_when(
      Dominance_Probability >= CONFIDENCE_HI       ~ "Solid (Pr >= 0.90)",
      Dominance_Probability >= DOMINANCE_THRESHOLD ~ "Mid (Pr >= 0.75)",
      TRUE                                          ~ "Light (Pr < 0.75)"
    ),
    Tier_Fill_Group = dplyr::case_when(
      consensus_tier == "strong"  ~ as.character(Detector_Label),
      consensus_tier == "partial" ~ as.character(Detector_Label),
      TRUE                        ~ "Contested"
    ),
    Tier_Linetype = dplyr::case_when(
      consensus_tier == "partial" ~ "dashed",
      TRUE                        ~ "solid"
    ),
    Tier_Linewidth = dplyr::case_when(
      consensus_tier == "strong"    ~ 0.25,
      consensus_tier == "partial"   ~ 0.55,
      consensus_tier == "lead_only" ~ 0.30,
      consensus_tier == "contested" ~ 0.30,
      TRUE                          ~ 0.30
    )
  )
map_df$Detector_Label <- factor(
  map_df$Detector_Label,
  levels = c("Constant TA", "Continuous TA", "Outbreak Threshold", "No consensus")
)
map_df$Confidence_Tier <- factor(
  map_df$Confidence_Tier,
  levels = c("Solid (Pr >= 0.90)", "Mid (Pr >= 0.75)", "Light (Pr < 0.75)")
)
map_df$Tier_Fill_Group <- factor(
  map_df$Tier_Fill_Group,
  levels = c("Constant TA", "Continuous TA", "Outbreak Threshold", "Contested")
)
stopifnot(nrow(map_df) >= 1)

# ------------------------------------------------------------------------------
# 24. FIGURE 5 PALETTES
# ------------------------------------------------------------------------------
# Okabe-Ito colourblind-safe palette. The previous scheme paired #2CA02C
# (green) with #D62728 (red), which deuteranopic and protanopic readers cannot
# separate; these three hues remain distinct under all common forms of colour
# vision deficiency and stay separable in greyscale.
# Six target detectors, all Okabe-Ito so the set stays colourblind-safe.
detector_palette <- c(
  "Constant TA"        = "#0072B2",  # blue
  "Continuous TA"      = "#E69F00",  # orange
  "Outbreak Threshold" = "#009E73",  # bluish green
  "Farrington"         = "#CC79A7",  # reddish purple
  "EARS"               = "#56B4E9",  # sky blue
  "EWARS"              = "#D55E00",  # vermillion
  "No consensus"       = "#9A9A9A"   # neutral grey
)
alpha_palette <- c(
  "Solid (Pr >= 0.90)" = 1.00,
  "Mid (Pr >= 0.75)"   = 0.82,
  "Light (Pr < 0.75)"  = 0.60
)
# Red used for the Continuous TA circles in Figure 5C only (see section 28).
PANEL_C_CONTINUOUS_TA_FILL <- "#D55E00"   # Okabe-Ito vermillion

tier_fill_palette <- c(
  "Constant TA"        = unname(detector_palette["Constant TA"]),
  "Continuous TA"      = unname(detector_palette["Continuous TA"]),
  "Outbreak Threshold" = unname(detector_palette["Outbreak Threshold"]),
  "Farrington"         = unname(detector_palette["Farrington"]),
  "EARS"               = unname(detector_palette["EARS"]),
  "EWARS"              = unname(detector_palette["EWARS"]),
  "Contested"          = "#9A9A9A"
)

# ------------------------------------------------------------------------------
# 25. FIGURE 5 PANEL A — CHOROPLETH DETECTOR MAP
# ------------------------------------------------------------------------------
# Geometry is projected to EPSG:3123 (PRS92 / Philippines Zone III) BEFORE
# centroids are computed. Two reasons this matters:
#   * st_point_on_surface() on unprojected lon/lat degrees places labels
#     slightly off for elongated islands; on a projected grid it is correct.
#   * a scale bar over geographic degrees is meaningless, because a degree of
#     longitude shrinks with latitude. On a projected CRS the bar represents a
#     true, constant ground distance.
map_df_proj <- sf::st_transform(map_df, CRS_PH_PROJECTED)

suppressWarnings({
  centroids_sf <- sf::st_point_on_surface(map_df_proj)

  # Plain region names only -- no detector line, no boxes. Dominant detector is
  # already encoded by the polygon fill and its legend, so the second label line
  # was redundant; removing it lets 17 single-line names sit legibly on the map.
  centroids_sf$label <- centroids_sf$REGION
})

p_map <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = map_df_proj,
    # alpha is not mapped here: fills are uniform, without
    # bootstrap-confidence shading, so there is nothing extra to explain in
    # the legend.
    ggplot2::aes(fill = Tier_Fill_Group),
    colour = "grey25",
    linetype = map_df_proj$Tier_Linetype,
    linewidth = pmax(map_df_proj$Tier_Linewidth, 0.30)  # keep boundaries above
                                                        # the press hairline floor
  ) +
  # Centroid markers, drawn from the sf geometry so they stay registered with
  # the polygons under any coord_sf transformation.
  ggplot2::geom_sf(
    data = centroids_sf, colour = "grey15", size = 0.6,
    show.legend = FALSE
  ) +
  # geom_TEXT_repel, not geom_label_repel: plain names drawn directly on the map
  # with no rectangle, fill or callout background. A soft white halo (bg.colour)
  # preserves contrast over dark choropleth fills without drawing a box.
  ggrepel::geom_text_repel(
    data = centroids_sf,
    ggplot2::aes(geometry = geometry, label = label),
    stat = "sf_coordinates",
    size = pub_text_size(PUB_ANNOT),   # points -> mm, matching the shared scale
    family = PUB_FAMILY,
    fontface = "bold", colour = "grey10",
    bg.colour = "white", bg.r = 0.14,
    box.padding = grid::unit(0.32, "lines"),
    point.padding = grid::unit(0.14, "lines"),
    force = 5, force_pull = 0.7, max.overlaps = Inf,
    min.segment.length = 0, segment.color = "grey45",
    segment.size = 0.28, seed = GLOBAL_SEED,
    show.legend = FALSE,
    # Confine repelled region names to the map's own bounding box so they
    # cannot be pushed off the canvas and cut.
    xlim = sf::st_bbox(map_df_proj)[c("xmin", "xmax")],
    ylim = sf::st_bbox(map_df_proj)[c("ymin", "ymax")]
  ) +
  # ONE legend only: the detector that dominates each region. The alpha
  # aesthetic still shades regions by bootstrap confidence, but its legend is
  # suppressed (guide = "none") -- the map then carries a single,
  # unambiguous key.
  ggplot2::scale_fill_manual(
    name   = "Detectors",
    values = tier_fill_palette,
    # Keys are the detectors that actually win at least one region, in the
    # canonical order. A detector with no region is dropped rather than shown
    # as an empty key, which is why "Outbreak Threshold" does not appear
    # unless it wins at least one region.
    limits = intersect(c("Constant TA", "Continuous TA", "Outbreak Threshold"),
                       unique(as.character(map_df_proj$Tier_Fill_Group))),
    drop   = TRUE, na.translate = FALSE
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(order = 1, nrow = 1, byrow = TRUE,
                                 override.aes = list(alpha = 0.9,
                                                     linewidth = 0.3))
  ) +
  ggplot2::labs(title = "a    Regional detector map",
               x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)") +
  theme_pub_map() +
  ggplot2::theme(
    # geom_sf()/coord_sf() default to literal "x" and "y" axis titles unless
    # overridden above; axis.title is re-enabled here (matching the bold,
    # PUB_AXIS_TIT house style used for the other Figure 5 panel) in case
    # theme_pub_map() blanks axis titles for other map panels.
    axis.title.x    = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                             colour = "grey25",
                                             margin = ggplot2::margin(t = 6)),
    axis.title.y    = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                             colour = "grey25",
                                             margin = ggplot2::margin(r = 6)),
    legend.position = "bottom",
    legend.box      = "horizontal",
    legend.box.just = "center",
    # EVEN SPACING between "Detectors", the first key and the second key.
    # The title-to-first-key gap is set by legend.title's right margin; the
    # key-to-key gap by legend.key.spacing.x plus the label's right margin.
    # They are matched at 22 pt so the three elements sit evenly apart.
    legend.title         = ggplot2::element_text(
                             size = PUB_LEG_TIT, face = "bold",
                             colour = "black", vjust = 0.5,
                             margin = ggplot2::margin(r = 22)),
    legend.text          = ggplot2::element_text(
                             size = PUB_LEG_TXT, colour = "black",
                             margin = ggplot2::margin(l = 6, r = 22)),
    legend.key.spacing.x = grid::unit(16, "pt"),
    legend.spacing.x     = grid::unit(16, "pt")
  )

# Cartographic furniture requested by the reviewer: scale bar (bottom-left),
# north arrow (top-right) and an explicit projected CRS. Added last so the
# annotations sit above the polygon layers.
p_map <- pub_map_frame(
  p_map,
  crs          = CRS_PH_PROJECTED,
  graticule    = TRUE,
  scalebar_loc = "bl",
  arrow_loc    = "tr"
)

# ------------------------------------------------------------------------------
# 26. FIGURE 5 PANEL B — METRIC TABLE WITH EMBEDDED PER-METRIC SIGNIFICANCE
# ------------------------------------------------------------------------------
fmt_int <- function(x) ifelse(is.na(x), "-",
                              formatC(round(x), format = "d", big.mark = ","))
fmt_1   <- function(x) ifelse(is.na(x), "-", formatC(x, format = "f", digits = 1))
fmt_2   <- function(x) ifelse(is.na(x), "-", formatC(x, format = "f", digits = 2))
fmt_pct <- function(x) ifelse(is.na(x), "-",
                              paste0(formatC(100 * x, format = "f", digits = 0), "%"))

detector_two_line <- c(
  "Constant Transmission Acceleration"   = "Constant\nTA",
  "Continuous Transmission Acceleration" = "Continuous\nTA",
  "Outbreak Threshold"                   = "Outbreak\nThreshold",
  "No consensus"                         = "No\nconsensus"
)

COLS_TABLE <- c("Region", "Detector",
                "Dominance\nprobability",
                "TAM",
                "N true\nalarms",
                "Sensitivity",
                "Mean lead\ntime (wk)",
                "WP\n(wk)")
N_COLS <- length(COLS_TABLE)

table_long <- winner_rows %>%
  dplyr::filter(REGION %in% CANONICAL_17, !is.na(Observed_Winner)) %>%
  dplyr::mutate(
    Detector_Label_short = dplyr::case_when(
      consensus_tier == "contested" ~ "Contested",
      Leader_Label_short %in% c("Constant TA", "Continuous TA",
                                "Outbreak Threshold") ~ Leader_Label_short,
      TRUE ~ "Contested"
    ),
    Detector_Label_two = dplyr::case_when(
      consensus_tier == "strong"    ~ unname(detector_two_line[Observed_Winner]),
      consensus_tier == "partial"   ~ paste0(Leader_Label_short, "\n(partial)"),
      consensus_tier == "lead_only" ~ paste0(Leader_Label_short, "\n(lead)"),
      consensus_tier == "contested" ~ "Contested",
      TRUE                          ~ "Contested"
    )
  ) %>%
  dplyr::arrange(REGION)

for (m in HH_METRICS_REGIONAL) {
  col_name <- paste0("sig_", m)
  if (!col_name %in% names(table_long)) table_long[[col_name]] <- "na"
  table_long[[col_name]][is.na(table_long[[col_name]])] <- "na"
}

N_ROWS_TABLE <- nrow(table_long)
HEADER_Y_TBL <- N_ROWS_TABLE + 1L
table_long <- table_long %>%
  dplyr::mutate(row_y = rev(seq_len(dplyr::n())))

fmt_sig_paren <- function(s) {
  s <- as.character(s)
  s[is.na(s) | s == "na" | s == ""] <- "-"
  paste0("(", s, ")")
}

cell_value_df <- table_long %>%
  dplyr::transmute(
    REGION, row_y, Detector_Label_short, Detector_Label_two,
    consensus_pass, consensus_tier,
    `1` = REGION,
    `2` = Detector_Label_two,
    `3` = fmt_2(Dominance_Probability),
    `4` = fmt_int(TAM),
    `5` = fmt_1(N_True_Alarms),
    `6` = fmt_pct(Sensitivity),
    `7` = fmt_1(Mean_Lead_Time),
    `8` = fmt_1(WP)
  ) %>%
  tidyr::pivot_longer(cols = `1`:`8`, names_to = "col_x",
                      values_to = "value_text") %>%
  dplyr::mutate(col_x = as.integer(col_x))

cell_sig_df <- table_long %>%
  dplyr::transmute(
    REGION, row_y, consensus_pass, consensus_tier,
    `3` = sig_star_dominance,
    `4` = sig_TAM,
    `5` = sig_N_True_Alarms,
    `6` = sig_Sensitivity,
    `7` = sig_Mean_Lead_Time,
    `8` = sig_WP
  ) %>%
  tidyr::pivot_longer(cols = `3`:`8`, names_to = "col_x",
                      values_to = "sig_star") %>%
  dplyr::mutate(col_x = as.integer(col_x),
                sig_label = fmt_sig_paren(sig_star))

.embedded_sig_cols <- c(3L, 4L, 5L, 6L, 7L, 8L)
.metric_data_cols  <- c(4L, 5L, 6L, 7L, 8L)

.blend_to_white <- function(hex, frac = 0.40) {
  m <- grDevices::col2rgb(hex) / 255
  w <- c(1, 1, 1)
  out <- (1 - frac) * m + frac * w
  grDevices::rgb(out[1, ], out[2, ], out[3, ])
}

.detector_cell_fill <- function(detector_short, tier) {
  base_col <- unname(detector_palette[as.character(detector_short)])
  base_col[is.na(base_col)] <- "#6B7280"
  out <- character(length(tier))
  for (i in seq_along(tier)) {
    out[i] <- switch(
      as.character(tier[i]),
      strong    = base_col[i],
      partial   = .blend_to_white(base_col[i], 0.40),
      lead_only = "#9CA3AF",
      contested = "#6B7280",
      "#6B7280"
    )
  }
  out
}

bg_df <- cell_value_df %>%
  dplyr::mutate(
    is_strong_tier = consensus_tier == "strong",
    fill = dplyr::case_when(
      col_x == 2 ~ .detector_cell_fill(Detector_Label_short, consensus_tier),
      row_y %% 2 == 0 ~ "#F2F2F2",
      TRUE ~ "#FFFFFF"
    ),
    fill = ifelse(is.na(fill), "grey90", fill),
    text_colour_value = dplyr::case_when(
      col_x == 2 ~ "white",
      !is_strong_tier & col_x %in% .metric_data_cols ~ "grey55",
      TRUE ~ "grey15"
    ),
    text_colour_value = ifelse(is.na(text_colour_value), "grey15",
                               text_colour_value),
    fontface_value = dplyr::case_when(
      col_x %in% c(1L, 2L) ~ "bold",
      !is_strong_tier & col_x %in% .metric_data_cols ~ "italic",
      TRUE ~ "plain"
    )
  )

sig_with_bg <- dplyr::left_join(
  cell_sig_df,
  bg_df %>% dplyr::select(row_y, col_x, fill, Detector_Label_short),
  by = c("row_y", "col_x")
) %>%
  dplyr::mutate(
    text_colour_sig = dplyr::case_when(
      sig_star == "***" ~ "#0B2447",
      sig_star == "**"  ~ "#27496D",
      sig_star == "*"   ~ "#3673B6",
      sig_star == "ns"  ~ "grey55",
      TRUE              ~ "grey75"
    ),
    fontface_sig = ifelse(sig_star %in% c("*", "**", "***"), "bold", "plain")
  )

bg_header <- data.frame(
  col_x = seq_len(N_COLS), row_y = HEADER_Y_TBL,
  fill = "#0B2447", text_colour = "white",
  fontface = "bold", label = COLS_TABLE, stringsAsFactors = FALSE
)

y_min_tbl <- 0.5
y_max_tbl <- HEADER_Y_TBL + 0.5
VAL_Y_OFFSET <- 0.16
SIG_Y_OFFSET <- -0.22

# NOTE: p_table (the per-region metrics table panel) is not part of
# Figure 5. The object is retained here because the underlying per-region
# metrics are still written to CSV, and keeping the builder makes it easy to
# reinstate the panel if required. It is not exported to any figure.
p_table <- ggplot2::ggplot() +
  ggplot2::geom_tile(
    data = bg_header,
    ggplot2::aes(x = col_x, y = row_y, fill = I(fill)),
    colour = "white", width = 0.98, height = 0.98
  ) +
  ggplot2::geom_text(
    data = bg_header,
    ggplot2::aes(x = col_x, y = row_y, label = label),
    colour = "white", fontface = "bold",
    size = 2.7, lineheight = 0.85
  ) +
  ggplot2::geom_tile(
    data = bg_df,
    ggplot2::aes(x = col_x, y = row_y, fill = I(fill)),
    colour = "white", width = 0.98, height = 0.98
  ) +
  ggplot2::geom_text(
    data = bg_df %>%
      dplyr::mutate(y_text = ifelse(col_x %in% .embedded_sig_cols,
                                    row_y + VAL_Y_OFFSET, row_y)),
    ggplot2::aes(x = col_x, y = y_text, label = value_text,
                 colour = I(text_colour_value),
                 fontface = I(fontface_value)),
    size = 2.6, lineheight = 0.85
  ) +
  ggplot2::geom_text(
    data = sig_with_bg %>%
      dplyr::mutate(y_text = row_y + SIG_Y_OFFSET),
    ggplot2::aes(x = col_x, y = y_text, label = sig_label,
                 colour = I(text_colour_sig),
                 fontface = I(fontface_sig)),
    size = 2.0, lineheight = 0.85
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq_len(N_COLS),
    limits = c(0.5, N_COLS + 0.5),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(y_min_tbl, y_max_tbl),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  # Header and edge-cell text can extend slightly beyond the cell grid.
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = paste0("b    Bootstrap-winning detector and its per-metric ",
                   "significance, by region (paired bootstrap, Bonf k = ",
                   BONF_LEVEL1_FACTOR, ")")
  ) +
  ggplot2::theme_void(base_size = PUB_BASE, base_family = PUB_FAMILY) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = PUB_TITLE + 1, hjust = 0,
                                       margin = ggplot2::margin(b = 6)),
    plot.margin = ggplot2::margin(10, 12, 8, 8)
  )

# ------------------------------------------------------------------------------
# 27. FIGURE 5 PANEL C — PER-DETECTOR DOT PLOT
# ------------------------------------------------------------------------------
# All six target detectors, pivoted from TARGET_P_COLS rather than a hard-coded
# set of three, so the dominance-probability figure covers every benchmarked
# method and cannot fall out of step with TARGET_DETECTORS.
# Maps a full detector name (as stored in Observed_Winner) to its short label.
WINNER_TO_SHORT <- c(
  "Constant Transmission Acceleration"   = "Constant TA",
  "Continuous Transmission Acceleration" = "Continuous TA",
  "Outbreak Threshold"                   = "Outbreak Threshold",
  "Farrington"                           = "Farrington",
  "EARS"                                 = "EARS",
  "EWARS"                                = "EWARS"
)

DETECTOR_SHORT_LABEL <- c(
  P_ConstantTA        = "Constant TA",
  P_ContinuousTA      = "Continuous TA",
  P_OutbreakThreshold = "Outbreak Threshold",
  P_Farrington        = "Farrington",
  P_EARS              = "EARS",
  P_EWARS             = "EWARS"
)

panel_c_long <- dominance_df_canonical %>%
  dplyr::select(REGION, dplyr::all_of(unname(TARGET_P_COLS))) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(unname(TARGET_P_COLS)),
                      names_to = "Detector_Code", values_to = "Dominance_P") %>%
  dplyr::mutate(
    Detector_Label = unname(DETECTOR_SHORT_LABEL[Detector_Code]),
    Detector_Label = factor(Detector_Label,
                            levels = unname(DETECTOR_SHORT_LABEL[
                              unname(TARGET_P_COLS)]))
  ) %>%
  dplyr::left_join(
    winner_rows %>%
      dplyr::select(REGION, consensus_tier, Observed_Winner, Leader_Label_short),
    by = "REGION"
  ) %>%
  dplyr::mutate(
    # Generalised: a dot is the leader when its detector IS the observed
    # winner, for any of the six target detectors. The previous case_when
    # listed only three and would have marked no leader for Farrington,
    # EARS or EWARS regions.
    is_leader_dot = !is.na(Observed_Winner) &
      Detector_Label == unname(WINNER_TO_SHORT[Observed_Winner]),
    border_colour = dplyr::case_when(
      is_leader_dot & consensus_tier == "strong"    ~ "#1F1F1F",
      is_leader_dot & consensus_tier == "partial"   ~ "#1F1F1F",
      is_leader_dot & consensus_tier == "lead_only" ~ "grey45",
      TRUE                                          ~ "grey55"
    ),
    border_stroke = dplyr::case_when(
      is_leader_dot & consensus_tier == "strong"    ~ 0.95,
      is_leader_dot & consensus_tier == "partial"   ~ 0.65,
      is_leader_dot & consensus_tier == "lead_only" ~ 0.55,
      TRUE                                          ~ 0.30
    )
  ) %>%
  dplyr::arrange(REGION, Detector_Label)

pj_c <- ggplot2::position_jitter(width = 0.18, height = 0, seed = 12345)
chance_p <- 1 / 3

# Significance brackets for the dominance-probability panel: Constant TA versus
# each comparator, paired by region and tested with the same paired
# Wilcoxon signed-rank used everywhere else. ONLY SIGNIFICANT comparisons get a
# bracket and a p-value; non-significant pairs are simply not drawn. Every
# comparison remains in the exported Wilcoxon CSV.
# Full pairwise test table FIRST: every comparison, significant or not, with
# V, p, CI and a Drawn_On_Figure flag. This is the file that supports the
# significance brackets on the figure -- exported below so every p-value shown
# on the panel can be traced to a row.
dom_tests <- dominance_pairwise_tests(
  d         = as.data.frame(panel_c_long),
  unit_col  = "REGION",
  level_col = "Detector_Label",
  value_col = "Dominance_P",
  focal     = "Constant TA")

if (!is.null(dom_tests)) {
  utils::write.csv(dom_tests, file.path(OUTPUT_DIR, "Figure5_DominanceProbability_Wilcoxon_Results.csv"),
                   row.names = FALSE)
  cat("  Saved: ", file.path(OUTPUT_DIR, "Figure5_DominanceProbability_Wilcoxon_Results.csv"), "\n", sep = "")
}

# Geometry is built FROM that same table, so figure and CSV cannot disagree.
dom_brackets <- build_sig_brackets(
  d         = as.data.frame(panel_c_long),
  unit_col  = "REGION",
  level_col = "Detector_Label",
  value_col = "Dominance_P",
  focal     = "Constant TA",
  tests     = dom_tests)
dom_vrange <- diff(range(panel_c_long$Dominance_P, na.rm = TRUE))
if (!is.finite(dom_vrange) || dom_vrange <= 0) dom_vrange <- 1
# build_sig_brackets() returns NULL when nothing is significant; nrow(NULL) is
# NULL, not 0, so the count is normalised here rather than inline.
n_dom_brackets <- if (is.null(dom_brackets)) 0L else nrow(dom_brackets)

p_dot <- ggplot2::ggplot(panel_c_long,
                         ggplot2::aes(x = Detector_Label, y = Dominance_P)) +
  # Decisive / chance reference lines, their labels and the shaded decisive
  # band removed. DOMINANCE_THRESHOLD still drives the tier
  # classification and the "Dominance (score >= x)" column; it is simply no
  # longer drawn on this panel.
  ggplot2::geom_point(
    ggplot2::aes(fill = Detector_Label,
                 stroke = border_stroke,
                 colour = border_colour),
    shape = 21, position = pj_c,
    size = 3.0, alpha = 0.92, na.rm = TRUE, show.legend = FALSE
  ) +
  ggplot2::stat_summary(
    fun = stats::median, geom = "crossbar",
    width = 0.55, linewidth = 0.45, colour = "grey15",
    fatten = 0, fill = NA, na.rm = TRUE
  ) +
  # Dot-plot-only override: Continuous TA circles are drawn red, as requested.
  # PANEL_C_CONTINUOUS_TA_FILL is Okabe-Ito vermillion -- a true red that stays
  # separable from the blue and green under deuteranopia and protanopia, so the
  # panel keeps the colourblind-safe property of the rest of the figure.
  #
  # NOTE: this makes Continuous TA red in panel b (the dot plot) but orange in
  # panel a of the same figure. To use red everywhere instead, set
  # detector_palette["Continuous TA"] <- PANEL_C_CONTINUOUS_TA_FILL where the
  # palette is defined (section 24) and delete this override.
  # All six target detectors, built from detector_palette so the scale cannot
  # fall out of step with what is plotted, then the panel-specific red override
  # for Continuous TA is applied on top.
  ggplot2::scale_fill_manual(
    values = {
      .v <- detector_palette[unname(DETECTOR_SHORT_LABEL[unname(TARGET_P_COLS)])]
      .v["Continuous TA"] <- PANEL_C_CONTINUOUS_TA_FILL
      .v
    }
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.05 + 0.16 * n_dom_brackets),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.0"),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_x_discrete(
    # All six detectors get an explicit two-line label; without this,
    # Farrington/EARS/EWARS fall back to their raw level names.
    labels = c("Constant TA"        = "Constant\nTA",
               "Continuous TA"      = "Continuous\nTA",
               "Outbreak Threshold" = "Outbreak\nThreshold",
               "Farrington"         = "Farrington",
               "EARS"               = "EARS",
               "EWARS"              = "EWARS"),
    expand = ggplot2::expansion(add = c(0.6, 0.6))
  ) +
  ggplot2::labs(
    title = "Per-detector regional dominance probabilities",
    x = NULL, y = "Dominance probability"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal(base_size = PUB_BASE, base_family = PUB_FAMILY) +
  theme_bold_axes() +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor   = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(size = PUB_AXIS_TXT, face = "bold",
                                         lineheight = 0.85,
                                         hjust = 0.5, vjust = 1,
                                         margin = ggplot2::margin(t = 4)),
    axis.text.y  = ggplot2::element_text(size = PUB_AXIS_TXT),
    axis.title.y = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                         colour = "grey25",
                                         margin = ggplot2::margin(r = 6)),
    plot.title   = ggplot2::element_text(face = "bold", size = 13,
                                         margin = ggplot2::margin(b = 8)),
    plot.margin  = ggplot2::margin(20, 30, 14, 18)
  )

# ------------------------------------------------------------------------------
# 28. PER-REGION PER-METRIC HEATMAP DATA (saved as CSV; not a figure panel)
# ------------------------------------------------------------------------------
PANEL_D_SUBGRIDS <- list(
  list(label = "Constant TA vs Outbreak Threshold",
       pair_code = "const_vs_ot"),
  list(label = "Continuous TA vs Outbreak Threshold",
       pair_code = "cont_vs_ot"),
  list(label = "Constant TA vs Continuous TA",
       pair_code = "const_vs_cont"),
  list(label = "Outbreak vs Constant TA",
       pair_code = "ot_vs_const"),
  list(label = "Outbreak vs Continuous TA",
       pair_code = "ot_vs_cont"),
  list(label = "Continuous TA vs Constant TA",
       pair_code = "cont_vs_const")
)

build_panel_d_data <- function() {
  rows <- list()
  for (sg in PANEL_D_SUBGRIDS) {
    sub <- level1_results_ordered %>%
      dplyr::filter(Pair_Code == sg$pair_code,
                    REGION %in% CANONICAL_17)
    if (nrow(sub) == 0L) next
    sub <- sub %>%
      dplyr::mutate(
        direction = dplyr::case_when(
          is.na(p_bonf)                      ~ "na",
          sig_star %in% c("*", "**", "***")  ~ "win",
          TRUE                               ~ "ns"
        ),
        fill_key = dplyr::case_when(
          direction == "na"  ~ "na",
          direction == "ns"  ~ "ns",
          direction == "win" ~ paste0("win_", sig_star),
          TRUE               ~ "na"
        ),
        sig_label = ifelse(direction == "win", sig_star, ""),
        Subgrid_Label = sg$label,
        Subgrid_Pair_Code = sg$pair_code
      ) %>%
      dplyr::select(REGION, Metric,
                    sig_star, p_bonf, p_one_sided, median_diff_XY,
                    direction, fill_key, sig_label,
                    Subgrid_Label, Subgrid_Pair_Code)
    rows[[length(rows) + 1L]] <- sub
  }
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}
panel_d_data <- build_panel_d_data()

# ------------------------------------------------------------------------------
# 29. COMPOSE AND SAVE FIGURE 5
# ------------------------------------------------------------------------------
# Panel B (the per-region metrics table) has been removed. Figure 5 is now two
# panels: the map (a) above the dominance-probability plot (b). Former panel c
# is therefore relabelled b throughout -- titles, filenames and captions.
fig5 <- (p_map / p_dot) +
  patchwork::plot_layout(heights = c(1.35, 1.00)) &
  # Padding between the map and dot panels so their titles, legends and axis
  # labels cannot collide in the composite.
  ggplot2::theme(plot.margin = ggplot2::margin(10, 12, 8, 8))

fig5 <- fig5 +
  patchwork::plot_annotation(
    title = "Bootstrap-supported regional dominance of three target detectors",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = PUB_TITLE + 1,
                                         family = base_family_global,
                                         margin = ggplot2::margin(b = 6)),
      plot.tag = ggplot2::element_text(size = PUB_TAG, face = "bold",
                                       family = base_family_global,
                                       hjust = 0, vjust = 1),
      # Same left edge as the title: tags anchor to the plot, and
      # plot.title.position = "plot" moves the title off the panel edge to match.
      plot.tag.position     = "topleft",
      plot.title.position   = "plot"
    )
  )

# Figure 5 (map + dominance probability), authored at final print size.
FIG5_W <- NC_W_DOUBLE
FIG5_H <- 8.40
# Attach the significance brackets (none are drawn if nothing is significant).
p_dot <- add_sig_brackets(p_dot, dom_brackets, dom_vrange)

save_plot_file("Figure5_Combined_RegionalDominance", fig5, FIG5_W, FIG5_H)

# ------------------------------------------------------------------------------
# 29b. SEPARATE PANEL EXPORTS (Figures 5A, 5B, 5C)
# ------------------------------------------------------------------------------
# Each panel is also exported on its own at full double-column width. Because a
# standalone panel has the whole width to itself rather than a share of the
# composite, titles are restated and the panels are given generous margins so
# no label, legend or annotation can collide.
panel_standalone_theme <- ggplot2::theme(
  plot.margin   = ggplot2::margin(10, 12, 8, 10),
  plot.title    = ggplot2::element_text(size = PUB_TITLE, face = "bold",
                                        family = base_family_global,
                                        margin = ggplot2::margin(b = 5)),
  legend.margin = ggplot2::margin(3, 3, 3, 3)
)

# Standalone titles carry no figure number, so each panel reads independently
# when reproduced on its own.
#
# NOTE ON WORDING: the requested examples were "National / Provincial /
# Regional Detector". Those do not describe these panels -- Stage 4 is entirely
# a regional analysis (17 Philippine administrative regions), and each
# panel differs by *what is shown* (map / dominance probability), not by
# administrative tier. There is no national or provincial stratification
# anywhere in this stage.
# Nature Communications style: lowercase letter, no trailing period.
# The per-region metrics table is not shown as a panel in Figure 5.
# The dominance-probability plot is panel b.
fig5a <- p_map + panel_standalone_theme +
  ggplot2::labs(title = "a    Regional detector map")
fig5b <- p_dot + panel_standalone_theme +
  ggplot2::labs(title = "b    Per-detector regional dominance probabilities")

# Heights are per-panel: the map is near-square, the table is tall and
# row-driven, the dot plot is short and wide.
save_plot_file("Figure5_panel_a_DetectorMap",          fig5a, NC_W_DOUBLE, 5.60)
save_plot_file("Figure5_panel_b_DominanceProbability", fig5b, NC_W_DOUBLE, 4.40)

# ------------------------------------------------------------------------------
# 30. WRITE FIGURE 5 CSV TABLES
# ------------------------------------------------------------------------------
readr::write_csv(winner_rows,
                 file.path(OUTPUT_DIR, "Stage4_Detector_Map_RegionTable.csv"))
readr::write_csv(level1_results_ordered,
                 file.path(OUTPUT_DIR, "Stage4_Detector_Map_PerMetricSignificance.csv"))
readr::write_csv(level2_results,
                 file.path(OUTPUT_DIR, "Stage4_Detector_Map_CrossRegion_Wilcoxon.csv"))
readr::write_csv(
  panel_c_long %>%
    dplyr::transmute(REGION, Detector = as.character(Detector_Label),
                     Dominance_P, is_leader_dot, consensus_tier,
                     Observed_Winner),
  file.path(OUTPUT_DIR, "Stage4_Detector_Map_PanelC_DotPlot.csv")
)
readr::write_csv(panel_d_data,
                 file.path(OUTPUT_DIR, "Stage4_Detector_Map_PanelD_HeatmapData.csv"))
readr::write_csv(audit_df,
                 file.path(OUTPUT_DIR, "Stage4_Detector_Map_JoinAudit.csv"))

# ==============================================================================
# 31. LEGENDS
# ==============================================================================
# Two text files written to OUTPUT_DIR:
#   Stage4_Figure4_Legend.txt   (one block per metric + method summary)
#   Stage4_Figure5_Legend.txt    (single block for the detector map figure)
# ------------------------------------------------------------------------------

build_figure4_metric_legend <- function(metric_id) {
  metric_label <- hh_metric_full_name[[metric_id]]
  res <- all_wilcoxon_results %>%
    dplyr::filter(Metric == metric_id) %>%
    dplyr::arrange(Comparison)
  n_pairs <- max(res$N_Regions_Paired, na.rm = TRUE)
  sig_pairs <- res %>%
    dplyr::filter(!is.na(Significant_005), Significant_005) %>%
    dplyr::pull(Comparison) %>% as.character()
  outcome <- if (length(sig_pairs) == 0L) {
    paste0("None of the three detector-pair comparisons reached the per-pairwise ",
           "significance threshold (p < ", sprintf("%.2f", HH_ALPHA),
           ") on the available evidence base of ", n_pairs,
           " paired regional observations.")
  } else if (length(sig_pairs) == 3L) {
    "All three detector-pair comparisons reached the per-pairwise significance threshold."
  } else {
    paste0(length(sig_pairs), " of 3 detector-pair comparisons reached the ",
           "per-pairwise significance threshold: ",
           paste(sig_pairs, collapse = "; "), ".")
  }
  paste0(
    "Figure 4 (", metric_label, ") | Regional generalisability of ",
    metric_label, " across the ", length(region_ppv_order),
    " evaluable Philippine regions, 2016-2024 (excluding 2020, 2021, ",
    "2025; n = ", length(EVALUABLE_YEARS), " evaluable years per region; ",
    "regional inclusion floor = ", MIN_PEAK_CASES_PER_YEAR,
    " annual peak cases and >= ", MIN_EVALUABLE_YEARS_PER_REGION,
    " evaluable years). This figure is one of the five per-metric panels ",
    "in the Figure 4 set (TAM, Number of True Alarms, Sensitivity, Mean Lead ",
    "Time, Warning Persistence). The early-warning timeliness metrics (Mean ",
    "Lead Time, Warning Persistence) use the same-denominator zero-coerced ",
    "scheme: years where the metric is not computable due to no qualifying ",
    "triggers contribute zero rather than being excluded. Sensitivity is A1-",
    "restricted (proportion of evaluable seasons with at least one True Alarm ",
    "in the Actionable Window).\n",
    "(a) Regional dominance matrix (method-centric layout). Rows = 11 ",
    "outbreak-detection methods, ordered top-to-bottom by sweep count ",
    "(descending; ties broken by descending mean score). Columns = ",
    length(region_ppv_order), " evaluable regions, ordered left-to-right by ",
    "region mean dominance score across the 11 methods (descending). Cell ",
    "shading encodes the normalized dominance score on this metric, computed ",
    "by min-max normalization within each region across detectors: 1.00 = ",
    "within-region leader on ", metric_label, "; 0.00 = within-region laggard. ",
    "Single-color blue ramp from white (0) through mid-blue (0.50) to dark ",
    "navy (1.00); a horizontal blue-intensity legend strip below the matrix ",
    "maps shade to score with reference labels at 0.00 (weak), 0.50 ",
    "(moderate), ", sprintf("%.2f", DOMINANCE_THRESHOLD),
    " (dominant), and 1.00 (sweep). The right-most annotation reports the ",
    "per-METHOD count of regions swept by that method on this metric (regions ",
    "where the method achieves dominance score >= ",
    sprintf("%.2f", DOMINANCE_THRESHOLD), "), out of ",
    length(region_ppv_order), " total. ",
    "(b) Per-detector regional values on ", metric_label,
    ". One dot per region per detector, jittered horizontally; mean across ",
    "regions shown as a horizontal crossbar; vertical lines indicate year-",
    "cluster bootstrap 95% CIs (B = ", BOOT_N_CI, " replicates) where ",
    "available. The pale-blue band marks the operationally favourable zone ",
    "where defined. ",
    "(c) Wilcoxon detector-paired significance bars. Three pairwise detector ",
    "comparisons: Constant TA vs Outbreak Threshold; Continuous TA vs ",
    "Outbreak Threshold; Constant TA vs Continuous TA. Wilcoxon signed-rank, ",
    "paired by region (n = ", n_pairs,
    "), two-sided, evaluated PER PAIRWISE COMPARISON at alpha = ",
    sprintf("%.2f", HH_ALPHA),
    " (no across-pair multiplicity correction). Bars show |median ",
    "difference|; navy bars (with adjacent '*') indicate p < ",
    sprintf("%.2f", HH_ALPHA), "; pale grey otherwise. ",
    "Outcome on this metric: ", outcome, "\n"
  )
}

build_figure4_summary_legend <- function() {
  TARGETS_LOCAL <- c("Constant Transmission Acceleration",
                     "Continuous Transmission Acceleration",
                     "Outbreak Threshold")
  PAIRS_LOCAL <- list(
    list(a = "Constant Transmission Acceleration",
         b = "Continuous Transmission Acceleration",
         label = "Constant TA vs Continuous TA"),
    list(a = "Constant Transmission Acceleration",
         b = "Outbreak Threshold",
         label = "Constant TA vs Outbreak Threshold"),
    list(a = "Continuous Transmission Acceleration",
         b = "Outbreak Threshold",
         label = "Continuous TA vs Outbreak Threshold")
  )
  pivot_score_lg <- regional_dominance_long %>%
    dplyr::filter(as.character(Method) %in% TARGETS_LOCAL,
                  as.character(Metric) %in% HH_METRICS_REGIONAL) %>%
    dplyr::select(REGION, Metric, Method, Dominance_Score) %>%
    dplyr::mutate(Method = as.character(Method)) %>%
    tidyr::pivot_wider(names_from = Method, values_from = Dominance_Score)
  sig_pairs_summary <- character(0)
  n_pairs_typical <- NA_integer_
  for (pr in PAIRS_LOCAL) {
    a_sc <- pivot_score_lg[[pr$a]]; b_sc <- pivot_score_lg[[pr$b]]
    keep <- !is.na(a_sc) & !is.na(b_sc)
    a_sc <- a_sc[keep]; b_sc <- b_sc[keep]
    if (length(a_sc) >= 2L && any(a_sc != b_sc)) {
      pv <- suppressWarnings(stats::wilcox.test(a_sc, b_sc, paired = TRUE,
                                                exact = FALSE))$p.value
      if (!is.na(pv) && pv < 0.05) {
        sig_pairs_summary <- c(sig_pairs_summary, pr$label)
      }
    }
    n_pairs_typical <- length(a_sc)
  }
  outcome <- if (length(sig_pairs_summary) == 0L) {
    paste0("None of the three detector-pair comparisons reached p < 0.05 on ",
           "the available evidence base of approximately ",
           n_pairs_typical, " region\u00D7metric paired observations.")
  } else if (length(sig_pairs_summary) == 3L) {
    "All three detector-pair comparisons reached p < 0.05."
  } else {
    paste0(length(sig_pairs_summary), " of 3 detector-pair comparisons ",
           "reached p < 0.05: ",
           paste(sig_pairs_summary, collapse = "; "), ".")
  }

  paste0(
    "Figure 4 (Method-Across-Metrics Summary) | Cross-metric summary of ",
    "detector performance across the 17 Philippine regions and the five ",
    "reported operational metrics (TAM, Number of True Alarms, Sensitivity, ",
    "Mean Lead Time, Warning Persistence), 2016-2024 (excluding 2020, 2021, ",
    "2025; n = ", length(EVALUABLE_YEARS),
    " evaluable years per region; regional inclusion floor = ",
    MIN_PEAK_CASES_PER_YEAR, " annual peak cases and >= ",
    MIN_EVALUABLE_YEARS_PER_REGION, " evaluable years). This composite ",
    "figure complements the five per-metric multipanels by aggregating ",
    "detector performance across all five metrics simultaneously, ",
    "answering: which detector achieves the best operational performance ",
    "across the broadest combination of metrics and regions?\n",
    "(a) Regional dominance matrix on five operational metrics. Rows = 11 ",
    "outbreak-detection methods, ordered top-to-bottom by mean sweep count ",
    "across the 5 metrics (descending; best-overall method at top). Columns ",
    "= the five reported metrics, ordered left-to-right as TAM, Number of ",
    "True Alarms, Sensitivity, Mean Lead Time, Warning Persistence. Cell ",
    "shading encodes the per-(method, metric) sweep count: the number of ",
    "regions (out of ", MAX_REGIONS_POSSIBLE,
    " evaluable) where that method achieves within-region dominance score ",
    ">= ", sprintf("%.2f", DOMINANCE_THRESHOLD),
    " on that metric. Single-color blue ramp from white (0 regions) through ",
    "mid-blue to dark navy (all ", MAX_REGIONS_POSSIBLE, " regions). ",
    "Numeric annotation inside each cell shows the raw sweep count (auto-",
    "contrast: dark text on light fills, white text on dark fills). Column ",
    "headers (TAM, Number of True Alarms, Sensitivity, Mean Lead Time, ",
    "Warning Persistence; and the rightmost Mean Sweep / z-score header) ",
    "are rendered in-cell as text inside neutral-grey header tiles ",
    "immediately above the top method row. The right-most data column ",
    "reports per-method aggregate values: the mean sweep count across the ",
    "5 metrics, followed by a parenthetical z-score computed across the 11 ",
    "methods. Classification thresholds: z >= +1.0 = 'Sweeper'; z <= -1.0 = ",
    "'Non-sweeper'; otherwise 'Average'. ",
    "(b) Per-method dot plot. One dot per metric per method (5 dots per row, ",
    "jittered horizontally for visibility) showing the raw sweep count on ",
    "that metric; method-mean across the 5 metrics shown as a horizontal ",
    "crossbar. ",
    "(c) Wilcoxon detector-paired significance bars (DUAL-STATISTIC). Three ",
    "pairwise detector contrasts: Constant TA vs Continuous TA; Constant TA ",
    "vs Outbreak Threshold; Continuous TA vs Outbreak Threshold. The bar ",
    "length per contrast is |median difference in regions-swept count across ",
    "the 5 metrics| (units: regions; computed from the 5 paired metric-level ",
    "counts per detector). The p-value per contrast is computed from a ",
    "separate properly-powered test: paired Wilcoxon signed-rank on within-",
    "region Dominance_Score, paired by (REGION x METRIC) cells, n_pairs ~ ",
    n_pairs_typical, " region\u00D7metric paired observations. Bars are ",
    "navy (with adjacent '*') if p < 0.05; pale grey otherwise. The dual ",
    "specification answers two distinct questions in one panel: is the ",
    "detector difference statistically real (n_pairs ~ ", n_pairs_typical,
    " powered test), and how big is the practical difference in operationally ",
    "interpretable units (n=5 sweep-count median across the 5 metrics)? ",
    "Outcome: ", outcome, "\n"
  )
}

# Combined Figure 4 legend (one block per metric + Method Summary).
fig4_legend_lines <- c(
  paste0("================================================================"),
  paste0("FIGURE 4  |  LEGENDS  (per-metric and method summary)"),
  paste0("================================================================")
)
for (m in HH_METRICS_REGIONAL) {
  fig4_legend_lines <- c(
    fig4_legend_lines, "",
    paste0("--- Figure4_", hh_metric_slug[[m]], " ---"),
    build_figure4_metric_legend(m)
  )
}
fig4_legend_lines <- c(
  fig4_legend_lines, "",
  "--- Figure4_Method_Summary ---",
  build_figure4_summary_legend()
)
fig4_legend_text <- paste(fig4_legend_lines, collapse = "\n")
writeLines(fig4_legend_text,
           file.path(OUTPUT_DIR, "Stage4_Figure4_Legend.txt"))
cat("\nSaved: ", file.path(OUTPUT_DIR, "Stage4_Figure4_Legend.txt"),
    "\n", sep = "")

# Figure 5 legend
build_figure5_legend <- function() {
  bs <- N_BOOTS
  n_total       <- nrow(map_df)
  n_decisive    <- sum(map_df$Dominance_Probability >= DOMINANCE_THRESHOLD)
  n_contested_d <- sum(map_df$Dominance_Probability <  DOMINANCE_THRESHOLD)
  n_total_b     <- nrow(winner_rows)
  n_signif_b    <- sum(winner_rows$sig_star_dominance %in% c("*", "**", "***"),
                       na.rm = TRUE)

  cross_region_lines <- if (nrow(level2_results) > 0L) {
    apply(level2_results, 1, function(r) {
      sprintf("        %-25s | %-30s | n=%s | median diff = %s | Bonf p = %s  %s",
              r["Metric"], r["Comparison"], r["n_pairs"],
              ifelse(is.na(r["median_diff"]), "NA",
                     sprintf("%.3f", as.numeric(r["median_diff"]))),
              ifelse(is.na(r["p_bonf"]), "NA",
                     ifelse(as.numeric(r["p_bonf"]) < 0.001, "< 0.001",
                            sprintf("%.3f", as.numeric(r["p_bonf"])))),
              r["sig_star"])
    })
  } else {
    "        (Cross-region Wilcoxon table not available.)"
  }
  cross_region_block <- paste(cross_region_lines, collapse = "\n")

  paste0(
    "================================================================\n",
    "FIGURE 5  |  LEGEND\n",
    "================================================================\n\n",

    "Figure 5 | Bootstrap-supported regional dominance of three target ",
    "outbreak detectors across the ", n_total, " administrative regions of ",
    "the Philippines. The figure has two panels: (a) a choropleth map ",
    "encoding the four-tier regional consensus, with plain region names ",
    "shown directly on the map face, and (b) a per-detector dot plot of ",
    "regional dominance probabilities. The per-region metric table ",
    "is not shown as a panel in the figure; those ",
    "values, together with the per-region per-metric significance results, ",
    "remain available as standalone CSV exports ",
    "(RegionTable.csv and Stage4_Detector_Map_PanelD_HeatmapData.csv).\n\n",

    "FOUR-TIER CONSENSUS CLASSIFICATION.\n",
    "Each region is classified into ONE of four tiers based on (i) the ",
    "Round 1 Observed_Winner - the detector with the lowest mean rank ",
    "across the five reported metrics in the point estimate - and (ii) ",
    "the Round 2 all-pairs head-to-head test on bootstrap dominance ",
    "win-counts (chi-squared, Fisher exact fallback when expected < 5, ",
    "Bonferroni-adjusted p < 0.05).\n",
    "    strong    : both pairs sig in the Observed_Winner's favour.\n",
    "    partial   : exactly ONE pair sig in the Observed_Winner's favour\n",
    "                (the other ns), no rival is sig against it.\n",
    "    lead_only : Observed_Winner has zero sig wins, zero sig losses,\n",
    "                AND P_Observed_Winner >= 0.50.\n",
    "    contested : a rival is sig against the Observed_Winner, OR\n",
    "                P_Observed_Winner < 0.50.\n\n",

    "DUAL BONFERRONI REPORTING.\n",
    "  Primary  : factor 3 within-region * factor ", n_total_b,
    " between-region\n",
    "             = ", N_PAIRS * n_total_b,
    ". Used to assign the consensus tier and\n",
    "             render the figure.\n",
    "  Sensitivity (within-region only): factor 3. The within-region-\n",
    "             only adjustment is reported as a parallel column\n",
    "             ('p_*_bonf_within' suffixes) in RegionTable.csv. Under\n",
    "             the within-only correction, ", n_strong_within,
    " regions would meet 'strong' (vs ", n_strong, " under primary).\n\n",

    "Consensus rule for the displayed winner.\n",
    "    Round 1: bootstrap argmax of COMPOSITE MEAN-RANK across the 5 ",
    "reported metrics (TAM, N_True_Alarms, Sensitivity, Mean_Lead_Time, ",
    "WP). For each replicate, each detector is ranked 1-3 on each metric ",
    "(1 = best, ties averaged); the replicate winner is the detector with ",
    "the lowest mean rank across the 5 metrics. PPV is NOT used. ",
    "Observed_Winner is the same composite rule applied to the point-",
    "estimate metrics.\n",
    "    Round 2: all-pairs head-to-head dominance probability. Pairs ",
    "tested in fixed order: Constant TA vs OT; Continuous TA vs OT; ",
    "Constant TA vs Continuous TA. Each via one-sided chi-squared on the ",
    "2x2 win-counts table (Fisher exact fallback when min expected < 5; ",
    "N = ", bs, " bootstrap replicates per region).\n",
    sprintf("    Tier counts: strong=%d, partial=%d, lead_only=%d, contested=%d (of %d).\n\n",
            n_strong, n_partial, n_lead_only, n_contested, n_total_b),

    "(a) Choropleth map. The figure legend shows FOUR fill swatches: ",
    "Constant TA (blue), Continuous TA (green), Outbreak Threshold (red), ",
    "and Contested (grey). The first three swatches cover the 'strong' ",
    "and 'partial' consensus tiers (both rendered in the leader's detector ",
    "colour); the fourth grey swatch ('Contested') covers BOTH the ",
    "'lead_only' and 'contested' computational tiers. A second legend ",
    "strip ('Bootstrap confidence - color shade') shows the alpha gradient: ",
    "solid >= 0.90, mid >= 0.75, light < 0.75 (the Dominance_Probability ",
    "tiers).\n",
    "    Within the four-swatch fill scheme, the four computational ",
    "consensus tiers are visually distinguished by border style and by ",
    "the per-region label annotation:\n",
    "      strong    : detector colour, solid thin border, label =\n",
    "                  'REGION\\n(P=0.XX)'.\n",
    "      partial   : detector colour, DASHED thicker border, label =\n",
    "                  'REGION\\n(<Leader>, P=0.XX, partial)'.\n",
    "      lead_only : Contested grey fill, solid thin border, label =\n",
    "                  'REGION\\n(lead: <Leader>, P=0.XX)'.\n",
    "      contested : Contested grey fill, solid thin border, label =\n",
    "                  'REGION\\n(contested, P=0.XX)'.\n\n",

    "(b) Per-region metric table with EMBEDDED PER-METRIC SIGNIFICANCE. ",
    "Eight columns: Region, Detector, Dominance probability, TAM, N true ",
    "alarms, Sensitivity, Mean lead time (wk), WP (wk). Each metric cell ",
    "renders TWO stacked text lines: the metric value (top, larger font) ",
    "and a parenthetical significance label below (smaller font, color-",
    "coded). Star encoding: *** p<0.001 (dark navy); ** p<0.01 (mid navy); ",
    "* p<0.05 (mid blue); ns otherwise (faded grey); (-) when inference is ",
    "unavailable. Per-metric significance is the WEAKEST-LINK Bonferroni-",
    "adjusted one-sided p across the two ordered pairs the consensus winner ",
    "participates in (where the winner is X). For 'No consensus' rows ",
    "(lead_only and contested tiers) the fallback rule uses the smallest ",
    "p_bonf across all six ordered pairs (most informative signal regardless ",
    "of direction). The Detector cell is colored by detector palette for ",
    "strong (full saturation) and partial (washed-out); light grey for ",
    "lead_only; darker grey for contested. The Detector cell text annotates ",
    "the tier with a '(partial)' or '(lead)' suffix where applicable.\n\n",

    "(c) Per-detector dot plot. One dot per region per detector, jittered ",
    "horizontally. Median crossbar per detector. The pale band over [",
    sprintf("%.2f", DOMINANCE_THRESHOLD), ", 1.0] is the decisive zone. ",
    "Dominance threshold (",
    sprintf("%.2f", DOMINANCE_THRESHOLD), "). Border encoding distinguishes ",
    "the leader from the other detectors per region:\n",
    "      strong leader    : black border, thickest stroke\n",
    "      partial leader   : black border, medium stroke\n",
    "      lead_only leader : grey border, medium stroke\n",
    "      non-leader / contested: grey border, thinnest stroke\n\n",

    "AUXILIARY OUTPUTS.\n\n",

    "  Per-region per-metric significance heatmap data.\n",
    "  File: Stage4_Detector_Map_PanelD_HeatmapData.csv\n",
    "  Rows: one per (sub-grid label, region, metric) cell - up to 6 ordered ",
    "pairs * 5 metrics * 17 regions = 510 entries.\n",
    "  The 6 ordered pairs (each is its own one-sided test of H1: X > Y on ",
    "diffs = X - Y; the two ordered directions of a comparison give ",
    "DIFFERENT one-sided p-values that sum to ~1):\n",
    "      Constant TA vs Outbreak Threshold     (X - Y = Const - OT)\n",
    "      Continuous TA vs Outbreak Threshold   (X - Y = Cont - OT)\n",
    "      Constant TA vs Continuous TA          (X - Y = Const - Cont)\n",
    "      Outbreak vs Constant TA               (X - Y = OT - Const)\n",
    "      Outbreak vs Continuous TA             (X - Y = OT - Cont)\n",
    "      Continuous TA vs Constant TA          (X - Y = Cont - Const)\n",
    "  Columns include sig_star (***/**/*/ns/na), p_bonf (one-sided, ",
    "Bonferroni-adjusted within ordered pair), p_one_sided (raw), ",
    "median_diff_XY (signed difference in this direction), direction ",
    "(win/ns/na), and Subgrid_Label.\n\n",

    "  Per-region per-metric full results.\n",
    "  File: Stage4_Detector_Map_PerMetricSignificance.csv\n",
    "  Method: paired bootstrap on year-cluster replicates, 6 ordered pairs.\n",
    sprintf("  Bonferroni: factor %d within (REGION, ordered_pair).\n",
            BONF_LEVEL1_FACTOR),

    "\n  Cross-region per-metric significance.\n",
    "  File: Stage4_Detector_Map_CrossRegion_Wilcoxon.csv\n",
    "  Population-level Wilcoxon signed-rank test paired by region (n = ",
    if (nrow(level2_results) > 0L)
      max(level2_results$n_pairs, na.rm = TRUE) else "?",
    "). 5 metrics x 3 pairs = 15 tests.\n",
    "  Bonferroni factor 3 within metric.\n\n",
    "  Results (Bonferroni-adjusted, k=3 within metric):\n",
    cross_region_block, "\n\n",

    "STATISTICS SUMMARY:\n",
    sprintf("    Decisive regions  (Pr >= %.2f) : %d / %d\n",
            DOMINANCE_THRESHOLD, n_decisive, n_total),
    sprintf("    Contested by Pr threshold (Pr < %.2f) : %d / %d\n",
            DOMINANCE_THRESHOLD, n_contested_d, n_total),
    sprintf("    Tier counts: strong=%d, partial=%d, lead_only=%d, contested=%d (of %d).\n",
            n_strong, n_partial, n_lead_only, n_contested, n_total_b),
    sprintf("    Within-region-only Bonferroni sensitivity:\n"),
    sprintf("        regions meeting 'strong' under within-only: %d (vs %d primary).\n",
            n_strong_within, n_strong),
    sprintf("    Dominance Sig. p<0.05 (consensus_dominator): %d / %d.\n",
            n_signif_b, n_total_b),
    sprintf("    Per-region per-metric tests: %d.\n",
            nrow(level1_results_ordered)),
    sprintf("    Per-region per-metric significant (one-sided p_bonf<0.05): %d / %d.\n",
            sum(level1_results_ordered$sig_star %in% c("*", "**", "***"),
                na.rm = TRUE),
            sum(!is.na(level1_results_ordered$p_bonf))),
    sprintf("    Cross-region tests: %d. Significant (Bonf k=%d): %d / %d.\n\n",
            nrow(level2_results), BONF_LEVEL2_FACTOR,
            sum(level2_results$sig_star %in% c("*", "**", "***"), na.rm = TRUE),
            nrow(level2_results)),

    "METHODS NOTE.\n",
    "Year-cluster (Cameron-Gelbach-Miller) bootstrap, B = ", bs,
    " replicates per region. The Round 1 composite mean-rank winner, the ",
    "Round 2 all-pairs head-to-head dominance check, the per-region per-",
    "metric ordered-pair tests, and the cross-region per-metric Wilcoxon ",
    "are all computed in this script. Joins between geometry and metrics ",
    "use canonical region keys (NCR, CAR, MIMAROPA, Region I-XIII, BARMM); ",
    "a province-to-region lookup ensures GADM/Natural Earth province-level ",
    "sources dissolve cleanly to the 17 administrative regions.\n",
    "================================================================\n"
  )
}

fig5_legend_text <- build_figure5_legend()
writeLines(fig5_legend_text,
           file.path(OUTPUT_DIR, "Stage4_Figure5_Legend.txt"))
cat("Saved: ", file.path(OUTPUT_DIR, "Stage4_Figure5_Legend.txt"),
    "\n", sep = "")

# Echo legends to console.
cat("\n", fig4_legend_text, "\n", sep = "")
cat("\n", fig5_legend_text, "\n", sep = "")

# ==============================================================================
# 32. PARAMETER REPORT AND END-OF-RUN BANNER
# ==============================================================================
cat("\n============================================================\n")
cat("STAGE 4 - REGIONAL ANALYSIS  ", SCRIPT_VERSION, "\n", sep = "")
cat("============================================================\n")
cat("Anchor framework parameters (two-anchor; T = A1 union A2):\n")
cat("  A1 lead window         = [", A1_LEAD_MIN, ", ", A1_LEAD_MAX,
    "] weeks\n", sep = "")
cat("  A2 burden fraction     = ", A2_BURDEN_FRAC, "\n", sep = "")
cat("  True-Alarm rule        = t in A1 OR t in A2\n")
cat("\nDetector parameters:\n")
cat("  Continuous TA threshold = ", ETA_ON_CLASSIC, "\n", sep = "")
cat("  Constant TA   eta_ON    = ", ETA_ON,         "\n", sep = "")
cat("  Constant TA   eta_OFF   = ", ETA_OFF,        "\n", sep = "")
cat("  STA / LTA / GUARD       = ", STA_WIN, " / ", LTA_WIN, " / ",
    GUARD, "\n", sep = "")
cat("\nRegional inclusion:\n")
cat("  MIN_PEAK_CASES_PER_YEAR        = ", MIN_PEAK_CASES_PER_YEAR,
    "\n", sep = "")
cat("  MIN_EVALUABLE_YEARS_PER_REGION = ", MIN_EVALUABLE_YEARS_PER_REGION,
    "\n", sep = "")
cat("  Excluded years                 = ",
    paste(EXCLUDED_YEARS, collapse = ", "), "\n", sep = "")
cat("\nBootstrap configuration:\n")
cat("  B = ", BOOT_N_CI, " year-cluster replicates per region\n", sep = "")
cat("  Wilcoxon detector-paired alpha (per pairwise comparison) = ",
    sprintf("%.2f", HH_ALPHA), "\n", sep = "")
cat("  Per-region per-metric Bonferroni factor                  = ",
    BONF_LEVEL1_FACTOR, " (within REGION x ordered pair)\n", sep = "")
cat("  Cross-region per-metric Bonferroni factor                = ",
    BONF_LEVEL2_FACTOR, " (within metric)\n", sep = "")

cat("\n=== Saved files (", OUTPUT_DIR, ") ===\n", sep = "")
cat("Figure 4 (per-metric, 5 figures + 1 method summary):\n")
cat("  Figure4_TAM.{pdf,png}\n")
cat("  Figure4_N_True_Alarms.{pdf,png}\n")
cat("  Figure4_Sensitivity.{pdf,png}\n")
cat("  Figure4_Mean_Lead_Time.{pdf,png}\n")
cat("  Figure4_Warning_Persistence.{pdf,png}\n")
cat("  Figure4_Method_Summary.{pdf,png}\n")
cat("Figure 5 (composite detector map):\n")
cat("  Figure5_Combined_RegionalDominance.{pdf,png}\n")
cat("  Figure5_panel_a_DetectorMap.{pdf,png}\n")
cat("  Figure5_panel_b_DominanceProbability.{pdf,png}\n")
cat("Tables (Figure 4 supporting):\n")
cat("  Stage4_Regional_Framework_Metrics.csv\n")
cat("  Stage4_Regional_Framework_Metrics_with_CIs.csv\n")
cat("  Stage4_Regional_8Metric_Summary.csv\n")
cat("  Stage4_Regional_8Metric_Summary_with_CIs.csv\n")
cat("  Stage4_Regional_Dominance_Matrix.csv\n")
cat("  Stage4_Regional_Dominance_Probabilities.csv\n")
cat("  Stage4_Regional_Wilcoxon_PerMetric.csv\n")
cat("  Stage4_Regional_Wilcoxon_ConstantTA_vs_Comparators.csv\n")
cat("  Stage4_Method_Summary_Long.csv\n")
cat("  Stage4_Method_Summary_Aggregate.csv\n")
cat("  Stage4_Regional_Bootstrap_Replicates.csv\n")
cat("Tables (Figure 5 supporting):\n")
cat("  Stage4_Detector_Map_RegionTable.csv\n")
cat("  Stage4_Detector_Map_PerMetricSignificance.csv\n")
cat("  Stage4_Detector_Map_CrossRegion_Wilcoxon.csv\n")
cat("  Stage4_Detector_Map_PanelC_DotPlot.csv\n")
cat("  Stage4_Detector_Map_PanelD_HeatmapData.csv\n")
cat("  Stage4_Detector_Map_JoinAudit.csv\n")
cat("Legends:\n")
cat("  Stage4_Figure4_Legend.txt\n")
cat("  Stage4_Figure5_Legend.txt\n")

cat("\nDone.\n")

# ==============================================================================
# END OF SCRIPT — STAGE 4 REGIONAL ANALYSIS
# ==============================================================================
