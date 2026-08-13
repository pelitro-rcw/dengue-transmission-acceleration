# ==============================================================================
# STAGE 5 — COUNTRY ANALYSIS
# 8 endemic countries; Figure 6 and 7, plus tables
# ------------------------------------------------------------------------------
# Country-level generalisability of dengue outbreak detection methods across
# 8 dengue-endemic countries (Brazil, Colombia, Mexico, Peru, Philippines,
# Singapore, Sri Lanka, Taiwan) over 2016-2024, excluding the 2020-2021
# pandemic disruption and the truncated 2025 surveillance year.
#
# OUTPUTS
# -------
# Figure 6 (six multipanel figures, one PDF/PNG per file):
#   Per-metric country dominance for the five reported operational metrics
#   plus a cross-metric Method Summary figure:
#
#     Figure6_TAM
#     Figure6_N_True_Alarms
#     Figure6_Sensitivity
#     Figure6_Mean_Lead_Time
#     Figure6_Warning_Persistence
#     Figure6_Method_Summary
#
#   Each per-metric figure has three sub-panels:
#     a. Country dominance matrix  - 11 methods (rows) by 17 countries
#                                    (columns); cell shading by within-country
#                                    min-max-normalised dominance score; right
#                                    column reports the per-method count of
#                                    countries swept at score >= 0.75.
#     b. Per-detector dot plot      - one dot per country per detector;
#                                    bootstrap 95% CI bars; group-mean
#                                    crossbar.
#     c. Wilcoxon detector-paired   - three pairwise contrasts among
#                                    Constant TA, Continuous TA, Outbreak
#                                    Threshold; paired by country (n = 17).
#
# Figure 7 (single two-panel composite figure):
#   a. Per-country metric table with embedded per-metric significance.
#   b. Per-detector dot plot of country dominance probabilities, with a
#      per-country tier-legend overlay panel.
#
# Tables (CSV):
#   Stage5_Country_Framework_Metrics.csv
#   Stage5_Country_Framework_Metrics_with_CIs.csv
#   Stage5_Country_8Metric_Summary.csv
#   Stage5_Country_8Metric_Summary_with_CIs.csv
#   Stage5_Country_Dominance_Matrix.csv
#   Stage5_Country_Dominance_Probabilities.csv
#   Stage5_Country_Wilcoxon_PerMetric.csv
#   Stage5_Country_Wilcoxon_ConstantTA_vs_Comparators.csv
#   Stage5_Country_ConstantTA_4Comparator_Consensus_Long.csv
#   Stage5_Country_ConstantTA_4Comparator_Consensus.csv
#   Stage5_Method_Summary_Long.csv
#   Stage5_Method_Summary_Aggregate.csv
#   Stage5_Country_Bootstrap_Replicates.csv
#   Stage5_Country_Summary_CountryTable.csv
#   Stage5_Country_Summary_PerMetricSignificance.csv
#   Stage5_Country_Summary_CrossCountry_Wilcoxon.csv
#   Stage5_Country_Summary_PanelB_DotPlot.csv
#   Stage5_Country_Summary_PanelB_AuxWilcoxon.csv
#   Stage5_Country_Summary_PanelD_HeatmapData.csv
#   Stage5_Figure6_Legend.txt
#   Stage5_Figure7_Legend.txt
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
#   Mean_Lead_Time and WP per (country, detector) divide by the count of
#   evaluable years; years that do not produce a qualifying trigger
#   contribute zero rather than NA.
#
# A1-restricted Sensitivity:
#   years_with_A1_true / years_evaluable.
#
# Year-cluster bootstrap (Cameron, Gelbach & Miller 2008):
#   B = 1000 replicates per country; year is the cluster unit because
#   within-year weekly observations are not independent.
#
# Country inclusion criteria:
#   Annual peak >= 30 cases per (country, year);
#   >= 5 evaluable years per country.
#
# REQUIREMENTS
# ------------
# R packages (auto-installed if missing):
#   readxl, dplyr, tidyr, purrr, ggplot2, zoo, ISOweek, scales, tibble,
#   grid, patchwork, rlang, readr, stringr.
#
# Input: country weekly dengue dataset specified by `INPUT_PATH`,
# sheet `Country Data`, columns COUNTRY, YR, WN, DC_OPENDENGUE.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------------------------
SCRIPT_TITLE   <- "STAGE 5 - Country analysis (8 endemic countries; Figure 6 and 7, plus tables)"
SCRIPT_VERSION <- "1.0 (final)"

cat("\n=============================================================\n")
cat(SCRIPT_TITLE, "\n", sep = "")
cat("Version: ", SCRIPT_VERSION, "\n", sep = "")
cat("=============================================================\n\n")

REQUIRED_PACKAGES <- c(
  "readxl", "dplyr", "tidyr", "purrr", "ggplot2", "zoo", "ISOweek",
  "scales", "tibble", "grid", "patchwork", "rlang",
  "readr", "stringr", "ggrepel", "viridisLite"
)

# --- Project bootstrap -------------------------------------------------------
# Portable paths + shared publication theme. Replaces the previous inline
# install.packages() loop and the hard-coded Desktop path.
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

require_packages(REQUIRED_PACKAGES, purpose = "Stage 5")
invisible(lapply(REQUIRED_PACKAGES, function(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))))

source(file.path(DIR_R, "01_publication_theme.R"), local = TRUE)

# Shared detection framework: build_ears(), build_farrington(), build_ewars(),
# the anchor/compartment helpers and hh_paired_wilcoxon().
# local = TRUE is REQUIRED -- see the note in Stage 3. run_all.R runs each stage
# with sys.source(envir = env); sourcing without local = TRUE would place these
# functions in globalenv, where they cannot see this stage's own constants.
source(file.path(DIR_R, "02_detection_framework.R"), local = TRUE)


set.seed(GLOBAL_SEED)
options(scipen = 999)

# ------------------------------------------------------------------------------
# 1. USER PARAMETERS
# ------------------------------------------------------------------------------
# Path to the input country dengue dataset. Edit as appropriate, or use
# `file.choose()` interactively.
INPUT_PATH <- DATA_FILE
SHEET_NAME <- SHEET_COUNTRY

OUTPUT_DIR <- file.path(DIR_OUTPUT, "Stage5_Country_analysis")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

base_family_global <- PUB_FAMILY

# Delegates to save_pub(): vector PDF + 600 dpi PNG, clamped to print area.
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

# Country inclusion criteria
MIN_PEAK_CASES_PER_YEAR        <- 30L
MIN_EVALUABLE_YEARS_PER_COUNTRY <- 5L

# Bootstrap configuration
BOOT_N_CI <- 1000L

# Significance and Bonferroni configuration
SIG_LEVEL          <- 0.05
HH_ALPHA           <- 0.05
BONF_LEVEL1_FACTOR <- 5L  # 5 metrics within (country x ordered pair)
BONF_LEVEL2_FACTOR <- 3L  # 3 pairs within metric

# Within-country dominance score threshold
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

# The Country Data sheet uses DC_OPENDENGUE for case counts; the rest of the
# pipeline expects DC_DOH, so we rename at load. RF_NASA, FLAG_SINGLE_CELL_RF,
# and FLAG_TERMINAL_GAP may also be present but are not used here.
required_cols <- c("COUNTRY", "YR", "WN", "DC_OPENDENGUE")
missing_cols  <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns in sheet '", SHEET_NAME, "': ",
              paste(missing_cols, collapse = ", "),
              "\n  Expected Country Data columns: COUNTRY, YR, WN, DC_OPENDENGUE",
              "\n  (Optional auxiliary columns: RF_NASA, FLAG_SINGLE_CELL_RF,",
              " FLAG_TERMINAL_GAP)"))
}

df_raw <- df_raw %>% dplyr::rename(DC_DOH = DC_OPENDENGUE)

df_all <- df_raw %>%
  dplyr::mutate(
    COUNTRY = as.character(COUNTRY),
    YR     = suppressWarnings(as.integer(YR)),
    WN     = suppressWarnings(as.integer(WN)),
    DC_DOH = suppressWarnings(as.numeric(DC_DOH))
  ) %>%
  dplyr::filter(!is.na(COUNTRY), !is.na(YR), !is.na(WN),
                WN >= 1, WN <= 53) %>%
  dplyr::mutate(
    ISOweek = sprintf("%d-W%02d", YR, WN),
    Date    = ISOweek::ISOweek2date(paste0(ISOweek, "-1"))
  ) %>%
  dplyr::filter(!is.na(Date)) %>%
  dplyr::arrange(COUNTRY, Date) %>%
  dplyr::filter(YR %in% TARGET_YEARS)

if (nrow(df_all) == 0) {
  stop("No rows remain after filtering to target years.")
}

# Diagnostic banner: how many countries and which years per country.
cat("\n=== Country Data load summary ===\n")
.cy_summary <- df_all %>%
  dplyr::group_by(COUNTRY) %>%
  dplyr::summarise(
    n_years = dplyr::n_distinct(YR),
    years   = paste(sort(unique(YR)), collapse = ", "),
    n_weeks = dplyr::n(),
    .groups = "drop"
  )
print(as.data.frame(.cy_summary), row.names = FALSE)
cat("Countries loaded: ", nrow(.cy_summary), "\n\n", sep = "")

# Anchor functions reference a column named DC_QC. Add an alias so the
# country pipeline can use them verbatim.
df_all <- df_all %>% dplyr::mutate(DC_QC = DC_DOH)

# ------------------------------------------------------------------------------
# 5. SURGE GENERATION (PER COUNTRY)
# ------------------------------------------------------------------------------
# Per-country detector pipeline:
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

create_surges <- function(df_country) {
  df_country <- df_country %>% dplyr::arrange(Date)
  n <- nrow(df_country)

  # 5.1 Per-country rolling donor map -----------------------------------------
  country_target_years <- sort(unique(df_country$YR))
  donor_pool_country   <- country_target_years
  rolling_map <- list()
  for (y in country_target_years) {
    donors <- intersect(seq(y - 5L, y - 1L), donor_pool_country)
    if (length(donors) < 3L) {
      donors <- tail(donor_pool_country[donor_pool_country < y], 5L)
    }
    rolling_map[[as.character(y)]] <- sort(unique(donors))
  }

  # 5.2 Rolling weekly baselines (week-specific, donor-year based) -----------
  bl_rows <- list()
  for (y in names(rolling_map)) {
    donors <- rolling_map[[y]]; y_int <- as.integer(y)
    if (length(donors) == 0L) next
    weeks_y <- df_country %>% dplyr::filter(YR == y_int) %>%
      dplyr::pull(WN) %>% unique() %>% sort()
    if (length(weeks_y) == 0L) next
    for (w in weeks_y) {
      vals <- df_country %>% dplyr::filter(YR %in% donors, WN == w) %>%
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
  thresholds_country <- if (length(bl_rows) > 0L) {
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
  df_country <- df_country %>%
    dplyr::left_join(thresholds_country, by = c("YR", "WN"))

  # 5.3 Derived series -------------------------------------------------------
  dc_ma3  <- zoo::rollmean(df_country$DC_DOH, 3,  fill = NA, align = "right")
  dc_ma12 <- zoo::rollmean(df_country$DC_DOH, 12, fill = NA, align = "right")
  dc_diff <- c(NA_real_, diff(dc_ma3))

  var8  <- zoo::rollapply(df_country$DC_DOH, 8,
                          function(x) stats::var(x, na.rm = TRUE),
                          fill = NA, align = "right")
  var26 <- zoo::rollapply(df_country$DC_DOH, 26,
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
  df_country$rise_thr_rolling  <- rolling_quantile_by_target_local(
    dc_diff,   df_country$YR, rolling_map, 0.75)
  df_country$roc80_thr_rolling <- rolling_quantile_by_target_local(
    dc_diff,   df_country$YR, rolling_map, 0.80)
  df_country$ct80_thr_rolling  <- rolling_quantile_by_target_local(
    ratio_var, df_country$YR, rolling_map, 0.80)

  # 5.5 Walk-forward CUSUM (year-boundary reset + rolling baseline) ----------
  k_mult <- 0.5; h_mult <- 5.0
  cusum_vec <- rep(0, n)
  cusum_h   <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    is_year_boundary <- i > 1L &&
      !is.na(df_country$YR[i]) && !is.na(df_country$YR[i - 1L]) &&
      df_country$YR[i] != df_country$YR[i - 1L]
    mu0    <- df_country$bl_mean[i]
    sigma0 <- df_country$bl_sd[i]
    prev_val <- if (i == 1L || is_year_boundary) 0 else cusum_vec[i - 1L]
    if (!is.na(df_country$DC_DOH[i]) && !is.na(mu0) &&
        !is.na(sigma0) && sigma0 > 0) {
      cusum_vec[i] <- max(0, prev_val +
                            (df_country$DC_DOH[i] - mu0 - k_mult * sigma0))
    } else {
      cusum_vec[i] <- prev_val
    }
    cusum_h[i] <- if (!is.na(sigma0) && sigma0 > 0) h_mult * sigma0 else NA_real_
  }
  df_country$cusum_val <- cusum_vec
  df_country$cusum_h   <- cusum_h

  # 5.6 Primary detector signals (using rolling thresholds) ------------------
  df_country <- df_country %>%
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
  dc <- df_country$DC_DOH
  triggered_v <- rep(FALSE, n)
  R_vaezi_v   <- rep(NA_real_, n)
  is_on <- FALSE; frozen_lta <- NA_real_; consec_off <- 0L
  for (t in seq_len(n)) {
    if (t > 1L && !is.na(df_country$YR[t]) &&
        !is.na(df_country$YR[t - 1L]) &&
        df_country$YR[t] != df_country$YR[t - 1L]) {
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
  df_country$surge_sta_lta_vaezi <- as.integer(triggered_v)
  df_country$R_vaezi             <- R_vaezi_v

  # 5.8 Composite signal -----------------------------------------------------
  df_country$surge_composite <- as.integer(
    df_country$surge_mean_2sd == 1L &
      (df_country$surge_sta_lta_vaezi == 1L | df_country$surge_critical_trans == 1L)
  )

  # --- Contemporary comparators: Farrington, EARS, EWARS --------------------
  # Uses the same shared constructors as the regional analysis.
  # EWARS is an EWARS-style alarm-indicator model, not the WHO software.
  .cs <- ifelse(is.na(df_country$DC_DOH), 0, df_country$DC_DOH)
  .rn <- if ("RF_NASA" %in% names(df_country)) df_country$RF_NASA else rep(NA_real_, nrow(df_country))
  .ev_idx <- which(df_country$YR %in% EVALUABLE_YEARS)

  df_country$surge_ears <- build_ears(.cs)
  df_country$surge_farrington <- build_farrington(.cs, df_country$YR, .ev_idx)$alarm
  df_country$surge_ewars <- build_ewars(
    .cs, .rn, df_country$YR,
    eval_years = EVALUABLE_YEARS,
    excluded_years = EXCLUDED_YEARS
  )$alarm

  df_country
}

cat("Generating per-country surge signals...\n")
df_all <- df_all %>%
  dplyr::group_by(COUNTRY) %>%
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
# 8. COUNTRY INCLUSION CRITERIA
# ------------------------------------------------------------------------------
peak_per_country_year <- df_all %>%
  dplyr::filter(YR %in% EVALUABLE_YEARS) %>%
  dplyr::group_by(COUNTRY, YR) %>%
  dplyr::summarise(
    annual_peak = suppressWarnings(max(DC_DOH, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    annual_peak = ifelse(is.finite(annual_peak), annual_peak, 0)
  )

evaluable_country_year <- peak_per_country_year %>%
  dplyr::filter(annual_peak >= MIN_PEAK_CASES_PER_YEAR) %>%
  dplyr::select(COUNTRY, YR)

country_year_count <- evaluable_country_year %>%
  dplyr::count(COUNTRY, name = "n_evaluable_years")

countries_in <- country_year_count %>%
  dplyr::filter(n_evaluable_years >= MIN_EVALUABLE_YEARS_PER_COUNTRY) %>%
  dplyr::pull(COUNTRY)

countries_excluded <- setdiff(unique(df_all$COUNTRY), countries_in)

cat("Included countries: ", length(countries_in), "; excluded: ",
    length(countries_excluded), "\n", sep = "")

# ------------------------------------------------------------------------------
# 9. PER-COUNTRY FRAMEWORK METRICS
# ------------------------------------------------------------------------------
build_trigger_detail_for_country <- function(df_country, evaluable_years_for_country) {
  rows <- list()
  for (yr in evaluable_years_for_country) {
    anchors <- compute_anchors_for_year(df_country, yr)
    df_y <- df_country %>% dplyr::filter(YR == yr) %>% dplyr::arrange(WN)
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

compute_yearly_lead_for_country <- function(trig_aug, evaluable_years_for_country) {
  if (nrow(trig_aug) == 0)
    return(data.frame(
      Method = character(), Year = integer(),
      First_A1_True_Week = integer(), Lead_Time_Yr = numeric(),
      stringsAsFactors = FALSE
    ))
  trig_aug %>%
    dplyr::filter(Year %in% evaluable_years_for_country, InA1, IsTrue) %>%
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

compute_method_metrics_country <- function(trig_aug, yearly_lead,
                                          evaluable_years_for_country) {
  rows <- list()
  n_eval <- length(evaluable_years_for_country)
  for (m in names(surge_defs)) {
    trig_sub <- trig_aug %>%
      dplyr::filter(Method == m, Year %in% evaluable_years_for_country)
    lead_sub <- yearly_lead %>%
      dplyr::filter(Method == m, Year %in% evaluable_years_for_country)

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
      per_year_wp <- vapply(evaluable_years_for_country, function(y) {
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
    tam_per_year <- vapply(evaluable_years_for_country, function(y) {
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

# Per-country pipeline driver
country_metrics_list <- list()
for (ctry in countries_in) {
  yrs_ctry <- evaluable_country_year %>%
    dplyr::filter(COUNTRY == ctry) %>%
    dplyr::pull(YR) %>% sort()
  if (length(yrs_ctry) < MIN_EVALUABLE_YEARS_PER_COUNTRY) next
  df_ctry     <- df_all %>% dplyr::filter(COUNTRY == ctry)
  trig_aug   <- build_trigger_detail_for_country(df_ctry, yrs_ctry)
  yearly_ld  <- compute_yearly_lead_for_country(trig_aug, yrs_ctry)
  metrics    <- compute_method_metrics_country(trig_aug, yearly_ld, yrs_ctry)
  metrics$COUNTRY <- ctry
  country_metrics_list[[ctry]] <- metrics
}

country_metrics <- dplyr::bind_rows(country_metrics_list) %>%
  dplyr::left_join(method_type_map, by = "Method") %>%
  dplyr::mutate(
    Method   = factor(Method,   levels = method_order),
    Paradigm = factor(Paradigm, levels = type_order)
  )

if (nrow(country_metrics) == 0)
  stop("No country metrics were computed. Check inclusion criteria and data.")

# ------------------------------------------------------------------------------
# 10. YEAR-CLUSTER BOOTSTRAP (Cameron-Gelbach-Miller, B = 1000)
# ------------------------------------------------------------------------------
# For each country, resample evaluable years with replacement (year is the
# cluster unit because within-year weekly observations are not independent).
# On each replicate, recompute all framework metrics from the cached trigger
# detail subset. Trigger columns and anchors are deterministic given case
# data and are cached outside the bootstrap loop for speed.
# ------------------------------------------------------------------------------
cat("\n=== Year-cluster bootstrap (B =", BOOT_N_CI, "per country) ===\n")
cat("This computation may take a few minutes for ",
    length(countries_in), " countries x ", BOOT_N_CI, " replicates.\n", sep = "")

# Anchor cache (one entry per country/year).
anchor_cache <- list()
for (ctry in countries_in) {
  yrs_ctry <- evaluable_country_year %>%
    dplyr::filter(COUNTRY == ctry) %>%
    dplyr::pull(YR) %>% sort()
  if (length(yrs_ctry) < MIN_EVALUABLE_YEARS_PER_COUNTRY) next
  df_ctry <- df_all %>% dplyr::filter(COUNTRY == ctry)
  anchor_cache[[ctry]] <- list()
  for (yr in yrs_ctry) {
    anchor_cache[[ctry]][[as.character(yr)]] <- compute_anchors_for_year(df_ctry, yr)
  }
}

# Bootstrap-specific trigger detail builder using the anchor cache.
build_trigger_detail_cached <- function(df_country, evaluable_years_for_country,
                                        country_name) {
  rows <- list()
  for (yr in evaluable_years_for_country) {
    anchors <- anchor_cache[[country_name]][[as.character(yr)]]
    if (is.null(anchors)) next
    df_y <- df_country %>% dplyr::filter(YR == yr) %>% dplyr::arrange(WN)
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

# Pre-compute the FULL (unsampled) trigger detail per country.
full_trig_detail <- list()
for (ctry in names(anchor_cache)) {
  df_ctry  <- df_all %>% dplyr::filter(COUNTRY == ctry)
  yrs_ctry <- as.integer(names(anchor_cache[[ctry]]))
  full_trig_detail[[ctry]] <- build_trigger_detail_cached(df_ctry, yrs_ctry, ctry)
}

# Single-replicate metric computation.
bootstrap_metrics_one_replicate <- function(ctry, yrs_resampled) {
  trig_full <- full_trig_detail[[ctry]]
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
    trig_m_unique <- full_trig_detail[[ctry]] %>% dplyr::filter(Method == m)
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
for (ctry in names(anchor_cache)) {
  yrs_ctry <- as.integer(names(anchor_cache[[ctry]]))
  if (length(yrs_ctry) < 2L) {
    cat("  Skipping country '", ctry, "' (only ",
        length(yrs_ctry), " evaluable year(s)).\n", sep = "")
    next
  }
  reps <- vector("list", BOOT_N_CI)
  for (b in seq_len(BOOT_N_CI)) {
    yrs_resampled <- sample(yrs_ctry, size = length(yrs_ctry), replace = TRUE)
    reps[[b]] <- bootstrap_metrics_one_replicate(ctry, yrs_resampled)
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
    dplyr::mutate(COUNTRY = ctry)

  boot_results[[ctry]] <- list(replicates = reps_df, ci = ci_summary)
  cat("  Country '", ctry, "' bootstrap complete (",
      BOOT_N_CI, " replicates).\n", sep = "")
}

boot_ci <- dplyr::bind_rows(lapply(boot_results, function(x) x$ci))

# Merge CIs into the metrics table
country_metrics <- country_metrics %>%
  dplyr::mutate(COUNTRY = as.character(COUNTRY)) %>%
  dplyr::left_join(
    boot_ci %>% dplyr::mutate(Method = as.character(Method)),
    by = c("COUNTRY", "Method")
  ) %>%
  dplyr::mutate(COUNTRY = factor(COUNTRY, levels = unique(COUNTRY)))

# ------------------------------------------------------------------------------
# 11. DOMINANCE PROBABILITY (composite mean-rank across the 5 reported metrics)
# ------------------------------------------------------------------------------
# For each replicate, rank the three target detectors on each of the five
# reported metrics (rank 1 = best; ties get average rank; all metrics are
# higher-is-better), then compute the mean rank across the five metrics.
# The replicate winner is the detector with the lowest mean rank. The
# Dominance_Probability for a country is the proportion of replicates in
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
for (ctry in names(boot_results)) {
  reps_df     <- boot_results[[ctry]]$replicates
  reps_target <- reps_df %>% dplyr::filter(Method %in% TARGET_DETECTORS)

  rep_winners <- .composite_winner_per_rep(reps_target,
                                           COMPOSITE_DOMINANCE_METRICS)

  obs_pt <- country_metrics %>%
    dplyr::filter(COUNTRY == ctry, Method %in% TARGET_DETECTORS) %>%
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

  dominance_results[[ctry]] <- data.frame(
    COUNTRY = ctry,
    Observed_Winner = if (length(obs_winner) > 0L) obs_winner else NA_character_,
    Dominance_Probability = dom_prob,
    Bootstrap_N    = total_reps,
    Composite_Rule = "mean_rank_5metrics",
    stringsAsFactors = FALSE
  )
  # One P_* column per target detector, added programmatically.
  for (dn in TARGET_DETECTORS) {
    dominance_results[[ctry]][[unname(TARGET_P_COLS[dn])]] <-
      if (dn %in% names(win_tab)) as.numeric(win_tab[dn]) / total_reps else 0
  }
}
dominance_df <- dplyr::bind_rows(dominance_results) %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("P_"),
                              ~ ifelse(is.na(.), 0, .)))

cat("\n=== Bootstrap dominance probability (composite mean-rank, 5 metrics) ===\n")
print(as.data.frame(dominance_df), row.names = FALSE)
cat("\n")

# Join dominance probability into country_metrics (country-level property,
# inherited by every method row in that country).
country_metrics <- country_metrics %>%
  dplyr::mutate(COUNTRY = as.character(COUNTRY)) %>%
  dplyr::left_join(
    dominance_df %>%
      dplyr::select(COUNTRY, Observed_Winner, Dominance_Probability,
                    dplyr::all_of(unname(TARGET_P_COLS)), Bootstrap_N),
    by = "COUNTRY"
  ) %>%
  dplyr::mutate(COUNTRY = factor(COUNTRY, levels = unique(COUNTRY)))

# ------------------------------------------------------------------------------
# 12. ORDER COUNTRYS AND METHODS; ESTABLISH METRIC METADATA
# ------------------------------------------------------------------------------
country_ppv_order <- country_metrics %>%
  dplyr::group_by(COUNTRY) %>%
  dplyr::summarise(mean_ppv = mean(PPV, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(mean_ppv)) %>%
  dplyr::pull(COUNTRY)

country_metrics <- country_metrics %>%
  dplyr::mutate(COUNTRY = factor(COUNTRY, levels = rev(country_ppv_order))) %>%
  dplyr::mutate(
    Method   = factor(as.character(Method),   levels = method_order),
    Paradigm = factor(as.character(Paradigm), levels = type_order)
  )

# Five reported head-to-head metrics (Figure 6 figure set).
HH_METRICS_COUNTRY <- c(
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

# Filename slug per metric (Figure 6 outputs).
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

# Dedicated Constant TA head-to-head comparators for supplementary CSV outputs.
# This does not replace DETECTOR_PAIRS used by the existing Figure 6 panels.
CONSTANT_TA_COMPARATOR_PAIRS <- list(
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
# COUNTRY (n = 17). Significance is evaluated PER PAIRWISE COMPARISON at
# alpha = 0.05; we do NOT multiply p-values across the three pairs within
# a metric (each pair answers a distinct scientific question).
# ------------------------------------------------------------------------------
compute_detector_paired_wilcoxon <- function(
  metric_id,
  country_metrics_df,
  detector_pairs = DETECTOR_PAIRS
) {
  metric_long <- country_metrics_df %>%
    dplyr::transmute(
      COUNTRY = as.character(COUNTRY),
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
      N_Countries_Paired = n_pairs,
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
  lapply(HH_METRICS_COUNTRY, function(m) {
    compute_detector_paired_wilcoxon(m, country_metrics)
  })
)

cat("\n=== Wilcoxon detector-paired (5 metrics x 3 pairs = 15 rows) ===\n")
print(all_wilcoxon_results, row.names = FALSE, digits = 3)
cat("\n")

# Additional Constant TA versus four prespecified comparators.
# Five metrics x four comparisons = 20 rows.
constant_ta_comparator_wilcoxon <- dplyr::bind_rows(
  lapply(HH_METRICS_COUNTRY, function(m) {
    compute_detector_paired_wilcoxon(
      m,
      country_metrics,
      detector_pairs = CONSTANT_TA_COMPARATOR_PAIRS
    )
  })
)

cat("\n=== Constant TA head-to-head (5 metrics x 4 comparators = 20 rows) ===\n")
print(constant_ta_comparator_wilcoxon, row.names = FALSE, digits = 3)
cat("\n")

# ------------------------------------------------------------------------------
# 14. COUNTRY DOMINANCE MATRIX (LONG)
# ------------------------------------------------------------------------------
# For each metric, compute per-(country, detector) normalized dominance score
# in [0, 1] using min-max normalization within (country, metric) across the
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

compute_country_dominance_long <- function() {
  rows <- list()
  for (m in HH_METRICS_COUNTRY) {
    higher_better <- hh_metric_higher_better[[m]]
    sub <- country_metrics %>%
      dplyr::transmute(
        COUNTRY = as.character(COUNTRY),
        Method = as.character(Method),
        Paradigm = as.character(Paradigm),
        Value = .data[[m]]
      )
    sub <- sub %>%
      dplyr::group_by(COUNTRY) %>%
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
country_dominance_long <- compute_country_dominance_long()

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
country_dominance_long$cell_fill <- vapply(
  country_dominance_long$Dominance_Score, mix_blue_intensity, character(1)
)

# ------------------------------------------------------------------------------
# 15. METHOD-ACROSS-METRICS AGGREGATION
# ------------------------------------------------------------------------------
compute_method_summary_long <- function() {
  country_dominance_long %>%
    dplyr::filter(Metric %in% HH_METRICS_COUNTRY) %>%
    dplyr::group_by(Method, Metric) %>%
    dplyr::summarise(
      Sweep_Count = sum(Is_Dominant, na.rm = TRUE),
      N_Regions   = sum(!is.na(Dominance_Score)),
      .groups     = "drop"
    ) %>%
    dplyr::mutate(
      Method = factor(as.character(Method), levels = method_order),
      Metric = factor(as.character(Metric), levels = HH_METRICS_COUNTRY)
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

MAX_COUNTRYS_POSSIBLE <- {
  .mrp <- suppressWarnings(max(method_summary_long$N_Regions, na.rm = TRUE))
  if (!is.finite(.mrp) || .mrp <= 0) 17L else as.integer(.mrp)
}

mix_blue_count <- function(count, max_count = MAX_COUNTRYS_POSSIBLE) {
  if (is.na(count) || !is.finite(count)) return("#FFFFFF")
  frac <- pmax(0, pmin(1, count / max(1, max_count)))
  mix_blue_intensity(frac)
}
method_summary_long$cell_fill <- vapply(
  method_summary_long$Sweep_Count, mix_blue_count, character(1)
)

# ------------------------------------------------------------------------------
# 16. FIGURE 6 BUILDERS — PER-METRIC PANELS AND METHOD SUMMARY
# ------------------------------------------------------------------------------

# 16a. Country dominance matrix (per metric) ---------------------------------
build_country_dominance_matrix <- function(metric_id) {
  metric_label <- hh_metric_full_name[[metric_id]]
  this_metric_df <- country_dominance_long %>%
    dplyr::filter(Metric == metric_id)

  # Per-method count of countries swept at score >= threshold.
  per_method_count <- this_metric_df %>%
    dplyr::group_by(Method) %>%
    dplyr::summarise(
      N_Countries_Swept = sum(Is_Dominant, na.rm = TRUE),
      N_Countries_Total = sum(!is.na(Dominance_Score)),
      Mean_Score      = mean(Dominance_Score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Sweep_Label = paste0(N_Countries_Swept, "/", N_Countries_Total))

  method_sweep_order <- per_method_count %>%
    dplyr::arrange(dplyr::desc(N_Countries_Swept),
                   dplyr::desc(Mean_Score)) %>%
    dplyr::pull(Method) %>% as.character()

  country_score_order <- this_metric_df %>%
    dplyr::group_by(COUNTRY) %>%
    dplyr::summarise(mean_score = mean(Dominance_Score, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_score)) %>%
    dplyr::pull(COUNTRY) %>% as.character()

  this_metric_df <- this_metric_df %>%
    dplyr::mutate(
      COUNTRY = factor(as.character(COUNTRY), levels = country_score_order),
      Method = factor(as.character(Method), levels = rev(method_sweep_order))
    )
  per_method_count <- per_method_count %>%
    dplyr::mutate(Method = factor(as.character(Method),
                                  levels = rev(method_sweep_order)))

  N_COUNTRYS     <- length(country_score_order)
  SWEEP_COUNT_X <- N_COUNTRYS + 1L

  ggplot2::ggplot(this_metric_df) +
    ggplot2::geom_tile(
      ggplot2::aes(x = COUNTRY, y = Method),
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
      label = paste0("Countries Swept\n(score \u2265 ",
                     sprintf("%.2f", DOMINANCE_THRESHOLD), ")"),
      size = 2.2, family = base_family_global, fontface = "bold",
      colour = "#0B2447", lineheight = 0.92, vjust = 0
    ) +
    ggplot2::scale_x_discrete(
      limits = c(country_score_order, "__SWEEP_COUNT__"),
      labels = {
        raw <- as.character(c(country_score_order, " "))
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
                     "    Country dominance matrix on ", metric_label,
                     " (rows = methods, columns = countries)"),
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
      title = paste0("Dominance score (min-max normalized within country; ",
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
build_country_dot_panel <- function(metric_id) {
  metric_label  <- hh_metric_full_name[[metric_id]]
  ylab          <- hh_metric_y_label[[metric_id]]
  higher_better <- hh_metric_higher_better[[metric_id]]
  threshold     <- hh_metric_threshold[[metric_id]]
  is_proportion <- hh_metric_is_proportion[[metric_id]]
  ci_cols       <- hh_metric_ci_cols[[metric_id]]

  plot_df <- country_metrics %>%
    dplyr::transmute(
      Method, COUNTRY,
      Value = .data[[metric_id]]
    ) %>%
    dplyr::mutate(COUNTRY_chr = as.character(COUNTRY),
                  Method_chr = as.character(Method))

  if (all(ci_cols %in% names(boot_ci))) {
    plot_df <- plot_df %>%
      dplyr::left_join(
        boot_ci %>%
          dplyr::transmute(
            COUNTRY_chr = as.character(COUNTRY),
            Method_chr = as.character(Method),
            ci_lo = .data[[ci_cols[1]]],
            ci_hi = .data[[ci_cols[2]]]
          ),
        by = c("COUNTRY_chr", "Method_chr")
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
                     "    Per-detector country values on ", metric_label),
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
build_country_sig_panel <- function(metric_id, all_wilcoxon_results) {
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
                       n_max <- suppressWarnings(max(res$N_Countries_Paired,
                                                     na.rm = TRUE))
                       if (!is.finite(n_max)) "?" else as.character(n_max)
                     },
                     " countries)"),
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
build_country_multipanel <- function(metric_id) {
  top          <- build_country_dominance_matrix(metric_id)
  legend_strip <- build_dominance_legend()
  bot_l        <- build_country_dot_panel(metric_id)
  bot_r        <- build_country_sig_panel(metric_id, all_wilcoxon_results)
  bottom       <- (bot_l | bot_r) + patchwork::plot_layout(widths = c(7, 3))

  (top / legend_strip / bottom) +
    patchwork::plot_layout(heights = c(1.55, 0.18, 1.00)) +
    patchwork::plot_annotation(
      title = paste0("Country generalisability: ",
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
      Metric = factor(as.character(Metric), levels = HH_METRICS_COUNTRY)
    )
  agg <- method_summary_agg %>%
    dplyr::mutate(
      Method = factor(as.character(Method), levels = rev(method_summary_order)),
      Summary_Label = sprintf("%.1f (z=%+.2f)", Mean_Sweep, Z_Sweep)
    )

  N_METRICS <- length(HH_METRICS_COUNTRY)
  SUMMARY_X <- N_METRICS + 1L
  N_METHODS <- length(method_summary_order)
  HEADER_Y  <- N_METHODS + 1L

  metric_full_names <- vapply(HH_METRICS_COUNTRY,
                              function(m) hh_metric_full_name[[m]],
                              character(1))
  metric_cell_labels <- c(
    "TAM"            = "True-Alarm\nMagnitude",
    "N_True_Alarms"  = "Number of\nTrue Alarms",
    "Sensitivity"    = "Sensitivity",
    "Mean_Lead_Time" = "Mean\nLead Time",
    "WP"             = "Warning\nPersistence"
  )
  for (m in HH_METRICS_COUNTRY) {
    if (is.na(metric_cell_labels[m]) || is.null(metric_cell_labels[[m]])) {
      metric_cell_labels[m] <- metric_full_names[[m]]
    }
  }

  header_df <- data.frame(
    x_pos = c(seq_len(N_METRICS), SUMMARY_X),
    y_pos = HEADER_Y,
    label = c(unname(metric_cell_labels[HH_METRICS_COUNTRY]),
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
                   colour = Sweep_Count >= round(MAX_COUNTRYS_POSSIBLE * 0.55)),
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
      title = "Country dominance matrix on 5 operational metrics",
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
      limits = c(0, MAX_COUNTRYS_POSSIBLE + 0.5),
      breaks = pretty(c(0, MAX_COUNTRYS_POSSIBLE), n = 5),
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    ) +
    ggplot2::scale_y_discrete(
      limits = rev(method_summary_order),
      labels = unname(method_two_line[rev(method_summary_order)])
    ) +
    ggplot2::labs(
      title = "Per-method sweep counts across the 5 metrics",
      x = "Countries swept (count)", y = NULL
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
#   p-value      : paired Wilcoxon signed-rank on within-country Dominance_Score
#                  paired by (COUNTRY x METRIC). n_pairs ~ 17 x 5 = 85.
#                  Properly powered; p < 0.05 is reachable.
#   Bar length   : |median(diff in countries-swept count)| across the 5 reported
#                  metrics. Interpretable in operational units of "countries".
#
# The dual specification answers two distinct questions in the same panel:
# (i) is the difference statistically real (n=85 powered test); and (ii) how
# big is the practical difference in operational countries (n=5 sweep-count
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

  pivot_score <- country_dominance_long %>%
    dplyr::filter(as.character(Method) %in% TARGETS,
                  as.character(Metric) %in% HH_METRICS_COUNTRY) %>%
    dplyr::select(COUNTRY, Metric, Method, Dominance_Score) %>%
    dplyr::mutate(Method = as.character(Method)) %>%
    tidyr::pivot_wider(names_from = Method, values_from = Dominance_Score)

  pivot_count <- method_summary_long %>%
    dplyr::filter(as.character(Method) %in% TARGETS,
                  as.character(Metric) %in% HH_METRICS_COUNTRY) %>%
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
  cat("    Test: paired Wilcoxon on Dominance_Score (17 COUNTRYS x 5 METRICS, n ~ 85)\n")
  cat("    Bar:  |median(diff in sweep counts)| across the 5 metrics, units 'countries'\n")
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
                     n_typical, ", 17 countries \u00D7 5 metrics)"),
      x = "Median diff in countries swept (across 5 metrics)",
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
      title = paste0("Figure 6 | Method-across-metrics summary ",
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
# 17. SAVE FIGURE 6 (PER-METRIC + METHOD SUMMARY)
# ------------------------------------------------------------------------------
# Figure 6 set, authored at final print size.
fig6_w <- NC_W_DOUBLE
# Composite stacks a dominance matrix, a legend strip and two sub-panels;
# the taller supplementary canvas prevents the rows from compressing.
fig6_h <- NC_H_SUPP

cat("\n=== Building 5 per-metric multipanel figures (Figure 6 set) ===\n")
for (m in HH_METRICS_COUNTRY) {
  fig <- build_country_multipanel(m)
  slug <- hh_metric_slug[[m]]
  save_plot_file(paste0("Figure6_", slug), fig, fig6_w, fig6_h)
}

cat("\n=== Building Method-Across-Metrics Summary figure (Figure 6 set) ===\n")
{
  fig_summary <- build_method_summary_multipanel()
  save_plot_file("Figure6_Method_Summary", fig_summary, fig6_w, fig6_h)
}

# ------------------------------------------------------------------------------
# 18. WRITE FIGURE 6 CSV TABLES
# ------------------------------------------------------------------------------
# Wide country metrics table (carries the 5 reported metrics plus diagnostic
# columns). RFR / RFR_RX columns dropped (not used by this pipeline).
utils::write.csv(
  country_metrics %>% dplyr::mutate(COUNTRY = as.character(COUNTRY),
                                     Method = as.character(Method)),
  file.path(OUTPUT_DIR, "Stage5_Country_Framework_Metrics.csv"),
  row.names = FALSE
)

# 8-metric harmonized country CSV.
country_8metric_csv <- country_metrics %>%
  dplyr::mutate(COUNTRY = as.character(COUNTRY),
                Method = as.character(Method),
                Paradigm = as.character(Paradigm)) %>%
  dplyr::transmute(
    COUNTRY, Method, Paradigm,
    TAM, N_True_Alarms, PPV, Sensitivity,
    Mean_Lead_Time, WP, ALY,
    Mean_Lead_Time_conditional, WP_conditional,
    N_False_Alarms,
    Total_Triggers, True_Alarms, False_Alarms,
    n_TrueActionable, n_Reactive,
    N_Years_with_A1_True, N_Years_Evaluable
  )
utils::write.csv(
  country_8metric_csv,
  file.path(OUTPUT_DIR, "Stage5_Country_8Metric_Summary.csv"),
  row.names = FALSE
)

# 8-metric CSV with bootstrap 95% CIs joined per (COUNTRY, Method).
country_8metric_with_ci <- country_8metric_csv %>%
  dplyr::left_join(
    boot_ci %>% dplyr::mutate(COUNTRY = as.character(COUNTRY),
                              Method = as.character(Method)) %>%
      dplyr::select(COUNTRY, Method,
                    TAM_lo, TAM_hi,
                    N_True_Alarms_lo, N_True_Alarms_hi,
                    PPV_lo, PPV_hi,
                    Sens_lo, Sens_hi,
                    MLT_lo, MLT_hi,
                    WP_lo, WP_hi,
                    ALY_lo, ALY_hi,
                    N_False_Alarms_lo, N_False_Alarms_hi),
    by = c("COUNTRY", "Method")
  )
utils::write.csv(
  country_8metric_with_ci,
  file.path(OUTPUT_DIR, "Stage5_Country_8Metric_Summary_with_CIs.csv"),
  row.names = FALSE
)

# Wide country metrics table with bootstrap CIs.
metrics_with_ci <- country_metrics %>%
  dplyr::mutate(COUNTRY = as.character(COUNTRY),
                Method = as.character(Method)) %>%
  dplyr::left_join(
    boot_ci %>%
      dplyr::mutate(COUNTRY = as.character(COUNTRY),
                    Method = as.character(Method)),
    by = c("COUNTRY", "Method")
  )
utils::write.csv(
  metrics_with_ci,
  file.path(OUTPUT_DIR, "Stage5_Country_Framework_Metrics_with_CIs.csv"),
  row.names = FALSE
)

# Dominance probability (country-level).
utils::write.csv(
  dominance_df,
  file.path(OUTPUT_DIR, "Stage5_Country_Dominance_Probabilities.csv"),
  row.names = FALSE
)

# Country Dominance Matrix (long).
utils::write.csv(
  country_dominance_long %>%
    dplyr::transmute(
      COUNTRY = as.character(COUNTRY),
      Method = as.character(Method),
      Paradigm = as.character(Paradigm),
      Metric = Metric,
      Value = Value,
      Dominance_Score = Dominance_Score,
      Is_Dominant_at_threshold_075 = Is_Dominant
    ),
  file.path(OUTPUT_DIR, "Stage5_Country_Dominance_Matrix.csv"),
  row.names = FALSE
)

# Wilcoxon per-metric detector-paired results.
utils::write.csv(
  all_wilcoxon_results,
  file.path(OUTPUT_DIR, "Stage5_Country_Wilcoxon_PerMetric.csv"),
  row.names = FALSE
)

# Dedicated Constant TA versus OT, Farrington, EWARS and EARS.
utils::write.csv(
  constant_ta_comparator_wilcoxon,
  file.path(OUTPUT_DIR, "Stage5_Country_Wilcoxon_ConstantTA_vs_Comparators.csv"),
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
  file.path(OUTPUT_DIR, "Stage5_Method_Summary_Long.csv"),
  row.names = FALSE
)
utils::write.csv(
  method_summary_agg %>%
    dplyr::mutate(Method = as.character(Method)) %>%
    dplyr::select(Method, Mean_Sweep, SD_Sweep, Z_Sweep, Sweep_Class),
  file.path(OUTPUT_DIR, "Stage5_Method_Summary_Aggregate.csv"),
  row.names = FALSE
)

# Bootstrap replicates (long; one row per COUNTRY x .rep x Method).
all_reps_df <- dplyr::bind_rows(
  lapply(names(boot_results), function(ctry) {
    reps <- boot_results[[ctry]]$replicates
    reps$COUNTRY <- ctry
    reps
  })
) %>%
  dplyr::select(COUNTRY, .rep, Method,
                TAM, N_True_Alarms, Sensitivity,
                Mean_Lead_Time, WP, PPV, ALY, N_False_Alarms)
utils::write.csv(
  all_reps_df,
  file.path(OUTPUT_DIR, "Stage5_Country_Bootstrap_Replicates.csv"),
  row.names = FALSE
)
cat("Saved bootstrap replicates: ",
    file.path(OUTPUT_DIR, "Stage5_Country_Bootstrap_Replicates.csv"),
    "\n", sep = "")

# ==============================================================================
# FIGURE 7 — DETECTOR MAP, METRIC TABLE, AND DOT PLOT
# ==============================================================================
# Composes a single composite figure summarising country dominance:
#   a. Choropleth detector map of the 17 administrative countries, encoding the
#      consensus winner under a four-tier classification (strong, partial,
#      lead_only, contested) with Bonferroni-adjusted head-to-head testing.
#   b. Per-country metric table with embedded per-metric significance below
#      each cell value.
#   c. Per-detector dot plot showing the bootstrap dominance probabilities
#      for each country and detector.
#
# Per-country per-metric significance is also computed for six ORDERED
# detector pairs (each ordered direction is its own one-sided test of the
# explicit signed difference X - Y). The full 6 x 5 x 17 table is saved as
# Stage5_Country_Summary_PanelD_HeatmapData.csv for downstream reporting; it
# is NOT rendered as a figure panel.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 19. COUNTRY CANONICALISER
# ------------------------------------------------------------------------------
.normalize_key <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s <- gsub("[().,]", " ", s)
  s <- gsub("\\s+", " ", s)
  trimws(s)
}

CANONICAL_COUNTRIES <- c("BRAZIL", "COLOMBIA", "MEXICO", "PERU",
                         "PHILIPPINES", "SINGAPORE", "SRI LANKA", "TAIWAN")

# Recoder for common name variants and ISO-2 codes. Falls through to the
# normalised key when no variant matches.
country_recoder <- c(
  "BRAZIL"           = "BRAZIL",
  "BRASIL"           = "BRAZIL",
  "BR"               = "BRAZIL",
  "COLOMBIA"         = "COLOMBIA",
  "CO"               = "COLOMBIA",
  "MEXICO"           = "MEXICO",
  "MX"               = "MEXICO",
  "PERU"             = "PERU",
  "PE"               = "PERU",
  "PHILIPPINES"      = "PHILIPPINES",
  "PHILIPPINE"       = "PHILIPPINES",
  "PH"               = "PHILIPPINES",
  "SINGAPORE"        = "SINGAPORE",
  "SG"               = "SINGAPORE",
  "SRI LANKA"        = "SRI LANKA",
  "SRILANKA"         = "SRI LANKA",
  "LK"               = "SRI LANKA",
  "TAIWAN"           = "TAIWAN",
  "TW"               = "TAIWAN",
  "CHINESE TAIPEI"   = "TAIWAN"
)

canonical_country <- function(x) {
  s   <- .normalize_key(x)
  out <- unname(country_recoder[s])
  out[is.na(out)] <- s[is.na(out)]    # keep unrecognised values for inspection
  out
}

# Apply canonicaliser to dominance and metrics tables.
dominance_df_canonical <- dominance_df %>%
  dplyr::mutate(COUNTRY = canonical_country(COUNTRY))
country_metrics_canonical <- country_metrics %>%
  dplyr::mutate(COUNTRY = canonical_country(COUNTRY))
boot_ci_canonical <- boot_ci %>%
  dplyr::mutate(COUNTRY = canonical_country(COUNTRY))

# ------------------------------------------------------------------------------
# 20. WINNER-ROW TABLE AND HEAD-TO-HEAD CONSENSUS TIER CLASSIFICATION
# ------------------------------------------------------------------------------
# For each country, the all-pairs head-to-head test is computed on the 2x2
# bootstrap win-counts table (chi-squared, Fisher exact fallback when any
# expected cell < 5), one-sided p-value based on which detector has the
# higher empirical win count. Three pairs per country:
#   ConstantTA  vs OutbreakThr
#   ContinuousTA vs OutbreakThr
#   ConstantTA  vs ContinuousTA
#
# Dual Bonferroni reporting:
#   Primary (strict)        : factor 3 within-country x N_regions cross-country
#                             = total_corr_strict
#   Sensitivity (within-only): factor 3 within-country only
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
  dplyr::select(COUNTRY, Observed_Winner, Dominance_Probability,
                P_ConstantTA, P_ContinuousTA, P_OutbreakThreshold) %>%
  dplyr::left_join(
    country_metrics_canonical %>%
      dplyr::select(COUNTRY, Method, TAM, N_True_Alarms,
                    Sensitivity, Mean_Lead_Time, WP) %>%
      dplyr::mutate(Method = as.character(Method)),
    by = c("COUNTRY" = "COUNTRY", "Observed_Winner" = "Method")
  )

n_missing_metric <- sum(is.na(winner_rows$TAM))
if (n_missing_metric > 0) {
  warning("Could not match metric values for ", n_missing_metric, " country(s).")
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

N_BOOTS  <- BOOT_N_CI
N_PAIRS  <- 3L
n_countries <- nrow(winner_rows)
total_corr_strict <- N_PAIRS * n_countries
total_corr_within <- N_PAIRS

p_const_vs_ot_raw   <- numeric(n_countries)
p_cont_vs_ot_raw    <- numeric(n_countries)
p_const_vs_cont_raw <- numeric(n_countries)
w_const_vs_ot       <- character(n_countries)
w_cont_vs_ot        <- character(n_countries)
w_const_vs_cont     <- character(n_countries)
for (i in seq_len(n_countries)) {
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
    # `Sig.` column of the per-country metrics table (CSV export): one-sample binomial
    # test of the displayed Dominance_Probability against the chance baseline
    # of 1/3 (three competing detectors).  We Bonferroni-adjust across the
    # n_countries countries so the family-wise error is controlled at
    # SIG_LEVEL across the column.  This makes the Sig. star a direct test
    # of the value in column 3 — high Dominance_Probability earns a star,
    # low/near-chance Dominance_Probability does not — rather than the
    # all-pairs consensus test, which can show "ns" even at high Pr if one
    # head-to-head pair fails the strict cross-country Bonferroni.
    n_wins_binom = as.integer(round(Dominance_Probability * N_BOOTS)),
    p_dominance_raw = ifelse(
      is.na(n_wins_binom) | N_BOOTS <= 0L,
      NA_real_,
      stats::pbinom(pmax(n_wins_binom - 1L, 0L),
                    size = N_BOOTS, prob = 1 / 3,
                    lower.tail = FALSE)
    ),
    p_dominance_bonf = pmin(p_dominance_raw * n_countries, 1),
    sig_star_dominance = dplyr::case_when(
      is.na(p_dominance_bonf)  ~ "ns",
      p_dominance_bonf < 0.001 ~ "***",
      p_dominance_bonf < 0.01  ~ "**",
      p_dominance_bonf < 0.05  ~ "*",
      TRUE                     ~ "ns"
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
  n_strong, n_partial, n_lead_only, n_contested, n_countries))
cat(sprintf(
  "Within-country-only Bonferroni sensitivity: %d countries meet 'strong' (vs %d under strict).\n",
  n_strong_within, n_strong))

# ------------------------------------------------------------------------------
# 20b. CONSTANT TA FOUR-COMPARATOR CONSENSUS OUTPUTS
# ------------------------------------------------------------------------------
# Supplementary-table-specific extension of the regional four-comparator logic.
# Existing Figure 7 all-pairs consensus objects are retained unchanged.
#
# Comparisons:
#   Constant TA vs Outbreak Threshold
#   Constant TA vs Farrington
#   Constant TA vs EWARS
#   Constant TA vs EARS
#
# Strict cross-country Bonferroni:
#   4 comparisons x number of included countries.
# Within-country Bonferroni:
#   4 comparisons.
# ------------------------------------------------------------------------------

ST_COUNTRY_COMPARATORS <- c(
  "Outbreak Threshold" = "P_OutbreakThreshold",
  "Farrington"         = "P_Farrington",
  "EWARS"              = "P_EWARS",
  "EARS"               = "P_EARS"
)

ST_COUNTRY_K_WITHIN <- length(ST_COUNTRY_COMPARATORS)
ST_COUNTRY_K_STRICT <- ST_COUNTRY_K_WITHIN *
  length(unique(as.character(dominance_df_canonical$COUNTRY)))

.st_country_pair_test <- function(p_const, p_comp, B = N_BOOTS) {
  if (!is.finite(p_const) || !is.finite(p_comp)) {
    return(list(p_raw = NA_real_, winner = NA_character_))
  }

  n_const <- round(p_const * B)
  n_comp  <- round(p_comp  * B)
  n_pair  <- n_const + n_comp

  if (n_pair <= 0L || n_const == n_comp) {
    return(list(
      p_raw = if (n_pair <= 0L) NA_real_ else 1.0,
      winner = NA_character_
    ))
  }

  if (n_const > n_comp) {
    p_raw <- stats::binom.test(
      n_const, n_pair, p = 0.5, alternative = "greater"
    )$p.value
    winner <- "Constant Transmission Acceleration"
  } else {
    p_raw <- stats::binom.test(
      n_comp, n_pair, p = 0.5, alternative = "greater"
    )$p.value
    winner <- "Comparator"
  }

  list(p_raw = p_raw, winner = winner)
}

st_country_rows <- list()

for (i in seq_len(nrow(dominance_df_canonical))) {
  rr <- dominance_df_canonical[i, , drop = FALSE]
  ctry <- as.character(rr$COUNTRY)
  p_const <- as.numeric(rr$P_ConstantTA)

  one_country <- lapply(names(ST_COUNTRY_COMPARATORS), function(comp) {
    p_comp_col <- unname(ST_COUNTRY_COMPARATORS[[comp]])
    p_comp <- as.numeric(rr[[p_comp_col]])

    tt <- .st_country_pair_test(p_const, p_comp, B = N_BOOTS)

    winner <- if (is.na(tt$winner)) {
      NA_character_
    } else if (tt$winner == "Comparator") {
      comp
    } else {
      tt$winner
    }

    data.frame(
      COUNTRY = ctry,
      Comparator = comp,
      P_ConstantTA = p_const,
      P_Comparator = p_comp,
      Pair_Winner = winner,
      p_raw = tt$p_raw,
      p_bonf_strict = ifelse(
        is.na(tt$p_raw), NA_real_,
        min(1, tt$p_raw * ST_COUNTRY_K_STRICT)
      ),
      p_bonf_within = ifelse(
        is.na(tt$p_raw), NA_real_,
        min(1, tt$p_raw * ST_COUNTRY_K_WITHIN)
      ),
      stringsAsFactors = FALSE
    )
  })

  st_country_rows[[ctry]] <- dplyr::bind_rows(one_country)
}

st_country_long <- dplyr::bind_rows(st_country_rows) %>%
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

st_country_wide <- st_country_long %>%
  dplyr::group_by(COUNTRY) %>%
  dplyr::summarise(
    Pr_ConstantTA = dplyr::first(P_ConstantTA),

    Strict_ConstantTA_vs_OT =
      Strict_Result[Comparator == "Outbreak Threshold"][1],
    Strict_ConstantTA_vs_Farrington =
      Strict_Result[Comparator == "Farrington"][1],
    Strict_ConstantTA_vs_EWARS =
      Strict_Result[Comparator == "EWARS"][1],
    Strict_ConstantTA_vs_EARS =
      Strict_Result[Comparator == "EARS"][1],

    Within_ConstantTA_vs_OT =
      Within_Result[Comparator == "Outbreak Threshold"][1],
    Within_ConstantTA_vs_Farrington =
      Within_Result[Comparator == "Farrington"][1],
    Within_ConstantTA_vs_EWARS =
      Within_Result[Comparator == "EWARS"][1],
    Within_ConstantTA_vs_EARS =
      Within_Result[Comparator == "EARS"][1],

    n_sig_wins_strict = sum(
      Sig_Strict &
        Pair_Winner == "Constant Transmission Acceleration",
      na.rm = TRUE
    ),
    n_sig_losses_strict = sum(
      Sig_Strict &
        !is.na(Pair_Winner) &
        Pair_Winner != "Constant Transmission Acceleration",
      na.rm = TRUE
    ),

    n_sig_wins_within = sum(
      Sig_Within &
        Pair_Winner == "Constant Transmission Acceleration",
      na.rm = TRUE
    ),
    n_sig_losses_within = sum(
      Sig_Within &
        !is.na(Pair_Winner) &
        Pair_Winner != "Constant Transmission Acceleration",
      na.rm = TRUE
    ),

    Weakest_link_P_strict = ifelse(
      all(is.na(
        p_bonf_strict[
          Pair_Winner == "Constant Transmission Acceleration"
        ]
      )),
      NA_real_,
      max(
        p_bonf_strict[
          Pair_Winner == "Constant Transmission Acceleration"
        ],
        na.rm = TRUE
      )
    ),

    Consensus_winner_strict = dplyr::case_when(
      sum(
        Sig_Strict &
          Pair_Winner == "Constant Transmission Acceleration",
        na.rm = TRUE
      ) == ST_COUNTRY_K_WITHIN ~
        "Constant Transmission Acceleration",

      sum(
        Sig_Strict &
          !is.na(Pair_Winner) &
          Pair_Winner != "Constant Transmission Acceleration",
        na.rm = TRUE
      ) > 0 ~ "Contested",

      TRUE ~ "No strict consensus"
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(COUNTRY)

utils::write.csv(
  st_country_long,
  file.path(
    OUTPUT_DIR,
    "Stage5_Country_ConstantTA_4Comparator_Consensus_Long.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  st_country_wide,
  file.path(
    OUTPUT_DIR,
    "Stage5_Country_ConstantTA_4Comparator_Consensus.csv"
  ),
  row.names = FALSE
)

cat("\nSaved country Constant TA four-comparator consensus outputs:\n")
cat("  Stage5_Country_ConstantTA_4Comparator_Consensus_Long.csv\n")
cat("  Stage5_Country_ConstantTA_4Comparator_Consensus.csv\n")
cat(
  "  Strict Bonferroni k = ", ST_COUNTRY_K_STRICT,
  "; within-country Bonferroni k = ", ST_COUNTRY_K_WITHIN,
  "\n", sep = ""
)

# ------------------------------------------------------------------------------
# 21. PER-COUNTRY PER-METRIC SIGNIFICANCE (six ordered pairs)
# ------------------------------------------------------------------------------
# Each ordered pair X vs Y is its own one-sided test of H1: X > Y on the
# explicit signed difference X - Y. The two ordered directions of an
# unordered comparison give DIFFERENT p-values that sum to ~1.
#
# Test method: paired bootstrap on year-cluster replicates (rigorous).
# Bonferroni: factor 5 within (COUNTRY x ordered_pair). Each ordered pair is
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

# Build long-format replicates frame (with canonical country names) for the
# paired-bootstrap test.
all_reps_canonical <- dplyr::bind_rows(
  lapply(names(boot_results), function(ctry) {
    reps <- boot_results[[ctry]]$replicates
    reps$COUNTRY <- canonical_country(ctry)
    reps
  })
)

compute_level1_paired_bootstrap_ordered <- function(reps_df, countries,
                                                    ordered_pairs, metrics) {
  out <- list()
  for (ctry in countries) {
    sub <- reps_df %>% dplyr::filter(COUNTRY == ctry)
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
            COUNTRY = ctry, Ordered_Pair = op$label, Pair_Code = op$code,
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
          COUNTRY = ctry, Ordered_Pair = op$label, Pair_Code = op$code,
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

countries_in_order <- unique(winner_rows$COUNTRY)

cat("\nPer-country per-metric inference: paired bootstrap, 6 ordered pairs.\n")
level1_results_ordered <- compute_level1_paired_bootstrap_ordered(
  all_reps_canonical, countries_in_order, ORDERED_PAIRS, HH_METRICS_COUNTRY
)
level1_method <- "bootstrap_paired"

# Bonferroni adjustment within (COUNTRY x Pair_Code).
if (nrow(level1_results_ordered) > 0L) {
  level1_results_ordered <- level1_results_ordered %>%
    dplyr::group_by(COUNTRY, Pair_Code) %>%
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

cat("=== Per-country per-metric significance (6 ordered pairs x 5 metrics x ",
    n_countries, " countries = ", nrow(level1_results_ordered), " tests) ===\n",
    sep = "")
cat("Method:                     ", level1_method, "\n", sep = "")
cat("Test:                       one-sided H1: X > Y per ordered pair\n")
cat("Bonferroni factor:          ", BONF_LEVEL1_FACTOR,
    " (within COUNTRY x ordered pair)\n", sep = "")
cat("Significant at p<0.05:      ",
    sum(level1_results_ordered$sig_star %in% c("*", "**", "***"), na.rm = TRUE),
    " of ", sum(!is.na(level1_results_ordered$p_bonf)), "\n", sep = "")

# ------------------------------------------------------------------------------
# 22. CROSS-COUNTRY PER-METRIC WILCOXON (formatted from in-memory results)
# ------------------------------------------------------------------------------
level2_results <- all_wilcoxon_results %>%
  dplyr::transmute(
    Metric, Comparison,
    Detector_A, Detector_B,
    n_pairs     = N_Countries_Paired,
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
build_per_metric_sig_for_country <- function(ctry, consensus_winner) {
  sub <- level1_results_ordered %>% dplyr::filter(COUNTRY == ctry)
  if (nrow(sub) == 0L) {
    return(data.frame(COUNTRY = ctry, Metric = HH_METRICS_COUNTRY,
                      sig_star = "na", p_bonf = NA_real_,
                      reference_pair = NA_character_,
                      stringsAsFactors = FALSE))
  }
  rows <- list()
  for (m in HH_METRICS_COUNTRY) {
    sub_m <- sub %>% dplyr::filter(Metric == m)
    if (nrow(sub_m) == 0L) {
      rows[[length(rows) + 1L]] <- data.frame(
        COUNTRY = ctry, Metric = m, sig_star = "na",
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
        COUNTRY = ctry, Metric = m, sig_star = "na",
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
      COUNTRY = ctry, Metric = m,
      sig_star = pick$sig_star, p_bonf = pick$p_bonf,
      reference_pair = pick$Ordered_Pair, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

per_country_per_metric_sig <- do.call(rbind, lapply(seq_len(nrow(winner_rows)),
  function(i) {
    build_per_metric_sig_for_country(
      winner_rows$COUNTRY[i], winner_rows$consensus_winner[i]
    )
  }))

sig_wide <- per_country_per_metric_sig %>%
  dplyr::select(COUNTRY, Metric, sig_star) %>%
  tidyr::pivot_wider(names_from = Metric, values_from = sig_star,
                     names_prefix = "sig_") %>%
  as.data.frame()

winner_rows <- winner_rows %>% dplyr::left_join(sig_wide, by = "COUNTRY")

# ------------------------------------------------------------------------------
# 23. FIGURE 7 PALETTES
# ------------------------------------------------------------------------------
detector_palette <- c(
  "Constant TA"        = "#1F77B4",
  "Continuous TA"      = "#2CA02C",
  "Outbreak Threshold" = "#D62728",
  "No consensus"       = "#6B7280"
)

# Confirm the winner_rows table has every canonical country represented.
country_set_present  <- intersect(CANONICAL_COUNTRIES, unique(winner_rows$COUNTRY))
country_set_missing  <- setdiff(CANONICAL_COUNTRIES, country_set_present)
if (length(country_set_missing) > 0) {
  warning("Countries missing from winner_rows: ",
          paste(country_set_missing, collapse = ", "),
          " (these will not appear in Figure 7).")
}
n_countries <- nrow(winner_rows)
stopifnot(n_countries >= 1)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 24. FIGURE 7 PANEL A — METRIC TABLE WITH EMBEDDED PER-METRIC SIGNIFICANCE
# ------------------------------------------------------------------------------
# Nine-column table:
#   1 Country | 2 Detector (tier-aware fill) | 3 Dominance probability |
#   4 Sig. (binomial test of col 3 vs chance 1/3, Bonferroni across countries) |
#   5 TAM | 6 N true alarms | 7 Sensitivity | 8 Mean lead time | 9 WP
# Each metric cell renders the value on its own line, then a parenthetical
# significance label "(*/**/***/ns/-)" below in a smaller subscript line.
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

COLS_TABLE <- c("Country", "Detector",
                "Dominance\nprobability",
                "Sig.",
                "TAM",
                "N true\nalarms",
                "Sensitivity",
                "Mean lead\ntime (wk)",
                "WP\n(wk)")
N_COLS <- length(COLS_TABLE)

table_long <- winner_rows %>%
  dplyr::filter(COUNTRY %in% CANONICAL_COUNTRIES, !is.na(Observed_Winner)) %>%
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
  dplyr::arrange(COUNTRY)

# Ensure all per-metric significance columns exist (they may be missing if
# Level 1 inference was unavailable for a particular country).
for (m in HH_METRICS_COUNTRY) {
  col_name <- paste0("sig_", m)
  if (!col_name %in% names(table_long)) table_long[[col_name]] <- "na"
  table_long[[col_name]][is.na(table_long[[col_name]])] <- "na"
}

N_ROWS_TABLE <- nrow(table_long)
HEADER_Y_TBL <- N_ROWS_TABLE + 1L
table_long <- table_long %>%
  dplyr::mutate(row_y = rev(seq_len(dplyr::n())))

# Helper: format metric value with parenthetical sig label below.
.fmt_metric_with_sig <- function(value_str, sig_label) {
  sig_paren <- ifelse(is.na(sig_label) | sig_label == "na" | sig_label == "ns",
                      "(ns)",
                      paste0("(", sig_label, ")"))
  ifelse(value_str == "-", value_str,
         paste0(value_str, "\n", sig_paren))
}

cell_df <- table_long %>%
  dplyr::transmute(
    COUNTRY, row_y, Detector_Label_short, Detector_Label_two,
    Dominance_Probability, sig_star = sig_star_dominance,
    sig_TAM, sig_N_True_Alarms, sig_Sensitivity,
    sig_Mean_Lead_Time, sig_WP,
    consensus_pass, consensus_tier,
    `1` = COUNTRY,
    `2` = Detector_Label_two,
    `3` = fmt_2(Dominance_Probability),
    `4` = sig_star,
    `5` = .fmt_metric_with_sig(fmt_int(TAM),            sig_TAM),
    `6` = .fmt_metric_with_sig(fmt_1(N_True_Alarms),    sig_N_True_Alarms),
    `7` = .fmt_metric_with_sig(fmt_pct(Sensitivity),    sig_Sensitivity),
    `8` = .fmt_metric_with_sig(fmt_1(Mean_Lead_Time),   sig_Mean_Lead_Time),
    `9` = .fmt_metric_with_sig(fmt_1(WP),               sig_WP)
  ) %>%
  tidyr::pivot_longer(cols = `1`:`9`,
                      names_to = "col_x", values_to = "cell_text") %>%
  dplyr::mutate(col_x = as.integer(col_x))

.metric_cols <- c(3L, 5L, 6L, 7L, 8L, 9L)

# Helpers for cell colours.
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

bg_df <- cell_df %>%
  dplyr::mutate(
    is_strong_tier = consensus_tier == "strong",
    fill = dplyr::case_when(
      col_x == 2 ~ .detector_cell_fill(Detector_Label_short, consensus_tier),
      row_y %% 2 == 0 ~ "#F2F2F2",
      TRUE ~ "#FFFFFF"
    ),
    fill = ifelse(is.na(fill), "grey90", fill),
    text_colour = dplyr::case_when(
      col_x == 2 ~ "white",
      col_x == 4 & sig_star == "ns"             ~ "grey55",
      col_x == 4                                ~ "#0B2447",
      !is_strong_tier & col_x %in% .metric_cols ~ "grey55",
      TRUE ~ "grey15"
    ),
    fontface = dplyr::case_when(
      col_x %in% c(1L, 2L, 4L)                  ~ "bold",
      !is_strong_tier & col_x %in% .metric_cols ~ "italic",
      TRUE                                      ~ "plain"
    )
  )

bg_header <- data.frame(
  col_x = seq_len(N_COLS), row_y = HEADER_Y_TBL,
  fill = "#0B2447", text_colour = "white",
  fontface = "bold", label = COLS_TABLE, stringsAsFactors = FALSE
)

y_min_tbl <- 0.5
y_max_tbl <- HEADER_Y_TBL + 0.5

# NOTE: p_table (the per-country metrics table panel) is not part of
# Figure 7. The object is retained because the underlying per-country
# metrics are still written to CSV, and keeping the builder makes it easy to
# reinstate the panel if required. Not exported.
p_table <- ggplot2::ggplot() +
  ggplot2::geom_tile(
    data = bg_header,
    ggplot2::aes(x = col_x, y = row_y, fill = I(fill)),
    colour = "white", width = 0.98, height = 0.98
  ) +
  ggplot2::geom_tile(
    data = bg_df,
    ggplot2::aes(x = col_x, y = row_y, fill = I(fill)),
    colour = "white", width = 0.98, height = 0.98
  ) +
  ggplot2::geom_text(
    data = bg_header,
    ggplot2::aes(x = col_x, y = row_y, label = label),
    colour = "white", fontface = "bold",
    size = 3.0, lineheight = 0.85
  ) +
  ggplot2::geom_text(
    data = dplyr::left_join(
      cell_df,
      bg_df %>% dplyr::select(col_x, row_y, text_colour, fontface),
      by = c("col_x", "row_y")
    ),
    ggplot2::aes(x = col_x, y = row_y, label = cell_text,
                 colour = I(text_colour), fontface = I(fontface)),
    size = 2.7, lineheight = 0.85
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
    title = "Bootstrap-winning detector and its per-metric significance, by country"
  ) +
  ggplot2::theme_void(base_size = PUB_BASE, base_family = PUB_FAMILY) +
  ggplot2::theme(
    plot.title  = ggplot2::element_text(face = "bold", size = 13, hjust = 0,
                                        margin = ggplot2::margin(b = 6)),
    plot.margin = ggplot2::margin(10, 30, 12, 18)
  )

# ------------------------------------------------------------------------------
# 25. FIGURE 7 PANEL B — PER-DETECTOR DOT PLOT (with per-country tier legend)
# ------------------------------------------------------------------------------
# One jittered dot per country per detector (n_countries x 3 dots). Median
# crossbar per detector. Reference lines at 1/3 (chance among three detectors)
# and DOMINANCE_THRESHOLD. The right-side x-extension hosts a per-country
# tier-legend overlay box.
# ------------------------------------------------------------------------------
# Maps a full detector name (as stored in Observed_Winner) to its short label.
WINNER_TO_SHORT <- c(
  "Constant Transmission Acceleration"   = "Constant TA",
  "Continuous Transmission Acceleration" = "Continuous TA",
  "Outbreak Threshold"                   = "Outbreak Threshold",
  "Farrington"                           = "Farrington",
  "EARS"                                 = "EARS",
  "EWARS"                                = "EWARS"
)

# Maps each P_* column to the label shown on the axis.
DETECTOR_SHORT_LABEL <- c(
  P_ConstantTA        = "Constant TA",
  P_ContinuousTA      = "Continuous TA",
  P_OutbreakThreshold = "Outbreak Threshold",
  P_Farrington        = "Farrington",
  P_EARS              = "EARS",
  P_EWARS             = "EWARS"
)

# All six target detectors, driven by TARGET_P_COLS rather than a hard-coded
# set of three, so this figure cannot fall out of step with TARGET_DETECTORS.
panel_b_long <- dominance_df_canonical %>%
  dplyr::select(COUNTRY, dplyr::all_of(unname(TARGET_P_COLS))) %>%
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
      dplyr::select(COUNTRY, consensus_tier, Observed_Winner, Leader_Label_short),
    by = "COUNTRY"
  ) %>%
  dplyr::mutate(
    # Generalised to any of the six target detectors; the previous case_when
    # listed only three and would mark no leader for Farrington/EARS/EWARS.
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
  dplyr::arrange(COUNTRY, Detector_Label)

# Auxiliary cross-country paired Wilcoxon (saved alongside the dot data).
.run_paired_wilcox <- function(A_lab, B_lab) {
  a <- panel_b_long %>%
    dplyr::filter(Detector_Label == A_lab) %>%
    dplyr::arrange(COUNTRY) %>% dplyr::pull(Dominance_P)
  b <- panel_b_long %>%
    dplyr::filter(Detector_Label == B_lab) %>%
    dplyr::arrange(COUNTRY) %>% dplyr::pull(Dominance_P)
  diffs <- a - b
  test  <- suppressWarnings(stats::wilcox.test(a, b, paired = TRUE,
                                               exact = FALSE))
  data.frame(
    Comparison         = paste(A_lab, "vs", B_lab),
    N_Countries_Paired = sum(!is.na(diffs)),
    Median_Diff        = stats::median(diffs, na.rm = TRUE),
    V_statistic        = unname(test$statistic),
    p_raw              = test$p.value,
    stringsAsFactors   = FALSE
  )
}
panel_b_pairs_for_csv <- list(
  c("Constant TA",   "Outbreak Threshold"),
  c("Continuous TA", "Outbreak Threshold"),
  c("Constant TA",   "Continuous TA")
)
panel_b_wilcox_csv <- do.call(rbind,
                              lapply(panel_b_pairs_for_csv,
                                     function(p) .run_paired_wilcox(p[1], p[2])))
panel_b_wilcox_csv <- panel_b_wilcox_csv %>%
  dplyr::mutate(p_bonferroni = pmin(p_raw * dplyr::n(), 1))

# Persist the dot data + auxiliary Wilcoxon table.
panel_b_dots <- panel_b_long %>%
  dplyr::transmute(COUNTRY, Detector = as.character(Detector_Label),
                   Dominance_P, is_leader_dot, consensus_tier,
                   Observed_Winner)

# Diagnostic banner.
cat("\n=== Per-detector dot plot data (",
    nrow(panel_b_long), " dots = ", n_countries,
    " countries x ", length(TARGET_DETECTORS), " detectors) ===\n", sep = "")
cat("Median Pr per detector:\n")
print(panel_b_long %>%
        dplyr::group_by(Detector_Label) %>%
        dplyr::summarise(median_Pr = stats::median(Dominance_P, na.rm = TRUE),
                         mean_Pr   = mean(Dominance_P, na.rm = TRUE),
                         n_countries = dplyr::n(), .groups = "drop"),
      row.names = FALSE, digits = 3)

pj_b     <- ggplot2::position_jitter(width = 0.18, height = 0, seed = 12345)
chance_p <- 1 / 3

# Significance brackets for the dominance-probability panel: Constant TA versus
# each comparator, paired by country and tested with the same paired
# Wilcoxon signed-rank used everywhere else. ONLY SIGNIFICANT comparisons get a
# bracket and a p-value; non-significant pairs are simply not drawn. Every
# comparison remains in the exported Wilcoxon CSV.
# Full pairwise test table FIRST: every comparison, significant or not, with
# V, p, CI and a Drawn_On_Figure flag. This is the file that supports the
# significance brackets on the figure -- exported below so every p-value shown
# on the panel can be traced to a row.
dom_tests <- dominance_pairwise_tests(
  d         = as.data.frame(panel_b_long),
  unit_col  = "COUNTRY",
  level_col = "Detector_Label",
  value_col = "Dominance_P",
  focal     = "Constant TA")

if (!is.null(dom_tests)) {
  utils::write.csv(dom_tests, file.path(OUTPUT_DIR, "Figure7_DominanceProbability_Wilcoxon_Results.csv"),
                   row.names = FALSE)
  cat("  Saved: ", file.path(OUTPUT_DIR, "Figure7_DominanceProbability_Wilcoxon_Results.csv"), "\n", sep = "")
}

# Geometry is built FROM that same table, so figure and CSV cannot disagree.
dom_brackets <- build_sig_brackets(
  d         = as.data.frame(panel_b_long),
  unit_col  = "COUNTRY",
  level_col = "Detector_Label",
  value_col = "Dominance_P",
  focal     = "Constant TA",
  tests     = dom_tests)
dom_vrange <- diff(range(panel_b_long$Dominance_P, na.rm = TRUE))
if (!is.finite(dom_vrange) || dom_vrange <= 0) dom_vrange <- 1
# build_sig_brackets() returns NULL when nothing is significant; nrow(NULL) is
# NULL, not 0, so the count is normalised here rather than inline.
n_dom_brackets <- if (is.null(dom_brackets)) 0L else nrow(dom_brackets)

p_dot <- ggplot2::ggplot(panel_b_long,
                         ggplot2::aes(x = Detector_Label, y = Dominance_P)) +
  # Decisive / chance reference lines, their labels and the shaded decisive
  # band removed. DOMINANCE_THRESHOLD still drives the tier
  # classification and the "Dominance (score >= x)" column; it is simply no
  # longer drawn on this panel.
  ggplot2::geom_point(
    ggplot2::aes(fill = Detector_Label,
                 stroke = border_stroke,
                 colour = border_colour),
    shape = 21, position = pj_b,
    size = 3.0, alpha = 0.92, na.rm = TRUE, show.legend = FALSE
  ) +
  ggplot2::stat_summary(
    fun = stats::median, geom = "crossbar",
    width = 0.55, linewidth = 0.45, colour = "grey15",
    fatten = 0, fill = NA, na.rm = TRUE
  ) +
  # Per-country tier-legend overlay REMOVED. It was an in-plot
  # annotation box (white rect + four example dots + explanatory text) pinned to
  # x = 3.55-4.70, which was positioned for three detectors and would in any
  # case now sit on top of the Farrington/EARS/EWARS columns. The tier encoding
  # it described is documented in the figure caption and the exported CSVs.

  # All six target detectors, taken straight from detector_palette so the fill
  # scale cannot fall out of step with the detectors actually plotted.
  ggplot2::scale_fill_manual(
    values = detector_palette[unname(DETECTOR_SHORT_LABEL[unname(TARGET_P_COLS)])]
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.05 + 0.16 * n_dom_brackets),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.0"),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_x_discrete(
    # All six detectors get a two-line label; without this,
    # Farrington/EARS/EWARS fall back to their raw level names.
    labels = c("Constant TA"        = "Constant\nTA",
               "Continuous TA"      = "Continuous\nTA",
               "Outbreak Threshold" = "Outbreak\nThreshold",
               "Farrington"         = "Farrington",
               "EARS"               = "EARS",
               "EWARS"              = "EWARS"),
    # The old right-hand expansion of 1.8 reserved a whole category's width for
    # the in-plot tier legend, which has been removed. That empty gap was what
    # stretched the axis. Symmetric 0.6 now, so the six categories fill the
    # panel instead of being pushed left.
    expand = ggplot2::expansion(add = c(0.6, 0.6))
  ) +
  ggplot2::labs(
    title = "Per-detector country dominance probabilities",
    x     = NULL,
    y     = "Dominance probability"
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
# 26. PER-COUNTRY PER-METRIC HEATMAP DATA (saved as CSV; not a figure panel)
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
                    COUNTRY %in% CANONICAL_COUNTRIES)
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
      dplyr::select(COUNTRY, Metric,
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
# 27. COMPOSE AND SAVE FIGURE 7
# ------------------------------------------------------------------------------
# Figure 7 is now a SINGLE standalone figure.
#
# The country metrics table (former panel a) has been removed, leaving only the
# bootstrap dominance-probability plot. Because one panel remains, it carries NO
# panel letter and no subtitle: panel lettering exists to disambiguate multiple
# panels, so a lone figure should not be labelled "b" (or anything else).
# A combined two-panel composite is not produced.
# 1.5-column width (120 mm). Six discrete categories spread across the full
# 180 mm double column left the points thinly scattered; this fills the panel.
FIG7_W <- NC_W_MEDIUM + 0.5

# Attach the significance brackets (none are drawn if nothing is significant).
p_dot <- add_sig_brackets(p_dot, dom_brackets, dom_vrange)

fig7 <- p_dot +
  ggplot2::labs(title = "Bootstrap dominance probability by detector") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = PUB_TITLE + 1,
                                       hjust = 0.5,   # centred: standalone figure
                                       family = base_family_global,
                                       margin = ggplot2::margin(b = 8)),
    # Symmetrical left/right margins so the single panel sits centred on the
    # canvas rather than offset as it was inside the two-panel stack.
    plot.margin = ggplot2::margin(10, 14, 9, 14),
    legend.position = "bottom",
    legend.justification = "center"
  )

# Named to mirror Figure5_panel_b_DominanceProbability at the regional scale, so
# the two scales' dominance-probability figures sit together in a file listing.
#
# DIAGNOSTIC GUARD. This figure has been reported as "not appearing". Every
# object it needs is verified here and the export is isolated in tryCatch, so a
# failure names its cause instead of the file silently not existing.
cat("\n[fig7] building country dominance-probability figure\n")
cat("  panel_b_long rows        : ", nrow(panel_b_long), "\n", sep = "")
cat("  detectors plotted        : ",
    paste(levels(panel_b_long$Detector_Label), collapse = ", "), "\n", sep = "")
cat("  significant brackets     : ", n_dom_brackets, "\n", sep = "")
cat("  output directory exists  : ", dir.exists(OUTPUT_DIR), "\n", sep = "")

if (nrow(panel_b_long) == 0L) {
  warning("[fig7] panel_b_long has zero rows - the dominance-probability ",
          "figure cannot be drawn. Check that dominance_df carries the ",
          paste(unname(TARGET_P_COLS), collapse = ", "), " columns.",
          call. = FALSE)
} else {
  .ok7 <- tryCatch({
    save_plot_file("Figure7_CountryDominanceProbability", fig7, FIG7_W, 5.00)
    TRUE
  }, error = function(e) {
    message("[fig7] EXPORT FAILED: ", conditionMessage(e))
    FALSE
  })
  if (isTRUE(.ok7)) {
    .f7 <- file.path(OUTPUT_DIR, "Figure7_CountryDominanceProbability.pdf")
    cat("  written                  : ", file.exists(.f7), "  ", .f7, "\n", sep = "")
  }
}

# ------------------------------------------------------------------------------
# 28. WRITE FIGURE 7 CSV TABLES
# ------------------------------------------------------------------------------
readr::write_csv(winner_rows,
                 file.path(OUTPUT_DIR, "Stage5_Country_Summary_CountryTable.csv"))
readr::write_csv(level1_results_ordered,
                 file.path(OUTPUT_DIR, "Stage5_Country_Summary_PerMetricSignificance.csv"))
readr::write_csv(level2_results,
                 file.path(OUTPUT_DIR, "Stage5_Country_Summary_CrossCountry_Wilcoxon.csv"))
readr::write_csv(panel_b_dots,
                 file.path(OUTPUT_DIR, "Stage5_Country_Summary_PanelB_DotPlot.csv"))
readr::write_csv(panel_b_wilcox_csv,
                 file.path(OUTPUT_DIR, "Stage5_Country_Summary_PanelB_AuxWilcoxon.csv"))
readr::write_csv(panel_d_data,
                 file.path(OUTPUT_DIR, "Stage5_Country_Summary_PanelD_HeatmapData.csv"))

# ==============================================================================
# 29. LEGENDS
# ==============================================================================
# Two text files written to OUTPUT_DIR:
#   Stage5_Figure6_Legend.txt   (one block per metric + method summary)
#   Stage5_Figure7_Legend.txt    (single block for the country summary figure)
# ------------------------------------------------------------------------------

build_figure6_metric_legend <- function(metric_id) {
  metric_label <- hh_metric_full_name[[metric_id]]
  res <- all_wilcoxon_results %>%
    dplyr::filter(Metric == metric_id) %>%
    dplyr::arrange(Comparison)
  n_pairs <- max(res$N_Countries_Paired, na.rm = TRUE)
  sig_pairs <- res %>%
    dplyr::filter(!is.na(Significant_005), Significant_005) %>%
    dplyr::pull(Comparison) %>% as.character()
  outcome <- if (length(sig_pairs) == 0L) {
    paste0("None of the three detector-pair comparisons reached the per-pairwise ",
           "significance threshold (p < ", sprintf("%.2f", HH_ALPHA),
           ") on the available evidence base of ", n_pairs,
           " paired country observations.")
  } else if (length(sig_pairs) == 3L) {
    "All three detector-pair comparisons reached the per-pairwise significance threshold."
  } else {
    paste0(length(sig_pairs), " of 3 detector-pair comparisons reached the ",
           "per-pairwise significance threshold: ",
           paste(sig_pairs, collapse = "; "), ".")
  }
  paste0(
    "Figure 6 (", metric_label, ") | Country generalisability of ",
    metric_label, " across the ", length(country_ppv_order),
    " evaluable dengue-endemic countries, 2016-2024 (excluding 2020, 2021, ",
    "2025; n = ", length(EVALUABLE_YEARS), " evaluable years per country; ",
    "country inclusion floor = ", MIN_PEAK_CASES_PER_YEAR,
    " annual peak cases and >= ", MIN_EVALUABLE_YEARS_PER_COUNTRY,
    " evaluable years). This figure is one of the five per-metric panels ",
    "in the Figure 6 set (TAM, Number of True Alarms, Sensitivity, Mean Lead ",
    "Time, Warning Persistence). The early-warning timeliness metrics (Mean ",
    "Lead Time, Warning Persistence) use the same-denominator zero-coerced ",
    "scheme: years where the metric is not computable due to no qualifying ",
    "triggers contribute zero rather than being excluded. Sensitivity is A1-",
    "restricted (proportion of evaluable seasons with at least one True Alarm ",
    "in the Actionable Window).\n",
    "(a) Country dominance matrix (method-centric layout). Rows = 11 ",
    "outbreak-detection methods, ordered top-to-bottom by sweep count ",
    "(descending; ties broken by descending mean score). Columns = ",
    length(country_ppv_order), " evaluable countries, ordered left-to-right by ",
    "country mean dominance score across the 11 methods (descending). Cell ",
    "shading encodes the normalized dominance score on this metric, computed ",
    "by min-max normalization within each country across detectors: 1.00 = ",
    "within-country leader on ", metric_label, "; 0.00 = within-country laggard. ",
    "Single-color blue ramp from white (0) through mid-blue (0.50) to dark ",
    "navy (1.00); a horizontal blue-intensity legend strip below the matrix ",
    "maps shade to score with reference labels at 0.00 (weak), 0.50 ",
    "(moderate), ", sprintf("%.2f", DOMINANCE_THRESHOLD),
    " (dominant), and 1.00 (sweep). The right-most annotation reports the ",
    "per-METHOD count of countries swept by that method on this metric (countries ",
    "where the method achieves dominance score >= ",
    sprintf("%.2f", DOMINANCE_THRESHOLD), "), out of ",
    length(country_ppv_order), " total. ",
    "(b) Per-detector country values on ", metric_label,
    ". One dot per country per detector, jittered horizontally; mean across ",
    "countries shown as a horizontal crossbar; vertical lines indicate year-",
    "cluster bootstrap 95% CIs (B = ", BOOT_N_CI, " replicates) where ",
    "available. The pale-blue band marks the operationally favourable zone ",
    "where defined. ",
    "(c) Wilcoxon detector-paired significance bars. Three pairwise detector ",
    "comparisons: Constant TA vs Outbreak Threshold; Continuous TA vs ",
    "Outbreak Threshold; Constant TA vs Continuous TA. Wilcoxon signed-rank, ",
    "paired by country (n = ", n_pairs,
    "), two-sided, evaluated PER PAIRWISE COMPARISON at alpha = ",
    sprintf("%.2f", HH_ALPHA),
    " (no across-pair multiplicity correction). Bars show |median ",
    "difference|; navy bars (with adjacent '*') indicate p < ",
    sprintf("%.2f", HH_ALPHA), "; pale grey otherwise. ",
    "Outcome on this metric: ", outcome, "\n"
  )
}

build_figure6_summary_legend <- function() {
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
  pivot_score_lg <- country_dominance_long %>%
    dplyr::filter(as.character(Method) %in% TARGETS_LOCAL,
                  as.character(Metric) %in% HH_METRICS_COUNTRY) %>%
    dplyr::select(COUNTRY, Metric, Method, Dominance_Score) %>%
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
           n_pairs_typical, " country\u00D7metric paired observations.")
  } else if (length(sig_pairs_summary) == 3L) {
    "All three detector-pair comparisons reached p < 0.05."
  } else {
    paste0(length(sig_pairs_summary), " of 3 detector-pair comparisons ",
           "reached p < 0.05: ",
           paste(sig_pairs_summary, collapse = "; "), ".")
  }

  paste0(
    "Figure 6 (Method-Across-Metrics Summary) | Cross-metric summary of ",
    "detector performance across the 8 endemic countries and the five ",
    "reported operational metrics (TAM, Number of True Alarms, Sensitivity, ",
    "Mean Lead Time, Warning Persistence), 2016-2024 (excluding 2020, 2021, ",
    "2025; n = ", length(EVALUABLE_YEARS),
    " evaluable years per country; country inclusion floor = ",
    MIN_PEAK_CASES_PER_YEAR, " annual peak cases and >= ",
    MIN_EVALUABLE_YEARS_PER_COUNTRY, " evaluable years). This composite ",
    "figure complements the five per-metric multipanels by aggregating ",
    "detector performance across all five metrics simultaneously, ",
    "answering: which detector achieves the best operational performance ",
    "across the broadest combination of metrics and countries?\n",
    "(a) Country dominance matrix on five operational metrics. Rows = 11 ",
    "outbreak-detection methods, ordered top-to-bottom by mean sweep count ",
    "across the 5 metrics (descending; best-overall method at top). Columns ",
    "= the five reported metrics, ordered left-to-right as TAM, Number of ",
    "True Alarms, Sensitivity, Mean Lead Time, Warning Persistence. Cell ",
    "shading encodes the per-(method, metric) sweep count: the number of ",
    "countries (out of ", MAX_COUNTRYS_POSSIBLE,
    " evaluable) where that method achieves within-country dominance score ",
    ">= ", sprintf("%.2f", DOMINANCE_THRESHOLD),
    " on that metric. Single-color blue ramp from white (0 countries) through ",
    "mid-blue to dark navy (all ", MAX_COUNTRYS_POSSIBLE, " countries). ",
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
    "length per contrast is |median difference in countries-swept count across ",
    "the 5 metrics| (units: countries; computed from the 5 paired metric-level ",
    "counts per detector). The p-value per contrast is computed from a ",
    "separate properly-powered test: paired Wilcoxon signed-rank on within-",
    "country Dominance_Score, paired by (COUNTRY x METRIC) cells, n_pairs ~ ",
    n_pairs_typical, " country\u00D7metric paired observations. Bars are ",
    "navy (with adjacent '*') if p < 0.05; pale grey otherwise. The dual ",
    "specification answers two distinct questions in one panel: is the ",
    "detector difference statistically real (n_pairs ~ ", n_pairs_typical,
    " powered test), and how big is the practical difference in operationally ",
    "interpretable units (n=5 sweep-count median across the 5 metrics)? ",
    "Outcome: ", outcome, "\n"
  )
}

# Combined Figure 6 legend (one block per metric + Method Summary).
fig6_legend_lines <- c(
  paste0("================================================================"),
  paste0("FIGURE 6  |  LEGENDS  (per-metric and method summary)"),
  paste0("================================================================")
)
for (m in HH_METRICS_COUNTRY) {
  fig6_legend_lines <- c(
    fig6_legend_lines, "",
    paste0("--- Figure6_", hh_metric_slug[[m]], " ---"),
    build_figure6_metric_legend(m)
  )
}
fig6_legend_lines <- c(
  fig6_legend_lines, "",
  "--- Figure6_Method_Summary ---",
  build_figure6_summary_legend()
)
fig6_legend_text <- paste(fig6_legend_lines, collapse = "\n")
writeLines(fig6_legend_text,
           file.path(OUTPUT_DIR, "Stage5_Figure6_Legend.txt"))
cat("\nSaved: ", file.path(OUTPUT_DIR, "Stage5_Figure6_Legend.txt"),
    "\n", sep = "")

# Figure 7 legend
build_figure7_legend <- function() {
  bs <- N_BOOTS
  n_total_b  <- nrow(winner_rows)
  n_decisive <- sum(winner_rows$Dominance_Probability >= DOMINANCE_THRESHOLD,
                    na.rm = TRUE)
  n_contested_d <- sum(winner_rows$Dominance_Probability < DOMINANCE_THRESHOLD,
                       na.rm = TRUE)
  n_signif_b <- sum(winner_rows$sig_star_dominance %in% c("*", "**", "***"),
                    na.rm = TRUE)

  cross_country_lines <- if (nrow(level2_results) > 0L) {
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
    "        (Cross-country Wilcoxon table not available.)"
  }
  cross_country_block <- paste(cross_country_lines, collapse = "\n")

  paste0(
    "================================================================\n",
    "FIGURE 7  |  LEGEND\n",
    "================================================================\n\n",

    "Figure 7 | Bootstrap-supported country dominance of three target ",
    "outbreak detectors across the ", n_total_b, " dengue-endemic countries ",
    "(Brazil, Colombia, Mexico, Peru, Philippines, Singapore, Sri Lanka, ",
    "Taiwan). The figure is a single panel and therefore carries no panel ",
    "letter: a per-detector dot plot of country dominance probabilities with ",
    "a per-country tier-legend overlay. The per-country metric table ",
    "is not shown as a panel in the figure; those ",
    "values, together with the per-country per-metric significance results, ",
    "remain available as standalone CSV exports ",
    "(CountryTable.csv and Stage5_Country_Summary_PanelD_HeatmapData.csv).\n\n",

    "FOUR-TIER CONSENSUS CLASSIFICATION.\n",
    "Each country is classified into ONE of four tiers based on (i) the ",
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
    "  Primary  : factor 3 within-country * factor ", n_total_b,
    " between-country\n",
    "             = ", N_PAIRS * n_total_b,
    ". Used to assign the consensus tier and\n",
    "             render the figure.\n",
    "  Sensitivity (within-country only): factor 3. The within-country-\n",
    "             only adjustment is reported as a parallel column\n",
    "             ('p_*_bonf_within' suffixes) in CountryTable.csv. Under\n",
    "             the within-only correction, ", n_strong_within,
    " countries would meet 'strong' (vs ", n_strong, " under primary).\n\n",

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
    "N = ", bs, " bootstrap replicates per country).\n",
    sprintf("    Tier counts: strong=%d, partial=%d, lead_only=%d, contested=%d (of %d).\n\n",
            n_strong, n_partial, n_lead_only, n_contested, n_total_b),

    "(a) Per-country metric table with EMBEDDED PER-METRIC SIGNIFICANCE. ",
    "Nine columns: Country, Detector, Dominance probability, Sig., TAM, ",
    "N true alarms, Sensitivity, Mean lead time (wk), WP (wk). The Sig. ",
    "column reports a one-sample binomial test of the dominance probability ",
    "in column 3 against the chance baseline of 1/3 (three competing ",
    "detectors), Bonferroni-adjusted across the ", n_total_b,
    " countries; this directly tests the value displayed next to it, so a ",
    "high dominance probability earns a star regardless of whether the ",
    "all-pairs head-to-head test reaches significance under the much ",
    "stricter cross-country correction. Star encoding: *** p<0.001 (dark ",
    "navy); ** p<0.01 (mid navy); * p<0.05 (mid blue); ns otherwise (faded ",
    "grey). Each metric cell renders TWO stacked text lines: the metric ",
    "value (top, larger font) and a parenthetical significance label below ",
    "(smaller font). Per-metric significance is the WEAKEST-LINK ",
    "Bonferroni-adjusted one-sided p across the two ordered pairs the ",
    "consensus winner participates in (where the winner is X). For 'No ",
    "consensus' rows (lead_only and contested tiers) the fallback rule ",
    "uses the smallest p_bonf across all six ordered pairs (most ",
    "informative signal regardless of direction). The Detector cell is ",
    "colored by detector palette for strong (full saturation) and partial ",
    "(washed-out); light grey for lead_only; darker grey for contested. ",
    "These four tiers come from the all-pairs head-to-head test (a ",
    "separate, more conservative check) and answer a different question ",
    "than the Sig. column: does the leader beat BOTH rivals significantly ",
    "under the cross-country Bonferroni? The Detector cell text annotates ",
    "the tier with a '(partial)' or '(lead)' suffix where applicable.\n\n",

    "(b) Per-detector dot plot. One dot per country per detector, jittered ",
    "horizontally. Median crossbar per detector. The pale band over [",
    sprintf("%.2f", DOMINANCE_THRESHOLD), ", 1.0] is the decisive zone. ",
    "Dominance threshold (",
    sprintf("%.2f", DOMINANCE_THRESHOLD), "). The right side of the panel ",
    "hosts a per-country tier legend showing how each dot border encodes ",
    "the country's consensus tier:\n",
    "      strong leader    : black border, thickest stroke\n",
    "      partial leader   : black border, medium stroke\n",
    "      lead_only leader : grey border, medium stroke\n",
    "      non-leader / contested: grey border, thinnest stroke\n\n",

    "AUXILIARY OUTPUTS.\n\n",

    "  Per-country per-metric significance heatmap data.\n",
    "  File: Stage5_Country_Summary_PanelD_HeatmapData.csv\n",
    "  Rows: one per (sub-grid label, country, metric) cell - up to 6 ordered ",
    "pairs * 5 metrics * ", n_total_b, " countries = ", 6 * 5 * n_total_b,
    " entries.\n",
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

    "  Per-country per-metric full results.\n",
    "  File: Stage5_Country_Summary_PerMetricSignificance.csv\n",
    "  Method: paired bootstrap on year-cluster replicates, 6 ordered pairs.\n",
    sprintf("  Bonferroni: factor %d within (COUNTRY, ordered_pair).\n",
            BONF_LEVEL1_FACTOR),

    "\n  Cross-country per-metric significance.\n",
    "  File: Stage5_Country_Summary_CrossCountry_Wilcoxon.csv\n",
    "  Population-level Wilcoxon signed-rank test paired by country (n = ",
    if (nrow(level2_results) > 0L)
      max(level2_results$n_pairs, na.rm = TRUE) else "?",
    "). 5 metrics x 3 pairs = 15 tests.\n",
    "  Bonferroni factor 3 within metric.\n\n",
    "  Results (Bonferroni-adjusted, k=3 within metric):\n",
    cross_country_block, "\n\n",

    "STATISTICS SUMMARY:\n",
    sprintf("    Decisive countries  (Pr >= %.2f) : %d / %d\n",
            DOMINANCE_THRESHOLD, n_decisive, n_total_b),
    sprintf("    Contested by Pr threshold (Pr < %.2f) : %d / %d\n",
            DOMINANCE_THRESHOLD, n_contested_d, n_total_b),
    sprintf("    Tier counts: strong=%d, partial=%d, lead_only=%d, contested=%d (of %d).\n",
            n_strong, n_partial, n_lead_only, n_contested, n_total_b),
    sprintf("    Within-country-only Bonferroni sensitivity:\n"),
    sprintf("        countries meeting 'strong' under within-only: %d (vs %d primary).\n",
            n_strong_within, n_strong),
    sprintf("    Dominance Sig. p<0.05 (binomial vs chance 1/3, Bonf k=%d): %d / %d.\n",
            n_total_b, n_signif_b, n_total_b),
    sprintf("    Per-country per-metric tests: %d.\n",
            nrow(level1_results_ordered)),
    sprintf("    Per-country per-metric significant (one-sided p_bonf<0.05): %d / %d.\n",
            sum(level1_results_ordered$sig_star %in% c("*", "**", "***"),
                na.rm = TRUE),
            sum(!is.na(level1_results_ordered$p_bonf))),
    sprintf("    Cross-country tests: %d. Significant (Bonf k=%d): %d / %d.\n\n",
            nrow(level2_results), BONF_LEVEL2_FACTOR,
            sum(level2_results$sig_star %in% c("*", "**", "***"), na.rm = TRUE),
            nrow(level2_results)),

    "METHODS NOTE.\n",
    "Year-cluster (Cameron-Gelbach-Miller) bootstrap, B = ", bs,
    " replicates per country. The Round 1 composite mean-rank winner, the ",
    "Round 2 all-pairs head-to-head dominance check, the per-country per-",
    "metric ordered-pair tests, and the cross-country per-metric Wilcoxon ",
    "are all computed in this script.\n",
    "================================================================\n"
  )
}

fig7_legend_text <- build_figure7_legend()
writeLines(fig7_legend_text,
           file.path(OUTPUT_DIR, "Stage5_Figure7_Legend.txt"))
cat("Saved: ", file.path(OUTPUT_DIR, "Stage5_Figure7_Legend.txt"),
    "\n", sep = "")

# Echo legends to console.
cat("\n", fig6_legend_text, "\n", sep = "")
cat("\n", fig7_legend_text, "\n", sep = "")

# ==============================================================================
# 30. PARAMETER REPORT AND END-OF-RUN BANNER
# ==============================================================================
cat("\n============================================================\n")
cat("STAGE 5 - COUNTRY ANALYSIS  ", SCRIPT_VERSION, "\n", sep = "")
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
cat("\nCountry inclusion:\n")
cat("  MIN_PEAK_CASES_PER_YEAR        = ", MIN_PEAK_CASES_PER_YEAR,
    "\n", sep = "")
cat("  MIN_EVALUABLE_YEARS_PER_COUNTRY = ", MIN_EVALUABLE_YEARS_PER_COUNTRY,
    "\n", sep = "")
cat("  Excluded years                 = ",
    paste(EXCLUDED_YEARS, collapse = ", "), "\n", sep = "")
cat("\nBootstrap configuration:\n")
cat("  B = ", BOOT_N_CI, " year-cluster replicates per country\n", sep = "")
cat("  Wilcoxon detector-paired alpha (per pairwise comparison) = ",
    sprintf("%.2f", HH_ALPHA), "\n", sep = "")
cat("  Per-country per-metric Bonferroni factor                  = ",
    BONF_LEVEL1_FACTOR, " (within COUNTRY x ordered pair)\n", sep = "")
cat("  Cross-country per-metric Bonferroni factor                = ",
    BONF_LEVEL2_FACTOR, " (within metric)\n", sep = "")

cat("\n=== Saved files (", OUTPUT_DIR, ") ===\n", sep = "")
cat("Figure 6 (per-metric, 5 figures + 1 method summary):\n")
cat("  Figure6_TAM.{pdf,png}\n")
cat("  Figure6_N_True_Alarms.{pdf,png}\n")
cat("  Figure6_Sensitivity.{pdf,png}\n")
cat("  Figure6_Mean_Lead_Time.{pdf,png}\n")
cat("  Figure6_Warning_Persistence.{pdf,png}\n")
cat("  Figure6_Method_Summary.{pdf,png}\n")
cat("Figure 7 (two-panel country summary):\n")
cat("  Figure7_CountryDominanceProbability.{pdf,png}\n")
cat("Tables (Figure 6 supporting):\n")
cat("  Stage5_Country_Framework_Metrics.csv\n")
cat("  Stage5_Country_Framework_Metrics_with_CIs.csv\n")
cat("  Stage5_Country_8Metric_Summary.csv\n")
cat("  Stage5_Country_8Metric_Summary_with_CIs.csv\n")
cat("  Stage5_Country_Dominance_Matrix.csv\n")
cat("  Stage5_Country_Dominance_Probabilities.csv\n")
cat("  Stage5_Country_Wilcoxon_PerMetric.csv\n")
cat("  Stage5_Country_Wilcoxon_ConstantTA_vs_Comparators.csv\n")
cat("  Stage5_Method_Summary_Long.csv\n")
cat("  Stage5_Method_Summary_Aggregate.csv\n")
cat("  Stage5_Country_Bootstrap_Replicates.csv\n")
cat("Tables (Figure 7 supporting):\n")
cat("  Stage5_Country_Summary_CountryTable.csv\n")
cat("  Stage5_Country_Summary_PerMetricSignificance.csv\n")
cat("  Stage5_Country_Summary_CrossCountry_Wilcoxon.csv\n")
cat("  Stage5_Country_Summary_PanelB_DotPlot.csv\n")
cat("  Stage5_Country_Summary_PanelB_AuxWilcoxon.csv\n")
cat("  Stage5_Country_Summary_PanelD_HeatmapData.csv\n")
cat("  Stage5_Country_ConstantTA_4Comparator_Consensus_Long.csv\n")
cat("  Stage5_Country_ConstantTA_4Comparator_Consensus.csv\n")
cat("Legends:\n")
cat("  Stage5_Figure6_Legend.txt\n")
cat("  Stage5_Figure7_Legend.txt\n")

cat("\nDone.\n")

# ==============================================================================
# END OF SCRIPT — STAGE 5 COUNTRY ANALYSIS
# ==============================================================================
