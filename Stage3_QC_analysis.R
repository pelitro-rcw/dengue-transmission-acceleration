# =============================================================================
# STAGE 3 - Quezon City analysis (Figures 2 and 3, plus tables)
# -----------------------------------------------------------------------------
# Part A : Setup (packages, paths, shared helpers)
# Part B : Figure 2 - Comparative dashboard of 11 outbreak-detection methods
#                     Panel A - dominance matrix across 8 performance metrics
#                     Panel B - True-Alarm Magnitude vs Mean Lead Time scatter
# Part C : Figure 3 - Per-year detection timing for Constant Transmission
#                     Acceleration (Figure 3A) and Continuous Transmission
#                     Acceleration (Supplementary Figure 3); three head-to-head
#                     Wilcoxon comparison figures
# Tables : Tables 1, 1b, 2, 2A, 2B; sensitivity tables S1-S3, S5;
#          head-to-head Wilcoxon results table
#
# TWO-ANCHOR TRUE/FALSE-ALARM FRAMEWORK
#   A1 = Pre-peak Actionable Window  (peak - 8 ... peak - 4 weeks)
#   A2 = Epidemic Burden block       (smallest contiguous block containing
#                                     the peak whose cumulative cases sum
#                                     to >= 70% of the annual total)
#   A trigger at week t is a True Alarm iff t in A1 OR t in A2.
#   True alarms partition into two compartments:
#     Actionable : lead in [4, 8] weeks (i.e., trigger lands in A1)
#     Reactive   : any other true alarm (lead < 4 weeks; A2-only hits)
#   Triggers with lead >= 9 weeks are False Alarms (outside both anchors).
#
# SAME-DENOMINATOR TIMELINESS METRICS
#   Mean Lead Time, Warning Persistence (WP), and Actionable Lead-Time
#   Yield (ALY) are computed with a uniform per-detector denominator
#   (count of evaluable years). Years where the metric is not computable
#   due to no qualifying triggers contribute zero rather than being
#   excluded. Conditional-on-firing diagnostic versions are preserved
#   as *_conditional columns in the CSV outputs.
#
# DETECTOR PARAMETERS
#   Constant Transmission Acceleration (Constant TA):
#     short window = 4 wk, long window = 26 wk, guard = 2 wk,
#     eta_on = 1.33, eta_off = 0.73
#   Continuous Transmission Acceleration (Continuous TA):
#     short window = 3 wk, long window = 12 wk, threshold = 1.33
#   Outbreak Threshold:
#     rolling week-specific Mean + 2 SD baseline
#
# YEAR EXCLUSIONS
#   2020, 2021 (COVID-19 surveillance disruption)
#   2025       (out-of-distribution / truncated reporting)
#
# Inputs   : Dengue-Rainfall_Dataset.xlsx, sheet "QC Data"
# Outputs  : Stage3_QC_analysis/                       (Figure 2 / 3 PDF + PNG)
#            Stage3_QC_analysis/evaluation_framework/  (CSV tables)
# Repro    : set.seed(12345); R >= 4.1
# =============================================================================


# =============================================================================
# PART A - PACKAGES, PATHS, AND SHARED HELPERS
# =============================================================================

required_packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "rlang",
  "ggplot2", "cowplot", "zoo", "ISOweek", "scales",
  "tibble", "grid", "patchwork", "ggrepel"
)
optional_packages <- c("irr")

installed_pkgs <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed_pkgs) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
for (pkg in optional_packages) {
  if (!pkg %in% installed_pkgs) {
    try(install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org"),
        silent = TRUE)
  }
  if (pkg %in% rownames(installed.packages())) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

set.seed(12345)
options(scipen = 999)

cat("\n=============================================================\n")
cat("STAGE 3 - Quezon City analysis (Figures 2 and 3, plus tables)\n")
cat("Two-anchor framework (T = A1 union A2)\n")
cat("Same-denominator early-warning timeliness metrics\n")
cat("=============================================================\n\n")


# -----------------------------------------------------------------------------
# Shared paths
# -----------------------------------------------------------------------------
PATH       <- "C:/Users/User/Desktop/Dengue-Rainfall_Dataset.xlsx"
SHEET_NAME <- "QC Data"

OUT_DIR <- file.path(getwd(), "Stage3_QC_analysis")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

EVAL_OUT_DIR <- file.path(OUT_DIR, "evaluation_framework")
if (!dir.exists(EVAL_OUT_DIR)) dir.create(EVAL_OUT_DIR, recursive = TRUE)


# -----------------------------------------------------------------------------
# Shared graphics helpers
# -----------------------------------------------------------------------------
available_fonts <- names(grDevices::pdfFonts())
base_family_global <- if ("Arial"     %in% available_fonts) "Arial"     else
  if ("Helvetica" %in% available_fonts) "Helvetica" else
    "sans"

safe_pdf_device <- function() {
  cap <- tryCatch(capabilities("cairo"), error = function(e) FALSE)
  if (isTRUE(cap["cairo"])) grDevices::cairo_pdf else grDevices::pdf
}

safe_quantile <- function(x, probs, na.rm = TRUE) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  suppressWarnings(as.numeric(stats::quantile(
    x, probs = probs, na.rm = na.rm, names = FALSE
  )))
}

safe_df_print <- function(x, title = NULL, round_cols = NULL, digits = 3) {
  if (!is.null(title)) cat("\n", title, "\n", sep = "")
  out <- as.data.frame(x)
  if (!is.null(round_cols)) {
    keep_cols <- intersect(round_cols, names(out))
    for (nm in keep_cols) out[[nm]] <- round(out[[nm]], digits)
  }
  print(out, row.names = FALSE, na.print = "NA")
}

rescale_01_safe <- function(x) {
  finite_x <- x[is.finite(x)]
  if (length(finite_x) < 2) return(rep(0.5, length(x)))
  xmin <- min(finite_x); xmax <- max(finite_x)
  if (!is.finite(xmin) || !is.finite(xmax) || xmax == xmin) return(rep(0.5, length(x)))
  (x - xmin) / (xmax - xmin)
}

mix_with_white <- function(color, strength) {
  strength <- pmin(pmax(strength, 0), 1)
  rgb_base <- grDevices::col2rgb(color) / 255
  rgb_mix  <- (1 - strength) * 1 + strength * rgb_base
  grDevices::rgb(rgb_mix[1], rgb_mix[2], rgb_mix[3])
}


# =============================================================================
# PART B - FIGURE 2: COMPARATIVE DASHBOARD OF 11 DETECTION METHODS
# =============================================================================

cat("\n-------------------------------------------------------------\n")
cat("PART B  Figure 2  Comparative dashboard\n")
cat("-------------------------------------------------------------\n")


# -----------------------------------------------------------------------------
# B.1 SETTINGS
# -----------------------------------------------------------------------------
panelA_width_in  <- 14.4
panelA_height_in <- 9.8

panelB_width_in  <- 11.2
panelB_height_in <- 8.6

BOOT_N_CI     <- 1000
SUSTAINED_RUN <- 2L

target_years    <- 2013:2025
EXCLUDED_YEARS  <- c(2020L, 2021L, 2025L)
EVALUABLE_YEARS <- setdiff(target_years, EXCLUDED_YEARS)   # 2013-19, 2022-24

A1_LEAD_MIN     <- 4L
A1_LEAD_MAX     <- 8L
A2_BURDEN_FRAC  <- 0.70

# Lead-time compartments (defined for True alarms only).
#   Actionable: lead in [A1_LEAD_MIN, A1_LEAD_MAX]   (matches A1 exactly)
#   Reactive  : any other true alarm (lead < 4 wk; A2-only hits)
# Pre-A1 too-early triggers (lead >= 9) are FALSE alarms and have no
# compartment label.
COMPARTMENT_ACTIONABLE_MIN <- A1_LEAD_MIN
COMPARTMENT_ACTIONABLE_MAX <- A1_LEAD_MAX

# STA/LTA detector thresholds
ETA_ON_CLASSIC <- 1.33
ETA_ON         <- 1.33
ETA_OFF        <- 0.73


save_plot_pair <- function(filename, plot, width, height, dpi = 600, device = NULL) {
  ggplot2::ggsave(
    filename  = file.path(OUT_DIR, filename),
    plot      = plot,
    device    = device,
    width     = width,
    height    = height,
    units     = "in",
    dpi       = dpi,
    limitsize = FALSE,
    bg        = "white"
  )
  cat("Saved: ", file.path(OUT_DIR, filename), "\n", sep = "")
}


# -----------------------------------------------------------------------------
# B.2 COLOR SYSTEM
# -----------------------------------------------------------------------------
TYPE_COLORS <- c(
  "Surveillance-Guideline Percentile Thresholds"        = "#F8766D",
  "Retrospective Thresholds" = "#00BA38",
  "Acceleration Measures"    = "#619CFF"
)

type_order <- c(
  "Surveillance-Guideline Percentile Thresholds",
  "Retrospective Thresholds",
  "Acceleration Measures"
)

method_group_order <- c(
  "Constant Transmission Acceleration",
  "Continuous Transmission Acceleration",
  "Composite Outbreak Signal",
  "Critical Transition Indicator",
  "Incidence Gradient",
  "Hydrological Inflection Measure",
  "Cumulative Sum Control",
  "Outbreak Threshold",
  "Alarm Threshold",
  "WHO 90th Percentile Threshold",
  "WHO 75th Percentile Threshold"
)

method_label_map <- c(
  "Constant Transmission Acceleration"   = "Constant Transmission\nAcceleration",
  "Continuous Transmission Acceleration" = "Continuous Transmission\nAcceleration",
  "Composite Outbreak Signal"            = "Composite Outbreak\nSignal",
  "Critical Transition Indicator"        = "Critical Transition\nIndicator",
  "Incidence Gradient"                   = "Incidence\nGradient",
  "Hydrological Inflection Measure"      = "Hydrological Inflection\nMeasure",
  "Cumulative Sum Control"               = "Cumulative Sum\nControl",
  "Outbreak Threshold"                   = "Outbreak\nThreshold",
  "Alarm Threshold"                      = "Alarm\nThreshold",
  "WHO 90th Percentile Threshold"        = "WHO 90th Percentile\nThreshold",
  "WHO 75th Percentile Threshold"        = "WHO 75th Percentile\nThreshold"
)


# -----------------------------------------------------------------------------
# B.3 LOAD DATA
# -----------------------------------------------------------------------------
if (!file.exists(PATH)) stop("Data file not found at:\n", PATH)

df_raw <- readxl::read_excel(PATH, sheet = SHEET_NAME)

required_cols <- c("YR", "WN", "DC_QC")
missing_cols  <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sheet '", SHEET_NAME, "': ",
       paste(missing_cols, collapse = ", "))
}

df <- df_raw %>%
  dplyr::mutate(
    YR    = suppressWarnings(as.integer(YR)),
    WN    = suppressWarnings(as.integer(WN)),
    DC_QC = suppressWarnings(as.numeric(DC_QC))
  ) %>%
  dplyr::filter(!is.na(YR), !is.na(WN), WN >= 1, WN <= 53) %>%
  dplyr::mutate(
    ISOweek = sprintf("%d-W%02d", YR, WN),
    Date    = ISOweek::ISOweek2date(paste0(ISOweek, "-1"))
  ) %>%
  dplyr::filter(!is.na(Date)) %>%
  dplyr::arrange(Date)

if (nrow(df) == 0) stop("No valid rows remain after cleaning.")
if (all(is.na(df$DC_QC))) stop("Column 'DC_QC' contains only NA values.")


# -----------------------------------------------------------------------------
# B.4 ROLLING WEEK-SPECIFIC BASELINES
# -----------------------------------------------------------------------------
build_rolling_map <- function(targets, donor_pool_years, exclude_years = integer(0)) {
  out <- list()
  donor_pool <- sort(setdiff(unique(donor_pool_years), exclude_years))
  for (y in targets) {
    donors <- intersect(seq(y - 5, y - 1), donor_pool)
    donors <- setdiff(donors, exclude_years)
    if (length(donors) < 3) donors <- tail(donor_pool[donor_pool < y], 5)
    out[[as.character(y)]] <- sort(unique(donors))
  }
  out
}

donor_pool_all   <- sort(unique(df$YR))
rolling_map_full <- build_rolling_map(target_years, donor_pool_all)

compute_rolling_weekly_thresholds <- function(df_input, rolling_map) {
  rows <- list()
  for (y in names(rolling_map)) {
    donors <- rolling_map[[y]]; y_int <- as.integer(y)
    weeks_y <- df_input %>% dplyr::filter(YR == y_int) %>%
      dplyr::pull(WN) %>% unique() %>% sort()
    if (length(weeks_y) == 0) next
    for (w in weeks_y) {
      vals <- df_input %>% dplyr::filter(YR %in% donors, WN == w) %>% dplyr::pull(DC_QC)
      vals <- vals[is.finite(vals)]; n_vals <- length(vals)
      if (n_vals == 0) {
        mean_val <- NA_real_; sd_val <- NA_real_; p75 <- NA_real_; p90 <- NA_real_
      } else if (n_vals == 1) {
        mean_val <- vals[1]; sd_val <- NA_real_; p75 <- vals[1]; p90 <- vals[1]
      } else {
        mean_val <- mean(vals, na.rm = TRUE)
        sd_val   <- stats::sd(vals, na.rm = TRUE)
        p75      <- safe_quantile(vals, 0.75)
        p90      <- safe_quantile(vals, 0.90)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        YR = y_int, WN = as.integer(w),
        bl_mean = mean_val, bl_sd = sd_val, bl_p75 = p75, bl_p90 = p90,
        bl_alarm    = if (is.na(sd_val)) mean_val else mean_val + 1 * sd_val,
        bl_outbreak = if (is.na(sd_val)) mean_val else mean_val + 2 * sd_val,
        bl_n = n_vals,
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(rows) %>% dplyr::arrange(YR, WN)
}

thresholds_full <- compute_rolling_weekly_thresholds(df, rolling_map_full)
df <- df %>% dplyr::left_join(thresholds_full, by = c("YR", "WN"))


# -----------------------------------------------------------------------------
# B.5 DERIVED SERIES
# -----------------------------------------------------------------------------
dc_ma3  <- zoo::rollmean(df$DC_QC,  3, fill = NA, align = "right")
dc_ma12 <- zoo::rollmean(df$DC_QC, 12, fill = NA, align = "right")
dc_diff <- c(NA_real_, diff(dc_ma3))

var8  <- zoo::rollapply(df$DC_QC,  8,
                        function(x) stats::var(x, na.rm = TRUE), fill = NA, align = "right")
var26 <- zoo::rollapply(df$DC_QC, 26,
                        function(x) stats::var(x, na.rm = TRUE), fill = NA, align = "right")

VAR26_GUARD <- 0.5
ratio_var <- ifelse(
  !is.na(var26) & is.finite(var26) & var26 > VAR26_GUARD,
  var8 / var26, NA_real_
)

rolling_quantile_by_target <- function(series, years_vec, rolling_map, prob) {
  out <- rep(NA_real_, length(series))
  for (y in names(rolling_map)) {
    donors <- rolling_map[[y]]; y_int <- as.integer(y)
    donor_idx <- which(years_vec %in% donors & is.finite(series))
    if (length(donor_idx) < 5) next
    out[which(years_vec == y_int)] <- safe_quantile(series[donor_idx], prob)
  }
  out
}

df$rise_thr_rolling  <- rolling_quantile_by_target(dc_diff,   df$YR, rolling_map_full, 0.75)
df$roc80_thr_rolling <- rolling_quantile_by_target(dc_diff,   df$YR, rolling_map_full, 0.80)
df$ct80_thr_rolling  <- rolling_quantile_by_target(ratio_var, df$YR, rolling_map_full, 0.80)


# -----------------------------------------------------------------------------
# B.6 WALK-FORWARD CUSUM
# -----------------------------------------------------------------------------
k_mult <- 0.5
h_mult <- 5.0

n_rows <- nrow(df)
cusum_vec <- rep(0, n_rows)
cusum_h   <- rep(NA_real_, n_rows)

for (i in seq_len(n_rows)) {
  is_year_boundary <- i > 1 &&
    !is.na(df$YR[i]) && !is.na(df$YR[i - 1]) && df$YR[i] != df$YR[i - 1]
  mu0    <- df$bl_mean[i]; sigma0 <- df$bl_sd[i]
  prev_val <- if (i == 1 || is_year_boundary) 0 else cusum_vec[i - 1]
  if (!is.na(df$DC_QC[i]) && !is.na(mu0) && !is.na(sigma0) && sigma0 > 0) {
    cusum_vec[i] <- max(0, prev_val + (df$DC_QC[i] - mu0 - k_mult * sigma0))
  } else {
    cusum_vec[i] <- prev_val
  }
  cusum_h[i] <- if (!is.na(sigma0) && sigma0 > 0) h_mult * sigma0 else NA_real_
}

df$cusum_val <- cusum_vec
df$cusum_h   <- cusum_h


# -----------------------------------------------------------------------------
# B.7 SURGE DEFINITIONS
# -----------------------------------------------------------------------------
df <- df %>%
  dplyr::mutate(
    surge_mean_2sd       = as.integer(!is.na(DC_QC) & !is.na(bl_outbreak) & DC_QC > bl_outbreak),
    surge_mean_1sd       = as.integer(!is.na(DC_QC) & !is.na(bl_alarm)    & DC_QC > bl_alarm),
    surge_who_75         = as.integer(!is.na(DC_QC) & !is.na(bl_p75)      & DC_QC > bl_p75),
    surge_who_90         = as.integer(!is.na(DC_QC) & !is.na(bl_p90)      & DC_QC > bl_p90),
    surge_rate_change    = as.integer(!is.na(dc_diff) & !is.na(roc80_thr_rolling) & dc_diff > roc80_thr_rolling),
    surge_sta_lta        = as.integer(!is.na(dc_ma3)  & !is.na(dc_ma12)   & dc_ma12 > 0 & (dc_ma3 / dc_ma12) > ETA_ON_CLASSIC),
    surge_cusum          = as.integer(!is.na(cusum_val) & !is.na(cusum_h) & cusum_val > cusum_h),
    surge_critical_trans = as.integer(!is.na(ratio_var) & !is.na(ct80_thr_rolling) & ratio_var > ct80_thr_rolling),
    surge_hydrology      = as.integer(
      !is.na(DC_QC) & !is.na(dc_diff) & !is.na(bl_p75) & !is.na(rise_thr_rolling) &
        DC_QC > bl_p75 & dc_diff > rise_thr_rolling
    )
  )


# -----------------------------------------------------------------------------
# B.8 CONSTANT TA (Vaezi-style STA/LTA with hysteresis)
# -----------------------------------------------------------------------------
STA_WIN <- 4L; LTA_WIN <- 26L; GUARD <- 2L
MIN_T   <- STA_WIN + GUARD + LTA_WIN
MIN_OFF_RESET <- 8L

n_rows <- nrow(df); dc <- df$DC_QC
R_vaezi_v   <- rep(NA_real_, n_rows)
triggered_v <- rep(FALSE, n_rows)
is_on <- FALSE; frozen_lta <- NA_real_; consec_off <- 0L

for (t in seq_len(n_rows)) {
  if (t > 1 && !is.na(df$YR[t]) && !is.na(df$YR[t - 1]) && df$YR[t] != df$YR[t - 1]) {
    is_on <- FALSE; frozen_lta <- NA_real_; consec_off <- 0L
  }
  if (!is_on) consec_off <- consec_off + 1L else consec_off <- 0L
  if (!is_on && consec_off >= MIN_OFF_RESET) frozen_lta <- NA_real_
  if (t < MIN_T) next
  sta_vals <- dc[(t - STA_WIN + 1):t]
  sta <- if (all(is.na(sta_vals))) NA_real_ else mean(sta_vals, na.rm = TRUE)
  if (!is_on || is.na(frozen_lta)) {
    lta_idx <- (t - STA_WIN - GUARD - LTA_WIN + 1):(t - STA_WIN - GUARD)
    if (length(lta_idx) == LTA_WIN && min(lta_idx) > 0) {
      lta_vals <- dc[lta_idx]
      frozen_lta <- if (all(is.na(lta_vals))) NA_real_ else mean(lta_vals, na.rm = TRUE)
    } else {
      frozen_lta <- NA_real_
    }
  }
  R_t <- if (!is.na(frozen_lta) && frozen_lta > 0 && !is.na(sta)) sta / frozen_lta else NA_real_
  R_vaezi_v[t] <- R_t
  if (!is_on && !is.na(R_t) && R_t >= ETA_ON) { is_on <- TRUE; consec_off <- 0L }
  if (is_on && !is.na(R_t) && R_t < ETA_OFF)  { is_on <- FALSE; frozen_lta <- NA_real_ }
  triggered_v[t] <- is_on
}

df$surge_sta_lta_vaezi <- as.integer(triggered_v)
df$R_vaezi <- R_vaezi_v


# -----------------------------------------------------------------------------
# B.9 COMPOSITE SIGNAL
# -----------------------------------------------------------------------------
df$surge_composite <- as.integer(
  df$surge_mean_2sd == 1L &
    (df$surge_sta_lta_vaezi == 1L | df$surge_critical_trans == 1L)
)


# -----------------------------------------------------------------------------
# B.10 METHOD MAP
# -----------------------------------------------------------------------------
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
  "WHO 75th Percentile Threshold"        = "surge_who_75"
)

missing_surge_cols <- setdiff(unname(surge_defs), names(df))
if (length(missing_surge_cols) > 0) {
  stop("Missing surge method columns: ", paste(missing_surge_cols, collapse = ", "))
}

method_type_map <- tibble::tibble(
  Method = names(surge_defs),
  Method_Type = c(
    "Acceleration Measures",
    "Acceleration Measures",
    "Acceleration Measures",
    "Acceleration Measures",
    "Retrospective Thresholds",
    "Acceleration Measures",
    "Retrospective Thresholds",
    "Surveillance-Guideline Percentile Thresholds",
    "Acceleration Measures",
    "Acceleration Measures",
    "Surveillance-Guideline Percentile Thresholds"
  )
)


# -----------------------------------------------------------------------------
# B.11 OPERATIONAL HELPERS
# -----------------------------------------------------------------------------
first_trigger_index <- function(trig_vec) {
  if (length(trig_vec) == 0) return(NA_integer_)
  hits <- which(!is.na(trig_vec) & trig_vec == 1L)
  if (length(hits) == 0) return(NA_integer_)
  as.integer(hits[1])
}

peak_index_whichmax <- function(dc_vec) {
  if (length(dc_vec) == 0) return(NA_integer_)
  valid_dc <- ifelse(is.na(dc_vec), -Inf, dc_vec)
  if (all(!is.finite(valid_dc))) return(NA_integer_)
  pk <- which.max(valid_dc)
  if (length(pk) == 0 || !is.finite(valid_dc[pk])) return(NA_integer_)
  as.integer(pk)
}


# -----------------------------------------------------------------------------
# B.12 TRUE / FALSE ALARM FRAMEWORK HELPERS
# -----------------------------------------------------------------------------
compute_A1 <- function(df_in, year, lead_min = A1_LEAD_MIN, lead_max = A1_LEAD_MAX) {
  df_y <- df_in %>% dplyr::filter(YR == year)
  if (nrow(df_y) == 0 || all(is.na(df_y$DC_QC)))
    return(list(peak_week = NA_integer_, A1_weeks = integer(0)))
  peak_idx <- which.max(df_y$DC_QC); peak_week <- df_y$WN[peak_idx]
  A1_weeks <- (peak_week - lead_max):(peak_week - lead_min)
  A1_weeks <- A1_weeks[A1_weeks >= 1L]
  list(peak_week = as.integer(peak_week), A1_weeks = as.integer(A1_weeks))
}

compute_A2 <- function(df_in, year, burden_frac = A2_BURDEN_FRAC) {
  df_y <- df_in %>% dplyr::filter(YR == year) %>% dplyr::arrange(WN)
  if (nrow(df_y) == 0 || all(is.na(df_y$DC_QC)))
    return(list(start_week = NA_integer_, end_week = NA_integer_, A2_weeks = integer(0)))
  cases <- ifelse(is.na(df_y$DC_QC), 0, df_y$DC_QC)
  weeks <- df_y$WN; n_w <- length(weeks); total <- sum(cases)
  if (total <= 0)
    return(list(start_week = NA_integer_, end_week = NA_integer_, A2_weeks = integer(0)))
  peak_idx <- which.max(cases); lo <- hi <- peak_idx; S <- cases[peak_idx]
  while (S / total < burden_frac && (lo > 1L || hi < n_w)) {
    L <- if (lo > 1L)  cases[lo - 1L] else -Inf
    R <- if (hi < n_w) cases[hi + 1L] else -Inf
    if (L >= R) { lo <- lo - 1L; S <- S + L } else { hi <- hi + 1L; S <- S + R }
  }
  list(start_week = as.integer(weeks[lo]), end_week = as.integer(weeks[hi]),
       A2_weeks = as.integer(weeks[lo:hi]))
}

compute_anchors_for_year <- function(df_in, year,
                                     lead_min = A1_LEAD_MIN, lead_max = A1_LEAD_MAX,
                                     burden_frac = A2_BURDEN_FRAC) {
  A1 <- compute_A1(df_in, year, lead_min, lead_max)
  A2 <- compute_A2(df_in, year, burden_frac)
  list(year = year, peak_week = A1$peak_week,
       A1_weeks = A1$A1_weeks, A2_weeks = A2$A2_weeks)
}

# Two-anchor classification: True Alarm iff t in A1 OR t in A2.
classify_trigger <- function(trigger_week, A1_weeks, A2_weeks) {
  in_A1 <- trigger_week %in% A1_weeks
  in_A2 <- trigger_week %in% A2_weeks
  is_true <- in_A1 || in_A2
  list(is_true = is_true, in_A1 = in_A1, in_A2 = in_A2)
}

# Two-bucket compartment scheme on True alarms only.
classify_compartment <- function(lead_time) {
  if (is.na(lead_time))                                return(NA_character_)
  if (lead_time >= COMPARTMENT_ACTIONABLE_MIN &&
      lead_time <= COMPARTMENT_ACTIONABLE_MAX)         return("Actionable")
  return("Reactive")
}


# -----------------------------------------------------------------------------
# B.13 PER-YEAR OPERATIONAL METRICS (TAM)
# -----------------------------------------------------------------------------
# TAM (True-Alarm Magnitude): per-year sum of weekly cases at weeks where
# the detector is alarm-on AND the week falls within T = A1 union A2 (the
# anchor-consistent epidemic-burden volume captured by true alarms).
compute_year_metrics <- function(data, year_label, surge_defs, df_full = NULL) {
  if (nrow(data) == 0) return(NULL)
  peak_idx <- peak_index_whichmax(data$DC_QC)
  if (is.na(peak_idx)) return(NULL)
  
  anchors_src <- if (is.null(df_full)) data else df_full
  anchors <- compute_anchors_for_year(anchors_src, as.integer(year_label))
  T_weeks <- union(anchors$A1_weeks, anchors$A2_weeks)
  
  out <- lapply(names(surge_defs), function(method_name) {
    col_name  <- surge_defs[[method_name]]
    trig_vec  <- data[[col_name]]
    triggered <- !is.na(trig_vec) & trig_vec == 1
    valid_dc  <- ifelse(is.na(data$DC_QC), 0, data$DC_QC)
    weeks_v   <- as.integer(data$WN)
    
    in_T_mask <- weeks_v %in% T_weeks
    on_in_T   <- triggered & in_T_mask
    tam_accum <- sum(valid_dc[on_in_T], na.rm = TRUE)
    n_on_in_T <- sum(on_in_T, na.rm = TRUE)
    
    trigger_rate <- if (length(triggered) > 0) mean(triggered, na.rm = TRUE) else NA_real_
    
    data.frame(
      Year = as.integer(year_label),
      Method = method_name,
      TAM_Accum    = tam_accum,
      N_Trig_Weeks = n_on_in_T,
      Trigger_Rate = trigger_rate,
      Is_Excluded  = as.integer(as.integer(year_label) %in% EXCLUDED_YEARS),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

yearly_metrics <- dplyr::bind_rows(
  lapply(target_years, function(yr) {
    year_data <- df %>% dplyr::filter(YR == yr)
    compute_year_metrics(year_data, yr, surge_defs, df_full = df)
  })
)
if (nrow(yearly_metrics) == 0) stop("No yearly metrics were generated.")


# -----------------------------------------------------------------------------
# B.14 BUILD TRIGGER DETAIL ACROSS ALL DETECTORS
# -----------------------------------------------------------------------------
build_trigger_detail <- function(df_in, surge_defs, evaluable_years) {
  rows <- list()
  for (yr in evaluable_years) {
    anchors <- compute_anchors_for_year(df_in, yr)
    df_y <- df_in %>% dplyr::filter(YR == yr) %>% dplyr::arrange(WN)
    for (method_name in names(surge_defs)) {
      col_name <- surge_defs[[method_name]]
      trig_idx <- which(df_y[[col_name]] == 1L)
      if (length(trig_idx) == 0L) next
      for (k in trig_idx) {
        wk  <- df_y$WN[k]
        cls <- classify_trigger(wk, anchors$A1_weeks, anchors$A2_weeks)
        lead_time <- if (!is.na(anchors$peak_week)) as.integer(anchors$peak_week - wk) else NA_integer_
        comp <- if (isTRUE(cls$is_true)) classify_compartment(lead_time) else NA_character_
        rows[[length(rows) + 1L]] <- data.frame(
          Year = yr, Method = method_name, Week = as.integer(wk),
          Peak_Week = anchors$peak_week, Lead_Time = lead_time,
          Compartment = comp,
          IsTrue = cls$is_true,
          InA1 = cls$in_A1, InA2 = cls$in_A2,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L)
    return(data.frame(Year = integer(), Method = character(), Week = integer(),
                      Peak_Week = integer(), Lead_Time = integer(),
                      Compartment = character(), IsTrue = logical(),
                      InA1 = logical(), InA2 = logical(),
                      stringsAsFactors = FALSE))
  dplyr::bind_rows(rows)
}

trigger_detail_aug <- build_trigger_detail(df, surge_defs, EVALUABLE_YEARS)


# -----------------------------------------------------------------------------
# B.15 PER-YEAR FIRST-A1-TRUE-ALARM LEAD TIME (Mean Lead Time metric)
# -----------------------------------------------------------------------------
# For each (Method, Year) pair, find the FIRST trigger that is both
# (a) classified as True, AND (b) lies in A1. The year's lead time =
# peak_week - first_a1_true_week. Bounded [4, 8].
compute_yearly_lead_data <- function(trig_aug, evaluable_years) {
  if (nrow(trig_aug) == 0)
    return(data.frame(Method = character(), Year = integer(),
                      First_A1_True_Week = integer(), Lead_Time_Yr = numeric(),
                      stringsAsFactors = FALSE))
  trig_aug %>%
    dplyr::filter(Year %in% evaluable_years, InA1, IsTrue) %>%
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

yearly_lead_data <- compute_yearly_lead_data(trigger_detail_aug, EVALUABLE_YEARS)


# -----------------------------------------------------------------------------
# B.16 AGGREGATE PER METHOD WITH BOOTSTRAP CIs (8 framework metrics)
# -----------------------------------------------------------------------------
# Metric directionality:
#   Epidemic Burden & Alarm Accuracy
#     TAM, N_True_Alarms, PPV, Sensitivity   - higher is better
#   Early Warning Timeliness
#     Mean_Lead_Time, WP, ALY                - higher is better
#   False Alarms
#     N_False_Alarms                         - lower is better

compute_method_point <- function(method_name, ym, trig_aug, yearly_lead, evaluable_years) {
  ym_sub    <- ym         %>% dplyr::filter(Method == method_name, Year %in% evaluable_years)
  trig_sub  <- trig_aug   %>% dplyr::filter(Method == method_name, Year %in% evaluable_years)
  lead_sub  <- yearly_lead%>% dplyr::filter(Method == method_name, Year %in% evaluable_years)
  
  tam       <- mean(ym_sub$TAM_Accum, na.rm = TRUE)
  trig_rate <- mean(ym_sub$Trigger_Rate, na.rm = TRUE)
  
  total       <- nrow(trig_sub)
  true_n      <- sum(trig_sub$IsTrue, na.rm = TRUE)
  false_n     <- total - true_n
  reactive_n  <- sum(trig_sub$Compartment == "Reactive",   na.rm = TRUE)
  truact      <- trig_sub %>% dplyr::filter(IsTrue, Compartment == "Actionable")
  truact_n    <- nrow(truact)
  
  n_eval_years <- length(evaluable_years)
  n_true_alarms_per_year  <- if (n_eval_years > 0) true_n  / n_eval_years else NA_real_
  n_false_alarms_per_year <- if (n_eval_years > 0) false_n / n_eval_years else NA_real_
  
  ppv  <- if (total > 0)  true_n / total                       else NA_real_
  
  # Sensitivity: A1-restricted (years with >= 1 True alarm in A1).
  years_with_a1_true <- length(unique(
    trig_sub$Year[trig_sub$IsTrue & trig_sub$Compartment == "Actionable"]
  ))
  sens <- if (n_eval_years > 0) years_with_a1_true / n_eval_years else NA_real_
  
  # Same-denominator timeliness: zero-coerce non-firing years; preserve
  # conditional-on-firing diagnostic versions.
  mean_lead_conditional <- if (nrow(lead_sub) > 0)
    mean(lead_sub$Lead_Time_Yr, na.rm = TRUE)
  else NA_real_
  mean_lead <- if (n_eval_years > 0)
    sum(lead_sub$Lead_Time_Yr, na.rm = TRUE) / n_eval_years
  else NA_real_
  
  wp_conditional <- if (truact_n > 0) mean(truact$Lead_Time, na.rm = TRUE)
  else NA_real_
  if (n_eval_years > 0L) {
    per_year_wp <- vapply(evaluable_years, function(y) {
      lt <- truact$Lead_Time[truact$Year == y]
      if (length(lt) > 0) mean(lt, na.rm = TRUE) else 0
    }, numeric(1))
    wp <- mean(per_year_wp, na.rm = TRUE)
    if (is.nan(wp)) wp <- NA_real_
  } else {
    wp <- NA_real_
  }
  
  aly_conditional <- if (true_n > 0) truact_n / true_n else NA_real_
  if (n_eval_years > 0L) {
    per_year_aly <- vapply(evaluable_years, function(y) {
      yr_trues <- sum(trig_sub$IsTrue[trig_sub$Year == y], na.rm = TRUE)
      yr_truact <- sum(trig_sub$IsTrue[trig_sub$Year == y] &
                         trig_sub$Compartment[trig_sub$Year == y] == "Actionable",
                       na.rm = TRUE)
      if (yr_trues > 0L) yr_truact / yr_trues else 0
    }, numeric(1))
    aly <- mean(per_year_aly, na.rm = TRUE)
    if (is.nan(aly)) aly <- NA_real_
  } else {
    aly <- NA_real_
  }
  
  list(
    TAM = tam, Trigger_Rate = trig_rate,
    N_True_Alarms_yr  = n_true_alarms_per_year,
    N_False_Alarms_yr = n_false_alarms_per_year,
    PPV = ppv, Sensitivity = sens,
    Mean_Lead_Time = mean_lead, WP = wp, ALY = aly,
    Mean_Lead_Time_conditional = mean_lead_conditional,
    WP_conditional = wp_conditional,
    ALY_conditional = aly_conditional,
    Total_Triggers = total, True_Alarms = true_n, False_Alarms = false_n,
    n_Reactive = reactive_n, n_TrueActionable = truact_n,
    N_Years_with_A1_True = nrow(lead_sub)
  )
}

bootstrap_method_metrics <- function(method_name, ym, trig_aug, yearly_lead,
                                     evaluable_years, B = BOOT_N_CI) {
  ym_method   <- ym         %>% dplyr::filter(Method == method_name, Year %in% evaluable_years)
  trig_method <- trig_aug   %>% dplyr::filter(Method == method_name, Year %in% evaluable_years)
  lead_method <- yearly_lead%>% dplyr::filter(Method == method_name, Year %in% evaluable_years)
  
  metric_names <- c("TAM", "N_True_yr", "N_False_yr", "PPV", "Sens",
                    "MeanLead", "WP", "ALY")
  
  uy <- intersect(unique(c(ym_method$Year, trig_method$Year)), evaluable_years)
  if (length(uy) < 3) {
    return(matrix(NA_real_, nrow = 2, ncol = length(metric_names),
                  dimnames = list(c("lo","hi"), metric_names)))
  }
  
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(metric_names),
                     dimnames = list(NULL, metric_names))
  
  for (b in seq_len(B)) {
    sy <- sample(uy, length(uy), replace = TRUE)
    
    yearly_vals <- vapply(sy, function(y) {
      v <- ym_method$TAM_Accum[ym_method$Year == y]
      if (length(v) == 0) NA_real_ else v[1]
    }, numeric(1))
    tam_b <- mean(yearly_vals, na.rm = TRUE)
    
    yearly_lead_vals <- vapply(sy, function(y) {
      v <- lead_method$Lead_Time_Yr[lead_method$Year == y]
      if (length(v) == 0) 0 else v[1]
    }, numeric(1))
    ml <- if (length(sy) > 0) mean(yearly_lead_vals) else NA_real_
    
    boot_idx <- unlist(lapply(sy, function(y) which(trig_method$Year == y)))
    if (length(boot_idx) == 0) {
      total <- 0; true_n <- 0; false_n <- 0
      truact_n <- 0
      wp_val <- 0
    } else {
      boot_sub   <- trig_method[boot_idx, , drop = FALSE]
      total      <- nrow(boot_sub)
      true_n     <- sum(boot_sub$IsTrue, na.rm = TRUE)
      false_n    <- total - true_n
      truact_mask <- boot_sub$IsTrue & boot_sub$Compartment == "Actionable"
      truact_n   <- sum(truact_mask, na.rm = TRUE)
      yearly_wp <- vapply(sy, function(y) {
        idx <- which(boot_sub$Year == y &
                       boot_sub$IsTrue &
                       boot_sub$Compartment == "Actionable")
        if (length(idx) > 0) mean(boot_sub$Lead_Time[idx], na.rm = TRUE) else 0
      }, numeric(1))
      wp_val <- mean(yearly_wp, na.rm = TRUE)
      if (is.nan(wp_val)) wp_val <- NA_real_
    }
    
    n_eval_b <- length(sy)
    
    sens_count <- sum(vapply(sy, function(y) {
      idx <- which(trig_method$Year == y)
      if (length(idx) == 0) return(FALSE)
      any(trig_method$IsTrue[idx] &
            trig_method$Compartment[idx] == "Actionable",
          na.rm = TRUE)
    }, logical(1)))
    
    if (length(boot_idx) == 0) {
      aly_val <- 0
    } else {
      yearly_aly <- vapply(sy, function(y) {
        yr_idx <- which(boot_sub$Year == y)
        if (length(yr_idx) == 0) return(0)
        yr_trues  <- sum(boot_sub$IsTrue[yr_idx], na.rm = TRUE)
        yr_truact <- sum(boot_sub$IsTrue[yr_idx] &
                           boot_sub$Compartment[yr_idx] == "Actionable",
                         na.rm = TRUE)
        if (yr_trues > 0L) yr_truact / yr_trues else 0
      }, numeric(1))
      aly_val <- mean(yearly_aly, na.rm = TRUE)
      if (is.nan(aly_val)) aly_val <- NA_real_
    }
    
    boot_mat[b, "TAM"]        <- tam_b
    boot_mat[b, "N_True_yr"]  <- if (n_eval_b > 0) true_n  / n_eval_b else NA_real_
    boot_mat[b, "N_False_yr"] <- if (n_eval_b > 0) false_n / n_eval_b else NA_real_
    boot_mat[b, "PPV"]        <- if (total > 0)  true_n / total      else NA_real_
    boot_mat[b, "Sens"]       <- if (n_eval_b > 0) sens_count / n_eval_b else NA_real_
    boot_mat[b, "MeanLead"]   <- ml
    boot_mat[b, "WP"]         <- wp_val
    boot_mat[b, "ALY"]        <- aly_val
  }
  
  apply(boot_mat, 2, function(v) safe_quantile(v, c(0.025, 0.975)))
}

aggregate_methods <- function(ym, trig_aug, yearly_lead, evaluable_years, B = BOOT_N_CI) {
  rows <- list()
  for (m in names(surge_defs)) {
    pt <- compute_method_point(m, ym, trig_aug, yearly_lead, evaluable_years)
    ci <- bootstrap_method_metrics(m, ym, trig_aug, yearly_lead, evaluable_years, B)
    rows[[m]] <- data.frame(
      Method               = m,
      TAM                  = pt$TAM,
      TAM_Lo  = ci[1, "TAM"],         TAM_Hi  = ci[2, "TAM"],
      N_True_Alarms_yr     = pt$N_True_Alarms_yr,
      NTrue_Lo  = ci[1, "N_True_yr"], NTrue_Hi  = ci[2, "N_True_yr"],
      PPV                  = pt$PPV,
      PPV_Lo  = ci[1, "PPV"],         PPV_Hi  = ci[2, "PPV"],
      Sensitivity          = pt$Sensitivity,
      Sens_Lo = ci[1, "Sens"],        Sens_Hi = ci[2, "Sens"],
      Mean_Lead_Time       = pt$Mean_Lead_Time,
      MeanLead_Lo = ci[1, "MeanLead"], MeanLead_Hi = ci[2, "MeanLead"],
      WP                   = pt$WP,
      WP_Lo  = ci[1, "WP"],            WP_Hi  = ci[2, "WP"],
      ALY                  = pt$ALY,
      ALY_Lo = ci[1, "ALY"],           ALY_Hi = ci[2, "ALY"],
      Mean_Lead_Time_conditional = pt$Mean_Lead_Time_conditional,
      WP_conditional             = pt$WP_conditional,
      ALY_conditional            = pt$ALY_conditional,
      N_False_Alarms_yr    = pt$N_False_Alarms_yr,
      NFalse_Lo = ci[1, "N_False_yr"], NFalse_Hi = ci[2, "N_False_yr"],
      Trigger_Rate         = pt$Trigger_Rate,
      Total_Triggers       = pt$Total_Triggers,
      True_Alarms          = pt$True_Alarms,
      False_Alarms         = pt$False_Alarms,
      n_Reactive           = pt$n_Reactive,
      n_TrueActionable     = pt$n_TrueActionable,
      N_Years_with_A1_True = pt$N_Years_with_A1_True,
      N_Years_Evaluable    = length(evaluable_years),
      stringsAsFactors     = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

# Primary specification: EVALUABLE_YEARS
summary_primary <- aggregate_methods(yearly_metrics, trigger_detail_aug,
                                     yearly_lead_data, EVALUABLE_YEARS, B = BOOT_N_CI)

# Sensitivity 1: + 2020/2021
sens1_years <- sort(c(EVALUABLE_YEARS, 2020L, 2021L))
trigger_detail_2021 <- build_trigger_detail(df, surge_defs, sens1_years)
yearly_lead_2021    <- compute_yearly_lead_data(trigger_detail_2021, sens1_years)
summary_with_covid  <- aggregate_methods(yearly_metrics, trigger_detail_2021,
                                         yearly_lead_2021, sens1_years, B = BOOT_N_CI)

# Sensitivity 2: + 2025
sens2_years <- sort(c(EVALUABLE_YEARS, 2025L))
trigger_detail_2025 <- build_trigger_detail(df, surge_defs, sens2_years)
yearly_lead_2025    <- compute_yearly_lead_data(trigger_detail_2025, sens2_years)
summary_with_2025   <- aggregate_methods(yearly_metrics, trigger_detail_2025,
                                         yearly_lead_2025, sens2_years, B = BOOT_N_CI)

# Sensitivity 3: all years
trigger_detail_all  <- build_trigger_detail(df, surge_defs, target_years)
yearly_lead_all     <- compute_yearly_lead_data(trigger_detail_all, target_years)
summary_all_years   <- aggregate_methods(yearly_metrics, trigger_detail_all,
                                         yearly_lead_all, target_years, B = BOOT_N_CI)

summary_df <- summary_primary %>%
  dplyr::left_join(method_type_map, by = "Method") %>%
  dplyr::mutate(
    Method_Type = factor(Method_Type, levels = type_order),
    Method      = factor(Method, levels = rev(method_group_order))
  )


# -----------------------------------------------------------------------------
# B.17 PANEL A DATA - DOMINANCE MATRIX
# -----------------------------------------------------------------------------
DOMINANCE_THRESHOLD <- 0.75

build_dominance_slice <- function(df_in, metric_key, val_col, lo_col, hi_col,
                                  direction, fmt_fn, n_col = NULL) {
  out <- data.frame(
    Method      = df_in$Method,
    Method_Type = df_in$Method_Type,
    Metric      = metric_key,
    Value       = df_in[[val_col]],
    Lo          = df_in[[lo_col]],
    Hi          = df_in[[hi_col]],
    Direction   = direction,
    stringsAsFactors = FALSE
  )
  out$Value_label <- fmt_fn(out$Value)
  out$CI_label    <- paste0("[", fmt_fn(out$Lo), ", ", fmt_fn(out$Hi), "]")
  out$N_label <- if (!is.null(n_col)) paste0("n=", df_in[[n_col]]) else NA_character_
  out
}

fmt_int   <- function(v) ifelse(is.na(v), "NA", format(round(v),     big.mark = ",", trim = TRUE))
fmt_int1  <- function(v) ifelse(is.na(v), "NA", format(round(v, 1),  nsmall = 1, trim = TRUE))
fmt_pct   <- function(v) ifelse(is.na(v), "NA", paste0(format(round(v * 100), trim = TRUE), "%"))
fmt_lead  <- function(v) ifelse(is.na(v), "NA", format(round(v, 1),  nsmall = 1, trim = TRUE))

heat_df <- dplyr::bind_rows(
  build_dominance_slice(summary_df, "True-Alarm Magnitude",
                        "TAM", "TAM_Lo", "TAM_Hi", "high", fmt_int),
  build_dominance_slice(summary_df, "Number of True Alarms",
                        "N_True_Alarms_yr", "NTrue_Lo", "NTrue_Hi", "high", fmt_int1,
                        n_col = "True_Alarms"),
  build_dominance_slice(summary_df, "Positive Predictive Value",
                        "PPV", "PPV_Lo", "PPV_Hi", "high", fmt_pct, n_col = "Total_Triggers"),
  build_dominance_slice(summary_df, "Sensitivity",
                        "Sensitivity", "Sens_Lo", "Sens_Hi", "high", fmt_pct,
                        n_col = "N_Years_Evaluable"),
  build_dominance_slice(summary_df, "Mean Lead Time",
                        "Mean_Lead_Time", "MeanLead_Lo", "MeanLead_Hi", "high", fmt_lead,
                        n_col = "N_Years_with_A1_True"),
  build_dominance_slice(summary_df, "Warning Persistence",
                        "WP", "WP_Lo", "WP_Hi", "high", fmt_lead, n_col = "n_TrueActionable"),
  build_dominance_slice(summary_df, "Actionable Lead-Time Yield",
                        "ALY", "ALY_Lo", "ALY_Hi", "high", fmt_pct, n_col = "True_Alarms"),
  build_dominance_slice(summary_df, "Number of False Alarms",
                        "N_False_Alarms_yr", "NFalse_Lo", "NFalse_Hi", "low", fmt_int1,
                        n_col = "False_Alarms")
)

metric_order <- c(
  "True-Alarm Magnitude",
  "Number of True Alarms",
  "Positive Predictive Value",
  "Sensitivity",
  "Mean Lead Time",
  "Warning Persistence",
  "Actionable Lead-Time Yield",
  "Number of False Alarms"
)
heat_df$Metric <- factor(heat_df$Metric, levels = metric_order)

heat_df$fill_norm <- NA_real_
for (mk in metric_order) {
  idx <- which(heat_df$Metric == mk)
  if (length(idx) == 0) next
  raw <- rescale_01_safe(heat_df$Value[idx])
  direction <- unique(heat_df$Direction[idx])[1]
  heat_df$fill_norm[idx] <- if (direction == "low") 1 - raw else raw
}

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
heat_df$cell_fill <- vapply(heat_df$fill_norm, mix_blue_intensity, character(1))

dominance_count <- heat_df %>%
  dplyr::group_by(Method) %>%
  dplyr::summarise(
    N_Dominant = sum(fill_norm >= DOMINANCE_THRESHOLD, na.rm = TRUE),
    N_Total    = sum(!is.na(fill_norm)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Dominance_Label = paste0(N_Dominant, "/", N_Total)
  )

separator_x_left  <- 4 + 0.5
separator_x_right <- 7 + 0.5

y_levels <- rev(method_group_order)
paradigms_in_y_order <- vapply(
  y_levels,
  function(m) {
    p <- method_type_map$Method_Type[method_type_map$Method == m]
    if (length(p) == 0) NA_character_ else as.character(p[1])
  },
  character(1)
)
y_breaks_between <- which(
  paradigms_in_y_order[-length(paradigms_in_y_order)] != paradigms_in_y_order[-1]
)
h_separator_y <- y_breaks_between + 0.5


# -----------------------------------------------------------------------------
# B.18 THEMES
# -----------------------------------------------------------------------------
theme_dashboard <- function(base_size = 9, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      panel.background   = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      plot.background    = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major   = ggplot2::element_line(colour = "#CCCCCC", linetype = "44", linewidth = 0.25),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(colour = "grey75", fill = NA, linewidth = 0.35),
      axis.title         = ggplot2::element_text(face = "bold", size = base_size - 0.2),
      axis.text          = ggplot2::element_text(colour = "grey25", size = base_size - 1.1),
      axis.ticks         = ggplot2::element_line(colour = "grey60", linewidth = 0.25),
      plot.title         = ggplot2::element_text(face = "bold", size = base_size + 0.2,
                                                 hjust = 0, lineheight = 0.96, margin = ggplot2::margin(t = 2, b = 8)),
      plot.subtitle      = ggplot2::element_blank(),
      legend.position    = "bottom",
      legend.direction   = "horizontal",
      legend.background  = ggplot2::element_rect(fill = "white", colour = "grey80", linewidth = 0.25),
      legend.key         = ggplot2::element_rect(fill = NA),
      legend.title       = ggplot2::element_text(face = "bold", size = base_size - 0.7),
      legend.text        = ggplot2::element_text(size = base_size - 1.2),
      legend.margin      = ggplot2::margin(4, 5, 4, 5),
      legend.spacing.x   = grid::unit(0.18, "cm"),
      plot.margin        = ggplot2::margin(10, 10, 10, 10)
    )
}

theme_heat_dashboard <- function(base_size = 8.9, base_family = "sans") {
  theme_dashboard(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black", size = base_size - 1.2,
                                          lineheight = 0.95, margin = ggplot2::margin(b = 6)),
      axis.text.y = ggplot2::element_text(color = "black", size = base_size - 1.1, lineheight = 0.92),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 6.9, lineheight = 0.92),
      legend.title = ggplot2::element_text(size = 7.8)
    )
}


# -----------------------------------------------------------------------------
# B.19 PANEL A - DOMINANCE MATRIX PLOT
# -----------------------------------------------------------------------------
dom_count_df <- dominance_count %>%
  dplyr::mutate(Method = factor(Method, levels = rev(method_group_order)))

DOM_COUNT_X <- length(metric_order) + 1L

metric_display_labels <- c(
  "True-Alarm Magnitude"        = "True-Alarm\nMagnitude",
  "Number of True Alarms"       = "Number of\nTrue Alarms",
  "Positive Predictive Value"   = "Positive\nPredictive Value",
  "Sensitivity"                 = "Sensitivity",
  "Mean Lead Time"              = "Mean Lead\nTime",
  "Warning Persistence"         = "Warning\nPersistence",
  "Actionable Lead-Time Yield"  = "Actionable\nLead-Time Yield",
  "Number of False Alarms"      = "Number of\nFalse Alarms"
)

panel_a <- ggplot2::ggplot(heat_df) +
  ggplot2::geom_tile(
    ggplot2::aes(x = Metric, y = Method),
    fill = heat_df$cell_fill,
    colour = "white", linewidth = 0.6, show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = dom_count_df,
    ggplot2::aes(x = DOM_COUNT_X, y = Method, label = Dominance_Label),
    size = 2.9, fontface = "bold", family = base_family_global,
    colour = "#0B2447", inherit.aes = FALSE, show.legend = FALSE
  ) +
  ggplot2::geom_vline(xintercept = separator_x_left, linetype = "dashed",
                      linewidth = 0.40, colour = "grey45") +
  ggplot2::geom_vline(xintercept = separator_x_right, linetype = "dashed",
                      linewidth = 0.40, colour = "grey45") +
  ggplot2::geom_vline(xintercept = DOM_COUNT_X - 0.5, linetype = "solid",
                      linewidth = 0.50, colour = "grey25") +
  {
    if (length(h_separator_y) > 0)
      ggplot2::geom_hline(yintercept = h_separator_y, colour = "grey60", linewidth = 0.35)
    else NULL
  } +
  ggplot2::annotate(
    "text", x = DOM_COUNT_X, y = length(method_group_order) + 0.55,
    label = paste0("Dominance\n(score \u2265 ",
                   sprintf("%.2f", DOMINANCE_THRESHOLD), ")"),
    size = 2.4, family = base_family_global, fontface = "bold",
    colour = "#0B2447", lineheight = 0.92
  ) +
  ggplot2::scale_x_discrete(
    limits = c(metric_order, "__DOM_COUNT__"),
    labels = c(unname(metric_display_labels[metric_order]), ""),
    position = "top",
    expand = ggplot2::expansion(add = c(0.04, 0.50))
  ) +
  ggplot2::scale_y_discrete(
    limits = rev(method_group_order),
    labels = method_label_map,
    expand = ggplot2::expansion(add = c(0.25, 0.30))
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::labs(
    title = "a    Dominance matrix across 8 performance metrics",
    x = NULL, y = NULL
  ) +
  theme_heat_dashboard(base_size = 8.9, base_family = base_family_global)

legend_score_seq <- seq(0, 1, length.out = 100L)
dominance_legend_df <- data.frame(
  x = legend_score_seq,
  y = 1L,
  fill_col = vapply(legend_score_seq, mix_blue_intensity, character(1))
)
dominance_legend <- ggplot2::ggplot(dominance_legend_df) +
  ggplot2::geom_tile(ggplot2::aes(x = x, y = y),
                     fill = dominance_legend_df$fill_col,
                     width = 1 / nrow(dominance_legend_df),
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
    title = paste0("Dominance score (min-max normalized within metric; ",
                   "darker = stronger dominance)"),
    x = NULL, y = NULL
  ) +
  ggplot2::theme_void(base_family = base_family_global) +
  ggplot2::theme(
    plot.title  = ggplot2::element_text(face = "bold", size = 7.5, hjust = 0,
                                        margin = ggplot2::margin(b = 4)),
    axis.text.x = ggplot2::element_text(size = 6.8, lineheight = 0.95,
                                        margin = ggplot2::margin(t = 2)),
    plot.margin = ggplot2::margin(2, 4, 2, 4)
  )


# -----------------------------------------------------------------------------
# B.20 PANEL B - TAM vs MEAN LEAD TIME SCATTER
# -----------------------------------------------------------------------------
scatter_df <- summary_df %>%
  dplyr::mutate(
    Method_chr = as.character(Method),
    Method_Type_f = factor(Method_Type, levels = type_order)
  )

x_limits <- range(c(scatter_df$TAM_Lo, scatter_df$TAM_Hi), na.rm = TRUE)
y_limits <- range(
  c(scatter_df$MeanLead_Lo, scatter_df$MeanLead_Hi, scatter_df$Mean_Lead_Time),
  na.rm = TRUE
)

x_pad <- diff(x_limits)
if (!is.finite(x_pad) || x_pad <= 0) x_pad <- max(1, abs(x_limits[1]) * 0.20, 25)

y_pad <- diff(y_limits)
if (!is.finite(y_pad) || y_pad <= 0) y_pad <- 1

x_limits <- c(x_limits[1] - 0.08 * x_pad, x_limits[2] + 0.12 * x_pad)
y_limits <- c(max(0, y_limits[1] - 0.20 * y_pad),
              min(10, y_limits[2] + 0.20 * y_pad))

tam_breaks <- pretty(scatter_df$TAM[is.finite(scatter_df$TAM)], n = 3)
tam_breaks <- tam_breaks[tam_breaks > 0]
if (length(tam_breaks) < 2) tam_breaks <- c(1000, 3000, 5000)

panel_b <- ggplot2::ggplot(
  scatter_df,
  ggplot2::aes(x = TAM, y = Mean_Lead_Time)
) +
  ggplot2::annotate("rect",
                    xmin = -Inf, xmax = Inf, ymin = 4, ymax = 8,
                    fill = "#cfe2f3", alpha = 0.30) +
  ggplot2::geom_point(
    ggplot2::aes(fill = Method_Type_f, size = TAM),
    shape = 21, color = "grey20", stroke = 0.35, alpha = 0.92, na.rm = TRUE
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = Method_chr),
    family = base_family_global, size = 2.5, color = "black",
    box.padding = 0.42, point.padding = 0.32,
    segment.color = "grey60", segment.size = 0.28,
    min.segment.length = 0, force = 1.5, force_pull = 0.30,
    max.overlaps = Inf, seed = 12345, direction = "both"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = TAM_Lo, xmax = TAM_Hi),
    height = 0, linewidth = 0.28, colour = "grey60", alpha = 0.55,
    show.legend = FALSE, na.rm = TRUE
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = MeanLead_Lo, ymax = MeanLead_Hi),
    width = 0, linewidth = 0.28, colour = "grey60", alpha = 0.55,
    show.legend = FALSE, na.rm = TRUE
  ) +
  ggplot2::scale_fill_manual(
    values = TYPE_COLORS, breaks = type_order, labels = type_order,
    drop = FALSE, name = "Paradigm"
  ) +
  ggplot2::scale_size_area(
    name = "True-Alarm Magnitude\n(cases/year)",
    breaks = tam_breaks,
    labels = scales::comma,
    max_size = 12
  ) +
  ggplot2::scale_x_continuous(
    limits = x_limits, breaks = pretty(x_limits, n = 6),
    labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_y_continuous(
    limits = y_limits, breaks = pretty(y_limits, n = 6),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "b    Tradeoff between True-Alarm Magnitude and Mean Lead Time",
    x = "True-Alarm Magnitude (cases/year captured during T = Actionable Window \u222A Epidemic Burden)",
    y = "Mean Lead Time (weeks before annual peak; bounded [4, 8])"
  ) +
  theme_dashboard(base_size = 8.8, base_family = base_family_global) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    size = ggplot2::guide_legend(
      order = 1, title.position = "left", label.position = "right",
      direction = "horizontal", nrow = 1, byrow = TRUE,
      override.aes = list(fill = NA, color = "grey20", shape = 21, alpha = 1, stroke = 0.6),
      keywidth = grid::unit(0.55, "cm"), keyheight = grid::unit(0.40, "cm")
    ),
    fill = ggplot2::guide_legend(
      order = 2, nrow = 1, byrow = TRUE,
      override.aes = list(shape = 21, size = 4.5, color = "grey20", alpha = 1),
      keywidth = grid::unit(0.90, "cm"), keyheight = grid::unit(0.36, "cm")
    )
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "center",
    legend.text = ggplot2::element_text(size = 6.9, lineheight = 0.92),
    legend.title = ggplot2::element_text(size = 7.8),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_rect(fill = NA, colour = NA),
    legend.margin = ggplot2::margin(2, 2, 2, 2),
    legend.spacing.y = grid::unit(0.18, "cm"),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0)
  )


# -----------------------------------------------------------------------------
# B.21 PRINT AND SAVE FIGURE 2
# -----------------------------------------------------------------------------
panel_a_with_legend <- panel_a / dominance_legend +
  patchwork::plot_layout(heights = c(10, 1))

print(panel_a_with_legend)
print(panel_b)

save_plot_pair("Figure2_PanelA_DominanceMatrix.pdf", panel_a_with_legend,
               panelA_width_in, panelA_height_in,
               dpi = 600, device = safe_pdf_device())
save_plot_pair("Figure2_PanelA_DominanceMatrix.png", panel_a_with_legend,
               panelA_width_in, panelA_height_in,
               dpi = 600, device = "png")

save_plot_pair("Figure2_PanelB_Scatter.pdf", panel_b, panelB_width_in, panelB_height_in,
               dpi = 600, device = safe_pdf_device())
save_plot_pair("Figure2_PanelB_Scatter.png", panel_b, panelB_width_in, panelB_height_in,
               dpi = 600, device = "png")


# -----------------------------------------------------------------------------
# B.22 FIGURE 2 CSV OUTPUTS
# -----------------------------------------------------------------------------
utils::write.csv(thresholds_full,
                 file.path(OUT_DIR, "Figure2_rolling_weekly_thresholds.csv"), row.names = FALSE)
utils::write.csv(yearly_metrics,
                 file.path(OUT_DIR, "Figure2_yearly_metrics.csv"), row.names = FALSE)
utils::write.csv(trigger_detail_aug,
                 file.path(OUT_DIR, "Figure2_trigger_detail_with_compartments.csv"), row.names = FALSE)
utils::write.csv(yearly_lead_data,
                 file.path(OUT_DIR, "Figure2_yearly_lead_data.csv"), row.names = FALSE)

utils::write.csv(summary_primary,
                 file.path(OUT_DIR, "Figure2_method_summary_primary.csv"), row.names = FALSE)
utils::write.csv(summary_with_covid,
                 file.path(OUT_DIR, "Figure2_method_summary_sensitivity_add_2020_2021.csv"), row.names = FALSE)
utils::write.csv(summary_with_2025,
                 file.path(OUT_DIR, "Figure2_method_summary_sensitivity_add_2025.csv"), row.names = FALSE)
utils::write.csv(summary_all_years,
                 file.path(OUT_DIR, "Figure2_method_summary_sensitivity_all_years.csv"), row.names = FALSE)

tam_diag <- yearly_metrics %>%
  dplyr::select(Year, Method, TAM_Accum) %>%
  tidyr::pivot_wider(names_from = Year, values_from = TAM_Accum)
utils::write.csv(tam_diag,
                 file.path(OUT_DIR, "Figure2_TAM_per_year.csv"), row.names = FALSE)

leadtopeak_diag <- yearly_lead_data %>%
  dplyr::select(Year, Method, Lead_Time_Yr) %>%
  tidyr::pivot_wider(names_from = Year, values_from = Lead_Time_Yr)
utils::write.csv(leadtopeak_diag,
                 file.path(OUT_DIR, "Figure2_mean_lead_time_per_year.csv"), row.names = FALSE)

utils::write.csv(heat_df,
                 file.path(OUT_DIR, "Figure2_dominance_matrix_long.csv"), row.names = FALSE)
utils::write.csv(dominance_count,
                 file.path(OUT_DIR, "Figure2_dominance_counts.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# B.23 PRINT FIGURE 2 TABLES
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("FIGURE 2  COMPARATIVE DASHBOARD SUMMARY\n")
cat("Primary specification: EVALUABLE_YEARS = ",
    paste(EVALUABLE_YEARS, collapse = ", "), " (n=", length(EVALUABLE_YEARS), ")\n", sep = "")
cat("Excluded years: ", paste(EXCLUDED_YEARS, collapse = ", "), "\n", sep = "")
cat("============================================================\n")

cat("\nSTA/LTA parameters:\n")
cat("  Continuous TA threshold = ", ETA_ON_CLASSIC, "\n", sep = "")
cat("  Constant TA   eta_ON    = ", ETA_ON,         "\n", sep = "")
cat("  Constant TA   eta_OFF   = ", ETA_OFF,        "\n", sep = "")

cat("\nTwo-anchor framework parameters:\n")
cat("  A1 lead window         = [", A1_LEAD_MIN, ", ", A1_LEAD_MAX, "] weeks\n", sep = "")
cat("  A2 burden fraction     = ", A2_BURDEN_FRAC, "\n", sep = "")
cat("  TAM integration window = T = A1 union A2\n")
cat("  Compartments           : Actionable (4-8 wks); Reactive (other true alarms)\n")

safe_df_print(
  summary_df %>% dplyr::arrange(Method_Type, Method),
  title = "Primary specification - 8 performance metrics:",
  round_cols = c(
    "TAM","TAM_Lo","TAM_Hi",
    "N_True_Alarms_yr","NTrue_Lo","NTrue_Hi",
    "PPV","PPV_Lo","PPV_Hi",
    "Sensitivity","Sens_Lo","Sens_Hi",
    "Mean_Lead_Time","MeanLead_Lo","MeanLead_Hi",
    "WP","WP_Lo","WP_Hi",
    "ALY","ALY_Lo","ALY_Hi",
    "N_False_Alarms_yr","NFalse_Lo","NFalse_Hi",
    "Trigger_Rate"
  ),
  digits = 3
)

cat("\n--- Sensitivity 1: + 2020/2021 (COVID-19 years) ---\n")
safe_df_print(
  summary_with_covid %>% dplyr::arrange(Method),
  round_cols = c("TAM","N_True_Alarms_yr","PPV","Sensitivity","Mean_Lead_Time",
                 "WP","ALY","N_False_Alarms_yr"),
  digits = 3
)

cat("\n--- Sensitivity 2: + 2025 (out-of-distribution year) ---\n")
safe_df_print(
  summary_with_2025 %>% dplyr::arrange(Method),
  round_cols = c("TAM","N_True_Alarms_yr","PPV","Sensitivity","Mean_Lead_Time",
                 "WP","ALY","N_False_Alarms_yr"),
  digits = 3
)

cat("\n--- Sensitivity 3: All years (2013-2025) ---\n")
safe_df_print(
  summary_all_years %>% dplyr::arrange(Method),
  round_cols = c("TAM","N_True_Alarms_yr","PPV","Sensitivity","Mean_Lead_Time",
                 "WP","ALY","N_False_Alarms_yr"),
  digits = 3
)

cat("\n--- Dominance counts (X / 8 metrics with normalized score >= ",
    sprintf("%.2f", DOMINANCE_THRESHOLD), ") ---\n", sep = "")
safe_df_print(dominance_count, digits = 3)


# -----------------------------------------------------------------------------
# B.24 FIGURE 2 LEGEND DRAFTS
# -----------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("FIGURE 2 - DRAFT FIGURE LEGENDS\n")
cat("=================================================================\n")

cat("\nDraft legend for Panel A (Dominance Matrix):\n")
cat(paste0(
  "Fig. 2a | Dominance matrix of 11 dengue outbreak-detection methods across 8 performance ",
  "metrics in Quezon City, Philippines, 2013-2024 (excluding 2020, 2021, and 2025; n = ",
  length(EVALUABLE_YEARS), " evaluable years). Each cell encodes the detector's normalized ",
  "dominance score on the corresponding metric, computed by min-max normalization within the ",
  "metric across all detectors: for higher-is-better metrics, score = (value - min) / (max - min); ",
  "for the lower-is-better Number of False Alarms metric, score = (max - value) / (max - min). ",
  "A score of 1.00 indicates the field leader on that metric; 0.00 indicates the field laggard. ",
  "Cell shading uses a single-color (blue) ramp from white (score 0; not dominant) through ",
  "mid-blue (0.50; moderate) to dark navy (1.00; total sweep). The 8 metrics are organized into ",
  "three categories (vertical dashed separators): Epidemic Burden & Alarm Accuracy (True-Alarm ",
  "Magnitude, Number of True Alarms, Positive Predictive Value, Sensitivity); Early Warning ",
  "Timeliness (Mean Lead Time, Warning Persistence, Actionable Lead-Time Yield); False Alarms ",
  "(Number of False Alarms). Sensitivity is computed as the proportion of evaluable seasons ",
  "with at least one True Alarm in the Actionable Window (A1). The three Early Warning ",
  "Timeliness metrics share the same denominator per detector (count of evaluable years): ",
  "years where the metric is not computable due to no qualifying triggers contribute zero ",
  "rather than being excluded. Conditional-on-firing diagnostic versions ",
  "(Mean_Lead_Time_conditional, WP_conditional, ALY_conditional) are preserved in the ",
  "companion CSV. The right-most column reports the per-detector dominance count: the number ",
  "of metrics on which the detector's normalized score is >= ",
  sprintf("%.2f", DOMINANCE_THRESHOLD),
  " (the dominance threshold), out of 8 total. Raw metric values and 95% year-cluster ",
  "bootstrap confidence intervals (", BOOT_N_CI, " replicates) are reported in the companion CSV.\n"
))

cat("\nDraft legend for Panel B (TAM Trade-off Scatter):\n")
cat(paste0(
  "Fig. 2b | Tradeoff between True-Alarm Magnitude (TAM) and Mean Lead Time among 11 dengue ",
  "outbreak-detection methods in Quezon City, Philippines, 2013-2024 (excluding 2020, 2021, ",
  "and 2025). The x-axis shows TAM as the mean across years of weekly case counts summed at ",
  "alarm-on weeks within T = A1 union A2 (the anchor-consistent epidemic-burden volume captured ",
  "by true alarms). The y-axis shows Mean Lead Time in weeks under the same-denominator scheme: ",
  "per-year first-A1-true-alarm lead time (zero for years where the detector did not produce ",
  "an A1 true alarm), averaged across all evaluable years. The light blue band marks the ",
  "Actionable Window (4 to 8 weeks before peak). POINT SIZE = TAM, with area scaling linearly ",
  "with TAM, so larger circles indicate detectors that capture more epidemic burden through ",
  "their true alarms. Horizontal and vertical error bars denote year-cluster bootstrap 95% ",
  "confidence intervals (", BOOT_N_CI, " replicates). Labels identify individual methods. ",
  "Colors denote outbreak-detection paradigms. Background gridlines are intentionally omitted ",
  "to reduce visual clutter.\n"
))

cat("\nDraft combined figure legend:\n")
cat(paste0(
  "Fig. 2 | Operational comparison of 11 dengue outbreak-detection methods in Quezon City, ",
  "Philippines, 2013-2024.\n\n",
  "a, Dominance matrix across 8 performance metrics organized in three categories (Epidemic ",
  "Burden & Alarm Accuracy: True-Alarm Magnitude, Number of True Alarms, Positive Predictive ",
  "Value, Sensitivity; Early Warning Timeliness: Mean Lead Time, Warning Persistence, ",
  "Actionable Lead-Time Yield; False Alarms: Number of False Alarms). Cell shading encodes ",
  "the normalized dominance score (min-max within each metric across detectors) on a single-",
  "color blue ramp from white (score 0) to dark navy (score 1). Right-most column reports the ",
  "dominance count: number of metrics where the detector's score is >= ",
  sprintf("%.2f", DOMINANCE_THRESHOLD),
  ", out of 8. Two anchors are used: the Actionable Window (A1, peak-8 to peak-4 weeks) and ",
  "the Epidemic Burden block (A2, smallest contiguous block containing the peak whose ",
  "cumulative cases sum to >= 70% of the annual total). A trigger is a True Alarm iff it ",
  "falls in A1 OR A2. Years 2020, 2021, and 2025 are excluded.\n\n",
  "b, Tradeoff between True-Alarm Magnitude (x-axis; cases/year captured by true alarms) and ",
  "Mean Lead Time (y-axis; weeks of warning per year, same-denominator). Point size reflects ",
  "True-Alarm Magnitude (area scales linearly with TAM). Horizontal and vertical error bars ",
  "denote year-cluster bootstrap 95% confidence intervals (", BOOT_N_CI, " replicates). Light ",
  "blue band marks the Actionable Window (4 to 8 weeks before peak). Colors denote outbreak-",
  "detection paradigms.\n"
))

cat("\nFigure 2 outputs saved to:\n  ", OUT_DIR, "\n", sep = "")


# =============================================================================
# END OF PART B (Figure 2)
# =============================================================================

# =============================================================================
# PART C - FIGURE 3: PER-YEAR DETECTION TIMING + HEAD-TO-HEAD COMPARISONS
# =============================================================================

cat("\n-------------------------------------------------------------\n")
cat("PART C  Figure 3  Per-year detection timing + head-to-head\n")
cat("-------------------------------------------------------------\n")


# -----------------------------------------------------------------------------
# C.1 SETTINGS
# -----------------------------------------------------------------------------
fig_width_in_main    <- 12
fig_height_in_main   <- 18
fig_width_in_single  <- 12
fig_height_in_single <- 9
Y_MAX                <- 800

# Head-to-head bootstrap constants
HH_BOOT_N <- 1000L           # year-cluster bootstrap replicates
HH_ALPHA  <- 0.05            # nominal significance per pairwise comparison

# Color palette for True/False alarm classification (markers)
COL_TRUE_ALARM     <- "#2ca02c"
COL_TRUE_OUTLINE   <- "#1b6b1b"
COL_FALSE_ALARM    <- "#e6550d"
COL_FALSE_OUTLINE  <- "#a83706"
COL_EXCLUDED_FILL  <- "white"
COL_EXCLUDED_LINE  <- "grey50"

# Color palette for lead-time compartments
COL_REFRACTORY <- "#a1d99b"
COL_ACTIONABLE <- "#2ca02c"
COL_PRE_PEAK   <- "#006d2c"
COL_REACTIVE   <- "#e31a1c"

# Bar fill highlight colors at first true-alarm week
COL_BAR_OT_FIRSTTRUE <- "#2c3e50"
COL_BAR_DEFAULT      <- "grey70"

# Y-positions for lead/lag arrows (one row per detector)
ARROW_Y_TA <- Y_MAX * 0.94
ARROW_Y_OT <- Y_MAX * 0.78

SUB_ARROW_OFFSET_Y <- 0.030 * Y_MAX
LABEL_OFFSET_Y     <- 0.020 * Y_MAX


# -----------------------------------------------------------------------------
# C.2 DETECTOR COLUMN ALIASES
# -----------------------------------------------------------------------------
# Map Part B detector columns onto the names Figure 3's framework expects.
# All algorithms and parameter values are identical, so aliasing produces
# byte-identical numerical results.
df$trig_classic            <- df$surge_sta_lta
df$trig_vaezi              <- df$surge_sta_lta_vaezi
df$surge_outbreak_threshold <- df$surge_mean_2sd
df$outbreak_threshold      <- df$bl_outbreak

# R(t) ratio series for Continuous TA (used for the per-year curve overlay).
.dc_sta_classic <- zoo::rollmean(df$DC_QC,  3, fill = NA, align = "right")
.dc_lta_classic <- zoo::rollmean(df$DC_QC, 12, fill = NA, align = "right")
df$R_classic <- ifelse(
  !is.na(.dc_sta_classic) & !is.na(.dc_lta_classic) & .dc_lta_classic > 0,
  .dc_sta_classic / .dc_lta_classic, NA_real_
)


# -----------------------------------------------------------------------------
# C.3 PANEL THEME (per-year detection panels)
# -----------------------------------------------------------------------------
theme_dashboard_panel <- function(base_size = 8, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = "#D0D0D0", linetype = "44",
                                               linewidth = 0.22),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "grey75", fill = NA,
                                           linewidth = 0.35),
      axis.title = ggplot2::element_text(face = "bold", size = base_size - 0.2),
      axis.text  = ggplot2::element_text(colour = "grey20", size = base_size - 1.2),
      axis.ticks = ggplot2::element_line(colour = "grey60", linewidth = 0.25),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5,
                                         size = base_size + 0.5,
                                         margin = ggplot2::margin(b = 3)),
      plot.margin = ggplot2::margin(8, 18, 8, 18)
    )
}


# -----------------------------------------------------------------------------
# C.4 PER-DETECTOR EVALUATION
# -----------------------------------------------------------------------------
evaluate_detector <- function(df_in, detector_col, evaluable_years,
                              lead_min = A1_LEAD_MIN, lead_max = A1_LEAD_MAX,
                              burden_frac = A2_BURDEN_FRAC) {
  per_trigger_rows <- list(); per_year_rows <- list()
  for (yr in evaluable_years) {
    anchors <- compute_anchors_for_year(df_in, yr,
                                        lead_min = lead_min, lead_max = lead_max,
                                        burden_frac = burden_frac)
    df_y <- df_in %>% dplyr::filter(YR == yr) %>% dplyr::arrange(WN)
    trig_idx <- which(df_y[[detector_col]] == 1L)
    year_true_count <- 0L; year_false_count <- 0L
    year_first_a1_true_week <- NA_integer_
    if (length(trig_idx) > 0L) {
      for (k in trig_idx) {
        wk <- df_y$WN[k]
        cls <- classify_trigger(wk, anchors$A1_weeks, anchors$A2_weeks)
        # Route attribute (route classifier preserved as a diagnostic field)
        route <- if (!cls$is_true) "FalseAlarm"
        else if (cls$in_A1 && cls$in_A2) "Both"
        else if (cls$in_A1) "Route1_EarlyWarning"
        else "Route2_WithinEpidemic"
        per_trigger_rows[[length(per_trigger_rows) + 1L]] <- data.frame(
          Year = yr, Week = as.integer(wk), Detector = detector_col,
          IsTrue = cls$is_true, Route = route,
          InA1 = cls$in_A1, InA2 = cls$in_A2,
          stringsAsFactors = FALSE)
        if (cls$is_true) {
          year_true_count <- year_true_count + 1L
          if (cls$in_A1 && (is.na(year_first_a1_true_week) || wk < year_first_a1_true_week))
            year_first_a1_true_week <- as.integer(wk)
        } else { year_false_count <- year_false_count + 1L }
      }
    }
    lead_time <- if (!is.na(year_first_a1_true_week) && !is.na(anchors$peak_week))
      as.integer(anchors$peak_week - year_first_a1_true_week) else NA_integer_
    per_year_rows[[length(per_year_rows) + 1L]] <- data.frame(
      Year = yr, Detector = detector_col, Peak_Week = anchors$peak_week,
      A1_start = if (length(anchors$A1_weeks)) min(anchors$A1_weeks) else NA_integer_,
      A1_end   = if (length(anchors$A1_weeks)) max(anchors$A1_weeks) else NA_integer_,
      A2_start = if (length(anchors$A2_weeks)) min(anchors$A2_weeks) else NA_integer_,
      A2_end   = if (length(anchors$A2_weeks)) max(anchors$A2_weeks) else NA_integer_,
      n_triggers = as.integer(length(trig_idx)),
      n_True = year_true_count, n_False = year_false_count,
      First_A1_True_Week = year_first_a1_true_week,
      Lead_Time_Weeks = lead_time, stringsAsFactors = FALSE)
  }
  per_trigger_df <- if (length(per_trigger_rows)) dplyr::bind_rows(per_trigger_rows) else
    data.frame(Year=integer(),Week=integer(),Detector=character(),IsTrue=logical(),
               Route=character(),InA1=logical(),InA2=logical(),
               stringsAsFactors=FALSE)
  per_year_df <- dplyr::bind_rows(per_year_rows)
  total_triggers <- sum(per_year_df$n_triggers, na.rm = TRUE)
  total_true     <- sum(per_year_df$n_True,     na.rm = TRUE)
  ppv <- if (total_triggers > 0) total_true / total_triggers else NA_real_
  evaluable_n  <- length(evaluable_years)
  # Sensitivity: A1-restricted (years with at least one True alarm in A1).
  years_with_A1 <- sum(!is.na(per_year_df$First_A1_True_Week))
  sensitivity   <- if (evaluable_n > 0) years_with_A1 / evaluable_n else NA_real_
  years_with_T  <- sum(per_year_df$n_True > 0L)
  # Same-denominator MLT: zero-coerce non-firing years.
  lead_zerocoerced <- ifelse(is.na(per_year_df$Lead_Time_Weeks), 0,
                             per_year_df$Lead_Time_Weeks)
  mean_lead             <- mean(lead_zerocoerced)
  mean_lead_conditional <- if (any(!is.na(per_year_df$Lead_Time_Weeks)))
    mean(per_year_df$Lead_Time_Weeks, na.rm = TRUE)
  else NA_real_
  if (is.nan(mean_lead))             mean_lead             <- NA_real_
  if (is.nan(mean_lead_conditional)) mean_lead_conditional <- NA_real_
  summary_row <- data.frame(
    Detector = detector_col, Total_Triggers = total_triggers,
    True_Alarms = total_true, False_Alarms = total_triggers - total_true,
    PPV = round(ppv, 3),
    Years_with_A1_True = years_with_A1,
    Years_with_True_Any = years_with_T,
    Years_Evaluable = evaluable_n,
    Sensitivity = round(sensitivity, 3),
    Mean_Lead_Time_wks = round(mean_lead, 2),
    Mean_Lead_Time_wks_conditional = round(mean_lead_conditional, 2),
    stringsAsFactors = FALSE)
  list(per_trigger = per_trigger_df, per_year = per_year_df, summary = summary_row)
}


# -----------------------------------------------------------------------------
# C.5 LEAD-TIME COMPARTMENT HELPERS (Figure 3 framework)
# -----------------------------------------------------------------------------
compartment_color <- function(compartment) {
  if (is.na(compartment)) return(NA_character_)
  switch(compartment,
         "Actionable" = COL_ACTIONABLE,
         "Reactive"   = COL_REACTIVE,
         NA_character_)
}

augment_compartments <- function(trigger_df, per_year_df) {
  if (nrow(trigger_df) == 0) {
    trigger_df$Lead_Time   <- integer(0)
    trigger_df$Compartment <- character(0)
    return(trigger_df)
  }
  peak_lookup <- per_year_df %>%
    dplyr::select(Year, Detector, Peak_Week) %>% dplyr::distinct()
  trigger_df %>%
    dplyr::left_join(peak_lookup, by = c("Year", "Detector")) %>%
    dplyr::mutate(
      Lead_Time   = as.integer(Peak_Week - Week),
      Compartment = ifelse(
        !is.na(IsTrue) & IsTrue,
        vapply(Lead_Time, classify_compartment, character(1)),
        NA_character_
      )
    )
}

compute_compartment_metrics <- function(trigger_df_aug) {
  if (nrow(trigger_df_aug) == 0) {
    return(data.frame(
      Detector = character(0), Total_Triggers = integer(0),
      True_Alarms = integer(0), False_Alarms = integer(0),
      n_Actionable = integer(0), n_Reactive = integer(0),
      n_TrueActionable = integer(0),
      ALY = numeric(0),    ALY_conditional = numeric(0),
      WP_wks = numeric(0), WP_wks_conditional = numeric(0)
    ))
  }
  # Same-denominator headline ALY and WP: per-year zero-coerced values
  # averaged across evaluable years.
  per_year_metrics <- trigger_df_aug %>%
    dplyr::group_by(Detector, Year) %>%
    dplyr::summarise(
      yr_n_true   = sum(IsTrue, na.rm = TRUE),
      yr_truact_n = sum(IsTrue & Compartment == "Actionable", na.rm = TRUE),
      yr_aly      = ifelse(yr_n_true > 0, yr_truact_n / yr_n_true, 0),
      yr_wp_wks   = ifelse(yr_truact_n > 0,
                           mean(Lead_Time[IsTrue & Compartment == "Actionable"],
                                na.rm = TRUE),
                           0),
      .groups = "drop"
    )
  per_detector_headline <- per_year_metrics %>%
    dplyr::group_by(Detector) %>%
    dplyr::summarise(
      ALY_headline    = mean(yr_aly,    na.rm = TRUE),
      WP_wks_headline = mean(yr_wp_wks, na.rm = TRUE),
      .groups = "drop"
    )
  
  trigger_df_aug %>%
    dplyr::group_by(Detector) %>%
    dplyr::summarise(
      Total_Triggers   = dplyr::n(),
      True_Alarms      = sum(IsTrue),
      False_Alarms     = sum(!IsTrue),
      n_Actionable     = sum(Compartment == "Actionable", na.rm = TRUE),
      n_Reactive       = sum(Compartment == "Reactive",   na.rm = TRUE),
      n_TrueActionable = sum(IsTrue & Compartment == "Actionable", na.rm = TRUE),
      ALY_conditional    = ifelse(True_Alarms > 0,
                                  n_TrueActionable / True_Alarms, NA_real_),
      WP_wks_conditional = ifelse(n_TrueActionable > 0,
                                  mean(Lead_Time[IsTrue & Compartment == "Actionable"], na.rm = TRUE),
                                  NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::left_join(per_detector_headline, by = "Detector") %>%
    dplyr::mutate(
      ALY                = round(ALY_headline,      3),
      WP_wks             = round(WP_wks_headline,   2),
      ALY_conditional    = round(ALY_conditional,   3),
      WP_wks_conditional = round(WP_wks_conditional, 2)
    ) %>%
    dplyr::select(-ALY_headline, -WP_wks_headline) %>%
    dplyr::select(
      Detector, Total_Triggers, True_Alarms, False_Alarms,
      n_Actionable, n_Reactive, n_TrueActionable,
      ALY, ALY_conditional,
      WP_wks, WP_wks_conditional
    )
}


# -----------------------------------------------------------------------------
# C.6 RUN PRIMARY SPECIFICATION (3 detectors)
# -----------------------------------------------------------------------------
DETECTOR_COLS <- c(
  "Constant_TA"        = "trig_vaezi",
  "Continuous_TA"      = "trig_classic",
  "Outbreak_Threshold" = "surge_outbreak_threshold"
)

primary_results <- lapply(names(DETECTOR_COLS), function(label) {
  res <- evaluate_detector(df, DETECTOR_COLS[[label]], EVALUABLE_YEARS)
  res$summary$Detector     <- label
  res$per_year$Detector    <- label
  res$per_trigger$Detector <- label
  res
})
names(primary_results) <- names(DETECTOR_COLS)

table1_primary  <- do.call(rbind, lapply(primary_results, function(r) r$summary))
rownames(table1_primary) <- NULL
table2_per_year <- do.call(rbind, lapply(primary_results, function(r) r$per_year))
rownames(table2_per_year) <- NULL
trigger_detail  <- do.call(rbind, lapply(primary_results, function(r) r$per_trigger))
rownames(trigger_detail) <- NULL

trigger_detail_aug_fig3 <- augment_compartments(trigger_detail, table2_per_year)
table1b_compartments    <- compute_compartment_metrics(trigger_detail_aug_fig3)


# -----------------------------------------------------------------------------
# C.7 PER-YEAR COMPARTMENT-AWARE METRICS (Table 2A)
# -----------------------------------------------------------------------------
compute_per_year_compartments <- function(trigger_df_aug, per_year_df) {
  if (nrow(trigger_df_aug) == 0) {
    return(per_year_df %>% dplyr::mutate(
      PPV_yr = NA_real_,
      n_Actionable_yr = 0L,
      n_Reactive_yr = 0L,
      n_TrueActionable_yr = 0L,
      Lead_Compartment = NA_character_,
      ALY_yr = 0,                ALY_yr_conditional = NA_real_,
      WP_yr_wks = 0,             WP_yr_wks_conditional = NA_real_,
      MLT_yr_wks = 0,            MLT_yr_wks_conditional = NA_real_,
      TAM_yr = NA_real_
    ))
  }
  per_year_compartments <- trigger_df_aug %>%
    dplyr::group_by(Year, Detector) %>%
    dplyr::summarise(
      n_Actionable_yr     = sum(Compartment == "Actionable", na.rm = TRUE),
      n_Reactive_yr       = sum(Compartment == "Reactive",   na.rm = TRUE),
      n_TrueActionable_yr = sum(IsTrue & Compartment == "Actionable", na.rm = TRUE),
      WP_yr_wks_conditional = ifelse(
        sum(IsTrue & Compartment == "Actionable", na.rm = TRUE) > 0,
        mean(Lead_Time[IsTrue & Compartment == "Actionable"], na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      WP_yr_wks_conditional = round(WP_yr_wks_conditional, 2),
      WP_yr_wks = ifelse(is.na(WP_yr_wks_conditional), 0, WP_yr_wks_conditional)
    )
  
  trigger_df_with_dc <- trigger_df_aug %>%
    dplyr::left_join(
      df %>% dplyr::transmute(Year = YR, Week = WN, DC = DC_QC),
      by = c("Year", "Week")
    )
  tam_per_year <- trigger_df_with_dc %>%
    dplyr::filter(IsTrue) %>%
    dplyr::group_by(Year, Detector) %>%
    dplyr::summarise(TAM_yr = sum(DC, na.rm = TRUE), .groups = "drop")
  
  per_year_df %>%
    dplyr::mutate(
      PPV_yr = ifelse(n_triggers > 0, round(n_True / n_triggers, 3), NA_real_),
      Lead_Compartment = vapply(Lead_Time_Weeks, classify_compartment, character(1))
    ) %>%
    dplyr::left_join(per_year_compartments, by = c("Year", "Detector")) %>%
    dplyr::left_join(tam_per_year,           by = c("Year", "Detector")) %>%
    dplyr::mutate(
      n_Actionable_yr     = ifelse(is.na(n_Actionable_yr),     0L, n_Actionable_yr),
      n_Reactive_yr       = ifelse(is.na(n_Reactive_yr),       0L, n_Reactive_yr),
      n_TrueActionable_yr = ifelse(is.na(n_TrueActionable_yr), 0L, n_TrueActionable_yr),
      TAM_yr = ifelse(is.na(TAM_yr), 0, TAM_yr),
      WP_yr_wks             = ifelse(is.na(WP_yr_wks), 0, WP_yr_wks),
      WP_yr_wks_conditional = WP_yr_wks_conditional,
      ALY_yr_conditional = ifelse(n_True > 0,
                                  round(n_TrueActionable_yr / n_True, 3),
                                  NA_real_),
      ALY_yr             = ifelse(n_True > 0,
                                  round(n_TrueActionable_yr / n_True, 3),
                                  0),
      MLT_yr_wks_conditional = Lead_Time_Weeks,
      MLT_yr_wks             = ifelse(is.na(Lead_Time_Weeks), 0, Lead_Time_Weeks)
    )
}

table2a_per_year_compartments <- compute_per_year_compartments(
  trigger_detail_aug_fig3, table2_per_year)


# -----------------------------------------------------------------------------
# C.8 TABLE 2B - PER-YEAR TA vs OT SIDE-BY-SIDE COMPARISON
# -----------------------------------------------------------------------------
build_table2b_side_by_side <- function(per_year_compartments_df) {
  if (nrow(per_year_compartments_df) == 0L)
    return(data.frame(Year = integer(0)))
  
  base_cols <- per_year_compartments_df %>%
    dplyr::filter(Detector == "Constant_TA") %>%
    dplyr::select(Year, Peak_Week, A1_start, A1_end, A2_start, A2_end)
  
  reshape_for_detector <- function(det_label, prefix) {
    per_year_compartments_df %>%
      dplyr::filter(Detector == det_label) %>%
      dplyr::transmute(
        Year,
        !!paste0(prefix, "_n_triggers")             := n_triggers,
        !!paste0(prefix, "_n_True")                 := n_True,
        !!paste0(prefix, "_n_False")                := n_False,
        !!paste0(prefix, "_PPV_yr")                 := PPV_yr,
        !!paste0(prefix, "_TAM_yr")                 := TAM_yr,
        !!paste0(prefix, "_First_True_Wk")          := First_A1_True_Week,
        !!paste0(prefix, "_Lead_Time_wks")          := Lead_Time_Weeks,
        !!paste0(prefix, "_Lead_Compartment")       := Lead_Compartment,
        !!paste0(prefix, "_n_Actionable_yr")        := n_Actionable_yr,
        !!paste0(prefix, "_n_Reactive_yr")          := n_Reactive_yr,
        !!paste0(prefix, "_n_TrueActionable_yr")    := n_TrueActionable_yr,
        !!paste0(prefix, "_MLT_yr_wks")             := MLT_yr_wks,
        !!paste0(prefix, "_MLT_yr_wks_conditional") := MLT_yr_wks_conditional,
        !!paste0(prefix, "_WP_yr_wks")              := WP_yr_wks,
        !!paste0(prefix, "_WP_yr_wks_conditional")  := WP_yr_wks_conditional,
        !!paste0(prefix, "_ALY_yr")                 := ALY_yr,
        !!paste0(prefix, "_ALY_yr_conditional")     := ALY_yr_conditional
      )
  }
  
  ta_wide <- reshape_for_detector("Constant_TA",        "TA")
  ot_wide <- reshape_for_detector("Outbreak_Threshold", "OT")
  
  base_cols %>%
    dplyr::left_join(ta_wide, by = "Year") %>%
    dplyr::left_join(ot_wide, by = "Year") %>%
    dplyr::arrange(Year)
}

table2b_ta_vs_ot <- build_table2b_side_by_side(table2a_per_year_compartments)


# -----------------------------------------------------------------------------
# C.9 SENSITIVITY ANALYSES (S2: A1 window; S3: A2 burden; S5: year inclusion)
# -----------------------------------------------------------------------------
S2_WINDOWS <- list("3-6" = c(3L, 6L), "4-8" = c(4L, 8L), "5-10" = c(5L, 10L))
S3_BURDEN  <- c(0.60, 0.70, 0.80)
S5_INCLUSIONS <- list(
  "primary"  = EVALUABLE_YEARS,
  "add_2020" = sort(c(EVALUABLE_YEARS, 2020L)),
  "add_2021" = sort(c(EVALUABLE_YEARS, 2021L)),
  "add_2025" = sort(c(EVALUABLE_YEARS, 2025L)),
  "add_all"  = 2013:2025
)

S2_rows <- list()
for (label in names(S2_WINDOWS)) {
  win <- S2_WINDOWS[[label]]
  for (det_label in names(DETECTOR_COLS)) {
    res <- evaluate_detector(df, DETECTOR_COLS[[det_label]], EVALUABLE_YEARS,
                             lead_min = win[1], lead_max = win[2])
    s <- res$summary; s$Detector <- det_label; s$Window <- label
    S2_rows[[length(S2_rows) + 1L]] <- s
  }
}
S2_table <- dplyr::bind_rows(S2_rows) %>%
  dplyr::select(Window, Detector, Total_Triggers, True_Alarms, False_Alarms,
                PPV, Sensitivity, Mean_Lead_Time_wks)

S3_rows <- list()
for (bf in S3_BURDEN) {
  for (det_label in names(DETECTOR_COLS)) {
    res <- evaluate_detector(df, DETECTOR_COLS[[det_label]], EVALUABLE_YEARS,
                             burden_frac = bf)
    s <- res$summary; s$Detector <- det_label; s$Burden_Frac <- bf
    S3_rows[[length(S3_rows) + 1L]] <- s
  }
}
S3_table <- dplyr::bind_rows(S3_rows) %>%
  dplyr::select(Burden_Frac, Detector, Total_Triggers, True_Alarms, False_Alarms,
                PPV, Sensitivity, Mean_Lead_Time_wks)

S5_rows <- list()
for (incl_label in names(S5_INCLUSIONS)) {
  yrs <- S5_INCLUSIONS[[incl_label]]
  for (det_label in names(DETECTOR_COLS)) {
    res <- evaluate_detector(df, DETECTOR_COLS[[det_label]], yrs)
    s <- res$summary; s$Detector <- det_label
    s$Inclusion <- incl_label; s$N_Years <- length(yrs)
    S5_rows[[length(S5_rows) + 1L]] <- s
  }
}
S5_table <- dplyr::bind_rows(S5_rows) %>%
  dplyr::select(Inclusion, N_Years, Detector, Total_Triggers, True_Alarms, False_Alarms,
                PPV, Sensitivity, Mean_Lead_Time_wks)

# Inter-rater agreement (Fleiss kappa) is not used in this two-anchor
# specification.
kappa_results <- NULL


# -----------------------------------------------------------------------------
# C.10 YEARLY SUMMARY TABLES (Constant TA and Continuous TA)
# -----------------------------------------------------------------------------
build_year_summary <- function(data, variant = c("classic", "vaezi")) {
  variant   <- match.arg(variant)
  trig_col  <- if (variant == "classic") "trig_classic" else "trig_vaezi"
  det_label <- if (variant == "classic") "Continuous_TA" else "Constant_TA"
  
  do.call(rbind, lapply(2013:2025, function(yr) {
    df_year <- data %>% dplyr::filter(YR == yr)
    if (nrow(df_year) == 0) {
      return(data.frame(Year = yr, Peak_Week = NA_integer_, Peak_DC = NA_real_,
                        First_Method_Week = NA_integer_, First_Outbreak_Threshold_Week = NA_integer_,
                        Lead_Lag = NA_integer_, Evaluable = !(yr %in% EXCLUDED_YEARS),
                        n_True_Alarms = NA_integer_, n_False_Alarms = NA_integer_,
                        stringsAsFactors = FALSE))
    }
    method_all   <- df_year %>% dplyr::filter(.data[[trig_col]] == 1)
    outbreak_all <- df_year %>% dplyr::filter(surge_outbreak_threshold == 1)
    first_method_week   <- if (nrow(method_all)   > 0) method_all$WN[1]   else NA_integer_
    first_outbreak_week <- if (nrow(outbreak_all) > 0) outbreak_all$WN[1] else NA_integer_
    peak_row <- df_year %>% dplyr::filter(DC_QC == max(DC_QC, na.rm = TRUE)) %>% dplyr::slice(1)
    peak_week <- if (nrow(peak_row) > 0) peak_row$WN[1]    else NA_integer_
    peak_dc   <- if (nrow(peak_row) > 0) peak_row$DC_QC[1] else NA_real_
    is_eval <- yr %in% EVALUABLE_YEARS
    if (is_eval) {
      eval_row <- table2_per_year %>% dplyr::filter(Year == yr, Detector == det_label) %>% dplyr::slice(1)
      n_T <- if (nrow(eval_row)) eval_row$n_True  else NA_integer_
      n_F <- if (nrow(eval_row)) eval_row$n_False else NA_integer_
    } else {
      n_T <- NA_integer_; n_F <- NA_integer_
    }
    data.frame(Year = yr, Peak_Week = peak_week, Peak_DC = peak_dc,
               First_Method_Week = first_method_week,
               First_Outbreak_Threshold_Week = first_outbreak_week,
               Lead_Lag = ifelse(!is.na(first_method_week) && !is.na(peak_week),
                                 peak_week - first_method_week, NA_integer_),
               Evaluable = is_eval, n_True_Alarms = n_T, n_False_Alarms = n_F,
               stringsAsFactors = FALSE)
  }))
}

classic_summary <- build_year_summary(df, "classic")
vaezi_summary   <- build_year_summary(df, "vaezi")


# -----------------------------------------------------------------------------
# C.11 PER-YEAR DETECTION PANEL BUILDER
# -----------------------------------------------------------------------------
make_detection_panel <- function(data, year_input, variant = c("classic", "vaezi")) {
  variant <- match.arg(variant)
  df_year <- data %>% dplyr::filter(YR == year_input)
  if (nrow(df_year) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = paste("No data:", year_input),
                               family = base_family_global))
  }
  
  is_evaluable <- year_input %in% EVALUABLE_YEARS
  det_label <- if (variant == "classic") "Continuous_TA" else "Constant_TA"
  
  if (variant == "classic") {
    r_col <- "R_classic"; trig_col <- "trig_classic"
    line_color <- "steelblue"
    method_label <- "Continuous TA"
    method_label_short <- "TA"
    panel_title <- as.character(year_input)
    R_MAX <- 4; eta_line <- ETA_ON_CLASSIC
  } else {
    r_col <- "R_vaezi"; trig_col <- "trig_vaezi"
    line_color <- "#d62728"
    method_label <- "Constant TA"
    method_label_short <- "TA"
    panel_title <- as.character(year_input)
    r_obs <- suppressWarnings(max(df_year[[r_col]], na.rm = TRUE))
    if (!is.finite(r_obs)) r_obs <- ETA_ON
    R_MAX <- max(r_obs * 1.15, ETA_ON * 1.33)
    eta_line <- ETA_ON
  }
  
  if (year_input %in% EXCLUDED_YEARS) {
    panel_title <- paste0(year_input, " (excluded)")
  }
  
  SCALE <- Y_MAX / R_MAX
  
  if (is_evaluable) {
    anchors <- compute_anchors_for_year(data, year_input)
  } else {
    anchors <- list(peak_week = NA_integer_,
                    A1_weeks = integer(0),
                    A2_weeks = integer(0))
  }
  
  # Trigger markers and lead arrows render only for triggers occurring on
  # or before the annual peak week. Post-peak triggers in the burden-block
  # tail are CSV-preserved but figure-suppressed.
  pk_wk_for_filter <- if (!is.na(anchors$peak_week)) anchors$peak_week else Inf
  
  df_year <- df_year %>%
    dplyr::mutate(
      .ta_trig_full = .data[[trig_col]] == 1L,
      .ot_trig_full = surge_outbreak_threshold == 1L,
      .ta_trig = .ta_trig_full & (WN <= pk_wk_for_filter),
      .ot_trig = .ot_trig_full & (WN <= pk_wk_for_filter),
      in_A1 = WN %in% anchors$A1_weeks,
      in_A2 = WN %in% anchors$A2_weeks,
      ta_class = dplyr::case_when(
        !is_evaluable & .ta_trig                          ~ "Excluded",
        .ta_trig & ( in_A1 | in_A2 )                      ~ "TrueAlarm",
        .ta_trig                                          ~ "FalseAlarm",
        TRUE                                              ~ "NoTrigger"
      ),
      ot_class = dplyr::case_when(
        !is_evaluable & .ot_trig                          ~ "Excluded",
        .ot_trig & ( in_A1 | in_A2 )                      ~ "TrueAlarm",
        .ot_trig                                          ~ "FalseAlarm",
        TRUE                                              ~ "NoTrigger"
      )
    )
  
  method_all       <- df_year %>% dplyr::filter(.ta_trig)
  method_all_true  <- df_year %>% dplyr::filter(ta_class == "TrueAlarm")
  
  outbreak_all       <- df_year %>% dplyr::filter(.ot_trig)
  outbreak_all_true  <- df_year %>% dplyr::filter(ot_class == "TrueAlarm")
  
  peak_row <- df_year %>% dplyr::filter(DC_QC == max(DC_QC, na.rm = TRUE)) %>% dplyr::slice(1)
  peak_week <- if (nrow(peak_row) > 0) peak_row$WN[1]    else NA_integer_
  peak_dc   <- if (nrow(peak_row) > 0) peak_row$DC_QC[1] else NA_real_
  
  # Render only the FIRST trigger per (detector, year) that lands in the
  # Actionable Window (A1). All other triggers are visually suppressed.
  # The full trigger inventory remains in the CSV outputs.
  first_method_a1 <- if (!is.na(peak_week)) {
    method_all_true %>%
      dplyr::mutate(.lt = peak_week - WN) %>%
      dplyr::filter(.lt >= A1_LEAD_MIN, .lt <= A1_LEAD_MAX) %>%
      dplyr::slice(1)
  } else {
    method_all_true[0, , drop = FALSE]
  }
  first_outbreak_a1 <- if (!is.na(peak_week)) {
    outbreak_all_true %>%
      dplyr::mutate(.lt = peak_week - WN) %>%
      dplyr::filter(.lt >= A1_LEAD_MIN, .lt <= A1_LEAD_MAX) %>%
      dplyr::slice(1)
  } else {
    outbreak_all_true[0, , drop = FALSE]
  }
  
  first_method_true   <- first_method_a1
  first_outbreak_true <- first_outbreak_a1
  first_method        <- first_method_a1
  first_outbreak      <- first_outbreak_a1
  
  get_safe_top_y <- function(val, upper_frac = 0.92, add_frac = 0.16) {
    pmin(val + add_frac * Y_MAX, Y_MAX * upper_frac)
  }
  get_safe_label_x <- function(wk) {
    if (wk <= 5)        list(x = wk + 1.2, hjust = 0)
    else if (wk >= 49)  list(x = wk - 1.2, hjust = 1)
    else                list(x = wk,       hjust = 0.5)
  }
  
  method_label_y <- if (nrow(first_method) > 0)
    get_safe_top_y(first_method$DC_QC, upper_frac = 0.68, add_frac = 0.18) else NA_real_
  outbreak_label_y <- if (nrow(first_outbreak) > 0)
    get_safe_top_y(first_outbreak$DC_QC, upper_frac = 0.60, add_frac = 0.10) else NA_real_
  
  if (nrow(first_method) > 0 && nrow(first_outbreak) > 0) {
    if (abs(method_label_y - outbreak_label_y) < 0.07 * Y_MAX) {
      method_label_y   <- pmin(method_label_y   + 0.08 * Y_MAX, Y_MAX * 0.68)
      outbreak_label_y <- pmin(outbreak_label_y - 0.04 * Y_MAX, Y_MAX * 0.56)
    }
  }
  
  first_method_true_week <-
    if (nrow(first_method_true) > 0) first_method_true$WN[1] else NA_integer_
  first_outbreak_true_week <-
    if (nrow(first_outbreak_true) > 0) first_outbreak_true$WN[1] else NA_integer_
  
  df_year <- df_year %>%
    dplyr::mutate(
      bar_fill = dplyr::case_when(
        !is.na(first_outbreak_true_week) & WN == first_outbreak_true_week ~ COL_BAR_OT_FIRSTTRUE,
        !is.na(first_method_true_week)   & WN == first_method_true_week   ~ line_color,
        TRUE                                                              ~ COL_BAR_DEFAULT
      )
    )
  
  p <- ggplot2::ggplot(df_year, ggplot2::aes(x = WN))
  
  if (!is_evaluable) {
    p <- p + ggplot2::annotate("rect",
                               xmin = 0.4, xmax = 53.6, ymin = 0, ymax = Inf,
                               fill = "grey85", alpha = 0.35)
  }
  
  p <- p +
    ggplot2::geom_col(ggplot2::aes(y = DC_QC, fill = bar_fill),
                      width = 0.9, alpha = 0.90, na.rm = TRUE, show.legend = FALSE) +
    ggplot2::scale_fill_identity()
  
  if (any(is.finite(df_year$outbreak_threshold))) {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = outbreak_threshold),
                                color = "black", linewidth = 0.75, linetype = "dashed",
                                na.rm = TRUE)
  }
  
  p <- p + ggplot2::geom_hline(yintercept = eta_line * SCALE,
                               linetype = "dashed", linewidth = 0.7, color = line_color)
  
  p <- p + ggplot2::geom_line(ggplot2::aes(y = .data[[r_col]] * SCALE),
                              color = line_color, linewidth = 0.85, na.rm = TRUE)
  
  if (is_evaluable && length(anchors$A2_weeks) > 0) {
    p <- p +
      ggplot2::geom_vline(xintercept = min(anchors$A2_weeks) - 0.5,
                          linetype = "dotted", colour = "#666666", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = max(anchors$A2_weeks) + 0.5,
                          linetype = "dotted", colour = "#666666", linewidth = 0.4)
  }
  
  if (nrow(first_outbreak_true) > 0) {
    p <- p + ggplot2::geom_point(data = first_outbreak_true,
                                 mapping = ggplot2::aes(x = WN, y = pmin(DC_QC + 0.05 * Y_MAX, Y_MAX * 0.95)),
                                 shape = 22, size = 3.0,
                                 fill = COL_TRUE_ALARM, colour = COL_TRUE_OUTLINE, inherit.aes = FALSE)
  }
  
  if (nrow(first_method_true) > 0) {
    p <- p + ggplot2::geom_point(data = first_method_true,
                                 mapping = ggplot2::aes(x = WN, y = pmin(DC_QC + 0.08 * Y_MAX, Y_MAX * 0.95)),
                                 shape = 24, size = 3.0,
                                 fill = COL_TRUE_ALARM, colour = COL_TRUE_OUTLINE, inherit.aes = FALSE)
  }
  
  if (nrow(first_outbreak) > 0) {
    pos_x <- get_safe_label_x(first_outbreak$WN[1])
    p <- p + ggplot2::annotate("text",
                               x = pos_x$x, y = outbreak_label_y,
                               label = paste0("Outbreak Threshold: W", first_outbreak$WN[1]),
                               hjust = pos_x$hjust, vjust = 0,
                               size = 2.0, colour = "#2c3e50", fontface = "bold",
                               family = base_family_global)
  }
  
  if (nrow(first_method) > 0) {
    pos_x <- get_safe_label_x(first_method$WN[1])
    p <- p + ggplot2::annotate("text",
                               x = pos_x$x, y = method_label_y,
                               label = paste0(method_label, ": W", first_method$WN[1]),
                               hjust = pos_x$hjust, vjust = 0,
                               size = 2.05, colour = line_color, fontface = "bold",
                               family = base_family_global)
  }
  
  # Lead/lag arrows: one per detector for its first A1 trigger.
  draw_lead_lag_arrow <- function(plot_obj, trigger_wk, peak_wk, arrow_y,
                                  detector_prefix, base_family) {
    if (is.na(trigger_wk) || is.na(peak_wk)) return(plot_obj)
    lt <- peak_wk - trigger_wk
    comp <- classify_compartment(lt)
    if (is.na(comp)) return(plot_obj)
    
    acol <- compartment_color(comp)
    if (is.na(acol)) acol <- "#555555"
    
    if (lt > 0) {
      lbl <- paste0(detector_prefix, ": ", lt, "wk ", comp)
    } else if (lt < 0) {
      lbl <- paste0(detector_prefix, ": ", abs(lt), "wk Reactive")
    } else {
      lbl <- paste0(detector_prefix, ": peak wk")
    }
    
    plot_obj +
      ggplot2::annotate("segment",
                        x = trigger_wk, xend = peak_wk,
                        y = arrow_y, yend = arrow_y,
                        colour = acol, linewidth = 0.9,
                        arrow = grid::arrow(type = "closed", length = grid::unit(0.13, "cm"))) +
      ggplot2::annotate("text",
                        x = (trigger_wk + peak_wk) / 2,
                        y = pmin(arrow_y + LABEL_OFFSET_Y, Y_MAX * 0.99),
                        label = lbl,
                        color = acol, fontface = "bold",
                        size = 2.0, family = base_family)
  }
  
  if (is_evaluable && is.finite(peak_week) && is.finite(peak_dc)) {
    if (nrow(first_method_true) > 0) {
      p <- draw_lead_lag_arrow(p,
                               trigger_wk = first_method_true$WN[1],
                               peak_wk    = peak_week,
                               arrow_y    = ARROW_Y_TA,
                               detector_prefix = method_label_short,
                               base_family = base_family_global)
    }
    if (nrow(first_outbreak_true) > 0) {
      p <- draw_lead_lag_arrow(p,
                               trigger_wk = first_outbreak_true$WN[1],
                               peak_wk    = peak_week,
                               arrow_y    = ARROW_Y_OT,
                               detector_prefix = "OT",
                               base_family = base_family_global)
    }
  }
  
  if (is_evaluable) {
    eval_row <- table2_per_year %>%
      dplyr::filter(Year == year_input, Detector == det_label) %>% dplyr::slice(1)
    eval_row_ot <- table2_per_year %>%
      dplyr::filter(Year == year_input, Detector == "Outbreak_Threshold") %>% dplyr::slice(1)
    badge_lines <- c()
    if (nrow(eval_row) > 0)
      badge_lines <- c(badge_lines,
                       paste0(method_label_short, "  T:", eval_row$n_True,
                              "  F:", eval_row$n_False))
    if (nrow(eval_row_ot) > 0)
      badge_lines <- c(badge_lines,
                       paste0("OT  T:", eval_row_ot$n_True, "  F:", eval_row_ot$n_False))
    if (length(badge_lines) > 0) {
      p <- p + ggplot2::annotate("text",
                                 x = 1.5, y = Y_MAX * 0.99,
                                 label = paste(badge_lines, collapse = "\n"),
                                 hjust = 0, vjust = 1, size = 2.1, colour = "black",
                                 fontface = "bold", family = base_family_global)
    }
  }
  
  # Two anchor bands stacked at panel bottom: Epidemic Burden (A2) below,
  # Actionable Window (A1) above.
  if (is_evaluable) {
    A2_BAND_YMIN <- 0.005 * Y_MAX
    A2_BAND_YMAX <- 0.030 * Y_MAX
    A1_BAND_YMIN <- 0.034 * Y_MAX
    A1_BAND_YMAX <- 0.059 * Y_MAX
    
    A2_BAND_FILL <- "#EF9F27"
    A1_BAND_FILL <- "#5DCAA5"
    A2_BAND_TICK <- "#854F0B"
    A1_BAND_TICK <- "#0F6E56"
    
    if (length(anchors$A2_weeks) > 0L) {
      a2_lo <- min(anchors$A2_weeks); a2_hi <- max(anchors$A2_weeks)
      p <- p + ggplot2::annotate("rect",
                                 xmin = a2_lo - 0.5, xmax = a2_hi + 0.5,
                                 ymin = A2_BAND_YMIN, ymax = A2_BAND_YMAX,
                                 fill = A2_BAND_FILL, alpha = 0.92,
                                 colour = A2_BAND_TICK, linewidth = 0.4)
      p <- p + ggplot2::annotate("segment",
                                 x = a2_lo - 0.5, xend = a2_lo - 0.5,
                                 y = A2_BAND_YMIN, yend = A2_BAND_YMAX,
                                 colour = A2_BAND_TICK, linewidth = 0.7)
      p <- p + ggplot2::annotate("segment",
                                 x = a2_hi + 0.5, xend = a2_hi + 0.5,
                                 y = A2_BAND_YMIN, yend = A2_BAND_YMAX,
                                 colour = A2_BAND_TICK, linewidth = 0.7)
      if ((a2_hi - a2_lo) >= 7L) {
        p <- p + ggplot2::annotate("text",
                                   x = (a2_lo + a2_hi) / 2,
                                   y = (A2_BAND_YMIN + A2_BAND_YMAX) / 2,
                                   label = paste0("Epidemic Burden: W", a2_lo, "-", a2_hi),
                                   size = 2.0, fontface = "bold", colour = "#412402")
      }
    }
    
    if (length(anchors$A1_weeks) > 0L) {
      a1_lo <- min(anchors$A1_weeks); a1_hi <- max(anchors$A1_weeks)
      p <- p + ggplot2::annotate("rect",
                                 xmin = a1_lo - 0.5, xmax = a1_hi + 0.5,
                                 ymin = A1_BAND_YMIN, ymax = A1_BAND_YMAX,
                                 fill = A1_BAND_FILL, alpha = 0.95,
                                 colour = A1_BAND_TICK, linewidth = 0.4)
      p <- p + ggplot2::annotate("segment",
                                 x = a1_lo - 0.5, xend = a1_lo - 0.5,
                                 y = A1_BAND_YMIN, yend = A1_BAND_YMAX,
                                 colour = A1_BAND_TICK, linewidth = 0.7)
      p <- p + ggplot2::annotate("segment",
                                 x = a1_hi + 0.5, xend = a1_hi + 0.5,
                                 y = A1_BAND_YMIN, yend = A1_BAND_YMAX,
                                 colour = A1_BAND_TICK, linewidth = 0.7)
      p <- p + ggplot2::annotate("text",
                                 x = (a1_lo + a1_hi) / 2,
                                 y = (A1_BAND_YMIN + A1_BAND_YMAX) / 2,
                                 label = paste0("Actionable Window: W", a1_lo, "-", a1_hi),
                                 size = 1.9, fontface = "bold", colour = "#04342C")
    }
  }
  
  p <- p +
    ggplot2::coord_cartesian(xlim = c(0.4, 53.6), ylim = c(0, Y_MAX), clip = "off") +
    ggplot2::scale_y_continuous(
      name = "Dengue cases", limits = c(0, Y_MAX),
      breaks = scales::pretty_breaks(n = 4),
      sec.axis = ggplot2::sec_axis(~ . / SCALE, name = "R(t)",
                                   breaks = seq(0, floor(R_MAX), by = 1))
    ) +
    ggplot2::scale_x_continuous(name = "Week number",
                                breaks = seq(1, 53, 13), labels = paste0("W", seq(1, 53, 13)),
                                expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
    ggplot2::labs(title = panel_title) +
    theme_dashboard_panel(base_size = 8, base_family = base_family_global) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 6.0),
      axis.text.y = ggplot2::element_text(size = 6.0),
      axis.text.y.right  = ggplot2::element_text(size = 6.0, color = line_color),
      axis.title.x = ggplot2::element_text(size = 6.5),
      axis.title.y = ggplot2::element_text(size = 6.5),
      axis.title.y.right = ggplot2::element_text(size = 6.5, color = line_color),
      plot.margin = ggplot2::margin(10, 24, 10, 12)
    )
  p
}


# -----------------------------------------------------------------------------
# C.12 GENERATE PER-YEAR PANELS (2013-2025)
# -----------------------------------------------------------------------------
years_to_plot <- 2013:2025

vaezi_panels_individual <- lapply(years_to_plot, function(y)
  make_detection_panel(df, y, "vaezi"))

cat("\n===============================================================\n")
cat("FIGURE 3A - CONSTANT TRANSMISSION ACCELERATION (PER YEAR)\n")
cat("===============================================================\n")
for (i in seq_along(years_to_plot)) print(vaezi_panels_individual[[i]])

classic_panels_individual <- lapply(years_to_plot, function(y)
  make_detection_panel(df, y, "classic"))

cat("\n========================================================================\n")
cat("SUPPLEMENTARY FIGURE 3 - CONTINUOUS TRANSMISSION ACCELERATION (PER YEAR)\n")
cat("========================================================================\n")
for (i in seq_along(years_to_plot)) print(classic_panels_individual[[i]])


# -----------------------------------------------------------------------------
# C.13 BUILD MULTIPANEL OBJECTS
# -----------------------------------------------------------------------------
years_main <- 2013:2024; year_last <- 2025

vaezi_main_panels   <- lapply(years_main, function(y)
  make_detection_panel(df, y, "vaezi"))
classic_main_panels <- lapply(years_main, function(y)
  make_detection_panel(df, y, "classic"))

vaezi_2025   <- make_detection_panel(df, year_last, "vaezi")
classic_2025 <- make_detection_panel(df, year_last, "classic")

blank_panel_left  <- ggplot2::ggplot() + ggplot2::theme_void()
blank_panel_right <- ggplot2::ggplot() + ggplot2::theme_void()

panel_export_theme <- ggplot2::theme(
  plot.margin  = ggplot2::margin(t = 8, r = 16, b = 8, l = 14),
  axis.title.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 6)),
  axis.title.y = ggplot2::element_text(size = 8, margin = ggplot2::margin(r = 6)),
  axis.text.x  = ggplot2::element_text(size = 7),
  axis.text.y  = ggplot2::element_text(size = 7),
  plot.title   = ggplot2::element_text(size = 10, margin = ggplot2::margin(b = 6))
)

vaezi_main_panels_styled   <- lapply(vaezi_main_panels,   function(p) p + panel_export_theme)
classic_main_panels_styled <- lapply(classic_main_panels, function(p) p + panel_export_theme)
vaezi_2025_styled   <- vaezi_2025   + panel_export_theme
classic_2025_styled <- classic_2025 + panel_export_theme

figure3a_multipanel <- cowplot::plot_grid(
  plotlist = c(vaezi_main_panels_styled,
               list(blank_panel_left, vaezi_2025_styled, blank_panel_right)),
  ncol = 3, nrow = 5, align = "hv", axis = "tblr",
  rel_widths = c(1, 1, 1), rel_heights = c(1, 1, 1, 1, 1)
)

suppfig_multipanel <- cowplot::plot_grid(
  plotlist = c(classic_main_panels_styled,
               list(blank_panel_left, classic_2025_styled, blank_panel_right)),
  ncol = 3, nrow = 5, align = "hv", axis = "tblr",
  rel_widths = c(1, 1, 1), rel_heights = c(1, 1, 1, 1, 1)
)

print(figure3a_multipanel)
print(suppfig_multipanel)


# -----------------------------------------------------------------------------
# C.14 SAVE FIGURE 3 MULTIPANELS
# -----------------------------------------------------------------------------
pdf_file_3a  <- file.path(OUT_DIR, "Figure3A_Constant_Transmission_Acceleration_Multipanel.pdf")
pdf_file_sup <- file.path(OUT_DIR, "Supplementary_Figure3_Continuous_Transmission_Acceleration_Multipanel.pdf")
png_file_3a  <- file.path(OUT_DIR, "Figure3A_Constant_Transmission_Acceleration_Multipanel.png")
png_file_sup <- file.path(OUT_DIR, "Supplementary_Figure3_Continuous_Transmission_Acceleration_Multipanel.png")

ggplot2::ggsave(filename = pdf_file_3a,  plot = figure3a_multipanel,
                device = safe_pdf_device(),
                width = fig_width_in_main, height = fig_height_in_main,
                units = "in", dpi = 600, limitsize = FALSE, bg = "white")
ggplot2::ggsave(filename = pdf_file_sup, plot = suppfig_multipanel,
                device = safe_pdf_device(),
                width = fig_width_in_main, height = fig_height_in_main,
                units = "in", dpi = 600, limitsize = FALSE, bg = "white")
ggplot2::ggsave(filename = png_file_3a,  plot = figure3a_multipanel, device = "png",
                width = fig_width_in_main, height = fig_height_in_main,
                units = "in", dpi = 600, limitsize = FALSE, bg = "white")
ggplot2::ggsave(filename = png_file_sup, plot = suppfig_multipanel, device = "png",
                width = fig_width_in_main, height = fig_height_in_main,
                units = "in", dpi = 600, limitsize = FALSE, bg = "white")

cat("\nSaved Figure 3 multipanels:\n")
cat("PDF : ", pdf_file_3a,  "\n", sep = "")
cat("PDF : ", pdf_file_sup, "\n", sep = "")
cat("PNG : ", png_file_3a,  "\n", sep = "")
cat("PNG : ", png_file_sup, "\n", sep = "")


# -----------------------------------------------------------------------------
# C.15 HEAD-TO-HEAD TA vs OT STATISTICAL COMPARISON
# -----------------------------------------------------------------------------
# Three pairwise comparisons across 8 metrics:
#   (1) Constant TA   vs Outbreak Threshold
#   (2) Continuous TA vs Outbreak Threshold
#   (3) Constant TA   vs Continuous TA
# Test : Wilcoxon signed-rank, paired by year.
# Significance: each pair tested at alpha = HH_ALPHA on its own (no
#   across-pair multiplicity correction; each pair answers a distinct
#   scientific question).
# Effect size: median paired difference; 95% percentile-method bootstrap CI
#   (year-cluster, B = HH_BOOT_N replicates).
# -----------------------------------------------------------------------------

cat("\n=================================================================\n")
cat("HEAD-TO-HEAD STATISTICAL COMPARISON\n")
cat("Wilcoxon signed-rank, paired by year, B =", HH_BOOT_N, "bootstrap reps\n")
cat("=================================================================\n")

# Long-format per-year metrics with same-denominator timeliness columns.
hh_long <- table2a_per_year_compartments %>%
  dplyr::transmute(
    Year, Detector,
    N_True_Alarms  = as.numeric(n_True),
    N_False_Alarms = as.numeric(n_False),
    TAM         = TAM_yr,
    PPV         = PPV_yr,
    Sensitivity = as.integer(!is.na(First_A1_True_Week)),
    MLT         = MLT_yr_wks,
    WP          = WP_yr_wks,
    ALY         = ALY_yr
  )

HH_METRICS <- c("TAM", "N_True_Alarms", "PPV", "Sensitivity",
                "MLT", "WP", "ALY",
                "N_False_Alarms")
HH_METRIC_CATEGORY <- c(
  TAM            = "Epidemic_Burden_and_Alarm_Accuracy",
  N_True_Alarms  = "Epidemic_Burden_and_Alarm_Accuracy",
  PPV            = "Epidemic_Burden_and_Alarm_Accuracy",
  Sensitivity    = "Epidemic_Burden_and_Alarm_Accuracy",
  MLT            = "Timeliness",
  WP             = "Timeliness",
  ALY            = "Timeliness",
  N_False_Alarms = "False_Alarms"
)

HH_PAIRS <- list(
  list(label = "ConstantTA_vs_OT",    A = "Constant_TA",   B = "Outbreak_Threshold"),
  list(label = "ContinuousTA_vs_OT",  A = "Continuous_TA", B = "Outbreak_Threshold"),
  list(label = "ConstantTA_vs_ContinuousTA",
       A = "Constant_TA",   B = "Continuous_TA")
)

hh_better_when <- function(metric) {
  if (metric == "N_False_Alarms") "lower" else "higher"
}

hh_bootstrap_ci <- function(diff_vec, B = HH_BOOT_N, alpha = 0.05) {
  diff_vec <- diff_vec[!is.na(diff_vec)]
  if (length(diff_vec) < 2L) {
    return(c(lo = NA_real_, hi = NA_real_))
  }
  set.seed(20260101L)
  reps <- vapply(seq_len(B), function(b) {
    idx <- sample(seq_along(diff_vec), size = length(diff_vec), replace = TRUE)
    median(diff_vec[idx], na.rm = TRUE)
  }, numeric(1))
  c(lo = unname(stats::quantile(reps, alpha/2,     na.rm = TRUE, type = 7)),
    hi = unname(stats::quantile(reps, 1 - alpha/2, na.rm = TRUE, type = 7)))
}

hh_results_rows <- list()
for (m in HH_METRICS) {
  for (pair in HH_PAIRS) {
    a_vec <- hh_long %>% dplyr::filter(Detector == pair$A) %>%
      dplyr::arrange(Year) %>% dplyr::pull(!!m)
    b_vec <- hh_long %>% dplyr::filter(Detector == pair$B) %>%
      dplyr::arrange(Year) %>% dplyr::pull(!!m)
    a_yrs <- hh_long %>% dplyr::filter(Detector == pair$A) %>%
      dplyr::arrange(Year) %>% dplyr::pull(Year)
    b_yrs <- hh_long %>% dplyr::filter(Detector == pair$B) %>%
      dplyr::arrange(Year) %>% dplyr::pull(Year)
    common_yrs <- intersect(a_yrs, b_yrs)
    a_aligned <- a_vec[match(common_yrs, a_yrs)]
    b_aligned <- b_vec[match(common_yrs, b_yrs)]
    diff_vec <- a_aligned - b_aligned
    
    n_pairs <- sum(!is.na(diff_vec))
    median_diff <- if (n_pairs > 0L) median(diff_vec, na.rm = TRUE) else NA_real_
    
    wt <- tryCatch(
      suppressWarnings(stats::wilcox.test(a_aligned, b_aligned, paired = TRUE,
                                          exact = FALSE, alternative = "two.sided")),
      error = function(e) NULL
    )
    v_stat <- if (!is.null(wt)) unname(wt$statistic) else NA_real_
    p_val  <- if (!is.null(wt)) wt$p.value else NA_real_
    
    ci <- hh_bootstrap_ci(diff_vec)
    
    hh_results_rows[[length(hh_results_rows) + 1L]] <- data.frame(
      Metric              = m,
      Category            = unname(HH_METRIC_CATEGORY[m]),
      Comparison          = pair$label,
      Detector_A          = pair$A,
      Detector_B          = pair$B,
      Better_When         = hh_better_when(m),
      N_Years_Paired      = n_pairs,
      Median_Diff_AminusB = round(median_diff, 3),
      Bootstrap_CI_lo     = round(ci["lo"], 3),
      Bootstrap_CI_hi     = round(ci["hi"], 3),
      V_statistic         = round(v_stat, 2),
      p_value             = signif(p_val, 4),
      stringsAsFactors    = FALSE
    )
  }
}
hh_results <- dplyr::bind_rows(hh_results_rows)

# Per-pairwise alpha = HH_ALPHA: each pair stands on its own. Backward-
# compatibility alias `p_bonferroni` retained (equals raw p-value).
hh_results <- hh_results %>%
  dplyr::mutate(
    p_pairwise      = p_value,
    p_bonferroni    = p_value,
    Significant_005 = !is.na(p_pairwise) & p_pairwise < HH_ALPHA
  ) %>%
  dplyr::mutate(
    p_pairwise   = signif(p_pairwise,   4),
    p_bonferroni = signif(p_bonferroni, 4)
  )

hh_csv_path <- file.path(EVAL_OUT_DIR, "Figure3_HeadToHead_Wilcoxon_Results.csv")
utils::write.csv(hh_results, hh_csv_path, row.names = FALSE)
cat("Saved CSV: ", hh_csv_path, "\n", sep = "")

cat("\nHead-to-head results (Wilcoxon signed-rank, year-paired):\n")
safe_df_print(hh_results)


# -----------------------------------------------------------------------------
# C.16 HEAD-TO-HEAD MULTIPANEL FIGURES (one per pairwise comparison)
# -----------------------------------------------------------------------------
hh_metric_levels <- c(
  "TAM", "N_True_Alarms", "PPV", "Sensitivity",
  "MLT", "WP", "ALY",
  "N_False_Alarms"
)
hh_metric_full_name <- c(
  N_True_Alarms  = "Number of True Alarms",
  N_False_Alarms = "Number of False Alarms",
  TAM            = "True-Alarm Magnitude",
  PPV            = "Positive Predictive Value",
  Sensitivity    = "Sensitivity",
  MLT            = "Mean Lead Time",
  WP             = "Warning Persistence",
  ALY            = "Actionable Lead-Time Yield"
)
hh_metric_y_label <- c(
  N_True_Alarms  = "Number of True Alarms (per year)",
  N_False_Alarms = "Number of False Alarms (per year)",
  TAM            = "True-Alarm Magnitude (cases per year)",
  PPV            = "Positive Predictive Value",
  Sensitivity    = "Sensitivity",
  MLT            = "Mean Lead Time (weeks before peak)",
  WP             = "Warning Persistence (weeks before peak)",
  ALY            = "Actionable Lead-Time Yield"
)
hh_metric_higher_better <- c(
  N_True_Alarms = TRUE,  N_False_Alarms = FALSE,
  TAM = TRUE,  PPV = TRUE,  Sensitivity = TRUE,
  MLT = TRUE,  WP  = TRUE,  ALY         = TRUE
)
hh_metric_threshold <- c(
  N_True_Alarms  = NA_real_,
  N_False_Alarms = NA_real_,
  TAM            = NA_real_,
  PPV            = 0.50,
  Sensitivity    = 1.00,
  MLT            = 4,
  WP             = 4,
  ALY            = 0.50
)
HH_DETECTOR_FILL <- c(
  "Constant_TA"        = "#D33F3F",
  "Continuous_TA"      = "#3673B6",
  "Outbreak_Threshold" = "#3F4955"
)
HH_DETECTOR_LABEL <- c(
  "Constant_TA"        = "Constant TA",
  "Continuous_TA"      = "Continuous TA",
  "Outbreak_Threshold" = "Outbreak Threshold"
)

build_metric_subpanel <- function(pair_long, metric_id, det_levels, stat_row = NULL) {
  is_higher_better <- unname(hh_metric_higher_better[metric_id])
  threshold        <- unname(hh_metric_threshold[metric_id])
  y_label          <- unname(hh_metric_y_label[metric_id])
  metric_full      <- unname(hh_metric_full_name[metric_id])
  
  vmin <- suppressWarnings(min(pair_long$Value, na.rm = TRUE))
  vmax <- suppressWarnings(max(pair_long$Value, na.rm = TRUE))
  if (!is.finite(vmin) || !is.finite(vmax)) { vmin <- 0; vmax <- 1 }
  vrange <- vmax - vmin
  if (vrange == 0) vrange <- max(abs(vmax), 1)
  y_lo <- vmin - 0.18 * vrange
  y_hi <- vmax + 0.20 * vrange
  if (metric_id %in% c("PPV", "Sensitivity", "ALY")) {
    y_lo <- max(0, y_lo); y_hi <- min(1.05, y_hi)
  }
  
  pair_long <- pair_long %>%
    dplyr::mutate(Detector = factor(Detector, levels = det_levels))
  
  p <- ggplot2::ggplot(pair_long, ggplot2::aes(x = Detector, y = Value))
  
  if (!is.na(threshold)) {
    if (is_higher_better) {
      p <- p + ggplot2::annotate(
        "rect", xmin = -Inf, xmax = Inf, ymin = threshold, ymax = Inf,
        fill = "#d8ecd2", alpha = 0.30
      )
    } else {
      p <- p + ggplot2::annotate(
        "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = threshold,
        fill = "#d8ecd2", alpha = 0.30
      )
    }
    p <- p + ggplot2::geom_hline(
      yintercept = threshold, linetype = "dashed",
      linewidth = 0.4, colour = "grey40"
    )
  }
  
  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(fill = Detector),
      shape = 21, color = "grey25", size = 2.4, stroke = 0.30,
      alpha = 0.85, na.rm = TRUE,
      position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 42)
    ) +
    ggplot2::stat_summary(
      fun = mean, geom = "crossbar",
      width = 0.55,
      colour = "grey15", linewidth = 0.35, na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = HH_DETECTOR_FILL[det_levels],
      labels = HH_DETECTOR_LABEL[det_levels],
      name   = "Detector"
    ) +
    ggplot2::scale_x_discrete(labels = HH_DETECTOR_LABEL[det_levels]) +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi)) +
    ggplot2::labs(title = metric_full, x = NULL, y = y_label) +
    ggplot2::theme_bw(base_size = 9, base_family = base_family_global) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(face = "bold", size = 9.2,
                                                 margin = ggplot2::margin(b = 4)),
      axis.text.x        = ggplot2::element_text(size = 7.5),
      axis.text.y        = ggplot2::element_text(size = 7.0),
      axis.title.y       = ggplot2::element_text(size = 7.0,
                                                 margin = ggplot2::margin(r = 4)),
      legend.position    = "none"
    )
  
  if (metric_id %in% c("PPV", "Sensitivity", "ALY")) {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    )
  }
  
  # Embed pairwise p-value annotation in the lower-left of each subpanel.
  p_disp_val <- if ("p_pairwise" %in% names(stat_row) &&
                    nrow(stat_row) > 0L &&
                    !is.na(stat_row$p_pairwise[1])) {
    stat_row$p_pairwise[1]
  } else if (!is.null(stat_row) && nrow(stat_row) > 0L &&
             "p_bonferroni" %in% names(stat_row) &&
             !is.na(stat_row$p_bonferroni[1])) {
    stat_row$p_bonferroni[1]
  } else {
    NA_real_
  }
  if (!is.na(p_disp_val)) {
    p_label <- if (p_disp_val < 0.001) "p < 0.001" else
      paste0("p = ", sprintf("%.3f", p_disp_val))
    sig_star <- if (!is.na(stat_row$Significant_005[1]) && stat_row$Significant_005[1]) " *" else ""
    annot_label <- paste0(p_label, sig_star)
    annot_color <- if (!is.na(stat_row$Significant_005[1]) && stat_row$Significant_005[1])
      "#1F2D5C" else "grey35"
    p <- p + ggplot2::annotate(
      "text",
      x = 0.55,
      y = y_lo + 0.04 * (y_hi - y_lo),
      label = annot_label,
      hjust = 0, vjust = 0,
      family = base_family_global, size = 2.7, fontface = "bold",
      colour = annot_color
    )
  }
  
  p
}

build_category_header <- function(category_label, n_metrics_in_section) {
  fill_col <- switch(
    category_label,
    "Epidemic Burden & Alarm Accuracy" = "#2C5DAA",
    "Early Warning Timeliness"         = "#5C9D7A",
    "False Alarms"                     = "#C0392B",
    "#374151"
  )
  ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                      fill = fill_col, alpha = 0.92) +
    ggplot2::annotate("text", x = 0.02, y = 0.5,
                      label = category_label,
                      hjust = 0, vjust = 0.5,
                      family = base_family_global,
                      size = 4.4, fontface = "bold", color = "white") +
    ggplot2::annotate("text", x = 0.98, y = 0.5,
                      label = paste0("(", n_metrics_in_section,
                                     " metric", if (n_metrics_in_section != 1L) "s" else "",
                                     ")"),
                      hjust = 1, vjust = 0.5,
                      family = base_family_global,
                      size = 3.0, fontface = "italic", color = "white") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
}

build_metric_section <- function(metric_blocks, metric_ids_in_section) {
  blocks <- metric_blocks[metric_ids_in_section]
  n <- length(blocks)
  
  if (n %% 2L == 0L) {
    return(patchwork::wrap_plots(blocks, ncol = 2, byrow = TRUE))
  }
  
  n_full_rows <- (n - 1L) %/% 2L
  if (n_full_rows == 0L) {
    centered_row <- patchwork::wrap_plots(
      list(patchwork::plot_spacer(), blocks[[1L]], patchwork::plot_spacer()),
      ncol = 3
    ) +
      patchwork::plot_layout(widths = c(1, 2, 1))
    return(centered_row)
  }
  
  full_blocks  <- blocks[seq_len(n - 1L)]
  full_section <- patchwork::wrap_plots(full_blocks, ncol = 2, byrow = TRUE)
  
  centered_row <- patchwork::wrap_plots(
    list(patchwork::plot_spacer(), blocks[[n]], patchwork::plot_spacer()),
    ncol = 3
  ) +
    patchwork::plot_layout(widths = c(1, 2, 1))
  
  patchwork::wrap_plots(list(full_section, centered_row), ncol = 1) +
    patchwork::plot_layout(heights = c(n_full_rows, 1))
}

build_headtohead_figure <- function(pair, hh_long_df, hh_results_df) {
  pair_label <- pair$label
  det_A <- pair$A; det_B <- pair$B
  det_levels <- c(det_A, det_B)
  
  metric_blocks <- list()
  for (m in hh_metric_levels) {
    pair_long <- hh_long_df %>%
      dplyr::filter(Detector %in% det_levels) %>%
      dplyr::transmute(Year, Detector, Value = .data[[m]])
    stat_row <- hh_results_df %>%
      dplyr::filter(Comparison == pair_label, Metric == m)
    
    metric_blocks[[m]] <- build_metric_subpanel(pair_long, m, det_levels,
                                                stat_row = stat_row)
  }
  
  category_definitions <- list(
    list(label = "Epidemic Burden & Alarm Accuracy",
         metrics = c("TAM", "N_True_Alarms", "PPV", "Sensitivity")),
    list(label = "Early Warning Timeliness",
         metrics = c("MLT", "WP", "ALY")),
    list(label = "False Alarms",
         metrics = c("N_False_Alarms"))
  )
  
  vertical_components <- list()
  vertical_heights    <- numeric(0)
  for (catdef in category_definitions) {
    metrics_in_cat <- catdef$metrics
    n_metrics <- length(metrics_in_cat)
    n_rows_in_section <- ceiling(n_metrics / 2L)
    header  <- build_category_header(catdef$label, n_metrics)
    section <- build_metric_section(metric_blocks, metrics_in_cat)
    
    vertical_components[[length(vertical_components) + 1L]] <- header
    vertical_heights <- c(vertical_heights, 0.20)
    vertical_components[[length(vertical_components) + 1L]] <- section
    vertical_heights <- c(vertical_heights, n_rows_in_section * 1.00)
  }
  
  grid_plot <- patchwork::wrap_plots(vertical_components, ncol = 1) +
    patchwork::plot_layout(heights = vertical_heights)
  
  title_str <- paste0("Head-to-head: ",
                      HH_DETECTOR_LABEL[det_A], "  vs  ",
                      HH_DETECTOR_LABEL[det_B])
  
  grid_plot +
    patchwork::plot_annotation(
      title = title_str,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12,
                                           family = base_family_global,
                                           margin = ggplot2::margin(b = 6))
      )
    )
}

HH_FIGURE_OUT <- list(
  list(label = "ConstantTA_vs_OT",           filename_stem = "Figure3_HeadToHead_ConstantTA_vs_OT"),
  list(label = "ContinuousTA_vs_OT",         filename_stem = "Figure3_HeadToHead_ContinuousTA_vs_OT"),
  list(label = "ConstantTA_vs_ContinuousTA", filename_stem = "Figure3_HeadToHead_ConstantTA_vs_ContinuousTA")
)

build_headtohead_legend <- function(pair, results_for_pair, fig_id) {
  pair_label <- pair$label
  det_A <- pair$A; det_B <- pair$B
  pretty_A <- HH_DETECTOR_LABEL[det_A]
  pretty_B <- HH_DETECTOR_LABEL[det_B]
  
  res <- results_for_pair %>%
    dplyr::filter(Comparison == pair_label) %>%
    dplyr::arrange(match(Metric, hh_metric_levels))
  sig_metrics <- res %>%
    dplyr::filter(!is.na(Significant_005), Significant_005) %>%
    dplyr::pull(Metric)
  n_sig <- length(sig_metrics)
  
  n_total <- length(hh_metric_levels)
  outcome_str <- if (n_sig == 0L) {
    paste0("No metric reached the per-pairwise significance threshold ",
           "(p < ", sprintf("%.2f", HH_ALPHA), ") on the available evidence base ",
           "of n = ", length(EVALUABLE_YEARS), " evaluable years.")
  } else if (n_sig == n_total) {
    paste0("All ", n_total, " metrics reached the per-pairwise significance ",
           "threshold (p < ", sprintf("%.2f", HH_ALPHA), ").")
  } else {
    pretty_sig <- vapply(sig_metrics, function(mc) unname(hh_metric_full_name[mc]),
                         character(1))
    paste0(n_sig, " of ", n_total, " metrics reached the per-pairwise significance ",
           "threshold (p < ", sprintf("%.2f", HH_ALPHA), "): ",
           paste(pretty_sig, collapse = "; "), ".")
  }
  
  paste0(
    fig_id, " | ", pretty_A, " versus ", pretty_B,
    " across ", n_total, " performance metrics in three operational categories ",
    "(Epidemic Burden & Alarm Accuracy: True-Alarm Magnitude, Number of True ",
    "Alarms, Positive Predictive Value, Sensitivity; Early Warning Timeliness: ",
    "Mean Lead Time, Warning Persistence, Actionable Lead-Time Yield; False ",
    "Alarms: Number of False Alarms). The figure is organised by category: ",
    "each of the three categories is preceded by a colored banner header ",
    "(Epidemic Burden & Alarm Accuracy in blue, Timeliness in green, False ",
    "Alarms in red) followed by the category's metric subpanels in a 2-column ",
    "layout. Within each metric subpanel, every dot represents one evaluable ",
    "year (n = ", length(EVALUABLE_YEARS),
    "; excluded years 2020, 2021, 2025); dots are jittered horizontally to ",
    "reduce overplotting and coloured by detector identity. The across-year ",
    "mean for each detector is shown as a horizontal crossbar. Pale-green ",
    "shaded regions indicate operationally favourable zones for each metric: ",
    "above the dashed reference threshold for higher-is-better metrics and ",
    "below the threshold for the lower-is-better Number of False Alarms ",
    "metric. The Wilcoxon paired signed-rank p value is annotated in BOLD at ",
    "the lower-left of each metric subpanel as 'p = X.XXX' (with adjacent ",
    "'*' if p < ", sprintf("%.2f", HH_ALPHA), "; navy if significant, grey ",
    "otherwise). The Wilcoxon test is paired by year, two-sided, and non-",
    "parametric. Significance is evaluated PER PAIRWISE COMPARISON: each ",
    "pair (Constant TA vs OT; Continuous TA vs OT; Constant TA vs Continuous ",
    "TA) is tested at alpha = ", sprintf("%.2f", HH_ALPHA),
    " on its own; p-values are NOT multiplied across pairs within a metric. ",
    "Effect sizes (median paired differences) and 95% year-cluster bootstrap ",
    "confidence intervals (B = ", HH_BOOT_N,
    " replicates) are reported in the companion CSV ",
    "(Figure3_HeadToHead_Wilcoxon_Results.csv). ",
    "Sensitivity is computed as the proportion of evaluable seasons with at ",
    "least one True Alarm in the Actionable Window (A1). The three Early ",
    "Warning Timeliness metrics share the same denominator per detector ",
    "(count of evaluable years): years where the metric is not computable ",
    "due to no qualifying triggers contribute zero rather than being ",
    "excluded. ", "Outcome on this pair: ", outcome_str, "\n"
  )
}

hh_w <- 12.0; hh_h <- 17.0

for (i in seq_along(HH_FIGURE_OUT)) {
  cfg <- HH_FIGURE_OUT[[i]]
  pair_idx <- which(vapply(HH_PAIRS, function(p) p$label == cfg$label, logical(1)))
  pair <- HH_PAIRS[[pair_idx]]
  fig  <- build_headtohead_figure(pair, hh_long, hh_results)
  
  pdf_out <- file.path(OUT_DIR, paste0(cfg$filename_stem, ".pdf"))
  png_out <- file.path(OUT_DIR, paste0(cfg$filename_stem, ".png"))
  
  ggplot2::ggsave(filename = pdf_out, plot = fig, device = safe_pdf_device(),
                  width = hh_w, height = hh_h, units = "in", dpi = 600,
                  limitsize = FALSE, bg = "white")
  ggplot2::ggsave(filename = png_out, plot = fig, device = "png",
                  width = hh_w, height = hh_h, units = "in", dpi = 600,
                  limitsize = FALSE, bg = "white")
  
  cat("Saved figure: ", pdf_out, "\n", sep = "")
  cat("Saved figure: ", png_out, "\n", sep = "")
  
  fig_id <- paste0("Figure 3 (head-to-head, ", i, " of ",
                   length(HH_FIGURE_OUT), ")")
  legend_text <- build_headtohead_legend(pair, hh_results, fig_id)
  
  cat("\n=================================================================\n")
  cat("DRAFT LEGEND  ", cfg$filename_stem, "\n", sep = "")
  cat("=================================================================\n")
  cat(legend_text)
}


# -----------------------------------------------------------------------------
# C.17 OPTIONAL: SAVE EACH YEAR AS INDIVIDUAL PDF
# -----------------------------------------------------------------------------
for (y in years_to_plot) {
  p_vaezi <- make_detection_panel(df, y, "vaezi")
  ggplot2::ggsave(
    filename = file.path(OUT_DIR, paste0("Figure3A_Constant_Transmission_Acceleration_", y, ".pdf")),
    plot = p_vaezi, device = safe_pdf_device(),
    width = fig_width_in_single, height = fig_height_in_single,
    units = "in", dpi = 600, limitsize = FALSE, bg = "white")
  p_classic <- make_detection_panel(df, y, "classic")
  ggplot2::ggsave(
    filename = file.path(OUT_DIR, paste0("Supplementary_Figure3_Continuous_Transmission_Acceleration_", y, ".pdf")),
    plot = p_classic, device = safe_pdf_device(),
    width = fig_width_in_single, height = fig_height_in_single,
    units = "in", dpi = 600, limitsize = FALSE, bg = "white")
}


# -----------------------------------------------------------------------------
# C.18 SAVE FIGURE 3 + EVALUATION FRAMEWORK CSVs
# -----------------------------------------------------------------------------
utils::write.csv(thresholds_full,
                 file.path(OUT_DIR, "Figure3_Rolling_Weekly_Outbreak_Thresholds.csv"),
                 row.names = FALSE)
utils::write.csv(vaezi_summary,
                 file.path(OUT_DIR, "Figure3A_Constant_Transmission_Acceleration_Yearly_Summary.csv"),
                 row.names = FALSE)
utils::write.csv(classic_summary,
                 file.path(OUT_DIR, "Supplementary_Figure3_Continuous_Transmission_Acceleration_Yearly_Summary.csv"),
                 row.names = FALSE)

utils::write.csv(table1_primary,
                 file.path(EVAL_OUT_DIR, "Table1_Primary_Summary.csv"), row.names = FALSE)
utils::write.csv(table1b_compartments,
                 file.path(EVAL_OUT_DIR, "Table1b_Compartment_Metrics.csv"), row.names = FALSE)
utils::write.csv(table2_per_year,
                 file.path(EVAL_OUT_DIR, "Table2_Per_Year_Detail.csv"), row.names = FALSE)
utils::write.csv(table2a_per_year_compartments,
                 file.path(EVAL_OUT_DIR, "Table2A_Per_Year_Detail_With_Compartments.csv"),
                 row.names = FALSE)
utils::write.csv(table2b_ta_vs_ot,
                 file.path(EVAL_OUT_DIR, "Table2B_Per_Year_TA_vs_OT_Comparison.csv"),
                 row.names = FALSE)

utils::write.csv(trigger_detail,
                 file.path(EVAL_OUT_DIR, "Trigger_Level_Detail.csv"), row.names = FALSE)
utils::write.csv(trigger_detail_aug_fig3,
                 file.path(EVAL_OUT_DIR, "Trigger_Level_Detail_with_Compartments.csv"), row.names = FALSE)

supp_S1 <- table2_per_year %>%
  dplyr::distinct(Year, Peak_Week, A1_start, A1_end, A2_start, A2_end)
utils::write.csv(supp_S1,
                 file.path(EVAL_OUT_DIR, "Supp_S1_Anchor_Boundaries.csv"), row.names = FALSE)

utils::write.csv(S2_table, file.path(EVAL_OUT_DIR, "Supp_S2_Sensitivity_A1_Window.csv"), row.names = FALSE)
utils::write.csv(S3_table, file.path(EVAL_OUT_DIR, "Supp_S3_Sensitivity_A2_Burden.csv"), row.names = FALSE)
utils::write.csv(S5_table, file.path(EVAL_OUT_DIR, "Supp_S5_Sensitivity_Year_Inclusion.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# C.19 PRINT FIGURE 3 TABULAR OUTPUTS
# -----------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("FIGURE 3A - CONSTANT TRANSMISSION ACCELERATION YEARLY SUMMARY\n")
cat("=================================================================\n")
safe_df_print(vaezi_summary, round_cols = c("Peak_DC"), digits = 3)

cat("\n=================================================================\n")
cat("SUPPLEMENTARY FIGURE 3 - CONTINUOUS TRANSMISSION ACCELERATION YEARLY SUMMARY\n")
cat("=================================================================\n")
safe_df_print(classic_summary, round_cols = c("Peak_DC"), digits = 3)

cat("\n=================================================================\n")
cat("TABLE 1 - PRIMARY SPECIFICATION SUMMARY (TRUE / FALSE ALARM)\n")
cat("Evaluable years: ", paste(EVALUABLE_YEARS, collapse = ", "), "\n", sep = "")
cat("Excluded years:  ", paste(EXCLUDED_YEARS,  collapse = ", "), "\n", sep = "")
cat("Framework:       two-anchor (T = A1 union A2)\n")
cat("Mean_Lead_Time_wks            : same-denominator HEADLINE metric\n")
cat("                                (years with no A1 true alarm contribute 0)\n")
cat("Mean_Lead_Time_wks_conditional: legacy conditional-on-firing diagnostic\n")
cat("=================================================================\n")
safe_df_print(table1_primary)

cat("\n=================================================================\n")
cat("TABLE 1b - LEAD-TIME COMPARTMENT METRICS (same-denominator)\n")
cat("Compartments (defined for True alarms only):\n")
cat("  Actionable: 4-8 wk      (operational sweet spot; matches A1)\n")
cat("  Reactive  : any other   (True alarms with lead < 4 wk; A2 only)\n")
cat("Triggers with lead >= 9 are FALSE alarms (not a separate compartment).\n")
cat("ALY    = Actionable Lead-Time Yield (same-denominator headline)\n")
cat("WP_wks = Warning Persistence        (same-denominator headline)\n")
cat("ALY_conditional, WP_wks_conditional: legacy conditional aggregates\n")
cat("=================================================================\n")
safe_df_print(table1b_compartments)

cat("\n=================================================================\n")
cat("TABLE 2A - PER-YEAR DETAIL WITH COMPARTMENT METRICS\n")
cat("Per-year compartment counts and per-year MLT, WP, ALY (same-\n")
cat("denominator headline columns) plus _conditional diagnostic versions,\n")
cat("alongside per-year PPV, TAM, and Lead Compartment of each year's\n")
cat("first True Alarm. One row per (Year, Detector).\n")
cat("=================================================================\n")
safe_df_print(table2a_per_year_compartments)

cat("\n=================================================================\n")
cat("TABLE 2B - PER-YEAR TA vs OT SIDE-BY-SIDE COMPARISON\n")
cat("Wide-format table mirroring Figure 3's per-year visual layout.\n")
cat("One row per Year. TA and OT shown side by side.\n")
cat("=================================================================\n")
safe_df_print(table2b_ta_vs_ot)

cat("\n=================================================================\n")
cat("SENSITIVITY S2 - A1 WINDOW (3-6 / 4-8 / 5-10)\n")
cat("=================================================================\n")
safe_df_print(S2_table)

cat("\n=================================================================\n")
cat("SENSITIVITY S3 - A2 BURDEN (0.60 / 0.70 / 0.80)\n")
cat("=================================================================\n")
safe_df_print(S3_table)

cat("\n=================================================================\n")
cat("SENSITIVITY S5 - YEAR INCLUSION\n")
cat("=================================================================\n")
safe_df_print(S5_table)


# -----------------------------------------------------------------------------
# C.20 FIGURE 3 LEGEND DRAFTS
# -----------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("FIGURE 3 - DRAFT FIGURE LEGENDS\n")
cat("=================================================================\n")

cat("\nFigure 3A:\n")
cat(paste0(
  "Fig. 3A | Constant Transmission Acceleration detection timing in Quezon City, 2013-2025. ",
  "Per-year panels show weekly dengue cases (grey bars), the rolling Outbreak Threshold derived ",
  "from week-specific historical baselines (black dashed line), and the Constant Transmission ",
  "Acceleration ratio R(t) calculated using a 4-week short-term average and a 26-week long-term ",
  "average with a 2-week guard (red line; right axis). The red dashed horizontal line indicates ",
  "the trigger threshold (eta_on = ", ETA_ON, "); the alarm-off threshold is eta_off = ", ETA_OFF,
  ". Trigger markers are shown ONLY for the FIRST trigger per detector per year that falls within ",
  "the Actionable Window (A1, peak-8 to peak-4 weeks); subsequent triggers and triggers outside ",
  "A1 are visually suppressed to remove over-triggering clutter (the full trigger inventory is ",
  "preserved in the underlying CSV outputs). Markers are shape-coded by detector: triangles for ",
  "the Constant Transmission Acceleration method (up-triangle) and squares for the Outbreak ",
  "Threshold. Two anchor bands are overlaid at the bottom of each panel: the wider amber band ",
  "labels the Epidemic Burden anchor (A2, 70% case-burden block); the narrower teal band labels ",
  "the Actionable Window (A1, peak-8..peak-4). A trigger is a True Alarm if it lies in A1 or A2; ",
  "otherwise it is a False Alarm. Per-panel badges (TA T:_ F:_ / OT T:_ F:_) report True and ",
  "False Alarm counts. Compartmentalized performance metrics (ALY, WP) are reported in Table 1b. ",
  "Years 2020 and 2021 are excluded due to COVID-19 surveillance disruption; year 2025 is ",
  "excluded as a methodologically special case (out-of-distribution relative to the historical ",
  "reference period used to derive framework parameters).\n"
))

cat("\nSupplementary Figure 3:\n")
cat(paste0(
  "Supplementary Fig. 3 | Continuous Transmission Acceleration detection timing in Quezon City, ",
  "2013-2025. Per-year panels show weekly dengue cases (grey bars), the rolling Outbreak Threshold ",
  "(black dashed line), and the Continuous Transmission Acceleration ratio R(t) calculated as the ",
  "3-week short-term average divided by the 12-week long-term average (blue line; right axis). ",
  "The dashed horizontal line indicates the trigger threshold (eta_on = ", ETA_ON_CLASSIC, "). ",
  "Triggers are classified under the same two-anchor rule used in Figure 3A. Anchor overlays ",
  "(Actionable Window in teal, Epidemic Burden in amber as horizontal bands), first-A1-trigger-",
  "only marker rendering, and exclusion conventions are identical to Figure 3A.\n"
))

cat("\nMethods note (paste into manuscript):\n")
cat(paste0(
  "Each weekly trigger from the three candidate detectors (Constant Transmission Acceleration, ",
  "Continuous Transmission Acceleration, and Outbreak Threshold) was classified as a True Alarm ",
  "or False Alarm using two pre-specified anchors and an explicit combination rule. The ",
  "Actionable Window anchor (A1) defined a trigger as actionable if it arrived 4-8 weeks before ",
  "the annual peak. The Epidemic Burden anchor (A2) defined the smallest contiguous block of ",
  "weeks containing the peak whose cumulative cases summed to >= 70% of the annual total. A ",
  "trigger at week t was a True Alarm if t was in A1 OR t was in A2; all other triggers were ",
  "False Alarms. Each True Alarm was additionally classified into one of two lead-time ",
  "compartments: Actionable (4-8 weeks; operational sweet spot, equivalent to A1) or Reactive ",
  "(any other True alarm; lead < 4 weeks). Compartmentalized performance was summarized by ",
  "8 metrics in three categories: Epidemic Burden & Alarm Accuracy (TAM = True-Alarm ",
  "Magnitude, Number of True Alarms, PPV, Sensitivity); Early Warning Timeliness (Mean Lead ",
  "Time, WP = Warning Persistence, ALY); and False Alarms (Number of False Alarms). ",
  "Sensitivity is computed as the proportion of evaluable seasons with at least one True ",
  "Alarm in the Actionable Window (A1). The three Early Warning Timeliness metrics (Mean ",
  "Lead Time, Warning Persistence, Actionable Lead-Time Yield) are computed with the SAME ",
  "DENOMINATOR per detector (count of evaluable years) by zero-coercing values for years ",
  "where the metric is not computable due to no qualifying triggers. Conditional-on-firing ",
  "diagnostic columns (suffix _conditional) are preserved in the per-year CSVs for reviewers ",
  "who prefer the 'mean over firing years only' semantics. Per-year panels render only the ",
  "first A1 trigger per detector to remove visual over-triggering clutter; the full trigger ",
  "inventory is preserved in the underlying CSVs. Years 2020 and 2021 were excluded due to ",
  "COVID-19 surveillance disruption. Year 2025 was excluded as a methodologically special case. ",
  "Sensitivity analyses on the A1 window, A2 burden fraction, and year-inclusion choices are ",
  "reported in Supplementary Tables S2, S3, S5. Pairwise statistical comparisons of TA versus ",
  "OT across all 8 metrics use the Wilcoxon signed-rank test (paired by year) evaluated PER ",
  "PAIRWISE COMPARISON at alpha = 0.05; each pair (Constant TA vs OT; Continuous TA vs OT; ",
  "Constant TA vs Continuous TA) is tested on its own without across-pair multiplicity ",
  "correction. Effect sizes are reported as median paired differences with 95% year-cluster ",
  "bootstrap confidence intervals.\n"
))


# -----------------------------------------------------------------------------
# C.21 FINAL FILE LIST
# -----------------------------------------------------------------------------
cat("\n=================================================================\n")
cat("ALL OUTPUTS SAVED\n")
cat("=================================================================\n")
cat("Main figures and per-year summaries in:\n  ", OUT_DIR, "\n", sep = "")
cat("\nEvaluation framework tables in:\n  ", EVAL_OUT_DIR, "\n", sep = "")
cat("\nFiles produced (evaluation framework):\n")
cat("  Table1_Primary_Summary.csv\n")
cat("  Table1b_Compartment_Metrics.csv\n")
cat("  Table2_Per_Year_Detail.csv\n")
cat("  Table2A_Per_Year_Detail_With_Compartments.csv\n")
cat("  Table2B_Per_Year_TA_vs_OT_Comparison.csv\n")
cat("  Trigger_Level_Detail.csv\n")
cat("  Trigger_Level_Detail_with_Compartments.csv\n")
cat("  Supp_S1_Anchor_Boundaries.csv\n")
cat("  Supp_S2_Sensitivity_A1_Window.csv\n")
cat("  Supp_S3_Sensitivity_A2_Burden.csv\n")
cat("  Supp_S5_Sensitivity_Year_Inclusion.csv\n")
cat("  Figure3_HeadToHead_ConstantTA_vs_OT.pdf\n")
cat("  Figure3_HeadToHead_ConstantTA_vs_OT.png\n")
cat("  Figure3_HeadToHead_ContinuousTA_vs_OT.pdf\n")
cat("  Figure3_HeadToHead_ContinuousTA_vs_OT.png\n")
cat("  Figure3_HeadToHead_ConstantTA_vs_ContinuousTA.pdf\n")
cat("  Figure3_HeadToHead_ConstantTA_vs_ContinuousTA.png\n")
cat("  Figure3_HeadToHead_Wilcoxon_Results.csv\n")

cat("\n=================================================================\n")
cat("END OF STAGE 3 - Quezon City analysis (Figures 2 and 3, plus tables)\n")
cat("=================================================================\n")

# =============================================================================
# END OF STAGE 3
# =============================================================================