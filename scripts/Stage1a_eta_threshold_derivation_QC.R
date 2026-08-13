# =============================================================================
# STAGE 1 - Empirical eta_ON / eta_OFF derivation (Quezon City)
# -----------------------------------------------------------------------------
# Empirical, multi-method derivation of the STA/LTA hysteresis thresholds
# (eta_ON, eta_OFF) for a Vaezi-style outbreak detector applied to weekly
# notified dengue case counts in Quezon City, Philippines.
#
# Computation excludes 2020, 2021 (COVID-19 surveillance disruption) and the
# truncated year 2025 to avoid biasing the null distribution and the
# operating-point calibration. The final adopted pair is the median across
# six independent statistical derivations (M1-M6).
#
# Inputs   : Dengue-Rainfall_Dataset.xlsx, sheet "QC Data"
#            Required columns: YR (ISO year), WN (ISO week), DC_QC (cases)
# Outputs  : Stage1_eta_thresholds_QC/  (4-panel figure, 6 CSV tables,
#            sourceable eta_thresholds_for_figure2.R, full analysis dataset)
#
# Repro    : set.seed(12345); R >= 4.1; UTF-8 locale recommended
# =============================================================================


# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------
REQUIRED_PACKAGES <- c(
  "readxl", "dplyr", "tidyr", "purrr", "pROC", "ggplot2",
  "cowplot", "zoo", "ISOweek", "scales", "tibble", "grid",
  "patchwork", "MASS", "viridisLite"
)

# -----------------------------------------------------------------------------
# 0. PROJECT BOOTSTRAP
# -----------------------------------------------------------------------------
# Portable path resolution, shared publication theme and dependency checking.
# Replaces the previous inline install.packages() loop and hard-coded
# "C:/Users/User/Desktop/..." path.
.bootstrap <- function() {
  root <- Sys.getenv("TA_PROJECT_ROOT", unset = "")
  if (!nzchar(root) || !dir.exists(root)) {
    d <- normalizePath(getwd(), winslash = "/")
    for (i in seq_len(6)) {
      if (dir.exists(file.path(d, "R")) && dir.exists(file.path(d, "scripts"))) {
        root <- d; break
      }
      p <- dirname(d); if (identical(p, d)) break; d <- p
    }
  }
  if (!nzchar(root)) {
    stop("Set the working directory to the project root before sourcing.",
         call. = FALSE)
  }
  source(file.path(root, "R", "00_config.R"))
}
.bootstrap()

require_packages(REQUIRED_PACKAGES, purpose = "Stage 1")
invisible(lapply(REQUIRED_PACKAGES, function(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))))

source(file.path(DIR_R, "01_publication_theme.R"), local = TRUE)

set.seed(GLOBAL_SEED)
options(scipen = 999)


# -----------------------------------------------------------------------------
# 1. USER SETTINGS
# -----------------------------------------------------------------------------
PATH       <- DATA_FILE
SHEET_NAME <- SHEET_QC

OUT_DIR <- file.path(DIR_OUTPUT, "Stage1_eta_thresholds_QC")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Years excluded from all threshold computations
EXCLUDED_YEARS <- c(2020L, 2021L, 2025L)

# STA/LTA window configuration (weeks)
STA_WIN <- 4L     # short-term average window
LTA_WIN <- 26L    # long-term average window
GUARD   <- 2L     # guard band between STA and LTA windows

# Derivation parameters
BASELINE_FRAC <- 0.30                    # baseline = weeks below 30% of annual peak
ARL_TARGETS   <- c(6, 12, 26)            # target ARL_0 values (weeks)
GRID_ETA_ON   <- seq(1.10, 2.50, by = 0.05)
GRID_ETA_OFF  <- seq(0.40, 1.50, by = 0.05)
BOOT_N        <- 2000                    # parametric-bootstrap replicates
MIN_OFF_RESET <- 8L                      # consecutive OFF weeks before LTA reset

FINAL_ROUND_DIGITS <- 1
ADOPTION_RULE      <- "raw_median"       # "raw_median" or "rounded_median"


# -----------------------------------------------------------------------------
# 2. HELPERS
# -----------------------------------------------------------------------------
safe_quantile <- function(x, probs, na.rm = TRUE) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  suppressWarnings(as.numeric(stats::quantile(
    x, probs = probs, na.rm = na.rm, names = FALSE
  )))
}

# safe_pdf_device(), PUB_FAMILY and save_pub() now come from
# R/01_publication_theme.R so that device selection, font fallback and export
# dimensions are identical across all five stages.
base_family_global <- PUB_FAMILY

# Thin wrapper retained for call-site compatibility. Writes vector PDF +
# 600 dpi PNG in one call and clamps the canvas to the journal's print area.
save_figure <- function(stem, plot, width, height, dpi = 600) {
  save_pub(stem = sub("\\.(pdf|png)$", "", stem),
           plot = plot, width = width, height = height,
           dir = OUT_DIR, dpi = dpi)
}


# -----------------------------------------------------------------------------
# 3. LOAD AND PREPARE DATA
# -----------------------------------------------------------------------------
if (!file.exists(PATH)) stop("Data file not found at:\n", PATH)

df_raw <- readxl::read_excel(PATH, sheet = SHEET_NAME)

required_cols <- c("YR", "WN", "DC_QC")
missing_cols  <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df_all <- df_raw %>%
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

df <- df_all %>%
  dplyr::filter(!YR %in% EXCLUDED_YEARS) %>%
  dplyr::arrange(Date)

if (nrow(df) == 0) stop("No valid rows remain after excluding selected years.")

included_years <- sort(unique(df$YR))

cat("\n--- Computation year filter ---\n")
cat("Excluded years: ", paste(EXCLUDED_YEARS, collapse = ", "), "\n")
cat("Included years: ", paste(included_years, collapse = ", "), "\n")
cat("Rows retained:  ", nrow(df), "\n")


# -----------------------------------------------------------------------------
# 4. COMPUTE STA/LTA RATIO SERIES
# -----------------------------------------------------------------------------
dc         <- df$DC_QC
sta_series <- zoo::rollmean(dc, STA_WIN, fill = NA, align = "right")

n          <- length(dc)
lta_series <- rep(NA_real_, n)

for (t in seq_len(n)) {
  end_idx   <- t - STA_WIN - GUARD
  start_idx <- end_idx - LTA_WIN + 1
  if (start_idx < 1 || end_idx < 1) next
  vals <- dc[start_idx:end_idx]
  if (all(is.na(vals))) next
  lta_series[t] <- mean(vals, na.rm = TRUE)
}

ratio_series <- ifelse(
  !is.na(lta_series) & lta_series > 0,
  sta_series / lta_series,
  NA_real_
)

df$sta_val   <- sta_series
df$lta_val   <- lta_series
df$ratio_val <- ratio_series


# -----------------------------------------------------------------------------
# 5. DEFINE BASELINE WEEKS
# -----------------------------------------------------------------------------
annual_max <- df %>%
  dplyr::group_by(YR) %>%
  dplyr::summarise(ymax = max(DC_QC, na.rm = TRUE), .groups = "drop")

df <- df %>%
  dplyr::left_join(annual_max, by = "YR") %>%
  dplyr::mutate(is_baseline = is.finite(DC_QC) & DC_QC < BASELINE_FRAC * ymax)

baseline_ratios <- df$ratio_val[df$is_baseline & is.finite(df$ratio_val)]
baseline_counts <- df$DC_QC[df$is_baseline & is.finite(df$DC_QC)]

cat("\n--- Baseline-week summary ---\n")
cat("Total retained weeks:   ", nrow(df), "\n")
cat("Baseline weeks:         ", sum(df$is_baseline, na.rm = TRUE), "\n")
cat("Baseline ratios valid:  ", length(baseline_ratios), "\n")
cat("Mean baseline ratio:    ", round(mean(baseline_ratios, na.rm = TRUE), 3), "\n")
cat("SD baseline ratio:      ", round(sd(baseline_ratios,   na.rm = TRUE), 3), "\n")


# -----------------------------------------------------------------------------
# 6. METHOD 1 - BASELINE-RATIO ASYMMETRIC QUANTILES
# -----------------------------------------------------------------------------
cat("\n[M1] Baseline-ratio asymmetric quantiles\n")

m1_eta_on  <- safe_quantile(baseline_ratios, 0.90)
m1_eta_off <- safe_quantile(baseline_ratios, 0.50)

cat("  eta_ON  =", round(m1_eta_on,  3), "\n")
cat("  eta_OFF =", round(m1_eta_off, 3), "\n")


# -----------------------------------------------------------------------------
# 7. METHOD 2 - HIGH-SENSITIVITY ROC OPERATING POINTS
# -----------------------------------------------------------------------------
cat("\n[M2] High-sensitivity ROC operating points\n")

outcome_label <- df %>%
  dplyr::group_by(YR) %>%
  dplyr::mutate(
    thr_elev = stats::quantile(DC_QC, 0.75, na.rm = TRUE, names = FALSE),
    is_elev  = as.integer(DC_QC > thr_elev)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::pull(is_elev)

valid_ix <- is.finite(df$ratio_val) & !is.na(outcome_label)

m2_eta_on  <- NA_real_
m2_eta_off <- NA_real_
m2_sens    <- NA_real_
m2_spec    <- NA_real_
m2_auc     <- NA_real_
roc_m2     <- NULL

if (sum(valid_ix) > 30 && length(unique(outcome_label[valid_ix])) == 2) {
  roc_m2 <- tryCatch(
    pROC::roc(
      response  = outcome_label[valid_ix],
      predictor = df$ratio_val[valid_ix],
      ci        = FALSE,
      quiet     = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(roc_m2)) {
    on_c <- tryCatch(
      pROC::coords(
        roc_m2,
        x         = 0.95,
        input     = "sensitivity",
        ret       = c("threshold", "sensitivity", "specificity"),
        transpose = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(on_c)) {
      m2_eta_on <- as.numeric(on_c$threshold[1])
      m2_sens   <- as.numeric(on_c$sensitivity[1])
      m2_spec   <- as.numeric(on_c$specificity[1])
      m2_auc    <- as.numeric(roc_m2$auc)
    }
    
    off_c <- tryCatch(
      pROC::coords(
        roc_m2,
        x         = 0.99,
        input     = "sensitivity",
        ret       = "threshold",
        transpose = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(off_c)) m2_eta_off <- as.numeric(off_c$threshold[1])
  }
}

cat("  eta_ON  =", round(m2_eta_on,  3), "\n")
cat("  eta_OFF =", round(m2_eta_off, 3), "\n")
cat("  AUC     =", round(m2_auc,     3), "\n")


# -----------------------------------------------------------------------------
# 8. METHOD 3 - ARL CALIBRATION
# -----------------------------------------------------------------------------
cat("\n[M3] Sensitivity-oriented ARL calibration\n")

eta_grid_arl <- seq(1.00, 3.50, by = 0.01)
arl0_at_eta  <- sapply(eta_grid_arl, function(e) {
  p <- mean(baseline_ratios > e, na.rm = TRUE)
  if (p <= 0) Inf else 1 / p
})

m3_df <- dplyr::bind_rows(lapply(ARL_TARGETS, function(tgt) {
  idx <- which.min(abs(arl0_at_eta - tgt))
  data.frame(
    Target_ARL0   = tgt,
    eta_ON        = eta_grid_arl[idx],
    Achieved_ARL0 = arl0_at_eta[idx],
    stringsAsFactors = FALSE
  )
}))

m3_eta_on_primary <- m3_df$eta_ON[m3_df$Target_ARL0 == 12]

print(m3_df)
cat("  eta_ON primary =", round(m3_eta_on_primary, 3), "\n")


# -----------------------------------------------------------------------------
# 9. METHOD 4 - LOYO-CV ASYMMETRIC UTILITY
# -----------------------------------------------------------------------------
cat("\n[M4] LOYO-CV grid search with asymmetric utility\n")

df <- df %>%
  dplyr::group_by(YR) %>%
  dplyr::mutate(
    thr_outbreak_yr = stats::quantile(DC_QC, 0.90, na.rm = TRUE, names = FALSE),
    is_outbreak_yr  = as.integer(!is.na(DC_QC) & DC_QC > thr_outbreak_yr)
  ) %>%
  dplyr::ungroup()

run_vaezi <- function(dc_vec, eta_on, eta_off) {
  nn         <- length(dc_vec)
  trig       <- rep(FALSE, nn)
  is_on      <- FALSE
  frozen_lta <- NA_real_
  consec_off <- 0L
  min_t      <- STA_WIN + GUARD + LTA_WIN
  
  for (t in seq_len(nn)) {
    if (!is_on) consec_off <- consec_off + 1L else consec_off <- 0L
    if (!is_on && consec_off >= MIN_OFF_RESET) frozen_lta <- NA_real_
    if (t < min_t) next
    
    sta <- mean(dc_vec[(t - STA_WIN + 1):t], na.rm = TRUE)
    
    if (!is_on || is.na(frozen_lta)) {
      idx <- (t - STA_WIN - GUARD - LTA_WIN + 1):(t - STA_WIN - GUARD)
      if (length(idx) == LTA_WIN && min(idx) > 0) {
        frozen_lta <- mean(dc_vec[idx], na.rm = TRUE)
      } else {
        frozen_lta <- NA_real_
      }
    }
    
    R_t <- if (!is.na(frozen_lta) && frozen_lta > 0 && !is.na(sta)) {
      sta / frozen_lta
    } else {
      NA_real_
    }
    
    if (!is_on && !is.na(R_t) && R_t >= eta_on) {
      is_on <- TRUE
      consec_off <- 0L
    }
    if (is_on && !is.na(R_t) && R_t < eta_off) {
      is_on <- FALSE
      frozen_lta <- NA_real_
    }
    
    trig[t] <- is_on
  }
  
  as.integer(trig)
}

compute_asymmetric_utility <- function(pred, truth, dc_year) {
  valid <- !is.na(pred) & !is.na(truth)
  if (sum(valid) == 0) return(NA_real_)
  
  pred  <- pred[valid]
  truth <- truth[valid]
  
  outbreak_weeks <- which(truth == 1L)
  baseline_weeks <- which(truth == 0L)
  
  first_trigger_wk  <- which(pred == 1L)[1]
  first_outbreak_wk <- if (length(outbreak_weeks) > 0) outbreak_weeks[1] else NA_integer_
  
  if (is.na(first_trigger_wk) || is.na(first_outbreak_wk)) {
    lead_norm <- -1
  } else {
    lead_wks  <- first_outbreak_wk - first_trigger_wk
    lead_norm <- pmin(pmax(lead_wks / 12, -1), 1)
  }
  
  recall  <- if (length(outbreak_weeks) > 0) mean(pred[outbreak_weeks] == 1L) else 0
  if (is.na(recall)) recall <- 0
  fa_rate <- if (length(baseline_weeks) > 0) mean(pred[baseline_weeks] == 1L) else 0
  
  0.5 * lead_norm + 0.4 * recall - 0.1 * fa_rate
}

years_all <- sort(unique(df$YR))

grid_combos <- expand.grid(
  eta_on  = GRID_ETA_ON,
  eta_off = GRID_ETA_OFF,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(eta_off < eta_on)

cat("  Grid size:", nrow(grid_combos), "combinations x",
    length(years_all), "included years\n")

cv_results <- purrr::map_dfr(seq_len(nrow(grid_combos)), function(i) {
  e_on  <- grid_combos$eta_on[i]
  e_off <- grid_combos$eta_off[i]
  
  pred_full <- run_vaezi(df$DC_QC, e_on, e_off)
  
  per_year <- sapply(years_all, function(yr) {
    test_idx <- which(df$YR == yr)
    if (length(test_idx) == 0) return(NA_real_)
    compute_asymmetric_utility(
      pred    = pred_full[test_idx],
      truth   = df$is_outbreak_yr[test_idx],
      dc_year = df$DC_QC[test_idx]
    )
  })
  
  data.frame(
    eta_on       = e_on,
    eta_off      = e_off,
    mean_utility = mean(per_year, na.rm = TRUE),
    sd_utility   = sd(per_year,   na.rm = TRUE),
    n_years      = sum(!is.na(per_year)),
    stringsAsFactors = FALSE
  )
})

best_cv <- cv_results %>%
  dplyr::arrange(dplyr::desc(mean_utility)) %>%
  dplyr::slice(1)

m4_eta_on  <- best_cv$eta_on
m4_eta_off <- best_cv$eta_off
m4_utility <- best_cv$mean_utility

cat("  Best CV utility =", round(m4_utility, 3), "\n")
cat("  eta_ON          =", round(m4_eta_on,  3), "\n")
cat("  eta_OFF         =", round(m4_eta_off, 3), "\n")


# -----------------------------------------------------------------------------
# 10. METHOD 5 - NEGATIVE-BINOMIAL THEORETICAL BOUNDS
# -----------------------------------------------------------------------------
cat("\n[M5] Negative-binomial theoretical bounds\n")

m5_eta_on   <- NA_real_
m5_eta_off  <- NA_real_
m5_mu       <- NA_real_
m5_theta    <- NA_real_
m5_sd_ratio <- NA_real_

if (length(baseline_counts) > 30 && mean(baseline_counts, na.rm = TRUE) > 0) {
  nb_fit <- tryCatch(
    MASS::glm.nb(baseline_counts ~ 1),
    error = function(e) NULL
  )
  
  if (!is.null(nb_fit)) {
    m5_mu    <- as.numeric(exp(stats::coef(nb_fit)))
    m5_theta <- nb_fit$theta
    
    var_count <- m5_mu + m5_mu^2 / m5_theta
    var_sta   <- var_count / STA_WIN
    var_lta   <- var_count / LTA_WIN
    
    m5_sd_ratio <- sqrt((var_sta + var_lta) / m5_mu^2)
    m5_eta_on   <- 1 + m5_sd_ratio
    m5_eta_off  <- 1 - m5_sd_ratio
    if (m5_eta_off < 0) m5_eta_off <- 0.01
  }
}

cat("  mu       =", round(m5_mu,      3), "\n")
cat("  theta    =", round(m5_theta,   3), "\n")
cat("  eta_ON   =", round(m5_eta_on,  3), "\n")
cat("  eta_OFF  =", round(m5_eta_off, 3), "\n")


# -----------------------------------------------------------------------------
# 11. METHOD 6 - PARAMETRIC BOOTSTRAP OF THE NULL STA/LTA RATIO
# -----------------------------------------------------------------------------
cat("\n[M6] Parametric bootstrap of null STA/LTA ratio\n")

m6_eta_on  <- NA_real_
m6_eta_off <- NA_real_

if (!is.na(m5_mu) && !is.na(m5_theta)) {
  bootstrap_ratios <- numeric(BOOT_N)
  total_win        <- STA_WIN + GUARD + LTA_WIN
  
  for (b in seq_len(BOOT_N)) {
    sim                  <- MASS::rnegbin(total_win, mu = m5_mu, theta = m5_theta)
    lta_b                <- mean(sim[1:LTA_WIN])
    sta_b                <- mean(sim[(LTA_WIN + GUARD + 1):total_win])
    bootstrap_ratios[b]  <- if (lta_b > 0) sta_b / lta_b else NA_real_
  }
  
  bootstrap_ratios <- bootstrap_ratios[is.finite(bootstrap_ratios)]
  m6_eta_on        <- safe_quantile(bootstrap_ratios, 0.90)
  m6_eta_off       <- safe_quantile(bootstrap_ratios, 0.25)
}

cat("  eta_ON  =", round(m6_eta_on,  3), "\n")
cat("  eta_OFF =", round(m6_eta_off, 3), "\n")


# -----------------------------------------------------------------------------
# 12. CONSOLIDATE METHOD-LEVEL RESULTS
# -----------------------------------------------------------------------------
results_summary <- data.frame(
  Method = c(
    "M1: Baseline quantile (90/50 pct)",
    "M2: ROC high-sensitivity (95%/99%)",
    "M3: ARL0 = 12-week calibration",
    "M4: LOYO-CV asymmetric utility",
    "M5: NegBin theoretical (+/- 1 sigma)",
    "M6: Parametric bootstrap (90/25 pct)"
  ),
  eta_ON = c(
    m1_eta_on, m2_eta_on, m3_eta_on_primary,
    m4_eta_on, m5_eta_on, m6_eta_on
  ),
  eta_OFF = c(
    m1_eta_off, m2_eta_off, NA_real_,
    m4_eta_off, m5_eta_off, m6_eta_off
  ),
  stringsAsFactors = FALSE
)

results_summary$eta_ON  <- round(results_summary$eta_ON,  3)
results_summary$eta_OFF <- round(results_summary$eta_OFF, 3)

cat("\n============================================================\n")
cat("SUMMARY OF THRESHOLD DERIVATIONS\n")
cat("Years excluded from computation:", paste(EXCLUDED_YEARS, collapse = ", "), "\n")
cat("============================================================\n")
print(results_summary, row.names = FALSE, na.print = "NA")


# -----------------------------------------------------------------------------
# 13. FINAL DECISION (MEDIAN ACROSS METHODS)
# -----------------------------------------------------------------------------
derived_on_vec  <- c(m1_eta_on, m2_eta_on, m3_eta_on_primary,
                     m4_eta_on, m5_eta_on, m6_eta_on)
derived_off_vec <- c(m1_eta_off, m2_eta_off, m4_eta_off,
                     m5_eta_off, m6_eta_off)

ETA_ON_FINAL_RAW      <- stats::median(derived_on_vec,  na.rm = TRUE)
ETA_OFF_FINAL_RAW     <- stats::median(derived_off_vec, na.rm = TRUE)
ETA_ON_FINAL_ROUNDED  <- round(ETA_ON_FINAL_RAW,  FINAL_ROUND_DIGITS)
ETA_OFF_FINAL_ROUNDED <- round(ETA_OFF_FINAL_RAW, FINAL_ROUND_DIGITS)

if (ADOPTION_RULE == "raw_median") {
  ETA_ON_ADOPTED  <- ETA_ON_FINAL_RAW
  ETA_OFF_ADOPTED <- ETA_OFF_FINAL_RAW
  adoption_note   <- "Using raw median across methods"
} else if (ADOPTION_RULE == "rounded_median") {
  ETA_ON_ADOPTED  <- ETA_ON_FINAL_ROUNDED
  ETA_OFF_ADOPTED <- ETA_OFF_FINAL_ROUNDED
  adoption_note   <- "Using 1-decimal rounded median across methods"
} else {
  stop("ADOPTION_RULE must be 'raw_median' or 'rounded_median'")
}

ETA_ON_FINAL_RANGE  <- c(min(derived_on_vec,  na.rm = TRUE),
                         max(derived_on_vec,  na.rm = TRUE))
ETA_OFF_FINAL_RANGE <- c(min(derived_off_vec, na.rm = TRUE),
                         max(derived_off_vec, na.rm = TRUE))

cat("\n============================================================\n")
cat("FINAL STA/LTA HYSTERESIS THRESHOLD DECISION\n")
cat("============================================================\n")
cat("Excluded years: ", paste(EXCLUDED_YEARS, collapse = ", "), "\n")
cat("Included years: ", paste(included_years, collapse = ", "), "\n")
cat(sprintf("eta_ON_RAW       = %.4f\n", ETA_ON_FINAL_RAW))
cat(sprintf("eta_OFF_RAW      = %.4f\n", ETA_OFF_FINAL_RAW))
cat(sprintf("eta_ON_ROUNDED   = %.1f\n", ETA_ON_FINAL_ROUNDED))
cat(sprintf("eta_OFF_ROUNDED  = %.1f\n", ETA_OFF_FINAL_ROUNDED))
cat(sprintf("ADOPTION RULE    = %s\n",   ADOPTION_RULE))
cat(sprintf("ADOPTED eta_ON   = %.4f\n", ETA_ON_ADOPTED))
cat(sprintf("ADOPTED eta_OFF  = %.4f\n", ETA_OFF_ADOPTED))
cat("Rationale: ", adoption_note, "\n")

final_decision <- data.frame(
  excluded_years        = paste(EXCLUDED_YEARS, collapse = ", "),
  included_years        = paste(included_years, collapse = ", "),
  eta_ON_FINAL_raw      = ETA_ON_FINAL_RAW,
  eta_OFF_FINAL_raw     = ETA_OFF_FINAL_RAW,
  eta_ON_FINAL_rounded  = ETA_ON_FINAL_ROUNDED,
  eta_OFF_FINAL_rounded = ETA_OFF_FINAL_ROUNDED,
  adoption_rule         = ADOPTION_RULE,
  eta_ON_ADOPTED        = ETA_ON_ADOPTED,
  eta_OFF_ADOPTED       = ETA_OFF_ADOPTED,
  eta_ON_method_min     = ETA_ON_FINAL_RANGE[1],
  eta_ON_method_max     = ETA_ON_FINAL_RANGE[2],
  eta_OFF_method_min    = ETA_OFF_FINAL_RANGE[1],
  eta_OFF_method_max    = ETA_OFF_FINAL_RANGE[2],
  m1_eta_on  = m1_eta_on,  m1_eta_off = m1_eta_off,
  m2_eta_on  = m2_eta_on,  m2_eta_off = m2_eta_off,
  m3_eta_on  = m3_eta_on_primary,
  m4_eta_on  = m4_eta_on,  m4_eta_off = m4_eta_off,
  m5_eta_on  = m5_eta_on,  m5_eta_off = m5_eta_off,
  m6_eta_on  = m6_eta_on,  m6_eta_off = m6_eta_off,
  selection_rule = paste0(
    "Thresholds derived after excluding 2020, 2021, and 2025 from all ",
    "computations. Median across six derivation methods. Adopted pair ",
    "defined by ADOPTION_RULE = '", ADOPTION_RULE, "'."
  ),
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 14. WRITE SOURCEABLE R FILE FOR DOWNSTREAM USE
# -----------------------------------------------------------------------------
eta_source_path <- file.path(OUT_DIR, "eta_thresholds_for_figure2.R")

# Also written under the canonical name that R/00_config.R::load_eta_thresholds()
# looks for, rather than Stage 3 hard-coding eta_ON / eta_OFF directly.
# file, so re-deriving the thresholds here had no downstream effect.
eta_canonical_path <- file.path(OUT_DIR, "eta_thresholds_derived.R")

writeLines(c(
  "# ============================================================================",
  "# eta_thresholds_for_figure2.R",
  "# Auto-generated by Stage1_eta_threshold_derivation_QC.R",
  "# Years excluded from computation: 2020, 2021, 2025",
  "# ============================================================================",
  "",
  paste0("ETA_ON  <- ", sprintf("%.6f", ETA_ON_ADOPTED)),
  paste0("ETA_OFF <- ", sprintf("%.6f", ETA_OFF_ADOPTED)),
  "",
  paste0("ETA_ON_RAW      <- ", sprintf("%.6f", ETA_ON_FINAL_RAW)),
  paste0("ETA_OFF_RAW     <- ", sprintf("%.6f", ETA_OFF_FINAL_RAW)),
  paste0("ETA_ON_ROUNDED  <- ", sprintf("%.1f", ETA_ON_FINAL_ROUNDED)),
  paste0("ETA_OFF_ROUNDED <- ", sprintf("%.1f", ETA_OFF_FINAL_ROUNDED)),
  "",
  paste0("EXCLUDED_YEARS <- c(", paste(EXCLUDED_YEARS, collapse = ", "), ")"),
  paste0("INCLUDED_YEARS <- c(", paste(included_years, collapse = ", "), ")"),
  "",
  "# End of file"
), con = eta_source_path)

# Canonical copy consumed by Stage 3 via load_eta_thresholds(). Carries the
# *_ADOPTED names that the loader validates.
writeLines(c(
  "# Auto-generated by Stage1_eta_threshold_derivation_QC.R -- do not edit.",
  paste0("ETA_ON_ADOPTED  <- ", sprintf("%.6f", ETA_ON_ADOPTED)),
  paste0("ETA_OFF_ADOPTED <- ", sprintf("%.6f", ETA_OFF_ADOPTED))
), con = eta_canonical_path)

cat("\nSourceable R file written: ", eta_source_path,    "\n", sep = "")
cat("Canonical eta file written: ", eta_canonical_path, "\n", sep = "")


# -----------------------------------------------------------------------------
# 15. PUBLICATION THEME
# -----------------------------------------------------------------------------
# The shared theme_pub() version in
# R/01_publication_theme.R now supplies it, so every stage renders with
# identical typography, gridlines and legend styling.


# -----------------------------------------------------------------------------
# 16. PANEL A - NULL STA/LTA RATIO DISTRIBUTION
# -----------------------------------------------------------------------------
hist_df <- data.frame(ratio = baseline_ratios)

panel_a <- ggplot2::ggplot(hist_df, ggplot2::aes(x = ratio)) +
  ggplot2::geom_histogram(bins = 40, fill = "#6BAED6",
                          colour = "white", linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = m1_eta_on,  colour = "#D7301F",
                      linetype = "dashed", linewidth = 0.7) +
  ggplot2::geom_vline(xintercept = m1_eta_off, colour = "#FC8D59",
                      linetype = "dashed", linewidth = 0.7) +
  ggplot2::geom_vline(xintercept = ETA_ON_ADOPTED,  colour = "#084594",
                      linetype = "solid",  linewidth = 1.0) +
  ggplot2::geom_vline(xintercept = ETA_OFF_ADOPTED, colour = "#2171B5",
                      linetype = "solid",  linewidth = 1.0) +
  ggplot2::labs(
    title    = "a    Null STA/LTA ratio distribution",
    subtitle = "Computation excludes 2020, 2021, and 2025",
    x        = "STA/LTA ratio during baseline weeks",
    y        = "Count"
  ) +
  theme_pub() +
  ggplot2::theme(
    plot.subtitle = ggplot2::element_text(size = PUB_SUBTITLE, colour = "grey30",
                                          margin = ggplot2::margin(b = 4))
  )


# -----------------------------------------------------------------------------
# 17. PANEL B - ROC CURVE WITH HIGH-SENSITIVITY OPERATING POINTS
# -----------------------------------------------------------------------------
if (!is.null(roc_m2)) {
  roc_plot_df <- data.frame(
    fpr = 1 - roc_m2$specificities,
    tpr = roc_m2$sensitivities
  )
  
  op_on <- tryCatch(
    pROC::coords(roc_m2, x = 0.95, input = "sensitivity",
                 ret = c("sensitivity", "specificity"),
                 transpose = FALSE),
    error = function(e) NULL
  )
  
  op_off <- tryCatch(
    pROC::coords(roc_m2, x = 0.99, input = "sensitivity",
                 ret = c("sensitivity", "specificity"),
                 transpose = FALSE),
    error = function(e) NULL
  )
  
  panel_b <- ggplot2::ggplot(roc_plot_df, ggplot2::aes(fpr, tpr)) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey50",
                         linewidth = 0.4) +
    ggplot2::geom_line(colour = "#2171B5", linewidth = 1.0)
  
  if (!is.null(op_on)) {
    panel_b <- panel_b +
      ggplot2::annotate("point",
                        x = 1 - as.numeric(op_on$specificity[1]),
                        y = as.numeric(op_on$sensitivity[1]),
                        colour = "#D7301F", size = pub_text_size(PUB_ANNOT))
  }
  
  if (!is.null(op_off)) {
    panel_b <- panel_b +
      ggplot2::annotate("point",
                        x = 1 - as.numeric(op_off$specificity[1]),
                        y = as.numeric(op_off$sensitivity[1]),
                        colour = "#FC8D59", size = pub_text_size(PUB_ANNOT))
  }
  
  panel_b <- panel_b +
    ggplot2::labs(
      title    = sprintf("b    ROC curve of STA/LTA ratio (AUC = %.2f)", m2_auc),
      subtitle = "High-sensitivity operating points; computation excludes 2020, 2021, and 2025",
      x        = "False positive rate",
      y        = "True positive rate"
    ) +
    ggplot2::coord_equal() +
    theme_pub() +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = PUB_SUBTITLE, colour = "grey30",
                                            margin = ggplot2::margin(b = 4))
    )
} else {
  panel_b <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5,
                      label = "ROC unavailable", size = pub_text_size(PUB_ANNOT + 1.5)) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::labs(title = "ROC curve of STA/LTA ratio") +
    theme_pub()
}


# -----------------------------------------------------------------------------
# 18. PANEL C - LOYO-CV ASYMMETRIC-UTILITY SURFACE
# -----------------------------------------------------------------------------
cv_heat <- cv_results %>%
  dplyr::mutate(
    eta_on_f  = factor(round(eta_on,  2)),
    eta_off_f = factor(round(eta_off, 2))
  )

best_point <- data.frame(
  eta_on_f  = factor(round(m4_eta_on,  2), levels = levels(cv_heat$eta_on_f)),
  eta_off_f = factor(round(m4_eta_off, 2), levels = levels(cv_heat$eta_off_f))
)

final_point <- data.frame(
  eta_on_f  = factor(round(ETA_ON_ADOPTED,  2), levels = levels(cv_heat$eta_on_f)),
  eta_off_f = factor(round(ETA_OFF_ADOPTED, 2), levels = levels(cv_heat$eta_off_f))
)

panel_c <- ggplot2::ggplot(
  cv_heat,
  ggplot2::aes(x = eta_on_f, y = eta_off_f, fill = mean_utility)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
  ggplot2::scale_fill_gradient2(
    low  = "#F7F7F7", mid = "#DEEBF7", high = "#084594",
    midpoint = 0, na.value = "#F2F2F2", name = "Mean utility"
  ) +
  ggplot2::geom_point(
    data = best_point,
    ggplot2::aes(x = eta_on_f, y = eta_off_f),
    inherit.aes = FALSE, shape = 1, size = 5, stroke = 1.2,
    colour = "#D7301F"
  ) +
  ggplot2::geom_point(
    data = final_point,
    ggplot2::aes(x = eta_on_f, y = eta_off_f),
    inherit.aes = FALSE, shape = 8, size = 4, stroke = 1.3,
    colour = "#084594"
  ) +
  ggplot2::labs(
    title    = "c    LOYO-CV asymmetric-utility surface",
    subtitle = "Red = M4 optimum; blue star = adopted pair; excludes 2020, 2021, and 2025",
    x        = "eta_ON",
    y        = "eta_OFF"
  ) +
  theme_pub() +
  ggplot2::theme(
    plot.subtitle = ggplot2::element_text(size = PUB_SUBTITLE, colour = "grey30",
                                          margin = ggplot2::margin(b = 4)),
    axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y   = ggplot2::element_text(size = PUB_AXIS_TXT)
  )


# -----------------------------------------------------------------------------
# 19. PANEL D - METHOD-LEVEL THRESHOLD COMPARISON
# -----------------------------------------------------------------------------
methods_df <- results_summary %>%
  dplyr::mutate(
    Method_short = c(
      "M1 quantile", "M2 ROC-sens", "M3 ARL=12",
      "M4 CV-util",  "M5 NegBin",   "M6 Bootstrap"
    )
  )

final_row <- data.frame(
  Method       = "ADOPTED",
  eta_ON       = ETA_ON_ADOPTED,
  eta_OFF      = ETA_OFF_ADOPTED,
  Method_short = "ADOPTED",
  stringsAsFactors = FALSE
)

methods_df <- dplyr::bind_rows(methods_df, final_row)

methods_df$Method_short <- factor(
  methods_df$Method_short,
  levels = c("M1 quantile", "M2 ROC-sens", "M3 ARL=12",
             "M4 CV-util",  "M5 NegBin",   "M6 Bootstrap", "ADOPTED")
)

methods_long <- methods_df %>%
  dplyr::select(Method_short, eta_ON, eta_OFF) %>%
  tidyr::pivot_longer(cols = c(eta_ON, eta_OFF),
                      names_to = "Parameter", values_to = "Value") %>%
  dplyr::filter(!is.na(Value))

panel_d <- ggplot2::ggplot(methods_df, ggplot2::aes(y = Method_short)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = eta_OFF, xend = eta_ON,
                 y = Method_short, yend = Method_short),
    colour = "grey60", linewidth = 0.6, na.rm = TRUE
  ) +
  ggplot2::geom_point(
    data = methods_long,
    ggplot2::aes(x = Value, y = Method_short,
                 colour = Parameter, shape = Parameter),
    size = 3.2, na.rm = TRUE
  ) +
  ggplot2::geom_text(
    data = methods_long,
    ggplot2::aes(x = Value, y = Method_short,
                 label = sprintf("%.2f", Value)),
    vjust = -1.1, size = 2.7, na.rm = TRUE
  ) +
  ggplot2::geom_vline(xintercept = ETA_ON_ADOPTED,  linetype = "solid",
                      colour = "#084594", linewidth = 0.5) +
  ggplot2::geom_vline(xintercept = ETA_OFF_ADOPTED, linetype = "solid",
                      colour = "#2171B5", linewidth = 0.5) +
  ggplot2::scale_colour_manual(
    values = c("eta_ON" = "#D7301F", "eta_OFF" = "#2171B5"),
    labels = c("eta_ON", "eta_OFF"), name = NULL
  ) +
  ggplot2::scale_shape_manual(
    values = c("eta_ON" = 16, "eta_OFF" = 17),
    labels = c("eta_ON", "eta_OFF"), name = NULL
  ) +
  ggplot2::labs(
    title    = "d    Derived thresholds across six methods",
    subtitle = sprintf(
      "Adopted eta_ON = %.2f, eta_OFF = %.2f; excludes 2020, 2021, and 2025",
      ETA_ON_ADOPTED, ETA_OFF_ADOPTED
    ),
    x = "Threshold value", y = NULL
  ) +
  theme_pub() +
  ggplot2::theme(
    legend.position = "bottom",
    plot.subtitle   = ggplot2::element_text(size = PUB_SUBTITLE, colour = "grey30",
                                            margin = ggplot2::margin(b = 4))
  )


# -----------------------------------------------------------------------------
# 20. ASSEMBLE AND SAVE FIGURE
# -----------------------------------------------------------------------------
fig_top    <- panel_a | panel_b
fig_bottom <- panel_c | panel_d

# Panel letters are INLINE in each panel's title ("a    Null STA/LTA ratio
# distribution", etc.), matching the house style. pub_tag_layout() is therefore
# not used here -- it would stack a second letter above the title.
fig_combined <- fig_top / fig_bottom +
  patchwork::plot_layout(heights = c(1, 1))

# Exported at full double-column width (180 mm); at smaller export sizes,
# all type falls below the 5 pt minimum.
save_figure("Fig_eta_threshold_derivation_QC", fig_combined,
            width = NC_W_DOUBLE, height = 6.30)


# -----------------------------------------------------------------------------
# 21. CSV OUTPUTS
# -----------------------------------------------------------------------------
utils::write.csv(results_summary,
                 file.path(OUT_DIR, "Table_threshold_summary.csv"),
                 row.names = FALSE)

utils::write.csv(cv_results,
                 file.path(OUT_DIR, "Table_CV_utility_grid.csv"),
                 row.names = FALSE)

utils::write.csv(m3_df,
                 file.path(OUT_DIR, "Table_ARL_calibration.csv"),
                 row.names = FALSE)

utils::write.csv(data.frame(ratio = baseline_ratios),
                 file.path(OUT_DIR, "Data_baseline_ratios.csv"),
                 row.names = FALSE)

utils::write.csv(final_decision,
                 file.path(OUT_DIR, "Table_FINAL_eta_decision.csv"),
                 row.names = FALSE)

utils::write.csv(df,
                 file.path(OUT_DIR, "Data_analysis_dataset.csv"),
                 row.names = FALSE)


# -----------------------------------------------------------------------------
# 22. STAGE 1 RESULT BANNER
# -----------------------------------------------------------------------------
cat("\n\n############################################################\n")
cat("### STAGE 1 RESULT - eta_ON / eta_OFF for downstream use ###\n")
cat("############################################################\n")
cat("Excluded years: ", paste(EXCLUDED_YEARS, collapse = ", "), "\n")
cat("Included years: ", paste(included_years, collapse = ", "), "\n")
cat(sprintf("ADOPTION RULE   : %s\n", ADOPTION_RULE))
cat(sprintf("ETA_ON_ADOPTED  = %.4f  (raw: %.4f, rounded: %.1f)\n",
            ETA_ON_ADOPTED, ETA_ON_FINAL_RAW, ETA_ON_FINAL_ROUNDED))
cat(sprintf("ETA_OFF_ADOPTED = %.4f  (raw: %.4f, rounded: %.1f)\n",
            ETA_OFF_ADOPTED, ETA_OFF_FINAL_RAW, ETA_OFF_FINAL_ROUNDED))
cat(sprintf("Method-derived range: eta_ON in [%.2f, %.2f], eta_OFF in [%.2f, %.2f]\n",
            ETA_ON_FINAL_RANGE[1], ETA_ON_FINAL_RANGE[2],
            ETA_OFF_FINAL_RANGE[1], ETA_OFF_FINAL_RANGE[2]))
cat("\nDecision record : Table_FINAL_eta_decision.csv\n")
cat("Sourceable file : eta_thresholds_for_figure2.R\n")
cat("Usage downstream:\n")
cat(sprintf('  source("%s")\n', file.path(OUT_DIR, "eta_thresholds_for_figure2.R")))
cat("  # ETA_ON and ETA_OFF are then defined in the calling environment\n")


# -----------------------------------------------------------------------------
# 23. DRAFT FIGURE LEGEND
# -----------------------------------------------------------------------------
cat("\nFigure legend:\n")
cat(
  paste0(
    "Figure | Empirical derivation of the STA/LTA hysteresis thresholds for ",
    "the Vaezi-style outbreak detector, applied to weekly notified dengue ",
    "case counts in Quezon City, Philippines, after excluding 2020, 2021 and ",
    "2025 from all threshold computations. The derivation follows an ",
    "asymmetric operational objective: eta_ON is chosen to trigger early on ",
    "the upslope, while eta_OFF is chosen to keep the detector ON throughout ",
    "the decline until cases return to baseline. Six independent statistical ",
    "approaches were applied: (M1) baseline-ratio asymmetric quantiles (90th ",
    "and 50th percentiles of the null STA/LTA distribution); (M2) ",
    "high-sensitivity receiver-operating-characteristic operating points ",
    "(95% and 99% sensitivity); (M3) sensitivity-oriented average run-length ",
    "(ARL_0) calibration at 12 weeks; (M4) leave-one-year-out cross-validation ",
    "maximising an asymmetric utility (lead time + recall - false-alarm rate); ",
    "(M5) negative-binomial theoretical +/- 1 sigma bounds derived from the ",
    "over-dispersion of baseline counts; and (M6) parametric bootstrap of the ",
    "null STA/LTA ratio (90th and 25th percentiles, 2,000 replicates). The ",
    "median across methods was adopted as the final pair. Under the ",
    "pre-specified adoption rule ('", ADOPTION_RULE, "'), the values used in ",
    "the downstream analysis are eta_ON = ", sprintf("%.2f", ETA_ON_ADOPTED),
    " and eta_OFF = ", sprintf("%.2f", ETA_OFF_ADOPTED),
    ". (a) Empirical distribution of the STA/LTA ratio during baseline weeks, ",
    "with M1 quantiles (dashed lines) and adopted thresholds (solid lines). ",
    "(b) ROC curve of the STA/LTA ratio against elevated weekly incidence ",
    "(>75th percentile within year), with M2 high-sensitivity operating ",
    "points marked. (c) Leave-one-year-out cross-validated asymmetric-utility ",
    "surface over the (eta_ON, eta_OFF) grid; the M4 optimum (red circle) and ",
    "the adopted final pair (blue star) are highlighted. (d) Comparison of ",
    "all six method-level estimates with the adopted final pair; horizontal ",
    "bars span (eta_OFF, eta_ON) per method. STA, short-term average window ",
    "(4 weeks); LTA, long-term average window (26 weeks); GUARD, 2-week ",
    "guard band between STA and LTA windows.\n"
  )
)

cat("\nDone. All outputs saved to:\n  ", OUT_DIR, "\n")


# =============================================================================
# END OF STAGE 1 - Empirical eta_ON / eta_OFF derivation (Quezon City)
# =============================================================================