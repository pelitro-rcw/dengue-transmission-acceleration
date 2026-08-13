# =============================================================================
# STAGE 1B - ANCHOR-FOCUSED STA/LTA WINDOW SELECTION (QUEZON CITY)
# =============================================================================
#
# PURPOSE
#   Select / validate STA-LTA window configurations by asking whether the
#   STA/LTA ratio captures the TWO EPIDEMIOLOGIC ANCHORS used by the study:
#
#     A1 = Actionable Window (PRIMARY: 4-8 weeks before annual peak)
#     A2 = Epidemic-Burden block (PRIMARY: contiguous block containing 70%
#          of annual cases, grown outward from the annual maximum)
#
#   T = A1 union A2
#
#   A candidate window is therefore rewarded when it:
#     1) produces True alarms in A1 (early warning);
#     2) captures epidemic burden during A2 / T (TAM);
#     3) provides useful lead time and warning persistence in A1; and
#     4) avoids alarms outside T.
#
#   The window-selection engine is NOT an unrelated "sensitivity analysis".
#   AW and EB are the TARGET FRAMEWORKS used to score the STA/LTA signal.
#
# PRIMARY CANDIDATE SETS
#   Continuous TA: 2/8, 3/12, 4/16
#      - primary/prespecified: 3/12
#      - continuously updated right-aligned STA/LTA
#
#   Constant TA: 2/13, 3/20, 4/26
#      - primary/prespecified: 4/26
#      - mandatory 2-week guard
#      - historical LTA frozen at activation
#      - frozen baseline retained through brief OFF periods
#      - release after 8 consecutive OFF weeks
#
# WHY THESE ALTERNATIVES?
#   They are scale-matched neighborhoods around the prespecified architectures:
#     Continuous maintains LTA/STA = 4.
#     Constant maintains LTA/STA approximately 6.5, while guard/freeze stay fixed.
#
# PRIMARY ANCHOR SPECIFICATION
#   AW = 4-8 weeks before peak
#   EB = 70%
#
# ROBUSTNESS (NOT THE WINDOW-SELECTION OBJECTIVE)
#   AW perturbation: 3-6, 4-8, 5-10 weeks; EB fixed at 70%.
#   EB perturbation: 60%, 70%, 80%; AW fixed at 4-8.
#   Full 3x3 AW x EB grid is exported as supplementary robustness.
#
# FRAMEWORK METRICS (same definitions as Stage 3, without a Reactive category)
#   TAM                  higher is better
#   N_True_Alarms_yr     higher is better
#   N_False_Alarms_yr    lower is better
#   PPV                  higher is better
#   Sensitivity          higher is better
#   Mean_Lead_Time       higher is better
#   WP                   higher is better
#   ALY                  higher is better
#
# TRUE/FALSE CLASSIFICATION
#   True alarm  = alarm-on week in A1 OR A2.
#   False alarm = alarm-on week outside A1 AND A2.
#   There is NO "Reactive" category.
#   "Actionable true alarm" = True alarm that lies in A1; it is a subset used
#   only for sensitivity/timeliness.
#
# THRESHOLDS
#   eta_ON  = 1.33 (default)
#   eta_OFF = 0.73 (Constant TA default)
#   These may be overridden with environment variables:
#     TA_ETA_ON, TA_ETA_OFF
#
#   Stage 1B also calculates THRESHOLD-INDEPENDENT A1/A2 alignment statistics
#   from the raw ratio R(t), so the window evidence is not dependent solely on
#   the particular eta threshold.
#
# PRIMARY YEARS
#   2013-2019 and 2022-2024 (10 seasons)
#
# OUTPUT
#   outputs/Stage1b_anchor_focused_window_selection_QC/
#
# HEAD-TO-HEAD SIGNIFICANCE
#   Prespecified versus each same-family alternative:
#     Continuous 3/12 vs 2/8 and 4/16
#     Constant   4/26 vs 2/13 and 3/20
#
#   Inference is repeated under:
#     Primary AW 4-8 / EB 70%
#     AW sensitivity 3-6 / 4-8 / 5-10 with EB 70% fixed
#     EB sensitivity 60% / 70% / 80% with AW 4-8 fixed
#     Full 3x3 AW x EB supplementary robustness
#
#   Operational metrics:
#     exact paired season-label randomization (2^10 = 1024 assignments)
#     plus 2,000 season-cluster bootstrap BCa confidence intervals.
#
#   AUC metrics:
#     paired season-level AUC differences, exact sign-flip test,
#     plus 2,000 paired season bootstrap BCa confidence intervals.
#
#   Multiplicity:
#     Bonferroni within each metric (2 same-family comparators) and
#     stricter Bonferroni across all head-to-head metrics within each
#     detector x anchor specification.
#
# REPRODUCIBILITY
#   Base R + readxl only; seed = 12345; 1,000 selection bootstraps and
#   2,000 head-to-head season-cluster bootstraps.
# =============================================================================

STAGE1B_VERSION <- "STAGE1B_ANCHOR_FOCUSED_WINDOW_SELECTION"
message("[Stage1B] ", STAGE1B_VERSION)

GLOBAL_SEED <- 12345L
BOOT_N <- 1000L
H2H_BOOT_N <- 2000L
H2H_ALPHA <- 0.05
set.seed(GLOBAL_SEED)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Run install.packages('readxl') once.",
       call. = FALSE)
}

# -----------------------------------------------------------------------------
# 1. PORTABLE PATH RESOLUTION
# -----------------------------------------------------------------------------
DATA_BASENAME <- "Dengue-Rainfall_Dataset.xlsx"
SHEET_QC <- "QC Data"

.cli_value <- function(name) {
  args <- commandArgs(trailingOnly = FALSE)
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return("")
  sub(prefix, "", hit[length(hit)], fixed = TRUE)
}

.norm_file <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return("")
  x <- path.expand(x)
  if (!file.exists(x)) return("")
  normalizePath(x, winslash = "/", mustWork = TRUE)
}

.norm_dir <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return("")
  x <- path.expand(x)
  if (!dir.exists(x)) return("")
  normalizePath(x, winslash = "/", mustWork = TRUE)
}

.detect_script_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- args[startsWith(args, "--file=")]
  if (length(hit)) {
    p <- .norm_file(sub("--file=", "", hit[length(hit)], fixed = TRUE))
    if (nzchar(p)) return(p)
  }

  frames <- sys.frames()
  if (length(frames)) {
    for (i in rev(seq_along(frames))) {
      if (exists("ofile", envir = frames[[i]], inherits = FALSE)) {
        p <- get("ofile", envir = frames[[i]], inherits = FALSE)
        if (is.character(p) && length(p) == 1L) {
          p <- .norm_file(p)
          if (nzchar(p)) return(p)
        }
      }
    }
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ok <- tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)
    if (isTRUE(ok)) {
      p <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                    error = function(e) "")
      p <- .norm_file(p)
      if (nzchar(p)) return(p)
    }
  }
  ""
}

.parent_chain <- function(start, max_up = 8L) {
  start <- .norm_dir(start)
  if (!nzchar(start)) return(character())
  out <- character()
  d <- start
  for (i in 0:max_up) {
    out <- c(out, d)
    p <- dirname(d)
    if (identical(p, d)) break
    d <- p
  }
  unique(out)
}

.resolve_data_file <- function(script_file = "") {
  explicit <- .cli_value("data-file")
  if (!nzchar(explicit)) explicit <- Sys.getenv("TA_DATA_FILE", unset = "")
  if (nzchar(explicit)) {
    p <- .norm_file(explicit)
    if (!nzchar(p)) stop("Specified data file does not exist: ", explicit,
                         call. = FALSE)
    return(p)
  }

  starts <- unique(c(
    if (nzchar(script_file)) dirname(script_file) else character(),
    getwd()
  ))
  dirs <- unique(unlist(lapply(starts, .parent_chain), use.names = FALSE))

  candidates <- character()
  for (d in dirs) {
    candidates <- c(
      candidates,
      file.path(d, DATA_BASENAME),
      file.path(d, "data", DATA_BASENAME),
      file.path(d, "Transmission_Acceleration", "data", DATA_BASENAME)
    )
    kids <- tryCatch(list.dirs(d, recursive = FALSE, full.names = TRUE),
                     error = function(e) character())
    if (length(kids)) {
      candidates <- c(
        candidates,
        unlist(lapply(kids, function(k) c(
          file.path(k, DATA_BASENAME),
          file.path(k, "data", DATA_BASENAME)
        )), use.names = FALSE)
      )
    }
  }

  candidates <- unique(candidates[file.exists(candidates)])
  if (!length(candidates)) {
    stop(
      "Could not locate ", DATA_BASENAME, ". Set it explicitly with:\n",
      "Sys.setenv(TA_DATA_FILE='C:/path/Dengue-Rainfall_Dataset.xlsx')",
      call. = FALSE
    )
  }
  normalizePath(candidates[1], winslash = "/", mustWork = TRUE)
}

.resolve_output_dir <- function(data_file) {
  explicit <- .cli_value("output-dir")
  if (!nzchar(explicit)) explicit <- Sys.getenv("TA_OUTPUT_DIR", unset = "")

  if (nzchar(explicit)) {
    dir.create(explicit, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  }

  data_dir <- dirname(data_file)
  root <- if (tolower(basename(data_dir)) == "data") dirname(data_dir) else data_dir
  out <- file.path(root, "outputs", "Stage1b_anchor_focused_window_selection_QC")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  normalizePath(out, winslash = "/", mustWork = TRUE)
}

SCRIPT_FILE <- .detect_script_file()
DATA_FILE <- .resolve_data_file(SCRIPT_FILE)
OUT_DIR <- .resolve_output_dir(DATA_FILE)

message("[Stage1B] Input:  ", DATA_FILE)
message("[Stage1B] Output: ", OUT_DIR)
message("[Stage1B] No setwd() / project-root requirement.")

# -----------------------------------------------------------------------------
# 2. ANALYSIS SETTINGS
# -----------------------------------------------------------------------------
PRIMARY_YEARS <- c(2013L:2019L, 2022L:2024L)
EXCLUDED_YEARS <- c(2020L, 2021L, 2025L)

PRIMARY_AW_MIN <- 4L
PRIMARY_AW_MAX <- 8L
PRIMARY_EB <- 0.70

AW_LEVELS <- data.frame(
  AW = c("3-6", "4-8", "5-10"),
  LeadMin = c(3L, 4L, 5L),
  LeadMax = c(6L, 8L, 10L),
  Is_Primary = c(FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

EB_LEVELS <- data.frame(
  EB = c("60%", "70%", "80%"),
  Burden = c(0.60, 0.70, 0.80),
  Is_Primary = c(FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

CONT_CANDIDATES <- data.frame(
  Candidate_ID = c("Continuous_2_8", "Continuous_3_12", "Continuous_4_16"),
  Detector = "Continuous TA",
  Window = c("2/8", "3/12", "4/16"),
  STA = c(2L, 3L, 4L),
  LTA = c(8L, 12L, 16L),
  Guard = 0L,
  Prespecified = c(FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

CONST_CANDIDATES <- data.frame(
  Candidate_ID = c("Constant_2_13", "Constant_3_20", "Constant_4_26"),
  Detector = "Constant TA",
  Window = c("2/13", "3/20", "4/26"),
  STA = c(2L, 3L, 4L),
  LTA = c(13L, 20L, 26L),
  Guard = 2L,
  Prespecified = c(FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

CANDIDATES <- rbind(CONT_CANDIDATES, CONST_CANDIDATES)

CONST_GUARD <- 2L
CONST_MIN_OFF_RESET <- 8L

ETA_ON <- suppressWarnings(as.numeric(Sys.getenv("TA_ETA_ON", unset = "1.33")))
ETA_OFF <- suppressWarnings(as.numeric(Sys.getenv("TA_ETA_OFF", unset = "0.73")))

if (!is.finite(ETA_ON) || ETA_ON <= 0) {
  stop("ETA_ON must be finite and positive.", call. = FALSE)
}
if (!is.finite(ETA_OFF) || ETA_OFF <= 0 || ETA_OFF >= ETA_ON) {
  stop("ETA_OFF must be finite, positive, and lower than ETA_ON.", call. = FALSE)
}

METRICS <- c(
  "TAM",
  "N_True_Alarms_yr",
  "N_False_Alarms_yr",
  "PPV",
  "Sensitivity",
  "Mean_Lead_Time",
  "WP",
  "ALY"
)
DIRECTION <- c(
  TAM = "high",
  N_True_Alarms_yr = "high",
  N_False_Alarms_yr = "low",
  PPV = "high",
  Sensitivity = "high",
  Mean_Lead_Time = "high",
  WP = "high",
  ALY = "high"
)

H2H_AUC_METRICS <- c(
  "A1_AUC",
  "A2_AUC",
  "A2_CaseWeighted_AUC",
  "T_AUC"
)
H2H_ALL_METRICS <- c(METRICS, H2H_AUC_METRICS)
H2H_DIRECTION <- c(
  DIRECTION,
  A1_AUC = "high",
  A2_AUC = "high",
  A2_CaseWeighted_AUC = "high",
  T_AUC = "high"
)

H2H_COMPARISONS <- data.frame(
  Detector = c(
    "Continuous TA", "Continuous TA",
    "Constant TA", "Constant TA"
  ),
  Target_ID = c(
    "Continuous_3_12", "Continuous_3_12",
    "Constant_4_26", "Constant_4_26"
  ),
  Comparator_ID = c(
    "Continuous_2_8", "Continuous_4_16",
    "Constant_2_13", "Constant_3_20"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 3. LOAD / VALIDATE QC DATA
# -----------------------------------------------------------------------------
sheets <- readxl::excel_sheets(DATA_FILE)
if (!(SHEET_QC %in% sheets)) {
  stop("Workbook does not contain sheet '", SHEET_QC, "'.", call. = FALSE)
}

raw <- as.data.frame(readxl::read_excel(DATA_FILE, sheet = SHEET_QC))
required <- c("YR", "WN", "DC_QC")
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols)) {
  stop("Missing QC columns: ", paste(missing_cols, collapse = ", "),
       call. = FALSE)
}

df <- raw[, required]
df$YR <- as.integer(df$YR)
df$WN <- as.integer(df$WN)
df$DC_QC <- as.numeric(df$DC_QC)
df <- df[order(df$YR, df$WN), , drop = FALSE]
row.names(df) <- NULL

if (anyDuplicated(df[, c("YR", "WN")])) {
  stop("Duplicate YR/WN rows detected.", call. = FALSE)
}
if (any(df$DC_QC < 0, na.rm = TRUE)) {
  stop("Negative dengue counts detected.", call. = FALSE)
}
if (!all(PRIMARY_YEARS %in% unique(df$YR))) {
  stop("Not all primary years are available.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 4. GENERIC HELPERS
# -----------------------------------------------------------------------------
roll_mean_right <- function(x, w) {
  n <- length(x)
  out <- rep(NA_real_, n)
  if (w < 1L || n < w) return(out)
  for (i in seq.int(w, n)) {
    z <- x[(i - w + 1L):i]
    out[i] <- if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  }
  out
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

rank_metric <- function(x, direction) {
  if (direction == "high") {
    rank(-x, na.last = "keep", ties.method = "average")
  } else {
    rank(x, na.last = "keep", ties.method = "average")
  }
}

auc_rank <- function(y, score, weights = NULL) {
  keep <- !is.na(y) & is.finite(score)
  y <- as.integer(y[keep])
  score <- score[keep]

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  } else {
    weights <- as.numeric(weights[keep])
    weights[!is.finite(weights) | weights < 0] <- 0
  }

  pos <- which(y == 1L)
  neg <- which(y == 0L)
  if (!length(pos) || !length(neg)) return(NA_real_)

  # Weighted pairwise concordance. Only 10 primary seasons and 6 candidates,
  # so the O(n_pos*n_neg) calculation is acceptable and transparent.
  num <- 0
  den <- 0
  for (i in pos) {
    for (j in neg) {
      w <- weights[i] * weights[j]
      if (w <= 0) next
      den <- den + w
      if (score[i] > score[j]) num <- num + w
      else if (score[i] == score[j]) num <- num + 0.5 * w
    }
  }
  if (den <= 0) NA_real_ else num / den
}

# -----------------------------------------------------------------------------
# 5. A1 / A2 ANCHORS
# -----------------------------------------------------------------------------
compute_A1 <- function(df_in, year, lead_min, lead_max) {
  z <- df_in[df_in$YR == year, , drop = FALSE]
  z <- z[order(z$WN), , drop = FALSE]

  if (!nrow(z) || all(is.na(z$DC_QC))) {
    return(list(peak_week = NA_integer_, A1_weeks = integer()))
  }

  peak_week <- as.integer(z$WN[which.max(z$DC_QC)])
  a1 <- (peak_week - lead_max):(peak_week - lead_min)
  a1 <- as.integer(a1[a1 >= 1L])

  list(peak_week = peak_week, A1_weeks = a1)
}

compute_A2 <- function(df_in, year, burden_frac) {
  z <- df_in[df_in$YR == year, , drop = FALSE]
  z <- z[order(z$WN), , drop = FALSE]

  if (!nrow(z) || all(is.na(z$DC_QC))) {
    return(list(
      start_week = NA_integer_,
      end_week = NA_integer_,
      A2_weeks = integer(),
      achieved_fraction = NA_real_
    ))
  }

  cases <- ifelse(is.na(z$DC_QC), 0, z$DC_QC)
  weeks <- as.integer(z$WN)
  total <- sum(cases)

  if (!is.finite(total) || total <= 0) {
    return(list(
      start_week = NA_integer_,
      end_week = NA_integer_,
      A2_weeks = integer(),
      achieved_fraction = NA_real_
    ))
  }

  p <- which.max(cases)
  lo <- p
  hi <- p
  accum <- cases[p]

  while (accum / total < burden_frac && (lo > 1L || hi < length(cases))) {
    left <- if (lo > 1L) cases[lo - 1L] else -Inf
    right <- if (hi < length(cases)) cases[hi + 1L] else -Inf

    if (left >= right) {
      lo <- lo - 1L
      accum <- accum + left
    } else {
      hi <- hi + 1L
      accum <- accum + right
    }
  }

  a2 <- as.integer(weeks[lo:hi])
  list(
    start_week = min(a2),
    end_week = max(a2),
    A2_weeks = a2,
    achieved_fraction = accum / total
  )
}

compute_anchors <- function(year, lead_min, lead_max, burden_frac) {
  a1 <- compute_A1(df, year, lead_min, lead_max)
  a2 <- compute_A2(df, year, burden_frac)

  list(
    peak_week = a1$peak_week,
    A1_weeks = a1$A1_weeks,
    A2_weeks = a2$A2_weeks,
    T_weeks = sort(unique(c(a1$A1_weeks, a2$A2_weeks))),
    achieved_fraction = a2$achieved_fraction
  )
}

# Anchor contract self-test over every robustness specification.
for (yr in PRIMARY_YEARS) {
  for (a in seq_len(nrow(AW_LEVELS))) {
    for (b in seq_len(nrow(EB_LEVELS))) {
      an <- compute_anchors(
        yr,
        AW_LEVELS$LeadMin[a],
        AW_LEVELS$LeadMax[a],
        EB_LEVELS$Burden[b]
      )
      if (!is.integer(an$A1_weeks) || !is.integer(an$A2_weeks)) {
        stop("Anchor return contract failed for year ", yr, ".", call. = FALSE)
      }
      if (is.finite(an$achieved_fraction) &&
          an$achieved_fraction + 1e-12 < EB_LEVELS$Burden[b]) {
        stop("Epidemic-burden block did not attain requested fraction for year ",
             yr, ".", call. = FALSE)
      }
    }
  }
}
message("[Stage1B] All AW/EB anchor contracts: PASS")

# -----------------------------------------------------------------------------
# 6. CONTINUOUS TA
# -----------------------------------------------------------------------------
build_continuous_ta <- function(dc, sta_win, lta_win, eta_on = ETA_ON) {
  sta <- roll_mean_right(dc, sta_win)
  lta <- roll_mean_right(dc, lta_win)
  ratio <- ifelse(
    is.finite(sta) & is.finite(lta) & lta > 0,
    sta / lta,
    NA_real_
  )

  data.frame(
    STA = sta,
    LTA_live = lta,
    LTA_used = lta,
    R = ratio,
    trigger = as.integer(is.finite(ratio) & ratio > eta_on),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 7. CONSTANT TA - GUARDED + FROZEN BASELINE
# -----------------------------------------------------------------------------
build_constant_ta <- function(dc, yr, sta_win, lta_win,
                              guard = CONST_GUARD,
                              eta_on = ETA_ON,
                              eta_off = ETA_OFF,
                              min_off_reset = CONST_MIN_OFF_RESET,
                              reset_at_year_boundary = TRUE) {
  if (guard < 1L) {
    stop("Constant TA requires a positive guard.", call. = FALSE)
  }
  if (!is.finite(eta_on) || !is.finite(eta_off) || eta_off >= eta_on) {
    stop("Constant TA requires finite eta_OFF < eta_ON.", call. = FALSE)
  }

  n <- length(dc)
  sta_v <- rep(NA_real_, n)
  live_lta_v <- rep(NA_real_, n)
  used_lta_v <- rep(NA_real_, n)
  frozen_lta_v <- rep(NA_real_, n)
  ratio_v <- rep(NA_real_, n)
  trigger_v <- integer(n)
  freeze_id_v <- rep(NA_integer_, n)
  activation_v <- logical(n)
  release_v <- logical(n)

  is_on <- FALSE
  frozen <- NA_real_
  off_streak <- 0L
  freeze_counter <- 0L
  current_freeze_id <- NA_integer_

  for (t in seq_len(n)) {
    if (reset_at_year_boundary && t > 1L &&
        !is.na(yr[t]) && !is.na(yr[t - 1L]) &&
        yr[t] != yr[t - 1L]) {
      is_on <- FALSE
      frozen <- NA_real_
      off_streak <- 0L
      current_freeze_id <- NA_integer_
    }

    sta_start <- t - sta_win + 1L
    lta_end <- t - sta_win - guard
    lta_start <- lta_end - lta_win + 1L

    if (sta_start < 1L || lta_start < 1L || lta_end < 1L) next

    sta_values <- dc[sta_start:t]
    lta_values <- dc[lta_start:lta_end]

    if (length(sta_values) != sta_win || length(lta_values) != lta_win) {
      stop("Constant TA internal window geometry error.", call. = FALSE)
    }

    sta <- if (all(is.na(sta_values))) NA_real_ else mean(sta_values, na.rm = TRUE)
    live_lta <- if (all(is.na(lta_values))) NA_real_ else
      mean(lta_values, na.rm = TRUE)

    sta_v[t] <- sta
    live_lta_v[t] <- live_lta

    used_lta <- if (is.finite(frozen)) frozen else live_lta
    ratio <- if (is.finite(sta) && is.finite(used_lta) && used_lta > 0) {
      sta / used_lta
    } else NA_real_

    if (!is_on && !is.finite(frozen) &&
        is.finite(ratio) && ratio >= eta_on) {
      if (!is.finite(live_lta) || live_lta <= 0) {
        stop("Activation attempted without valid guarded LTA.", call. = FALSE)
      }

      frozen <- live_lta
      freeze_counter <- freeze_counter + 1L
      current_freeze_id <- freeze_counter
      is_on <- TRUE
      off_streak <- 0L
      activation_v[t] <- TRUE

      used_lta <- frozen
      ratio <- sta / frozen

    } else if (is.finite(frozen)) {
      used_lta <- frozen
      ratio <- if (is.finite(sta) && frozen > 0) sta / frozen else NA_real_

      if (is_on) {
        if (is.finite(ratio) && ratio < eta_off) {
          is_on <- FALSE
          off_streak <- 1L
        } else {
          off_streak <- 0L
        }
      } else {
        if (is.finite(ratio) && ratio >= eta_on) {
          is_on <- TRUE
          off_streak <- 0L
        } else {
          off_streak <- off_streak + 1L
        }
      }
    }

    used_lta_v[t] <- used_lta
    frozen_lta_v[t] <- frozen
    ratio_v[t] <- ratio
    trigger_v[t] <- as.integer(is_on)
    freeze_id_v[t] <- current_freeze_id

    if (!is_on && is.finite(frozen) && off_streak >= min_off_reset) {
      release_v[t] <- TRUE
      frozen <- NA_real_
      current_freeze_id <- NA_integer_
      off_streak <- 0L
    }
  }

  out <- data.frame(
    STA = sta_v,
    LTA_live = live_lta_v,
    LTA_used = used_lta_v,
    LTA_frozen = frozen_lta_v,
    R = ratio_v,
    trigger = trigger_v,
    freeze_id = freeze_id_v,
    activation = activation_v,
    release = release_v,
    stringsAsFactors = FALSE
  )

  # Audit frozen baseline within each episode.
  ids <- sort(unique(out$freeze_id[is.finite(out$freeze_id)]))
  for (id in ids) {
    ii <- which(out$freeze_id == id)
    vals <- out$LTA_frozen[ii]
    vals <- vals[is.finite(vals)]
    if (length(vals) > 1L && max(vals) - min(vals) > 1e-12) {
      stop("Frozen LTA drift detected within freeze_id=", id, ".", call. = FALSE)
    }

    aa <- ii[out$activation[ii]]
    if (length(aa) != 1L) {
      stop("Each Constant-TA freeze_id must have exactly one activation.",
           call. = FALSE)
    }
    if (!isTRUE(all.equal(
      out$LTA_frozen[aa], out$LTA_live[aa], tolerance = 1e-12
    ))) {
      stop("Frozen LTA does not equal guarded LTA at activation.",
           call. = FALSE)
    }
  }

  out
}

# -----------------------------------------------------------------------------
# 8. BUILD SIX WINDOW CANDIDATES
# -----------------------------------------------------------------------------
DETECTOR_OBJECTS <- vector("list", nrow(CANDIDATES))
names(DETECTOR_OBJECTS) <- CANDIDATES$Candidate_ID

for (i in seq_len(nrow(CANDIDATES))) {
  cnd <- CANDIDATES[i, ]

  if (cnd$Detector == "Continuous TA") {
    DETECTOR_OBJECTS[[cnd$Candidate_ID]] <- build_continuous_ta(
      df$DC_QC, cnd$STA, cnd$LTA, ETA_ON
    )
  } else {
    DETECTOR_OBJECTS[[cnd$Candidate_ID]] <- build_constant_ta(
      df$DC_QC, df$YR,
      sta_win = cnd$STA,
      lta_win = cnd$LTA,
      guard = CONST_GUARD,
      eta_on = ETA_ON,
      eta_off = ETA_OFF,
      min_off_reset = CONST_MIN_OFF_RESET,
      reset_at_year_boundary = TRUE
    )
  }
}

message("[Stage1B] All six candidate detector objects built: PASS")

# -----------------------------------------------------------------------------
# 9. EXACT STAGE-3-STYLE OPERATIONAL METRICS, ANCHORED TO A1/A2
# -----------------------------------------------------------------------------
compute_year_raw <- function(trigger, year, lead_min, lead_max, burden_frac) {
  ii <- which(df$YR == year)
  z <- df[ii, , drop = FALSE]
  z <- z[order(z$WN), , drop = FALSE]
  trig <- as.integer(trigger[ii][order(df$WN[ii])] == 1L)

  a <- compute_anchors(year, lead_min, lead_max, burden_frac)
  weeks <- as.integer(z$WN)
  cases <- ifelse(is.na(z$DC_QC), 0, z$DC_QC)

  in_a1 <- weeks %in% a$A1_weeks
  in_a2 <- weeks %in% a$A2_weeks
  in_t <- in_a1 | in_a2

  is_true <- trig == 1L & in_t
  is_false <- trig == 1L & !in_t
  is_actionable <- is_true & in_a1

  true_n <- sum(is_true)
  false_n <- sum(is_false)
  total_n <- true_n + false_n
  actionable_n <- sum(is_actionable)

  tam <- sum(cases[is_true], na.rm = TRUE)

  actionable_weeks <- weeks[is_actionable]
  actionable_leads <- if (length(actionable_weeks) && is.finite(a$peak_week)) {
    a$peak_week - actionable_weeks
  } else numeric()

  first_a1_week <- if (length(actionable_weeks)) min(actionable_weeks) else NA_integer_
  lead_time <- if (is.finite(first_a1_week) && is.finite(a$peak_week)) {
    as.numeric(a$peak_week - first_a1_week)
  } else 0

  wp <- if (length(actionable_leads)) mean(actionable_leads) else 0
  aly <- if (true_n > 0L) actionable_n / true_n else 0

  data.frame(
    Year = as.integer(year),
    TAM_yr = tam,
    Total_Triggers_yr = total_n,
    True_Alarms_yr = true_n,
    False_Alarms_yr = false_n,
    PPV_yr = if (total_n > 0L) true_n / total_n else NA_real_,
    Has_A1_True = as.integer(actionable_n > 0L),
    Lead_Time_yr = lead_time,
    WP_yr = wp,
    ALY_yr = aly,
    stringsAsFactors = FALSE
  )
}

aggregate_raw <- function(year_df) {
  n_years <- nrow(year_df)
  total_triggers <- sum(year_df$Total_Triggers_yr, na.rm = TRUE)
  total_true <- sum(year_df$True_Alarms_yr, na.rm = TRUE)

  data.frame(
    TAM = mean(year_df$TAM_yr, na.rm = TRUE),
    N_True_Alarms_yr = total_true / n_years,
    N_False_Alarms_yr = sum(year_df$False_Alarms_yr, na.rm = TRUE) / n_years,
    PPV = if (total_triggers > 0L) total_true / total_triggers else NA_real_,
    Sensitivity = mean(year_df$Has_A1_True, na.rm = TRUE),
    Mean_Lead_Time = mean(year_df$Lead_Time_yr, na.rm = TRUE),
    WP = mean(year_df$WP_yr, na.rm = TRUE),
    ALY = mean(year_df$ALY_yr, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

evaluate_candidate <- function(candidate_id,
                               lead_min = PRIMARY_AW_MIN,
                               lead_max = PRIMARY_AW_MAX,
                               burden_frac = PRIMARY_EB,
                               years = PRIMARY_YEARS) {
  trigger <- DETECTOR_OBJECTS[[candidate_id]]$trigger

  yr <- do.call(rbind, lapply(years, function(y) {
    compute_year_raw(trigger, y, lead_min, lead_max, burden_frac)
  }))

  list(yearly = yr, summary = aggregate_raw(yr))
}

# -----------------------------------------------------------------------------
# 10. THRESHOLD-INDEPENDENT RATIO-TO-ANCHOR ALIGNMENT
# -----------------------------------------------------------------------------
# These statistics ask directly whether high R(t) aligns with A1/A2/T,
# independently of eta_ON.
compute_ratio_alignment <- function(candidate_id,
                                    lead_min = PRIMARY_AW_MIN,
                                    lead_max = PRIMARY_AW_MAX,
                                    burden_frac = PRIMARY_EB,
                                    years = PRIMARY_YEARS) {
  obj <- DETECTOR_OBJECTS[[candidate_id]]
  ratio <- obj$R

  y_a1 <- integer()
  y_a2 <- integer()
  y_t <- integer()
  scores <- numeric()
  case_weights <- numeric()
  r_a1 <- numeric()
  r_a2 <- numeric()
  r_out <- numeric()

  for (yr in years) {
    ii <- which(df$YR == yr)
    a <- compute_anchors(yr, lead_min, lead_max, burden_frac)
    weeks <- df$WN[ii]
    cases <- ifelse(is.na(df$DC_QC[ii]), 0, df$DC_QC[ii])
    rr <- ratio[ii]

    in_a1 <- weeks %in% a$A1_weeks
    in_a2 <- weeks %in% a$A2_weeks
    in_t <- in_a1 | in_a2
    valid <- is.finite(rr)

    # Compare anchor weeks with weeks outside the complete T region.
    use_a1 <- valid & (in_a1 | !in_t)
    use_a2 <- valid & (in_a2 | !in_t)
    use_t <- valid

    y_a1 <- c(y_a1, as.integer(in_a1[use_a1]))
    y_a2 <- c(y_a2, as.integer(in_a2[use_a2]))
    y_t <- c(y_t, as.integer(in_t[use_t]))

    # Separate score vectors are reconstructed below from saved partitions.
    # For AUC we retain a common concatenated vector with masks.
    # Store raw partitions for enrichment statistics.
    r_a1 <- c(r_a1, rr[valid & in_a1])
    r_a2 <- c(r_a2, rr[valid & in_a2])
    r_out <- c(r_out, rr[valid & !in_t])

    scores <- c(scores, rr[use_t])
    case_weights <- c(case_weights, cases[use_t] + 1)
  }

  # Reconstruct dedicated AUC vectors to avoid mismatched lengths.
  a1_y <- integer(); a1_s <- numeric()
  a2_y <- integer(); a2_s <- numeric(); a2_w <- numeric()
  t_y <- integer(); t_s <- numeric()

  for (yr in years) {
    ii <- which(df$YR == yr)
    a <- compute_anchors(yr, lead_min, lead_max, burden_frac)
    weeks <- df$WN[ii]
    cases <- ifelse(is.na(df$DC_QC[ii]), 0, df$DC_QC[ii])
    rr <- ratio[ii]
    valid <- is.finite(rr)
    in_a1 <- weeks %in% a$A1_weeks
    in_a2 <- weeks %in% a$A2_weeks
    in_t <- in_a1 | in_a2

    k1 <- valid & (in_a1 | !in_t)
    a1_y <- c(a1_y, as.integer(in_a1[k1]))
    a1_s <- c(a1_s, rr[k1])

    k2 <- valid & (in_a2 | !in_t)
    a2_y <- c(a2_y, as.integer(in_a2[k2]))
    a2_s <- c(a2_s, rr[k2])
    a2_w <- c(a2_w, cases[k2] + 1)

    kt <- valid
    t_y <- c(t_y, as.integer(in_t[kt]))
    t_s <- c(t_s, rr[kt])
  }

  med_out <- stats::median(r_out, na.rm = TRUE)
  med_a1 <- stats::median(r_a1, na.rm = TRUE)
  med_a2 <- stats::median(r_a2, na.rm = TRUE)

  data.frame(
    A1_AUC = auc_rank(a1_y, a1_s),
    A2_AUC = auc_rank(a2_y, a2_s),
    A2_CaseWeighted_AUC = auc_rank(a2_y, a2_s, a2_w),
    T_AUC = auc_rank(t_y, t_s),
    A1_Median_R = med_a1,
    A2_Median_R = med_a2,
    OutsideT_Median_R = med_out,
    A1_Enrichment = if (is.finite(med_out) && med_out > 0) med_a1 / med_out else NA_real_,
    A2_Enrichment = if (is.finite(med_out) && med_out > 0) med_a2 / med_out else NA_real_,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 11. PRIMARY AW=4-8 / EB=70% WINDOW SCORING
# -----------------------------------------------------------------------------
PRIMARY_ROWS <- vector("list", nrow(CANDIDATES))
PRIMARY_YEAR_ROWS <- vector("list", nrow(CANDIDATES))

for (i in seq_len(nrow(CANDIDATES))) {
  cnd <- CANDIDATES[i, ]
  ev <- evaluate_candidate(cnd$Candidate_ID)
  al <- compute_ratio_alignment(cnd$Candidate_ID)

  PRIMARY_ROWS[[i]] <- cbind(
    cnd,
    ev$summary,
    al
  )

  y <- ev$yearly
  y$Candidate_ID <- cnd$Candidate_ID
  y$Detector <- cnd$Detector
  y$Window <- cnd$Window
  y$STA <- cnd$STA
  y$LTA <- cnd$LTA
  y$Prespecified <- cnd$Prespecified
  PRIMARY_YEAR_ROWS[[i]] <- y
}

PRIMARY_RESULTS <- do.call(rbind, PRIMARY_ROWS)
PRIMARY_YEARLY <- do.call(rbind, PRIMARY_YEAR_ROWS)
row.names(PRIMARY_RESULTS) <- NULL
row.names(PRIMARY_YEARLY) <- NULL

# Rank operational metrics within detector family.
for (m in METRICS) PRIMARY_RESULTS[[paste0(m, "_Rank")]] <- NA_real_
PRIMARY_RESULTS$Operational_Mean_Rank <- NA_real_

# Direct ratio-alignment ranks, all higher is better.
ALIGNMENT_METRICS <- c(
  "A1_AUC",
  "A2_CaseWeighted_AUC",
  "T_AUC",
  "A1_Enrichment",
  "A2_Enrichment"
)
for (m in ALIGNMENT_METRICS) PRIMARY_RESULTS[[paste0(m, "_Rank")]] <- NA_real_
PRIMARY_RESULTS$Alignment_Mean_Rank <- NA_real_
PRIMARY_RESULTS$Anchor_Focused_Mean_Rank <- NA_real_
PRIMARY_RESULTS$Anchor_Focused_Final_Rank <- NA_real_

for (det in unique(PRIMARY_RESULTS$Detector)) {
  ii <- which(PRIMARY_RESULTS$Detector == det)

  operational_rank_mat <- matrix(
    NA_real_, nrow = length(ii), ncol = length(METRICS)
  )
  for (j in seq_along(METRICS)) {
    m <- METRICS[j]
    rr <- rank_metric(PRIMARY_RESULTS[[m]][ii], DIRECTION[[m]])
    PRIMARY_RESULTS[[paste0(m, "_Rank")]][ii] <- rr
    operational_rank_mat[, j] <- rr
  }
  PRIMARY_RESULTS$Operational_Mean_Rank[ii] <-
    rowMeans(operational_rank_mat, na.rm = TRUE)

  alignment_rank_mat <- matrix(
    NA_real_, nrow = length(ii), ncol = length(ALIGNMENT_METRICS)
  )
  for (j in seq_along(ALIGNMENT_METRICS)) {
    m <- ALIGNMENT_METRICS[j]
    rr <- rank_metric(PRIMARY_RESULTS[[m]][ii], "high")
    PRIMARY_RESULTS[[paste0(m, "_Rank")]][ii] <- rr
    alignment_rank_mat[, j] <- rr
  }
  PRIMARY_RESULTS$Alignment_Mean_Rank[ii] <-
    rowMeans(alignment_rank_mat, na.rm = TRUE)

  # Equal rank-level contribution:
  # 50% operational A1/A2 alarm behavior + 50% raw ratio alignment to A1/A2.
  # This avoids letting the larger number of operational metrics mechanically
  # swamp direct ratio-to-anchor evidence.
  PRIMARY_RESULTS$Anchor_Focused_Mean_Rank[ii] <-
    (
      PRIMARY_RESULTS$Operational_Mean_Rank[ii] +
      PRIMARY_RESULTS$Alignment_Mean_Rank[ii]
    ) / 2

  PRIMARY_RESULTS$Anchor_Focused_Final_Rank[ii] <- rank(
    PRIMARY_RESULTS$Anchor_Focused_Mean_Rank[ii],
    ties.method = "average"
  )
}

# -----------------------------------------------------------------------------
# 12. AW ROBUSTNESS: vary A1 only; EB remains 70%
# -----------------------------------------------------------------------------
AW_ROBUST_ROWS <- list()
k <- 0L

for (a in seq_len(nrow(AW_LEVELS))) {
  for (i in seq_len(nrow(CANDIDATES))) {
    cnd <- CANDIDATES[i, ]
    ev <- evaluate_candidate(
      cnd$Candidate_ID,
      lead_min = AW_LEVELS$LeadMin[a],
      lead_max = AW_LEVELS$LeadMax[a],
      burden_frac = PRIMARY_EB
    )
    al <- compute_ratio_alignment(
      cnd$Candidate_ID,
      lead_min = AW_LEVELS$LeadMin[a],
      lead_max = AW_LEVELS$LeadMax[a],
      burden_frac = PRIMARY_EB
    )

    k <- k + 1L
    AW_ROBUST_ROWS[[k]] <- cbind(
      data.frame(
        AW = AW_LEVELS$AW[a],
        Is_Primary_AW = AW_LEVELS$Is_Primary[a],
        EB = PRIMARY_EB,
        stringsAsFactors = FALSE
      ),
      cnd,
      ev$summary,
      al
    )
  }
}

AW_ROBUST <- do.call(rbind, AW_ROBUST_ROWS)
row.names(AW_ROBUST) <- NULL

# -----------------------------------------------------------------------------
# 13. EB ROBUSTNESS: vary A2 only; AW remains 4-8
# -----------------------------------------------------------------------------
EB_ROBUST_ROWS <- list()
k <- 0L

for (b in seq_len(nrow(EB_LEVELS))) {
  for (i in seq_len(nrow(CANDIDATES))) {
    cnd <- CANDIDATES[i, ]
    ev <- evaluate_candidate(
      cnd$Candidate_ID,
      lead_min = PRIMARY_AW_MIN,
      lead_max = PRIMARY_AW_MAX,
      burden_frac = EB_LEVELS$Burden[b]
    )
    al <- compute_ratio_alignment(
      cnd$Candidate_ID,
      lead_min = PRIMARY_AW_MIN,
      lead_max = PRIMARY_AW_MAX,
      burden_frac = EB_LEVELS$Burden[b]
    )

    k <- k + 1L
    EB_ROBUST_ROWS[[k]] <- cbind(
      data.frame(
        AW = "4-8",
        EB = EB_LEVELS$EB[b],
        Burden = EB_LEVELS$Burden[b],
        Is_Primary_EB = EB_LEVELS$Is_Primary[b],
        stringsAsFactors = FALSE
      ),
      cnd,
      ev$summary,
      al
    )
  }
}

EB_ROBUST <- do.call(rbind, EB_ROBUST_ROWS)
row.names(EB_ROBUST) <- NULL

# -----------------------------------------------------------------------------
# 14. FULL 3x3 AW x EB ROBUSTNESS (SUPPLEMENTARY)
# -----------------------------------------------------------------------------
FULL_ROBUST_ROWS <- list()
k <- 0L

for (a in seq_len(nrow(AW_LEVELS))) {
  for (b in seq_len(nrow(EB_LEVELS))) {
    for (i in seq_len(nrow(CANDIDATES))) {
      cnd <- CANDIDATES[i, ]
      ev <- evaluate_candidate(
        cnd$Candidate_ID,
        lead_min = AW_LEVELS$LeadMin[a],
        lead_max = AW_LEVELS$LeadMax[a],
        burden_frac = EB_LEVELS$Burden[b]
      )
      al <- compute_ratio_alignment(
        cnd$Candidate_ID,
        lead_min = AW_LEVELS$LeadMin[a],
        lead_max = AW_LEVELS$LeadMax[a],
        burden_frac = EB_LEVELS$Burden[b]
      )

      k <- k + 1L
      FULL_ROBUST_ROWS[[k]] <- cbind(
        data.frame(
          AW = AW_LEVELS$AW[a],
          EB = EB_LEVELS$EB[b],
          LeadMin = AW_LEVELS$LeadMin[a],
          LeadMax = AW_LEVELS$LeadMax[a],
          Burden = EB_LEVELS$Burden[b],
          Is_Primary = AW_LEVELS$Is_Primary[a] & EB_LEVELS$Is_Primary[b],
          stringsAsFactors = FALSE
        ),
        cnd,
        ev$summary,
        al
      )
    }
  }
}

FULL_ROBUST <- do.call(rbind, FULL_ROBUST_ROWS)
row.names(FULL_ROBUST) <- NULL

# Rank candidates within detector x anchor specification using exactly the same
# anchor-focused ranking logic as the primary analysis.
rank_anchor_table <- function(tab, group_cols) {
  out <- tab
  for (m in METRICS) out[[paste0(m, "_Rank")]] <- NA_real_
  for (m in ALIGNMENT_METRICS) out[[paste0(m, "_Rank")]] <- NA_real_
  out$Operational_Mean_Rank <- NA_real_
  out$Alignment_Mean_Rank <- NA_real_
  out$Anchor_Focused_Mean_Rank <- NA_real_
  out$Anchor_Focused_Final_Rank <- NA_real_

  keys <- interaction(out[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)

  for (key in levels(keys)) {
    ii <- which(keys == key)
    # group_cols include Detector, so each group should contain three windows.
    if (length(ii) != 3L) {
      stop("Anchor robustness rank group must contain exactly 3 candidates.",
           call. = FALSE)
    }

    op <- matrix(NA_real_, nrow = 3L, ncol = length(METRICS))
    for (j in seq_along(METRICS)) {
      m <- METRICS[j]
      rr <- rank_metric(out[[m]][ii], DIRECTION[[m]])
      out[[paste0(m, "_Rank")]][ii] <- rr
      op[, j] <- rr
    }
    out$Operational_Mean_Rank[ii] <- rowMeans(op, na.rm = TRUE)

    al <- matrix(NA_real_, nrow = 3L, ncol = length(ALIGNMENT_METRICS))
    for (j in seq_along(ALIGNMENT_METRICS)) {
      m <- ALIGNMENT_METRICS[j]
      rr <- rank_metric(out[[m]][ii], "high")
      out[[paste0(m, "_Rank")]][ii] <- rr
      al[, j] <- rr
    }
    out$Alignment_Mean_Rank[ii] <- rowMeans(al, na.rm = TRUE)

    out$Anchor_Focused_Mean_Rank[ii] <-
      (out$Operational_Mean_Rank[ii] + out$Alignment_Mean_Rank[ii]) / 2

    out$Anchor_Focused_Final_Rank[ii] <- rank(
      out$Anchor_Focused_Mean_Rank[ii],
      ties.method = "average"
    )
  }

  out
}

AW_ROBUST_RANKED <- rank_anchor_table(
  AW_ROBUST,
  c("AW", "Detector")
)
EB_ROBUST_RANKED <- rank_anchor_table(
  EB_ROBUST,
  c("EB", "Detector")
)
FULL_ROBUST_RANKED <- rank_anchor_table(
  FULL_ROBUST,
  c("AW", "EB", "Detector")
)

# -----------------------------------------------------------------------------
# 15. HEAD-TO-HEAD SIGNIFICANCE: PRESPECIFIED VS SAME-FAMILY ALTERNATIVES
# -----------------------------------------------------------------------------
# Positive ORIENTED differences always favor the prespecified configuration.
# For false alarms, the raw difference is target - comparator but the oriented
# difference is comparator - target because lower false-alarm burden is better.
#
# Operational inference:
#   - effect estimate uses the SAME aggregate_raw() definitions as the headline
#     tables, including pooled PPV.
#   - exact paired randomization swaps the complete year record between the two
#     configurations within each season (2^10 = 1024 assignments).
#   - paired season-cluster bootstrap resamples the same year indices for both
#     configurations and returns a BCa 95% CI.
#
# AUC inference:
#   - pooled AUC estimates shown in the table match compute_ratio_alignment().
#   - inferential effect is the mean paired SEASON-LEVEL AUC difference.
#   - exact sign-flip and paired season bootstrap are applied to those 10
#     season-level AUC differences. This avoids treating weekly observations as
#     independent and therefore respects the surveillance-season clustering.

oriented_difference <- function(target_value, comparator_value, metric) {
  raw <- target_value - comparator_value
  if (H2H_DIRECTION[[metric]] == "high") raw else -raw
}

bca_limits <- function(boot, observed, jack, alpha = H2H_ALPHA) {
  boot <- boot[is.finite(boot)]
  jack <- jack[is.finite(jack)]

  if (length(boot) < 100L || !is.finite(observed)) {
    return(c(
      lower = NA_real_,
      upper = NA_real_,
      method = NA_character_
    ))
  }

  # Degenerate bootstrap distributions require no quantile transformation.
  if (max(boot) - min(boot) < .Machine$double.eps^0.5) {
    return(c(
      lower = boot[1],
      upper = boot[1],
      method = "degenerate"
    ))
  }

  # Bias correction.
  prop_less <- mean(boot < observed)
  eps <- 1 / (2 * length(boot))
  prop_less <- min(max(prop_less, eps), 1 - eps)
  z0 <- stats::qnorm(prop_less)

  # Acceleration from leave-one-season-out jackknife.
  if (length(jack) >= 3L && max(jack) - min(jack) > .Machine$double.eps^0.5) {
    jbar <- mean(jack)
    u <- jbar - jack
    den <- 6 * (sum(u^2)^(3/2))
    accel <- if (is.finite(den) && den > 0) sum(u^3) / den else 0
  } else {
    accel <- 0
  }

  zlo <- stats::qnorm(alpha / 2)
  zhi <- stats::qnorm(1 - alpha / 2)

  adj_prob <- function(zalpha) {
    denom <- 1 - accel * (z0 + zalpha)
    if (!is.finite(denom) || abs(denom) < 1e-12) return(NA_real_)
    stats::pnorm(z0 + (z0 + zalpha) / denom)
  }

  plo <- adj_prob(zlo)
  phi <- adj_prob(zhi)

  if (!is.finite(plo) || !is.finite(phi) ||
      plo <= 0 || phi >= 1 || plo >= phi) {
    q <- stats::quantile(
      boot,
      probs = c(alpha / 2, 1 - alpha / 2),
      na.rm = TRUE,
      names = FALSE,
      type = 8
    )
    return(c(
      lower = q[1],
      upper = q[2],
      method = "percentile_fallback"
    ))
  }

  q <- stats::quantile(
    boot,
    probs = c(plo, phi),
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )
  c(lower = q[1], upper = q[2], method = "BCa")
}

exact_signflip_mean <- function(d) {
  d <- as.numeric(d[is.finite(d)])
  n <- length(d)
  if (!n) {
    return(c(
      observed = NA_real_,
      p_two_sided = NA_real_,
      p_target_better = NA_real_,
      n_pairs = 0
    ))
  }

  observed <- mean(d)
  n_perm <- 2^n
  perm_stats <- numeric(n_perm)

  for (b in 0:(n_perm - 1L)) {
    bits <- as.integer(intToBits(b))[seq_len(n)]
    signs <- ifelse(bits == 1L, 1, -1)
    perm_stats[b + 1L] <- mean(signs * d)
  }

  tol <- 1e-12
  c(
    observed = observed,
    p_two_sided = mean(abs(perm_stats) >= abs(observed) - tol),
    p_target_better = mean(perm_stats >= observed - tol),
    n_pairs = n
  )
}

align_year_tables <- function(target_yearly, comparator_yearly) {
  ty <- target_yearly[match(PRIMARY_YEARS, target_yearly$Year), , drop = FALSE]
  cy <- comparator_yearly[match(PRIMARY_YEARS, comparator_yearly$Year), , drop = FALSE]

  if (any(is.na(ty$Year)) || any(is.na(cy$Year)) ||
      !identical(as.integer(ty$Year), as.integer(cy$Year))) {
    stop("Head-to-head year alignment failed.", call. = FALSE)
  }
  list(target = ty, comparator = cy)
}

operational_metric_value <- function(year_df, metric) {
  agg <- aggregate_raw(year_df)
  as.numeric(agg[[metric]][1])
}

operational_h2h_metric <- function(target_yearly,
                                   comparator_yearly,
                                   metric,
                                   seed_offset) {
  aligned <- align_year_tables(target_yearly, comparator_yearly)
  ty <- aligned$target
  cy <- aligned$comparator
  n <- nrow(ty)

  target_est <- operational_metric_value(ty, metric)
  comp_est <- operational_metric_value(cy, metric)
  raw_diff <- target_est - comp_est
  obs <- oriented_difference(target_est, comp_est, metric)

  # Exact paired season-label randomization.
  n_perm <- 2^n
  perm_stats <- numeric(n_perm)

  for (b in 0:(n_perm - 1L)) {
    bits <- as.integer(intToBits(b))[seq_len(n)]
    swap <- bits == 1L

    pt <- ty
    pc <- cy
    if (any(swap)) {
      pt[swap, ] <- cy[swap, ]
      pc[swap, ] <- ty[swap, ]
    }

    tv <- operational_metric_value(pt, metric)
    cv <- operational_metric_value(pc, metric)
    perm_stats[b + 1L] <- oriented_difference(tv, cv, metric)
  }

  tol <- 1e-12
  p_two <- mean(abs(perm_stats) >= abs(obs) - tol)
  p_better <- mean(perm_stats >= obs - tol)

  # Paired season-cluster bootstrap.
  set.seed(GLOBAL_SEED + as.integer(seed_offset))
  boot <- rep(NA_real_, H2H_BOOT_N)
  for (b in seq_len(H2H_BOOT_N)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    tv <- operational_metric_value(ty[idx, , drop = FALSE], metric)
    cv <- operational_metric_value(cy[idx, , drop = FALSE], metric)
    boot[b] <- oriented_difference(tv, cv, metric)
  }

  # Leave-one-season-out jackknife for BCa acceleration.
  jack <- rep(NA_real_, n)
  for (j in seq_len(n)) {
    keep <- setdiff(seq_len(n), j)
    tv <- operational_metric_value(ty[keep, , drop = FALSE], metric)
    cv <- operational_metric_value(cy[keep, , drop = FALSE], metric)
    jack[j] <- oriented_difference(tv, cv, metric)
  }

  ci <- bca_limits(boot, obs, jack, H2H_ALPHA)

  data.frame(
    Metric_Domain = "Operational",
    Metric = metric,
    Direction = H2H_DIRECTION[[metric]],
    Prespecified_Estimate = target_est,
    Comparator_Estimate = comp_est,
    Raw_Difference_TargetMinusComparator = raw_diff,
    Oriented_Difference_PositiveFavorsPrespecified = obs,
    Paired_Season_Effect = obs,
    Bootstrap_Lower95 = as.numeric(ci["lower"]),
    Bootstrap_Upper95 = as.numeric(ci["upper"]),
    CI_Method = as.character(ci["method"]),
    Exact_P_TwoSided = p_two,
    Exact_P_PrespecifiedBetter = p_better,
    N_Paired_Seasons = n,
    stringsAsFactors = FALSE
  )
}

season_auc_table <- function(candidate_id,
                             lead_min,
                             lead_max,
                             burden_frac) {
  rows <- vector("list", length(PRIMARY_YEARS))

  for (i in seq_along(PRIMARY_YEARS)) {
    yr <- PRIMARY_YEARS[i]
    al <- compute_ratio_alignment(
      candidate_id,
      lead_min = lead_min,
      lead_max = lead_max,
      burden_frac = burden_frac,
      years = yr
    )
    rows[[i]] <- data.frame(
      Year = yr,
      A1_AUC = al$A1_AUC,
      A2_AUC = al$A2_AUC,
      A2_CaseWeighted_AUC = al$A2_CaseWeighted_AUC,
      T_AUC = al$T_AUC,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

auc_h2h_metric <- function(target_id,
                           comparator_id,
                           metric,
                           lead_min,
                           lead_max,
                           burden_frac,
                           seed_offset) {
  target_pooled <- compute_ratio_alignment(
    target_id,
    lead_min = lead_min,
    lead_max = lead_max,
    burden_frac = burden_frac
  )
  comp_pooled <- compute_ratio_alignment(
    comparator_id,
    lead_min = lead_min,
    lead_max = lead_max,
    burden_frac = burden_frac
  )

  target_est <- as.numeric(target_pooled[[metric]][1])
  comp_est <- as.numeric(comp_pooled[[metric]][1])
  raw_diff <- target_est - comp_est

  ty <- season_auc_table(
    target_id, lead_min, lead_max, burden_frac
  )
  cy <- season_auc_table(
    comparator_id, lead_min, lead_max, burden_frac
  )
  cy <- cy[match(ty$Year, cy$Year), , drop = FALSE]

  tvals <- ty[[metric]]
  cvals <- cy[[metric]]
  keep <- is.finite(tvals) & is.finite(cvals)
  d <- tvals[keep] - cvals[keep]

  exact <- exact_signflip_mean(d)
  obs <- as.numeric(exact["observed"])
  n <- length(d)

  set.seed(GLOBAL_SEED + as.integer(seed_offset))
  boot <- rep(NA_real_, H2H_BOOT_N)
  for (b in seq_len(H2H_BOOT_N)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    boot[b] <- mean(d[idx])
  }

  jack <- rep(NA_real_, n)
  if (n >= 2L) {
    for (j in seq_len(n)) {
      jack[j] <- mean(d[-j])
    }
  }

  ci <- bca_limits(boot, obs, jack, H2H_ALPHA)

  data.frame(
    Metric_Domain = "AUC",
    Metric = metric,
    Direction = "high",
    Prespecified_Estimate = target_est,
    Comparator_Estimate = comp_est,
    Raw_Difference_TargetMinusComparator = raw_diff,
    Oriented_Difference_PositiveFavorsPrespecified = raw_diff,
    # Inference is paired mean within-season AUC difference.
    Paired_Season_Effect = obs,
    Bootstrap_Lower95 = as.numeric(ci["lower"]),
    Bootstrap_Upper95 = as.numeric(ci["upper"]),
    CI_Method = as.character(ci["method"]),
    Exact_P_TwoSided = as.numeric(exact["p_two_sided"]),
    Exact_P_PrespecifiedBetter = as.numeric(exact["p_target_better"]),
    N_Paired_Seasons = as.integer(exact["n_pairs"]),
    stringsAsFactors = FALSE
  )
}

run_h2h_spec <- function(analysis_set,
                         spec_label,
                         aw_label,
                         eb_label,
                         lead_min,
                         lead_max,
                         burden_frac,
                         seed_base) {
  out <- list()
  k <- 0L

  for (c in seq_len(nrow(H2H_COMPARISONS))) {
    cmp <- H2H_COMPARISONS[c, ]

    target_cand <- CANDIDATES[
      CANDIDATES$Candidate_ID == cmp$Target_ID,
      ,
      drop = FALSE
    ]
    comparator_cand <- CANDIDATES[
      CANDIDATES$Candidate_ID == cmp$Comparator_ID,
      ,
      drop = FALSE
    ]

    if (nrow(target_cand) != 1L || nrow(comparator_cand) != 1L) {
      stop("Head-to-head candidate lookup failed.", call. = FALSE)
    }

    target_ev <- evaluate_candidate(
      cmp$Target_ID,
      lead_min = lead_min,
      lead_max = lead_max,
      burden_frac = burden_frac
    )
    comparator_ev <- evaluate_candidate(
      cmp$Comparator_ID,
      lead_min = lead_min,
      lead_max = lead_max,
      burden_frac = burden_frac
    )

    for (m in METRICS) {
      k <- k + 1L
      res <- operational_h2h_metric(
        target_ev$yearly,
        comparator_ev$yearly,
        metric = m,
        seed_offset = seed_base + 1000L * c + k
      )

      out[[k]] <- cbind(
        data.frame(
          Analysis_Set = analysis_set,
          Specification = spec_label,
          AW = aw_label,
          EB = eb_label,
          LeadMin = lead_min,
          LeadMax = lead_max,
          Burden = burden_frac,
          Detector = cmp$Detector,
          Prespecified_ID = cmp$Target_ID,
          Prespecified_Window = target_cand$Window,
          Comparator_ID = cmp$Comparator_ID,
          Comparator_Window = comparator_cand$Window,
          stringsAsFactors = FALSE
        ),
        res
      )
    }

    for (m in H2H_AUC_METRICS) {
      k <- k + 1L
      res <- auc_h2h_metric(
        cmp$Target_ID,
        cmp$Comparator_ID,
        metric = m,
        lead_min = lead_min,
        lead_max = lead_max,
        burden_frac = burden_frac,
        seed_offset = seed_base + 20000L + 1000L * c + k
      )

      out[[k]] <- cbind(
        data.frame(
          Analysis_Set = analysis_set,
          Specification = spec_label,
          AW = aw_label,
          EB = eb_label,
          LeadMin = lead_min,
          LeadMax = lead_max,
          Burden = burden_frac,
          Detector = cmp$Detector,
          Prespecified_ID = cmp$Target_ID,
          Prespecified_Window = target_cand$Window,
          Comparator_ID = cmp$Comparator_ID,
          Comparator_Window = comparator_cand$Window,
          stringsAsFactors = FALSE
        ),
        res
      )
    }
  }

  ans <- do.call(rbind, out)
  row.names(ans) <- NULL
  ans
}

apply_h2h_multiplicity <- function(tab) {
  out <- tab
  out$P_Bonferroni_WithinMetric <- NA_real_
  out$P_Bonferroni_Family <- NA_real_
  out$Significant_WithinMetric_0_05 <- FALSE
  out$Significant_Familywise_0_05 <- FALSE
  out$CI_Favors <- NA_character_
  out$Conclusion_Familywise <- NA_character_

  # Two same-family alternatives for each metric.
  metric_keys <- interaction(
    out$Analysis_Set,
    out$Specification,
    out$Detector,
    out$Metric,
    drop = TRUE,
    lex.order = TRUE
  )

  for (key in levels(metric_keys)) {
    ii <- which(metric_keys == key)
    finite_n <- sum(is.finite(out$Exact_P_TwoSided[ii]))
    if (finite_n > 0L) {
      out$P_Bonferroni_WithinMetric[ii] <- pmin(
        1,
        out$Exact_P_TwoSided[ii] * finite_n
      )
    }
  }

  # Confirmatory family: all operational + AUC comparisons for one detector
  # under one anchor specification.
  family_keys <- interaction(
    out$Analysis_Set,
    out$Specification,
    out$Detector,
    drop = TRUE,
    lex.order = TRUE
  )

  for (key in levels(family_keys)) {
    ii <- which(family_keys == key)
    finite_n <- sum(is.finite(out$Exact_P_TwoSided[ii]))
    if (finite_n > 0L) {
      out$P_Bonferroni_Family[ii] <- pmin(
        1,
        out$Exact_P_TwoSided[ii] * finite_n
      )
    }
  }

  out$Significant_WithinMetric_0_05 <-
    is.finite(out$P_Bonferroni_WithinMetric) &
    out$P_Bonferroni_WithinMetric < H2H_ALPHA

  out$Significant_Familywise_0_05 <-
    is.finite(out$P_Bonferroni_Family) &
    out$P_Bonferroni_Family < H2H_ALPHA

  out$CI_Favors <- ifelse(
    is.finite(out$Bootstrap_Lower95) & out$Bootstrap_Lower95 > 0,
    "Prespecified",
    ifelse(
      is.finite(out$Bootstrap_Upper95) & out$Bootstrap_Upper95 < 0,
      "Comparator",
      "No clear difference"
    )
  )

  out$Conclusion_Familywise <- ifelse(
    out$Significant_Familywise_0_05 &
      out$Oriented_Difference_PositiveFavorsPrespecified > 0,
    "Prespecified significantly better",
    ifelse(
      out$Significant_Familywise_0_05 &
        out$Oriented_Difference_PositiveFavorsPrespecified < 0,
      "Comparator significantly better",
      "No familywise-significant difference"
    )
  )

  out
}

message("[Stage1B] Running primary head-to-head significance...")
H2H_PRIMARY <- apply_h2h_multiplicity(
  run_h2h_spec(
    analysis_set = "Primary",
    spec_label = "AW4-8_EB70",
    aw_label = "4-8",
    eb_label = "70%",
    lead_min = PRIMARY_AW_MIN,
    lead_max = PRIMARY_AW_MAX,
    burden_frac = PRIMARY_EB,
    seed_base = 100000L
  )
)

message("[Stage1B] Running AW-sensitivity head-to-head significance...")
H2H_AW_ROWS <- vector("list", nrow(AW_LEVELS))
for (a in seq_len(nrow(AW_LEVELS))) {
  H2H_AW_ROWS[[a]] <- run_h2h_spec(
    analysis_set = "AW_sensitivity",
    spec_label = paste0("AW", AW_LEVELS$AW[a], "_EB70"),
    aw_label = AW_LEVELS$AW[a],
    eb_label = "70%",
    lead_min = AW_LEVELS$LeadMin[a],
    lead_max = AW_LEVELS$LeadMax[a],
    burden_frac = PRIMARY_EB,
    seed_base = 200000L + 20000L * a
  )
}
H2H_AW <- apply_h2h_multiplicity(do.call(rbind, H2H_AW_ROWS))

message("[Stage1B] Running EB-sensitivity head-to-head significance...")
H2H_EB_ROWS <- vector("list", nrow(EB_LEVELS))
for (b in seq_len(nrow(EB_LEVELS))) {
  H2H_EB_ROWS[[b]] <- run_h2h_spec(
    analysis_set = "EB_sensitivity",
    spec_label = paste0("AW4-8_EB", EB_LEVELS$EB[b]),
    aw_label = "4-8",
    eb_label = EB_LEVELS$EB[b],
    lead_min = PRIMARY_AW_MIN,
    lead_max = PRIMARY_AW_MAX,
    burden_frac = EB_LEVELS$Burden[b],
    seed_base = 300000L + 20000L * b
  )
}
H2H_EB <- apply_h2h_multiplicity(do.call(rbind, H2H_EB_ROWS))

message("[Stage1B] Running full AW x EB supplementary head-to-head significance...")
H2H_FULL_ROWS <- list()
hk <- 0L
for (a in seq_len(nrow(AW_LEVELS))) {
  for (b in seq_len(nrow(EB_LEVELS))) {
    hk <- hk + 1L
    H2H_FULL_ROWS[[hk]] <- run_h2h_spec(
      analysis_set = "Full_AW_x_EB",
      spec_label = paste0("AW", AW_LEVELS$AW[a], "_EB", EB_LEVELS$EB[b]),
      aw_label = AW_LEVELS$AW[a],
      eb_label = EB_LEVELS$EB[b],
      lead_min = AW_LEVELS$LeadMin[a],
      lead_max = AW_LEVELS$LeadMax[a],
      burden_frac = EB_LEVELS$Burden[b],
      seed_base = 400000L + 30000L * hk
    )
  }
}
H2H_FULL <- apply_h2h_multiplicity(do.call(rbind, H2H_FULL_ROWS))

H2H_MASTER <- rbind(
  H2H_PRIMARY,
  H2H_AW,
  H2H_EB,
  H2H_FULL
)
row.names(H2H_MASTER) <- NULL

# Compact target-centric significance summary.
H2H_SUMMARY <- do.call(rbind, lapply(
  split(
    H2H_MASTER,
    interaction(
      H2H_MASTER$Analysis_Set,
      H2H_MASTER$Specification,
      H2H_MASTER$Detector,
      H2H_MASTER$Prespecified_Window,
      drop = TRUE,
      lex.order = TRUE
    )
  ),
  function(z) {
    data.frame(
      Analysis_Set = z$Analysis_Set[1],
      Specification = z$Specification[1],
      AW = z$AW[1],
      EB = z$EB[1],
      Detector = z$Detector[1],
      Prespecified_Window = z$Prespecified_Window[1],
      N_HeadToHead_Tests = nrow(z),
      N_Familywise_Significant_Better = sum(
        z$Conclusion_Familywise == "Prespecified significantly better",
        na.rm = TRUE
      ),
      N_Familywise_Significant_Worse = sum(
        z$Conclusion_Familywise == "Comparator significantly better",
        na.rm = TRUE
      ),
      N_No_Familywise_Difference = sum(
        z$Conclusion_Familywise == "No familywise-significant difference",
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }
))

message("[Stage1B] Head-to-head significance calculations: PASS")

# -----------------------------------------------------------------------------
# 16. SEASON-CLUSTER BOOTSTRAP WINDOW-SELECTION STABILITY
# -----------------------------------------------------------------------------
# Primary AW=4-8 / EB=70%; thresholds remain fixed. Each replicate resamples
# complete seasons and recomputes the operational rank. Ratio-alignment scores
# are fixed at the full primary-series estimate; the bootstrap therefore tests
# stability of the operational anchor capture while preserving direct ratio
# evidence in the combined score.

# Build per-candidate per-year tables under primary anchors.
YEAR_BY_CAND <- split(PRIMARY_YEARLY, PRIMARY_YEARLY$Candidate_ID)

bootstrap_detector <- function(detector, B = BOOT_N) {
  cands <- CANDIDATES[CANDIDATES$Detector == detector, , drop = FALSE]
  n_c <- nrow(cands)

  full_align <- PRIMARY_RESULTS[
    match(cands$Candidate_ID, PRIMARY_RESULTS$Candidate_ID),
    ,
    drop = FALSE
  ]

  # Alignment ranks from primary data.
  align_rank <- full_align$Alignment_Mean_Rank

  win1 <- setNames(integer(n_c), cands$Candidate_ID)
  top2 <- setNames(integer(n_c), cands$Candidate_ID)
  target_rank <- rep(NA_real_, B)

  set.seed(GLOBAL_SEED + ifelse(detector == "Continuous TA", 501L, 601L))

  for (b in seq_len(B)) {
    sampled_years <- sample(PRIMARY_YEARS, length(PRIMARY_YEARS), replace = TRUE)

    metric_rows <- vector("list", n_c)
    for (i in seq_len(n_c)) {
      z <- YEAR_BY_CAND[[cands$Candidate_ID[i]]]
      z <- z[match(sampled_years, z$Year), , drop = FALSE]
      metric_rows[[i]] <- aggregate_raw(z)
    }
    met <- do.call(rbind, metric_rows)

    op_rank_mat <- matrix(NA_real_, nrow = n_c, ncol = length(METRICS))
    for (j in seq_along(METRICS)) {
      m <- METRICS[j]
      op_rank_mat[, j] <- rank_metric(met[[m]], DIRECTION[[m]])
    }
    op_rank <- rowMeans(op_rank_mat, na.rm = TRUE)
    combined <- (op_rank + align_rank) / 2
    final_rank <- rank(combined, ties.method = "average")

    ord <- order(final_rank, cands$STA, cands$LTA)
    win1[cands$Candidate_ID[ord[1]]] <- win1[cands$Candidate_ID[ord[1]]] + 1L
    top_ids <- cands$Candidate_ID[ord[seq_len(min(2L, n_c))]]
    top2[top_ids] <- top2[top_ids] + 1L

    ti <- which(cands$Prespecified)
    target_rank[b] <- final_rank[ti]
  }

  data.frame(
    Candidate_ID = cands$Candidate_ID,
    Detector = detector,
    Window = cands$Window,
    STA = cands$STA,
    LTA = cands$LTA,
    Prespecified = cands$Prespecified,
    P_Rank1 = unname(win1[cands$Candidate_ID]) / B,
    P_Top2 = unname(top2[cands$Candidate_ID]) / B,
    Target_Rank_Median = ifelse(
      cands$Prespecified, stats::median(target_rank, na.rm = TRUE), NA_real_
    ),
    Target_Rank_Lower95 = ifelse(
      cands$Prespecified, safe_quantile(target_rank, 0.025), NA_real_
    ),
    Target_Rank_Upper95 = ifelse(
      cands$Prespecified, safe_quantile(target_rank, 0.975), NA_real_
    ),
    stringsAsFactors = FALSE
  )
}

BOOT_SELECTION <- rbind(
  bootstrap_detector("Continuous TA", BOOT_N),
  bootstrap_detector("Constant TA", BOOT_N)
)

# -----------------------------------------------------------------------------
# 17. LEAVE-ONE-YEAR-OUT WINDOW-SELECTION STABILITY
# -----------------------------------------------------------------------------
# No parameter is re-estimated here; this is a direct influence analysis:
# score the three candidate windows on 9 of 10 seasons and record the rank of
# the prespecified target. This asks whether one season drives the window choice.
LOYO_ROWS <- list()
k <- 0L

for (omit in PRIMARY_YEARS) {
  keep_years <- setdiff(PRIMARY_YEARS, omit)

  for (det in c("Continuous TA", "Constant TA")) {
    cands <- CANDIDATES[CANDIDATES$Detector == det, , drop = FALSE]
    rows <- vector("list", nrow(cands))

    for (i in seq_len(nrow(cands))) {
      cnd <- cands[i, ]
      z <- YEAR_BY_CAND[[cnd$Candidate_ID]]
      z <- z[z$Year %in% keep_years, , drop = FALSE]
      met <- aggregate_raw(z)

      # Ratio alignment is recomputed using the same anchor definitions but the
      # ratio evidence itself is not year-restricted in this compact LOYO.
      # The operational component is the influence diagnostic of interest.
      rows[[i]] <- cbind(cnd, met)
    }

    tab <- do.call(rbind, rows)
    op <- matrix(NA_real_, nrow = nrow(tab), ncol = length(METRICS))
    for (j in seq_along(METRICS)) {
      m <- METRICS[j]
      op[, j] <- rank_metric(tab[[m]], DIRECTION[[m]])
    }
    tab$Operational_Mean_Rank <- rowMeans(op, na.rm = TRUE)
    tab$Final_Rank <- rank(tab$Operational_Mean_Rank, ties.method = "average")

    ti <- which(tab$Prespecified)
    bi <- which.min(tab$Operational_Mean_Rank)

    k <- k + 1L
    LOYO_ROWS[[k]] <- data.frame(
      Omitted_Year = omit,
      Detector = det,
      Target_Window = tab$Window[ti],
      Target_Rank = tab$Final_Rank[ti],
      Best_Window = tab$Window[bi],
      Best_STA = tab$STA[bi],
      Best_LTA = tab$LTA[bi],
      Target_Is_Best = tab$Candidate_ID[ti] == tab$Candidate_ID[bi],
      stringsAsFactors = FALSE
    )
  }
}

LOYO <- do.call(rbind, LOYO_ROWS)

# -----------------------------------------------------------------------------
# 18. ROBUSTNESS SUMMARY OF TARGET WINDOWS
# -----------------------------------------------------------------------------
target_robustness <- function(detector, target_window) {
  z <- FULL_ROBUST_RANKED[
    FULL_ROBUST_RANKED$Detector == detector &
      FULL_ROBUST_RANKED$Window == target_window,
    ,
    drop = FALSE
  ]

  data.frame(
    Detector = detector,
    Target_Window = target_window,
    N_AW_EB_Specifications = nrow(z),
    Mean_Rank = mean(z$Anchor_Focused_Final_Rank, na.rm = TRUE),
    Median_Rank = stats::median(z$Anchor_Focused_Final_Rank, na.rm = TRUE),
    Best_Rank = min(z$Anchor_Focused_Final_Rank, na.rm = TRUE),
    Worst_Rank = max(z$Anchor_Focused_Final_Rank, na.rm = TRUE),
    Rank1_Fraction = mean(z$Anchor_Focused_Final_Rank == 1, na.rm = TRUE),
    Top2_Fraction = mean(z$Anchor_Focused_Final_Rank <= 2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

TARGET_ROBUSTNESS <- rbind(
  target_robustness("Continuous TA", "3/12"),
  target_robustness("Constant TA", "4/26")
)

TARGET_SUMMARY <- merge(
  PRIMARY_RESULTS[PRIMARY_RESULTS$Prespecified, c(
    "Detector", "Window", "STA", "LTA",
    METRICS,
    "A1_AUC", "A2_CaseWeighted_AUC", "T_AUC",
    "A1_Enrichment", "A2_Enrichment",
    "Operational_Mean_Rank", "Alignment_Mean_Rank",
    "Anchor_Focused_Mean_Rank", "Anchor_Focused_Final_Rank"
  )],
  TARGET_ROBUSTNESS,
  by.x = c("Detector", "Window"),
  by.y = c("Detector", "Target_Window"),
  all.x = TRUE
)

TARGET_SUMMARY <- merge(
  TARGET_SUMMARY,
  BOOT_SELECTION[
    BOOT_SELECTION$Prespecified,
    c("Detector", "Window", "P_Rank1", "P_Top2",
      "Target_Rank_Median", "Target_Rank_Lower95", "Target_Rank_Upper95")
  ],
  by = c("Detector", "Window"),
  all.x = TRUE
)

# -----------------------------------------------------------------------------
# 19. CONSTANT GUARD/FREEZE AUDIT
# -----------------------------------------------------------------------------
CONST_AUDIT_ROWS <- list()
k <- 0L

for (i in seq_len(nrow(CONST_CANDIDATES))) {
  cnd <- CONST_CANDIDATES[i, ]
  obj <- DETECTOR_OBJECTS[[cnd$Candidate_ID]]

  ids <- sort(unique(obj$freeze_id[is.finite(obj$freeze_id)]))
  invariant <- TRUE
  activation_match <- TRUE

  for (id in ids) {
    ii <- which(obj$freeze_id == id)
    vals <- obj$LTA_frozen[ii]
    vals <- vals[is.finite(vals)]
    if (length(vals) > 1L && max(vals) - min(vals) > 1e-12) invariant <- FALSE

    aa <- ii[obj$activation[ii]]
    if (length(aa) != 1L ||
        !isTRUE(all.equal(
          obj$LTA_frozen[aa], obj$LTA_live[aa], tolerance = 1e-12
        ))) {
      activation_match <- FALSE
    }
  }

  k <- k + 1L
  CONST_AUDIT_ROWS[[k]] <- data.frame(
    Window = cnd$Window,
    STA = cnd$STA,
    LTA = cnd$LTA,
    Guard = CONST_GUARD,
    Baseline_Freeze = TRUE,
    Min_OFF_Reset = CONST_MIN_OFF_RESET,
    N_Freeze_Episodes = length(ids),
    Frozen_LTA_Invariant = invariant,
    Freeze_Equals_Guarded_LTA_At_Activation = activation_match,
    PASS = invariant && activation_match,
    stringsAsFactors = FALSE
  )
}

CONST_AUDIT <- do.call(rbind, CONST_AUDIT_ROWS)
if (!all(CONST_AUDIT$PASS)) {
  stop("Constant-TA guard/freeze audit failed.", call. = FALSE)
}
message("[Stage1B] Constant guard/freeze audit: PASS")

# -----------------------------------------------------------------------------
# 20. FIGURE: PRIMARY ANCHOR-FOCUSED WINDOW RANKS
# -----------------------------------------------------------------------------
render_rank_figure <- function() {
  draw <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))

    for (det in c("Continuous TA", "Constant TA")) {
      z <- PRIMARY_RESULTS[PRIMARY_RESULTS$Detector == det, , drop = FALSE]
      z <- z[order(z$STA, z$LTA), , drop = FALSE]

      plot(
        seq_len(nrow(z)),
        z$Anchor_Focused_Final_Rank,
        type = "b",
        pch = 19,
        xaxt = "n",
        ylim = c(0.8, 3.2),
        xlab = "STA/LTA window",
        ylab = "Anchor-focused rank (lower is better)",
        main = det
      )
      axis(1, at = seq_len(nrow(z)), labels = z$Window)
      pri <- which(z$Prespecified)
      if (length(pri) == 1L) {
        points(pri, z$Anchor_Focused_Final_Rank[pri],
               pch = 8, cex = 1.8, lwd = 2)
        text(pri, z$Anchor_Focused_Final_Rank[pri],
             labels = " Primary", pos = 4, cex = 0.85)
      }
      grid()
    }
  }

  pdf_path <- file.path(OUT_DIR, "Figure_Stage1B_anchor_focused_window_rank.pdf")
  grDevices::pdf(pdf_path, width = 10, height = 4.8)
  tryCatch(draw(), finally = grDevices::dev.off())

  png_path <- file.path(OUT_DIR, "Figure_Stage1B_anchor_focused_window_rank.png")
  grDevices::png(png_path, width = 1800, height = 850, res = 180)
  tryCatch(draw(), finally = grDevices::dev.off())

  if (!file.exists(pdf_path) || !file.exists(png_path)) {
    stop("Anchor-focused figure did not render.", call. = FALSE)
  }
}

render_rank_figure()

# -----------------------------------------------------------------------------
# 21. WRITE OUTPUTS
# -----------------------------------------------------------------------------
utils::write.csv(
  CANDIDATES,
  file.path(OUT_DIR, "00_candidate_windows.csv"),
  row.names = FALSE
)

utils::write.csv(
  PRIMARY_RESULTS,
  file.path(OUT_DIR, "01_primary_AW4_8_EB70_anchor_focused_scores.csv"),
  row.names = FALSE
)

utils::write.csv(
  PRIMARY_YEARLY,
  file.path(OUT_DIR, "02_primary_AW4_8_EB70_per_year_metrics.csv"),
  row.names = FALSE
)

utils::write.csv(
  AW_ROBUST_RANKED,
  file.path(OUT_DIR, "03_AW_robustness_EB70_fixed.csv"),
  row.names = FALSE
)

utils::write.csv(
  EB_ROBUST_RANKED,
  file.path(OUT_DIR, "04_EB_robustness_AW4_8_fixed.csv"),
  row.names = FALSE
)

utils::write.csv(
  FULL_ROBUST_RANKED,
  file.path(OUT_DIR, "05_full_AW_x_EB_robustness.csv"),
  row.names = FALSE
)

utils::write.csv(
  BOOT_SELECTION,
  file.path(OUT_DIR, "06_bootstrap_anchor_focused_window_selection.csv"),
  row.names = FALSE
)

utils::write.csv(
  LOYO,
  file.path(OUT_DIR, "07_LOYO_window_influence.csv"),
  row.names = FALSE
)

utils::write.csv(
  TARGET_ROBUSTNESS,
  file.path(OUT_DIR, "08_target_AW_EB_robustness_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  TARGET_SUMMARY,
  file.path(OUT_DIR, "09_target_3_12_and_4_26_evidence_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  CONST_AUDIT,
  file.path(OUT_DIR, "10_constant_guard_freeze_audit.csv"),
  row.names = FALSE
)

# Head-to-head significance outputs.
utils::write.csv(
  H2H_PRIMARY[H2H_PRIMARY$Metric_Domain == "Operational", , drop = FALSE],
  file.path(OUT_DIR, "11_H2H_primary_operational_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  H2H_PRIMARY[H2H_PRIMARY$Metric_Domain == "AUC", , drop = FALSE],
  file.path(OUT_DIR, "12_H2H_primary_AUC_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  H2H_AW[H2H_AW$Metric_Domain == "Operational", , drop = FALSE],
  file.path(OUT_DIR, "13_H2H_AW_sensitivity_operational_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  H2H_AW[H2H_AW$Metric_Domain == "AUC", , drop = FALSE],
  file.path(OUT_DIR, "14_H2H_AW_sensitivity_AUC_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  H2H_EB[H2H_EB$Metric_Domain == "Operational", , drop = FALSE],
  file.path(OUT_DIR, "15_H2H_EB_sensitivity_operational_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  H2H_EB[H2H_EB$Metric_Domain == "AUC", , drop = FALSE],
  file.path(OUT_DIR, "16_H2H_EB_sensitivity_AUC_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  H2H_FULL,
  file.path(OUT_DIR, "17_H2H_full_AW_x_EB_all_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  H2H_MASTER,
  file.path(OUT_DIR, "18_H2H_all_significance_master.csv"),
  row.names = FALSE
)
utils::write.csv(
  H2H_SUMMARY,
  file.path(OUT_DIR, "19_H2H_significance_summary.csv"),
  row.names = FALSE
)

# Sourceable handoff: windows are PRESPECIFIED targets under empirical
# anchor-focused validation. Do not mislabel them as forced optimizer outputs.
handoff <- c(
  "# Generated by Stage1b_anchor_focused_window_selection_QC.R",
  sprintf("CONTINUOUS_STA_WIN_ADOPTED <- %dL", 3L),
  sprintf("CONTINUOUS_LTA_WIN_ADOPTED <- %dL", 12L),
  sprintf("CONSTANT_STA_WIN_ADOPTED <- %dL", 4L),
  sprintf("CONSTANT_LTA_WIN_ADOPTED <- %dL", 26L),
  sprintf("CONSTANT_GUARD_ADOPTED <- %dL", CONST_GUARD),
  sprintf("CONSTANT_MIN_OFF_RESET_ADOPTED <- %dL", CONST_MIN_OFF_RESET),
  "CONSTANT_BASELINE_FREEZE_REQUIRED <- TRUE",
  "WINDOW_VALIDATION_FRAMEWORK <- \"A1_ACTIONABLE_WINDOW_PLUS_A2_EPIDEMIC_BURDEN\"",
  "PRIMARY_ACTIONABLE_WINDOW <- \"4-8 weeks before peak\"",
  "PRIMARY_EPIDEMIC_BURDEN <- 0.70"
)
writeLines(
  handoff,
  file.path(OUT_DIR, "sta_lta_anchor_focused_handoff.R")
)

# -----------------------------------------------------------------------------
# 22. REPRODUCIBILITY / RUNTIME CONTRACTS
# -----------------------------------------------------------------------------
if (nrow(PRIMARY_RESULTS) != 6L) {
  stop("Primary result must contain exactly 6 window configurations.",
       call. = FALSE)
}
if (nrow(AW_ROBUST_RANKED) != 18L) {
  stop("AW robustness must contain 18 rows (3 AW x 6 windows).",
       call. = FALSE)
}
if (nrow(EB_ROBUST_RANKED) != 18L) {
  stop("EB robustness must contain 18 rows (3 EB x 6 windows).",
       call. = FALSE)
}
if (nrow(FULL_ROBUST_RANKED) != 54L) {
  stop("Full AWxEB robustness must contain 54 rows (3x3x6).",
       call. = FALSE)
}
if (nrow(BOOT_SELECTION) != 6L) {
  stop("Bootstrap selection table must contain 6 rows.", call. = FALSE)
}
if (nrow(LOYO) != 20L) {
  stop("LOYO influence table must contain 20 rows.", call. = FALSE)
}
if (nrow(TARGET_SUMMARY) != 2L) {
  stop("Target evidence summary must contain 2 rows.", call. = FALSE)
}

# Head-to-head deterministic row contracts:
# 4 comparisons x (8 operational + 4 AUC) = 48 rows per specification.
EXPECTED_H2H_PER_SPEC <- nrow(H2H_COMPARISONS) * length(H2H_ALL_METRICS)
if (nrow(H2H_PRIMARY) != EXPECTED_H2H_PER_SPEC) {
  stop("Primary H2H row count mismatch.", call. = FALSE)
}
if (nrow(H2H_AW) != nrow(AW_LEVELS) * EXPECTED_H2H_PER_SPEC) {
  stop("AW-sensitivity H2H row count mismatch.", call. = FALSE)
}
if (nrow(H2H_EB) != nrow(EB_LEVELS) * EXPECTED_H2H_PER_SPEC) {
  stop("EB-sensitivity H2H row count mismatch.", call. = FALSE)
}
if (nrow(H2H_FULL) !=
    nrow(AW_LEVELS) * nrow(EB_LEVELS) * EXPECTED_H2H_PER_SPEC) {
  stop("Full AWxEB H2H row count mismatch.", call. = FALSE)
}

# Every head-to-head row must have an exact p-value and finite comparison count.
if (any(!is.finite(H2H_MASTER$Exact_P_TwoSided))) {
  stop("At least one H2H exact two-sided p-value is non-finite.", call. = FALSE)
}
if (any(H2H_MASTER$N_Paired_Seasons < 1L)) {
  stop("At least one H2H comparison has no paired seasons.", call. = FALSE)
}
if (any(
  H2H_MASTER$P_Bonferroni_Family < 0 |
  H2H_MASTER$P_Bonferroni_Family > 1,
  na.rm = TRUE
)) {
  stop("H2H Bonferroni-adjusted p-value outside [0,1].", call. = FALSE)
}

message(
  "[Stage1B] H2H runtime contracts: PASS; primary=",
  nrow(H2H_PRIMARY), ", AW=", nrow(H2H_AW),
  ", EB=", nrow(H2H_EB), ", full=", nrow(H2H_FULL)
)

required_files <- c(
  "00_candidate_windows.csv",
  "01_primary_AW4_8_EB70_anchor_focused_scores.csv",
  "02_primary_AW4_8_EB70_per_year_metrics.csv",
  "03_AW_robustness_EB70_fixed.csv",
  "04_EB_robustness_AW4_8_fixed.csv",
  "05_full_AW_x_EB_robustness.csv",
  "06_bootstrap_anchor_focused_window_selection.csv",
  "07_LOYO_window_influence.csv",
  "08_target_AW_EB_robustness_summary.csv",
  "09_target_3_12_and_4_26_evidence_summary.csv",
  "10_constant_guard_freeze_audit.csv",
  "11_H2H_primary_operational_metrics.csv",
  "12_H2H_primary_AUC_metrics.csv",
  "13_H2H_AW_sensitivity_operational_metrics.csv",
  "14_H2H_AW_sensitivity_AUC_metrics.csv",
  "15_H2H_EB_sensitivity_operational_metrics.csv",
  "16_H2H_EB_sensitivity_AUC_metrics.csv",
  "17_H2H_full_AW_x_EB_all_metrics.csv",
  "18_H2H_all_significance_master.csv",
  "19_H2H_significance_summary.csv",
  "sta_lta_anchor_focused_handoff.R",
  "Figure_Stage1B_anchor_focused_window_rank.pdf",
  "Figure_Stage1B_anchor_focused_window_rank.png"
)

missing_outputs <- required_files[
  !file.exists(file.path(OUT_DIR, required_files))
]
if (length(missing_outputs)) {
  stop(
    "Final output contract failed. Missing: ",
    paste(missing_outputs, collapse = ", "),
    call. = FALSE
  )
}

manifest <- data.frame(
  Item = c(
    "Version",
    "Input",
    "Input_MD5",
    "Seed",
    "Bootstrap_N",
    "Primary_years",
    "Excluded_years",
    "Primary_AW",
    "Primary_EB",
    "AW_robustness",
    "EB_robustness",
    "Continuous_candidates",
    "Constant_candidates",
    "eta_ON",
    "eta_OFF",
    "Constant_guard",
    "Constant_min_off_reset",
    "Selection_framework",
    "H2H_bootstrap_N",
    "H2H_exact_test",
    "H2H_multiplicity",
    "H2H_AUC_metrics",
    "Output_dir"
  ),
  Value = c(
    STAGE1B_VERSION,
    DATA_FILE,
    unname(tools::md5sum(DATA_FILE)),
    GLOBAL_SEED,
    BOOT_N,
    paste(PRIMARY_YEARS, collapse = ","),
    paste(EXCLUDED_YEARS, collapse = ","),
    "4-8",
    "70%",
    "3-6 | 4-8 | 5-10",
    "60% | 70% | 80%",
    paste(CONT_CANDIDATES$Window, collapse = " | "),
    paste(CONST_CANDIDATES$Window, collapse = " | "),
    ETA_ON,
    ETA_OFF,
    CONST_GUARD,
    CONST_MIN_OFF_RESET,
    "A1/A2 anchor-focused operational metrics + threshold-independent ratio alignment",
    H2H_BOOT_N,
    "paired season exact 2^10 randomization/sign-flip",
    "Bonferroni within metric plus detector x specification family",
    paste(H2H_AUC_METRICS, collapse = " | "),
    OUT_DIR
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  manifest,
  file.path(OUT_DIR, "00_reproducibility_manifest.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "00_sessionInfo.txt")
)

if (nzchar(SCRIPT_FILE) && file.exists(SCRIPT_FILE)) {
  writeLines(
    paste0("MD5: ", unname(tools::md5sum(SCRIPT_FILE))),
    file.path(OUT_DIR, "00_script_md5.txt")
  )
}

message("[Stage1B] Final output contract: PASS")

cat("\n")
cat("====================================================================\n")
cat("STAGE 1B ANCHOR-FOCUSED STA/LTA WINDOW ANALYSIS COMPLETE\n")
cat("====================================================================\n")
cat("Window selection objective:\n")
cat("  capture A1 Actionable Window + A2 Epidemic-Burden framework.\n\n")
cat("Primary anchors:\n")
cat("  A1 = 4-8 weeks before annual peak\n")
cat("  A2 = contiguous block containing 70% of annual cases\n")
cat("  T  = A1 union A2\n\n")
cat("Continuous candidates: 2/8 | 3/12 PRIMARY | 4/16\n")
cat("Constant candidates:   2/13 | 3/20 | 4/26 PRIMARY\n")
cat("Constant architecture: 2-week guard + activation baseline freeze + 8 OFF reset\n\n")
cat("Selection evidence:\n")
cat("  Stage-3 operational metrics anchored to A1/A2\n")
cat("  threshold-independent A1/A2/T ratio alignment\n")
cat("  prespecified-vs-alternative H2H exact paired significance\n")
cat("  2,000 season-cluster BCa H2H confidence intervals\n")
cat("  Bonferroni within-metric and family-wise multiplicity control\n")
cat("  1,000 season-cluster bootstrap selection stability\n")
cat("  leave-one-year-out influence\n")
cat("  AW robustness: 3-6 / 4-8 / 5-10\n")
cat("  EB robustness: 60% / 70% / 80%\n")
cat("  full 3x3 AW x EB supplementary robustness\n\n")
cat("Output: ", OUT_DIR, "\n", sep = "")
cat("====================================================================\n")
