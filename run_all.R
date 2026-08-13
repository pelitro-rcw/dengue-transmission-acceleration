# =============================================================================
# RUN ALL  --  reproduce every manuscript figure and table
# -----------------------------------------------------------------------------
# Usage
#   From R      :  setwd("<project root>"); source("run_all.R")
#   From a shell:  Rscript run_all.R
#
# Prerequisite: R/00_config.R, R/01_publication_theme.R, and (Stages 4-5)
# R/02_detection_framework.R must be present in R/, alongside the dataset
# file, before this script is sourced -- each Stage script provisions any
# missing CRAN packages itself on first run via require_packages().
#
# Stages must run in order: Stage 1a derives the eta_ON / eta_OFF thresholds
# that Stage 3 consumes. Stage 1b is an independent window-selection robustness
# check (Supplementary Table 21) and does not feed any other stage. Stages 4
# and 5 are independent of 1-3 and of each other, but are run last so the
# console log follows manuscript figure order.
# =============================================================================

t_start <- Sys.time()

# --- Locate the project root and load configuration --------------------------
if (!dir.exists("R") || !dir.exists("scripts")) {
  stop("run_all.R must be sourced from the project root ",
       "(the folder containing R/ and scripts/).", call. = FALSE)
}
Sys.setenv(TA_PROJECT_ROOT = normalizePath(getwd(), winslash = "/"))
source(file.path("R", "00_config.R"))

# --- Stage registry ----------------------------------------------------------
STAGES <- list(
  list(id = "Stage 1a", file = "Stage1a_eta_threshold_derivation_QC.R",
       desc = "Empirical eta_ON / eta_OFF derivation (Quezon City)"),
  list(id = "Stage 1b", file = "Stage1b_anchor_focused_window_selection_QC.R",
       desc = "Window-selection robustness check (Supplementary Table 21)",
       optional = TRUE),
  list(id = "Stage 2", file = "Stage2_outbreak_threshold_drift_QC.R",
       desc = "Outbreak threshold drift (Figure 1)"),
  list(id = "Stage 3", file = "Stage3_QC_analysis.R",
       desc = "Quezon City analysis (Figures 2 and 3, Tables 1-2)"),
  list(id = "Stage 4", file = "Stage4_Regional_Analysis.R",
       desc = "Regional analysis (Figures 4 and 5, incl. the map)"),
  list(id = "Stage 5", file = "Stage5_Country_Analysis.R",
       desc = "Country analysis (Figures 6 and 7)")
)

# Stage 1b does not feed any downstream stage; set to FALSE to skip it.
RUN_OPTIONAL_STAGES <- TRUE

# --- Execution ---------------------------------------------------------------
# Each stage runs in its own environment so that variables cannot leak between
# stages -- important because several stages reuse names such as OUT_DIR, df and
# theme_dashboard with different meanings.
results <- list()

for (st in STAGES) {
  if (isTRUE(st$optional) && !isTRUE(RUN_OPTIONAL_STAGES)) {
    cat("\n[skip] ", st$id, " (optional)\n", sep = "")
    next
  }
  # DIR_SCRIPTS is defined by R/00_config.R as <project root>/scripts, which
  # is where the Stage*.R scripts live in this package.
  path <- file.path(DIR_SCRIPTS, st$file)
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat(st$id, " -- ", st$desc, "\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")

  if (!file.exists(path)) {
    stop("Missing stage script: ", path, call. = FALSE)
  }

  t0 <- Sys.time()
  env <- new.env(parent = globalenv())
  ok  <- TRUE
  err <- NULL

  tryCatch(
    sys.source(path, envir = env, keep.source = FALSE),
    error = function(e) {
      ok  <<- FALSE
      err <<- conditionMessage(e)
    }
  )

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  results[[st$id]] <- list(ok = ok, secs = elapsed, error = err)

  if (ok) {
    cat("\n[", st$id, "] completed in ", elapsed, " s\n", sep = "")
  } else {
    cat("\n[", st$id, "] FAILED after ", elapsed, " s\n", sep = "")
    cat("  ", err, "\n", sep = "")
    # Stage 3 depends on Stage 1; stopping early avoids cascading failures that
    # obscure the original error.
    if (identical(st$id, "Stage 1a")) {
      stop("Stage 1a failed; Stage 3 depends on its output. Fix and re-run.",
           call. = FALSE)
    }
  }
  rm(env); invisible(gc(verbose = FALSE))
}

# --- Summary -----------------------------------------------------------------
cat("\n", strrep("=", 78), "\n", sep = "")
cat("RUN SUMMARY\n")
cat(strrep("=", 78), "\n", sep = "")
for (nm in names(results)) {
  r <- results[[nm]]
  cat(sprintf("  %-9s %-8s %6.1f s\n", nm, if (r$ok) "OK" else "FAILED", r$secs))
}

figs <- list.files(DIR_OUTPUT, pattern = "\\.(pdf|png)$",
                   recursive = TRUE, full.names = FALSE)
csvs <- list.files(DIR_OUTPUT, pattern = "\\.csv$",
                   recursive = TRUE, full.names = FALSE)

cat("\n  Figures written : ", length(figs), " files (",
    sum(grepl("\\.pdf$", figs)), " PDF + ",
    sum(grepl("\\.png$", figs)), " PNG)\n", sep = "")
cat("  Tables written  : ", length(csvs), " CSV files\n", sep = "")
cat("  Output root     : ", DIR_OUTPUT, "\n", sep = "")
cat("  Total elapsed   : ",
    round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1),
    " min\n", sep = "")

failed <- names(results)[!vapply(results, `[[`, logical(1), "ok")]
if (length(failed) > 0) {
  cat("\n  Stages with errors: ", paste(failed, collapse = ", "), "\n", sep = "")
} else {
  cat("\n  All stages completed successfully.\n")
}

cat("\nSession information\n")
print(utils::sessionInfo())
