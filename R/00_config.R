# =============================================================================
# PROJECT CONFIGURATION
# -----------------------------------------------------------------------------
# Portable paths and shared constants. Sourced first by every stage script.
#
# Replaces the hard-coded "C:/Users/User/Desktop/..." paths that appeared in all
# five original scripts. The project root is discovered automatically, so the
# project runs unmodified on Windows, macOS and Linux from any location.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. PROJECT ROOT DISCOVERY
# -----------------------------------------------------------------------------
# Resolution order:
#   1. TA_PROJECT_ROOT environment variable (explicit override)
#   2. here::here() if the 'here' package is installed
#   3. Walk up from the working directory looking for the project marker
.find_project_root <- function() {
  env <- Sys.getenv("TA_PROJECT_ROOT", unset = "")
  if (nzchar(env) && dir.exists(env)) return(normalizePath(env, winslash = "/"))

  if (requireNamespace("here", quietly = TRUE)) {
    r <- try(here::here(), silent = TRUE)
    if (!inherits(r, "try-error") && dir.exists(file.path(r, "R"))) {
      return(normalizePath(r, winslash = "/"))
    }
  }

  # Walk upward for a directory containing both R/ and scripts/
  d <- normalizePath(getwd(), winslash = "/")
  for (i in seq_len(6)) {
    if (dir.exists(file.path(d, "R")) && dir.exists(file.path(d, "scripts"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }

  stop("Could not locate the project root.\n",
       "Set the working directory to the project folder (the one containing ",
       "R/ and scripts/), or set TA_PROJECT_ROOT.", call. = FALSE)
}

PROJECT_ROOT <- .find_project_root()


# -----------------------------------------------------------------------------
# 2. STANDARD DIRECTORIES
# -----------------------------------------------------------------------------
DIR_R       <- file.path(PROJECT_ROOT, "R")
DIR_SCRIPTS <- file.path(PROJECT_ROOT, "scripts")
DIR_DATA    <- file.path(PROJECT_ROOT, "data")
DIR_OUTPUT  <- file.path(PROJECT_ROOT, "outputs")

for (d in c(DIR_DATA, DIR_OUTPUT)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# -----------------------------------------------------------------------------
# 3. INPUT DATA
# -----------------------------------------------------------------------------
DATA_FILE <- file.path(DIR_DATA, "Dengue-Rainfall_Dataset.xlsx")

# Sheet names, verified against the supplied workbook.
SHEET_QC       <- "QC Data"        # YR, WN, DC_QC, RF_NASA, RF_PAGASA, flags
SHEET_REGIONAL <- "Regional Data"  # REGION, YR, WN, DC_DOH, RF_HDX, flags
SHEET_COUNTRY  <- "Country Data"   # COUNTRY, YR, WN, DC_OPENDENGUE, RF_NASA, flags

if (!file.exists(DATA_FILE)) {
  stop("Input dataset not found:\n  ", DATA_FILE,
       "\nPlace 'Dengue-Rainfall_Dataset.xlsx' in the data/ directory.",
       call. = FALSE)
}


# -----------------------------------------------------------------------------
# 4. OPTIONAL GEOMETRY FOR THE STAGE 4 MAP
# -----------------------------------------------------------------------------
# Stage 4 draws a Philippine regional choropleth. Geometry is resolved in this
# order (see scripts/Stage4_Regional_analysis.R):
#   1. A shapefile / GeoPackage placed in data/geometry/  (fully offline)
#   2. geodata::gadm()          (downloads, cached under data/geometry/)
#   3. rnaturalearth::ne_states() (downloads)
#
# Supplying local geometry is strongly recommended for a reproducible archive:
# it removes the only remaining network dependency and pins the boundary
# vintage, which matters because GADM revises Philippine regions periodically.
DIR_GEOMETRY <- file.path(DIR_DATA, "geometry")
if (!dir.exists(DIR_GEOMETRY)) {
  dir.create(DIR_GEOMETRY, recursive = TRUE, showWarnings = FALSE)
}

# Set to an explicit file path to force a specific boundary file.
PH_SHAPEFILE <- {
  cands <- list.files(DIR_GEOMETRY, pattern = "\\.(shp|gpkg|geojson)$",
                      full.names = TRUE, recursive = TRUE)
  if (length(cands) > 0) cands[1] else NULL
}


# -----------------------------------------------------------------------------
# 5. ANALYSIS CONSTANTS SHARED ACROSS STAGES
# -----------------------------------------------------------------------------
# Years excluded from threshold derivation and evaluation:
#   2020, 2021 - COVID-19 surveillance disruption
#   2025       - truncated (incomplete) year
EXCLUDED_YEARS <- c(2020L, 2021L, 2025L)

# Reproducibility
GLOBAL_SEED <- 12345L
set.seed(GLOBAL_SEED)
options(scipen = 999, stringsAsFactors = FALSE)


# -----------------------------------------------------------------------------
# 6. ETA THRESHOLDS (Stage 1 -> Stage 3 handoff)
# -----------------------------------------------------------------------------
# Stage 1 derives eta_ON / eta_OFF empirically and writes them to
# outputs/Stage1_eta_thresholds_QC/eta_thresholds_derived.R.
#
# In the original code Stage 3 hard-coded ETA_ON <- 1.33 / ETA_OFF <- 0.73 and
# never read Stage 1's output, so re-running Stage 1 could not change Stage 3 --
# a silent reproducibility break. load_eta_thresholds() now prefers the derived
# values and falls back to the published constants with a warning.
ETA_ON_PUBLISHED  <- 1.33
ETA_OFF_PUBLISHED <- 0.73

load_eta_thresholds <- function(verbose = TRUE) {
  derived <- file.path(DIR_OUTPUT, "Stage1_eta_thresholds_QC",
                       "eta_thresholds_derived.R")

  r <- read_value_file(derived,
                       required = c("ETA_ON_ADOPTED", "ETA_OFF_ADOPTED"),
                       label = "eta_thresholds_derived.R")

  if (r$ok) {
    on  <- r$values[["ETA_ON_ADOPTED"]]
    off <- r$values[["ETA_OFF_ADOPTED"]]
    if (is.numeric(on) && is.numeric(off) && is.finite(on) && is.finite(off)) {
      if (verbose) {
        message(sprintf(
          "[eta] Using Stage 1 derived thresholds: eta_ON = %.2f, eta_OFF = %.2f",
          on, off))
      }
      return(list(ETA_ON = on, ETA_OFF = off, source = "stage1"))
    }
  }

  # The fallback previously fired SILENTLY because of the emptyenv bug above,
  # so Stage 1's derived thresholds never actually reached Stage 3. The reason
  # is now surfaced in the warning so a masked failure cannot recur unnoticed.
  if (verbose) {
    warning("Falling back to published eta thresholds ",
            sprintf("(eta_ON = %.2f, eta_OFF = %.2f).\n  Reason: ",
                    ETA_ON_PUBLISHED, ETA_OFF_PUBLISHED),
            if (nzchar(r$message)) r$message else "unknown",
            "\n  Run Stage 1 for a fully reproducible chain.", call. = FALSE)
  }
  list(ETA_ON = ETA_ON_PUBLISHED, ETA_OFF = ETA_OFF_PUBLISHED,
       source = "published")
}

# -----------------------------------------------------------------------------
# 6a. READING AUTO-GENERATED VALUE FILES
# -----------------------------------------------------------------------------
# Shared reader for the small sourceable files that Stage 1 and Stage 3b emit.
#
# ROOT CAUSE OF THE STAGE 3-5 FAILURE (fixed here):
#   These files were previously sourced into new.env(parent = emptyenv()).
#   In R an assignment such as `X <- 0.74` is a call to the function `<-`,
#   resolved lexically through the environment's parent chain. With emptyenv()
#   as the parent that chain is empty, so `<-` itself cannot be found and
#   sourcing failed with "could not find function \"<-\"". The generated files
#   were always correct; the environment they were read into was not.
#   Parent is now baseenv(), which exposes base functions while still isolating
#   the file from the caller's workspace.
read_value_file <- function(path, required, label = basename(path)) {
  if (!file.exists(path)) {
    return(list(ok = FALSE, reason = "missing", values = NULL,
                message = paste0("File not found: ", path)))
  }
  e <- new.env(parent = baseenv())
  res <- tryCatch({ sys.source(path, envir = e); NULL },
                  error = function(err) conditionMessage(err))
  if (!is.null(res)) {
    return(list(ok = FALSE, reason = "source_error", values = NULL,
                message = paste0("Could not source ", label, ": ", res)))
  }
  missing <- required[!vapply(required, exists, logical(1),
                              envir = e, inherits = FALSE)]
  if (length(missing) > 0) {
    return(list(ok = FALSE, reason = "missing_object", values = NULL,
                message = paste0(label, " does not define: ",
                                 paste(missing, collapse = ", "),
                                 "\n  It defines: ",
                                 paste(ls(e), collapse = ", "))))
  }
  vals <- mget(required, envir = e)
  list(ok = TRUE, reason = "ok", values = vals, message = "")
}

#' Validate a value is a single finite proportion strictly inside (0, 1).
is_valid_proportion <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0 && x < 1
}

# -----------------------------------------------------------------------------
# 6b. EPIDEMIC BURDEN FRACTION  (fixed at 70%)
# -----------------------------------------------------------------------------
# The A2 anchor ("Epidemic Burden" block), and hence the Activation Threshold,
# uses a FIXED burden fraction of 0.70, matching the published manuscript.
#
# This value is defined once here and read by Stages 3, 4 and 5. It is a plain
# constant: no file is read, nothing is sourced, and no stage depends on the
# output of any other stage for it. That removes the cross-stage handoff which
# was the source of the Stage 3-5 startup failures.
#
# Stage 3b (scripts/Stage3b_burden_optimization.R) remains available as an
# OPTIONAL, standalone supporting analysis that estimates this fraction
# empirically. It is descriptive only: nothing in the pipeline reads its output,
# and it can be deleted without affecting any result.
A2_BURDEN_FRAC <- 0.70
A2_BURDEN_PCT  <- sprintf("%.0f%%", 100 * A2_BURDEN_FRAC)   # "70%" for captions

# -----------------------------------------------------------------------------
# 7. DEPENDENCY CHECK
# -----------------------------------------------------------------------------
# The original scripts called install.packages() at run time. That mutates the
# user's library without consent, fails on offline / HPC nodes, and breaks
# non-interactive execution. Dependencies are now verified, never installed;
# R/install_dependencies.R performs installation as a deliberate one-off step.
require_packages <- function(pkgs, purpose = NULL) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required package(s)",
         if (!is.null(purpose)) paste0(" for ", purpose) else "", ":\n  ",
         paste(missing, collapse = ", "),
         "\n\nInstall them with:\n  source(\"R/install_dependencies.R\")",
         call. = FALSE)
  }
  invisible(TRUE)
}

message("[config] Project root: ", PROJECT_ROOT)
