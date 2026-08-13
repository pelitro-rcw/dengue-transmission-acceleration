# =============================================================================
# STAGE 2 - Figure 1: Outbreak threshold drift (Quezon City)
# -----------------------------------------------------------------------------
# Weekly notified dengue case counts in Quezon City, Philippines, 2013-2025,
# overlaid with rolling week-specific (Mean + 2 SD) outbreak thresholds. Two
# vertically stacked panels compare threshold behaviour:
#   (a) with pandemic years included
#   (b) with pandemic years (2020-2021) excluded from the rolling baseline
#
# Inputs   : Dengue-Rainfall_Dataset.xlsx, sheet "QC Data"
#            Required columns: YR (ISO year), WN (ISO week), DC_QC (cases)
#            Optional column : RF_NASA (rainfall; not used in this figure)
#
# Outputs  : Stage2_outbreak_threshold_drift_QC/
#              fig1_dengue_thresholds_vertical.pdf  (vector, 600 dpi)
#              fig1_dengue_thresholds_vertical.png  (raster, 600 dpi)
#              fig1_thresholds_with_pandemic.csv
#              fig1_thresholds_without_pandemic.csv
#
# Style    : sans-serif text, 5-7 pt body, vector export at 600 dpi,
#            final width near 180 mm (two-column journal format)
# Repro    : R >= 4.1; UTF-8 locale recommended
# =============================================================================


# -----------------------------------------------------------------------------
# 1. PACKAGES
# -----------------------------------------------------------------------------
REQUIRED_PACKAGES <- c(
  "readxl", "dplyr", "ggplot2", "cowplot",
  "ISOweek", "scales", "grid", "patchwork", "viridisLite"
)

# --- Project bootstrap -------------------------------------------------------
# Portable path resolution + shared publication theme. Replaces the previous
# inline install.packages() call and the hard-coded Desktop path.
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

require_packages(REQUIRED_PACKAGES, purpose = "Stage 2")
invisible(lapply(REQUIRED_PACKAGES, function(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))))

source(file.path(DIR_R, "01_publication_theme.R"), local = TRUE)
# Keep a handle on the shared theme before Stage 2 shadows the name.
theme_pub_shared <- theme_pub

set.seed(GLOBAL_SEED)


# -----------------------------------------------------------------------------
# 2. USER SETTINGS
# -----------------------------------------------------------------------------
DATA_FILE_FIG1 <- DATA_FILE
SHEET_FIG1     <- SHEET_QC
FIGURES_DIR    <- file.path(DIR_OUTPUT, "Stage2_outbreak_threshold_drift_QC")

if (!dir.exists(FIGURES_DIR)) {
  dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
}


# -----------------------------------------------------------------------------
# 3. LINE STYLE SETTINGS
# -----------------------------------------------------------------------------
# Line widths are absolute (mm) and do NOT scale with canvas size. The previous
# values were set for an 11 in canvas that production reduced to ~7.09 in, so
# the printed widths were ~0.64x the nominal figures. Now that the figure is
# authored at final size, the values are restated as true printed widths --
# preserving the intended appearance while staying above the ~0.25 pt stroke
# floor below which journal presses drop hairlines.
LW_SOLID_MAIN   <- 0.90
LW_SOLID_MED    <- 0.72
LW_SOLID_LIGHT  <- 0.58
LW_DASH_MAIN    <- 0.68
LW_DASH_MED     <- 0.58
LW_DASH_LIGHT   <- 0.50
LW_VLINE        <- 0.45
LW_BORDER       <- 0.30

LT_LONG_DASH    <- "33"
LT_MED_DASH     <- "42"
LT_FINE_DASH    <- "44"
LT_DOTTED_LIGHT <- "13"
LT_DASH_ALARM   <- "13"
LT_DASH_OUT     <- "42"


# -----------------------------------------------------------------------------
# 4. COLOUR PALETTE
# -----------------------------------------------------------------------------
COL <- list(
  cases        = "#1a1a1a",
  alarm        = "#E69F00",
  threshold    = "#D55E00",
  danger       = "#C0392B",
  grey         = "#7F8C8D",
  dark_grey    = "#4D4D4D",
  light_grey   = "#E6E6E6",
  panel_bg     = "#FFFFFF",  # white: grey fills lose contrast in CMYK print
  panel_border = "grey70",
  pandemic     = "#FADBD8"
)


# -----------------------------------------------------------------------------
# 5. PUBLICATION THEME
# -----------------------------------------------------------------------------
# Stage 2 now inherits the shared theme from R/01_publication_theme.R and only
# overrides the legend placement it needs (horizontal, below the panels).
theme_pub <- function(base_size = PUB_BASE) {
  theme_pub_shared(base_size) %+replace%
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.spacing.x = grid::unit(0.25, "cm"),
      legend.box.margin = ggplot2::margin(2, 2, 2, 2)
    )
}


# -----------------------------------------------------------------------------
# 6. HELPERS
# -----------------------------------------------------------------------------
# Delegates to save_pub() so PDF + PNG export, dimension clamping and device
# selection are handled identically in every stage.
save_figure <- function(stem, plot, width, height, dpi = 600) {
  save_pub(stem = sub("\\.(pdf|png)$", "", stem), plot = plot,
           width = width, height = height, dir = FIGURES_DIR, dpi = dpi)
}

compute_week_specific_thresholds <- function(df_input, rolling_map) {
  out <- lapply(names(rolling_map), function(y) {
    yrs <- rolling_map[[y]]
    
    weeks_y <- df_input %>%
      filter(YR == as.integer(y)) %>%
      pull(WN) %>%
      unique() %>%
      sort()
    
    if (length(weeks_y) == 0) weeks_y <- 1:53
    
    per_week <- lapply(weeks_y, function(w) {
      vals <- df_input %>%
        filter(YR %in% yrs, WN == w) %>%
        pull(DC_QC)
      
      vals <- vals[is.finite(vals)]
      n_vals <- length(vals)
      
      if (n_vals == 0) {
        mean_val <- NA_real_; sd_val <- NA_real_
        alarm_th <- NA_real_; out_th <- NA_real_
      } else if (n_vals == 1) {
        mean_val <- mean(vals); sd_val <- NA_real_
        alarm_th <- mean_val;   out_th <- mean_val
      } else {
        mean_val <- mean(vals)
        sd_val   <- stats::sd(vals)
        alarm_th <- mean_val + 1 * sd_val
        out_th   <- mean_val + 2 * sd_val
      }
      
      data.frame(
        YR = as.integer(y),
        WN = as.integer(w),
        Mean = mean_val,
        SD = sd_val,
        AlarmThreshold = alarm_th,
        Threshold = out_th,
        PriorYears = paste(yrs, collapse = ", "),
        stringsAsFactors = FALSE
      )
    })
    
    bind_rows(per_week)
  })
  
  bind_rows(out) %>% arrange(YR, WN)
}

make_year_axis_df <- function(df_input) {
  df_input %>%
    group_by(YR) %>%
    summarise(
      x_pos = mean(WeekSeq, na.rm = TRUE),
      xmin  = min(WeekSeq,  na.rm = TRUE),
      xmax  = max(WeekSeq,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(YR)
}

# make_threshold_label_df() removed: the per-year numeric threshold values it
# printed between the alarm and outbreak lines cluttered the panels. The values
# are unchanged and remain available in df_thresh and the exported CSV tables.


# -----------------------------------------------------------------------------
# 7. LOAD DATA
# -----------------------------------------------------------------------------
if (!file.exists(DATA_FILE_FIG1)) {
  stop(paste0(
    "Data file not found at:\n", DATA_FILE_FIG1, "\n\n",
    "Please edit DATA_FILE_FIG1 or place the file in the working directory."
  ))
}

df_raw_fig1 <- readxl::read_excel(DATA_FILE_FIG1, sheet = SHEET_FIG1) %>%
  as.data.frame()

required_cols_fig1 <- c("YR", "WN", "DC_QC")
missing_cols_fig1  <- setdiff(required_cols_fig1, names(df_raw_fig1))
if (length(missing_cols_fig1) > 0) {
  stop("Missing required columns in '", SHEET_FIG1, "': ",
       paste(missing_cols_fig1, collapse = ", "))
}

df_fig1 <- df_raw_fig1 %>%
  mutate(
    YR    = suppressWarnings(as.integer(YR)),
    WN    = suppressWarnings(as.integer(WN)),
    DC_QC = suppressWarnings(as.numeric(DC_QC)),
    RF_NASA = if ("RF_NASA" %in% names(.))
      suppressWarnings(as.numeric(RF_NASA)) else NA_real_
  )

if (any(is.na(df_fig1$YR))) stop("Column 'YR' contains invalid or missing values.")
if (any(is.na(df_fig1$WN))) stop("Column 'WN' contains invalid or missing values.")
if (all(is.na(df_fig1$DC_QC))) stop("Column 'DC_QC' contains only NA values.")

df_fig1 <- df_fig1 %>%
  filter(WN >= 1, WN <= 53) %>%
  mutate(
    ISOweek = sprintf("%d-W%02d", YR, WN),
    Date    = ISOweek::ISOweek2date(paste0(ISOweek, "-1"))
  ) %>%
  arrange(Date, YR, WN)

if (any(is.na(df_fig1$Date))) {
  stop("Some ISO week/year combinations could not be converted to dates.")
}


# -----------------------------------------------------------------------------
# 8. ANALYSIS SETTINGS
# -----------------------------------------------------------------------------
pandemic_years <- c(2020, 2021)
target_years   <- 2013:2025

rolling_map_with_pandemic <- list(
  "2013" = 2008:2012, "2014" = 2009:2013, "2015" = 2010:2014,
  "2016" = 2011:2015, "2017" = 2012:2016, "2018" = 2013:2017,
  "2019" = 2014:2018, "2020" = 2015:2019, "2021" = 2016:2020,
  "2022" = 2017:2021, "2023" = 2018:2022, "2024" = 2019:2023,
  "2025" = 2020:2024
)

rolling_map_no_pandemic <- list(
  "2013" = 2008:2012,
  "2014" = 2009:2013,
  "2015" = 2010:2014,
  "2016" = 2011:2015,
  "2017" = 2012:2016,
  "2018" = 2013:2017,
  "2019" = 2014:2018,
  "2022" = 2015:2019,
  "2023" = c(2016:2019, 2022),
  "2024" = c(2017:2019, 2022:2023),
  "2025" = c(2018:2019, 2022:2024)
)


# -----------------------------------------------------------------------------
# 9. PREPARE PLOTTING DATA
# -----------------------------------------------------------------------------
df_plot_full <- df_fig1 %>%
  filter(YR %in% target_years) %>%
  arrange(Date, WN) %>%
  mutate(WeekSeq = row_number())

df_thresh_full_year <- compute_week_specific_thresholds(
  df_input    = df_fig1,
  rolling_map = rolling_map_with_pandemic
)

df_thresh_full <- df_plot_full %>%
  dplyr::select(WeekSeq, Date, YR, WN) %>%
  left_join(df_thresh_full_year, by = c("YR", "WN")) %>%
  arrange(WeekSeq)

df_plot_no_pandemic <- df_fig1 %>%
  filter(YR %in% target_years, !YR %in% pandemic_years) %>%
  arrange(Date, WN) %>%
  mutate(WeekSeq = row_number())

df_thresh_no_pandemic_year <- compute_week_specific_thresholds(
  df_input    = df_fig1 %>% filter(!YR %in% pandemic_years),
  rolling_map = rolling_map_no_pandemic
)

df_thresh_no_pandemic <- df_plot_no_pandemic %>%
  dplyr::select(WeekSeq, Date, YR, WN) %>%
  dplyr::left_join(df_thresh_no_pandemic_year, by = c("YR", "WN")) %>%
  dplyr::arrange(WeekSeq)


# -----------------------------------------------------------------------------
# 10. PANEL BUILDER
# -----------------------------------------------------------------------------
build_fig1_panel <- function(df_data,
                             df_thresh,
                             panel_title,
                             shade_pandemic = FALSE,
                             show_legend = FALSE) {
  axis_df  <- make_year_axis_df(df_data)
  
  ymax_cases  <- suppressWarnings(max(df_data$DC_QC,            na.rm = TRUE))
  ymax_alarm  <- suppressWarnings(max(df_thresh$AlarmThreshold, na.rm = TRUE))
  ymax_thresh <- suppressWarnings(max(df_thresh$Threshold,      na.rm = TRUE))
  ymax_total  <- max(c(ymax_cases, ymax_alarm, ymax_thresh), na.rm = TRUE)
  
  if (!is.finite(ymax_total) || ymax_total <= 0) ymax_total <- 1
  
  pandemic_mid <- if (all(c(2020, 2021) %in% axis_df$YR)) {
    mean(c(
      axis_df$xmin[axis_df$YR == 2020],
      axis_df$xmax[axis_df$YR == 2021]
    ))
  } else {
    NA_real_
  }
  
  p <- ggplot(df_data, aes(x = WeekSeq)) +
    {
      if (shade_pandemic && all(c(2020, 2021) %in% axis_df$YR)) {
        geom_rect(
          data = data.frame(
            xmin = axis_df$xmin[axis_df$YR == 2020] - 0.5,
            xmax = axis_df$xmax[axis_df$YR == 2021] + 0.5
          ),
          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
          inherit.aes = FALSE,
          fill  = COL$pandemic,
          alpha = 0.35
        )
      }
    } +
    # Mapped to `colour` (not set as a constant) so the case series appears in
    # the legend alongside the two thresholds.
    geom_line(
      aes(y = DC_QC, colour = "Dengue cases"),
      linewidth = LW_SOLID_MED,
      lineend   = "round",
      linejoin  = "round",
      na.rm = TRUE
    ) +
    geom_line(
      data = df_thresh,
      aes(y = AlarmThreshold, colour = "Alarm threshold"),
      linewidth = LW_DASH_LIGHT,
      linetype  = LT_DASH_ALARM,
      lineend   = "round",
      linejoin  = "round",
      na.rm = TRUE
    ) +
    geom_line(
      data = df_thresh,
      aes(y = Threshold, colour = "Outbreak threshold"),
      linewidth = LW_DASH_MED,
      linetype  = LT_DASH_OUT,
      lineend   = "round",
      linejoin  = "round",
      na.rm = TRUE
    ) +
    {
      if (shade_pandemic && is.finite(pandemic_mid)) {
        annotate(
          "text",
          x = pandemic_mid,
          y = ymax_total * 1.09,
          label = "Pandemic period\n2020-2021",
          hjust = 0.5, vjust = 1,
          size = pub_text_size(PUB_ANNOT),
          colour = COL$danger,
          fontface = "bold"
        )
      }
    } +
    scale_colour_manual(
      values = c(
        "Dengue cases"       = COL$cases,
        "Alarm threshold"    = COL$alarm,
        "Outbreak threshold" = COL$threshold
      ),
      # limits + drop = FALSE guarantee all three keys exist, so the
      # length-3 override.aes vectors below always match the key count. Without
      # this, a panel where a threshold series is entirely NA would produce
      # fewer keys and fail with "replacement has 3 rows, data has N".
      limits = c("Dengue cases", "Alarm threshold", "Outbreak threshold"),
      breaks = c("Dengue cases", "Alarm threshold", "Outbreak threshold"),
      labels = c("Dengue cases", "Alarm threshold", "Outbreak threshold"),
      drop   = FALSE,
      name = NULL
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          # Three entries now: solid case line first, then the two dashed
          # thresholds. The vectors must match the `breaks` length or the
          # override is silently recycled.
          linewidth = c(LW_SOLID_MED, LW_DASH_LIGHT, LW_DASH_MED),
          linetype  = c("solid", LT_DASH_ALARM, LT_DASH_OUT)
        )
      )
    ) +
    scale_x_continuous(
      breaks = axis_df$x_pos,
      labels = axis_df$YR,
      expand = expansion(mult = c(0.004, 0.004))
    ) +
    scale_y_continuous(
      breaks = pretty(c(0, ymax_total), n = 5),
      # Headroom raised from 1.03 to 1.12 so the pandemic annotation sits clear
      # of the case trace instead of overprinting its peak.
      limits = c(0, ymax_total * 1.12),
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(
      title = panel_title,
      x = "Year",
      y = "Weekly dengue cases"
    ) +
    theme_pub() +
    theme_bold_axes() +
    theme(
      plot.title   = element_text(size = PUB_TITLE, face = "bold",
                                  colour = "black"),
      axis.title.x = element_text(margin = margin(t = 7)),
      axis.title.y = element_text(margin = margin(r = 7)),
      axis.text.x  = element_text(size = PUB_AXIS_TXT),
      axis.text.y  = element_text(size = PUB_AXIS_TXT),
      legend.position      = if (show_legend) "bottom" else "none",
      legend.direction     = "horizontal",
      legend.justification = "center",
      legend.box           = "horizontal",
      # Alarm / Outbreak threshold entries set two points above the shared
      # legend size, with wider keys and spacing so the two
      # dashed patterns stay distinguishable at print size.
      legend.text          = element_text(size = PUB_LEG_TXT + 2.0,
                                          margin = margin(l = 3, r = 9)),
      legend.key.width     = grid::unit(22, "pt"),
      legend.key.height    = grid::unit(11, "pt"),
      legend.spacing.x     = grid::unit(9, "pt"),
      legend.box.spacing   = grid::unit(7, "pt"),
      legend.margin        = margin(4, 4, 2, 4),
      # Extra bottom margin when the legend is shown keeps it clear of the
      # x-axis title; extra top margin keeps the panel title clear of the data.
      plot.margin = if (show_legend) margin(9, 12, 14, 9)
      else              margin(9, 12, 7,  9)
    )
  
  p
}


# -----------------------------------------------------------------------------
# 11. BUILD PANELS
# -----------------------------------------------------------------------------
panel_a <- build_fig1_panel(
  df_data        = df_plot_full,
  df_thresh      = df_thresh_full,
  panel_title    = "a    With pandemic years",
  shade_pandemic = FALSE,   # pandemic band removed
  show_legend    = FALSE
)

panel_b <- build_fig1_panel(
  df_data        = df_plot_no_pandemic,
  df_thresh      = df_thresh_no_pandemic,
  panel_title    = "b    Without pandemic years",
  shade_pandemic = FALSE,
  show_legend    = TRUE
)


# -----------------------------------------------------------------------------
# 12. COMBINE FIGURE
# -----------------------------------------------------------------------------
# PANEL LETTERS are inline in the panel titles ("a    With pandemic years"),
# matching the house style used by Figures 2, 4, 5 and 6: lowercase letter, no
# trailing period, four spaces, then the title -- all on one line. cowplot's
# `labels` argument is deliberately NOT used, since it would stack a second
# letter above the title.
fig1_combined <- cowplot::plot_grid(
  panel_a,
  panel_b,
  ncol  = 1,
  align = "v",
  axis  = "lr",
  rel_heights = c(1, 1.08)
)

print(fig1_combined)


# -----------------------------------------------------------------------------
# 13. SAVE FIGURE
# -----------------------------------------------------------------------------
# Stage 2 is a tall stacked panel: full double-column width, full page height.
# Previously 11 x 12 in, which production reduced by ~0.56x.
save_figure("fig1_dengue_thresholds_vertical", fig1_combined,
            width = NC_W_DOUBLE, height = NC_H_MAX)


# -----------------------------------------------------------------------------
# 14. SAVE THRESHOLD TABLES
# -----------------------------------------------------------------------------
utils::write.csv(
  df_thresh_full_year,
  file.path(FIGURES_DIR, "fig1_thresholds_with_pandemic.csv"),
  row.names = FALSE
)

utils::write.csv(
  df_thresh_no_pandemic_year,
  file.path(FIGURES_DIR, "fig1_thresholds_without_pandemic.csv"),
  row.names = FALSE
)

cat("\nThreshold tables saved.\n")


# -----------------------------------------------------------------------------
# 15. PRINT TABULAR OUTPUTS
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("FIGURE 1 THRESHOLD OUTPUTS\n")
cat("============================================================\n")

cat("\nWith Pandemic Years Included:\n")
print(as.data.frame(df_thresh_full_year), row.names = FALSE, na.print = "NA")

cat("\nWithout Pandemic Years (Excluding 2020-2021):\n")
print(as.data.frame(df_thresh_no_pandemic_year), row.names = FALSE, na.print = "NA")

cat("\nRounded Threshold Table -- With Pandemic Years:\n")
print(
  as.data.frame(
    df_thresh_full_year %>%
      mutate(
        Mean           = round(Mean,           2),
        SD             = round(SD,             2),
        AlarmThreshold = round(AlarmThreshold, 2),
        Threshold      = round(Threshold,      2)
      )
  ),
  row.names = FALSE,
  na.print = "NA"
)

cat("\nRounded Threshold Table -- Without Pandemic Years:\n")
print(
  as.data.frame(
    df_thresh_no_pandemic_year %>%
      mutate(
        Mean           = round(Mean,           2),
        SD             = round(SD,             2),
        AlarmThreshold = round(AlarmThreshold, 2),
        Threshold      = round(Threshold,      2)
      )
  ),
  row.names = FALSE,
  na.print = "NA"
)


# -----------------------------------------------------------------------------
# 16. DRAFT FIGURE LEGEND
# -----------------------------------------------------------------------------
cat("\nFigure legend:\n")
cat(
  paste0(
    "Fig. 1 | Weekly notified dengue cases and shifting alarm and outbreak ",
    "thresholds in Quezon City, Philippines, 2013-2025. (a) Weekly dengue ",
    "cases with rolling week-specific thresholds (Mean + 1 SD = alarm ",
    "threshold; Mean + 2 SD = outbreak threshold) computed from prior-year ",
    "windows that include the pandemic years 2020 and 2021. The shaded band ",
    "and annotation indicate the pandemic period. (b) The same series with ",
    "rolling thresholds computed from prior-year windows that exclude 2020 ",
    "and 2021, illustrating the threshold drift attributable to pandemic-",
    "era reporting disruption. Threshold labels show the annual mean ",
    "outbreak-threshold value assigned to each year.\n"
  )
)


# =============================================================================
# END OF STAGE 2 - Figure 1: Outbreak threshold drift (Quezon City)
# =============================================================================