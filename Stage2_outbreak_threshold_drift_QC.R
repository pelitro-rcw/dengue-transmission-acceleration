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
required_fig1 <- c("readxl", "dplyr", "ggplot2", "cowplot",
                   "ISOweek", "scales", "grid")

new_pkgs_fig1 <- required_fig1[!required_fig1 %in% rownames(installed.packages())]
if (length(new_pkgs_fig1) > 0) {
  install.packages(new_pkgs_fig1, dependencies = TRUE,
                   repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(ISOweek)
  library(scales)
  library(grid)
})


# -----------------------------------------------------------------------------
# 2. USER SETTINGS
# -----------------------------------------------------------------------------
# Option A:
# DATA_FILE_FIG1 <- file.path(getwd(), "Dengue-Rainfall_Dataset.xlsx")

# Option B:
DATA_FILE_FIG1 <- "C:/Users/User/Desktop/Dengue-Rainfall_Dataset.xlsx"

# Option C:
# DATA_FILE_FIG1 <- file.choose()

SHEET_FIG1  <- "QC Data"
FIGURES_DIR <- file.path(getwd(), "Stage2_outbreak_threshold_drift_QC")

if (!dir.exists(FIGURES_DIR)) {
  dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
}


# -----------------------------------------------------------------------------
# 3. LINE STYLE SETTINGS
# -----------------------------------------------------------------------------
LW_SOLID_MAIN   <- 1.45
LW_SOLID_MED    <- 1.20
LW_SOLID_LIGHT  <- 0.95
LW_DASH_MAIN    <- 1.05
LW_DASH_MED     <- 0.90
LW_DASH_LIGHT   <- 0.80
LW_VLINE        <- 0.65
LW_BORDER       <- 0.35

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
  panel_bg     = "#FAFAFA",
  panel_border = "grey70",
  pandemic     = "#FADBD8"
)


# -----------------------------------------------------------------------------
# 5. PUBLICATION THEME
# -----------------------------------------------------------------------------
theme_pub <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.background = element_rect(fill = COL$panel_bg, colour = NA),
      panel.grid.major = element_line(colour = "#CCCCCC",
                                      linetype = LT_FINE_DASH,
                                      linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = COL$panel_border, fill = NA,
                                  linewidth = LW_BORDER),
      
      axis.title = element_text(face = "bold", size = base_size - 0.5),
      axis.text  = element_text(colour = "grey30", size = base_size - 1.5),
      axis.ticks = element_line(colour = "grey60", linewidth = 0.25),
      
      plot.title    = element_text(face = "bold", size = base_size + 0.1,
                                   hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "grey40", size = base_size - 1.4,
                                   hjust = 0, margin = margin(b = 6)),
      plot.caption  = element_text(colour = "grey50", size = base_size - 2,
                                   hjust = 0, margin = margin(t = 6)),
      
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.background = element_rect(fill = "#F2F2F2", colour = "grey80",
                                       linewidth = 0.25),
      legend.key        = element_rect(fill = "#F2F2F2", colour = "#F2F2F2"),
      legend.title      = element_text(face = "bold", size = base_size - 1),
      legend.text       = element_text(size = base_size - 2, colour = "grey20"),
      legend.margin     = margin(4, 6, 4, 6),
      legend.box.margin = margin(2, 2, 2, 2),
      legend.spacing.x  = unit(0.3, "cm"),
      
      plot.margin = margin(8, 10, 6, 8)
    )
}


# -----------------------------------------------------------------------------
# 6. HELPERS
# -----------------------------------------------------------------------------
save_figure <- function(filename, plot, width, height, dpi = 600, device = NULL) {
  ggsave(
    filename  = file.path(FIGURES_DIR, filename),
    plot      = plot,
    device    = device,
    width     = width,
    height    = height,
    dpi       = dpi,
    bg        = "white",
    units     = "in",
    limitsize = FALSE
  )
  cat("  Saved:", filename, "\n")
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

make_threshold_label_df <- function(df_plot, df_thresholds) {
  annual_thresholds <- df_thresholds %>%
    group_by(YR) %>%
    summarise(
      Threshold = mean(Threshold, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Threshold = ifelse(is.infinite(Threshold), NA_real_, Threshold)
    )
  
  df_plot %>%
    group_by(YR) %>%
    summarise(
      x_pos = mean(range(WeekSeq, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    left_join(annual_thresholds, by = "YR") %>%
    filter(is.finite(Threshold)) %>%
    mutate(threshold_lab = number(Threshold, accuracy = 0.1)) %>%
    arrange(YR)
}


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
  label_df <- make_threshold_label_df(df_data, df_thresh)
  
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
    geom_line(
      aes(y = DC_QC),
      colour    = COL$cases,
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
    geom_text(
      data = label_df,
      aes(x = x_pos, y = Threshold, label = threshold_lab),
      inherit.aes = FALSE,
      size   = 2.45,
      family = "sans",
      colour = COL$threshold,
      vjust  = -0.35,
      na.rm  = TRUE
    ) +
    {
      if (shade_pandemic && is.finite(pandemic_mid)) {
        annotate(
          "text",
          x = pandemic_mid,
          y = ymax_total * 0.96,
          label = "Pandemic period\n2020-2021",
          hjust = 0.5, vjust = 1,
          size = 3.0,
          colour = COL$danger,
          fontface = "bold"
        )
      }
    } +
    scale_colour_manual(
      values = c(
        "Alarm threshold"    = COL$alarm,
        "Outbreak threshold" = COL$threshold
      ),
      breaks = c("Alarm threshold", "Outbreak threshold"),
      labels = c("Alarm threshold", "Outbreak threshold"),
      name = NULL
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = c(LW_DASH_LIGHT, LW_DASH_MED),
          linetype  = c(LT_DASH_ALARM, LT_DASH_OUT)
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
      limits = c(0, ymax_total * 1.03),
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(
      title = panel_title,
      x = "Year",
      y = "Weekly dengue cases"
    ) +
    theme_pub(base_size = 10.5) +
    theme(
      plot.title   = element_text(size = 10.5, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 7)),
      axis.title.y = element_text(margin = margin(r = 7)),
      axis.text.x  = element_text(size = 8.0),
      axis.text.y  = element_text(size = 8.0),
      legend.position      = if (show_legend) "bottom" else "none",
      legend.direction     = "horizontal",
      legend.justification = "center",
      legend.box           = "horizontal",
      plot.margin = if (show_legend) margin(8, 10, 18, 8)
      else              margin(8, 10, 6,  8)
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
  shade_pandemic = TRUE,
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
save_figure(
  filename = "fig1_dengue_thresholds_vertical.pdf",
  plot     = fig1_combined,
  width    = 11.0,
  height   = 12.0,
  dpi      = 600,
  device   = cairo_pdf
)

save_figure(
  filename = "fig1_dengue_thresholds_vertical.png",
  plot     = fig1_combined,
  width    = 11.0,
  height   = 12.0,
  dpi      = 600,
  device   = "png"
)


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
cat("\nDraft figure legend:\n")
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