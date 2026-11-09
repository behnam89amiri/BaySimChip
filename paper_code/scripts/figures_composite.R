# =============================================================================
# figures_composite.R - composite Figs. 4-7 at print size
#
# Builds the multi-panel figures from the data written by the panel scripts
# into the same run folder (no simulations here):
#   Fig. 4: A fig4a_profiles.csv, B fig4b_exposure_grid.csv, C fig456c_bias_Evap.csv
#   Fig. 5: A fig5a_profiles.csv, B fig5b_exposure_grid.csv, C fig456c_bias_Samp.csv
#   Fig. 6: A fig6a_profiles.csv, B fig6b_exposure_grid.csv, C fig456c_bias_NSB.csv,
#           D fig6d_panel_*.csv
#   Fig. 7: A fig7a_fold_errors.csv, B fig7b_alldrugs_results.csv
# Layout: RSC double-column width (17.1 cm) with an 8 pt base font, so that the
# text is printed at 6-9 pt; one row per panel group (A profiles, B heatmaps
# with colour bars below, C bias maps; Fig. 6: C and D in one row; Fig. 7:
# A radars, B scatter plots).
# Outputs: Fig4.png/.pdf ... Fig7.png/.pdf (600 dpi PNG via ragg, PDF via cairo_pdf)
# Runtime: < 1 min
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---- parameters ----------------------------------------------------------------
FIG_W <- 6.73                  # in, RSC double-column width (17.1 cm)
BASE <- 8                      # base font size (pt) at print size
COMPARTMENTS <- c("1" = "Gut apical", "3" = "Liver")
TAGS <- theme(plot.tag = element_text(size = BASE + 3, face = "bold"))

read_run <- function(name) read.csv(pc_in_run(name))

# A: profiles of the two compartments side by side (legend inside panel 1, or
# collected to the right of the row)
row_profiles <- function(prof, ylab, cols, outside = FALSE, ...) {
  pa <- lapply(names(COMPARTMENTS), function(j) {
    plot_profiles(prof, as.numeric(j), COMPARTMENTS[[j]], ylab, cols,
                  legend_inside = if (!outside && j == "1") c(0.98, 0.97) else NULL,
                  legend_outside = outside, base_size = BASE, linewidth = 0.6, ...)
  })
  if (outside) (pa[[1]] | pa[[2]]) + plot_layout(guides = "collect") else pa[[1]] | pa[[2]]
}
# B: the four exposure heatmaps in one row, colour bars below
row_heatmaps <- function(res, xcol, ycol, xlab, ylab, x_n_breaks = NULL) {
  hm <- function(z, title) plot_heatmap(res, xcol, ycol, z, title, xlab, ylab, base_size = BASE, bar_height_mm = 28,
                                        bar_horizontal = TRUE, n_breaks = 3, x_n_breaks = x_n_breaks) +
    theme(plot.margin = margin(2, 5, 2, 2))
  hm("AUC", lab("AUC")) | hm("Cmax", lab("Cmax")) | hm("T12", lab("T12")) | hm("Tmax", lab("Tmax"))
}
# C: bias maps, wide (Figs. 4, 5) or narrow with colour bars below (Fig. 6)
row_bias <- function(res) {
  b <- function(z) plot_bias_map(res, z, base_size = BASE, bar_height_mm = 22) +
    theme(legend.key.width = unit(2.5, "mm"), legend.text = element_text(size = BASE * 0.75),
          plot.margin = margin(2, 8, 2, 2))
  b("ratio_CL") | b("ratio_PS")
}
pair_bias_narrow <- function(res) {
  b <- function(z) plot_bias_map(res, z, base_size = BASE, bar_height_mm = 24, bar_horizontal = TRUE, n_breaks = 4,
                                 ylab = lab("CLint2")) + theme(plot.margin = margin(2, 8, 2, 2))
  b("ratio_CL") | b("ratio_PS")
}

# ---- Fig. 4 --------------------------------------------------------------------
prof <- read_run("fig4a_profiles.csv")
prof$cond <- factor(prof$cond, levels = c("w/o evaporation", "w/ evaporation"))
A <- row_profiles(prof, lab("conc_au"), setNames(HUE2, levels(prof$cond)))
B <- row_heatmaps(read_run("fig4b_exposure_grid.csv"), "Kev1", "Kev3", lab("kev_apical"), lab("kev_liver"), x_n_breaks = 4)
C <- row_bias(read_run("fig456c_bias_Evap.csv"))
fig4 <- wrap_elements(full = A) / wrap_elements(full = B) / wrap_elements(full = C) +
  plot_layout(heights = c(1.9, 1.9, 2.1)) + plot_annotation(tag_levels = "A") & TAGS
save_figure(fig4, "Fig4", FIG_W, 6.0)

# ---- Fig. 5 --------------------------------------------------------------------
prof <- read_run("fig5a_profiles.csv")
prof$cond <- factor(prof$cond, levels = c("0 mL", "0.05 mL", "0.10 mL", "0.15 mL"))
A <- row_profiles(prof, lab("conc_au"), setNames(HUE4, levels(prof$cond)), colour_title = "Sample volume")
B <- row_heatmaps(read_run("fig5b_exposure_grid.csv"), "N_apical", "N_liver", lab("n_apical"), lab("n_liver"))
C <- row_bias(read_run("fig456c_bias_Samp.csv"))
fig5 <- wrap_elements(full = A) / wrap_elements(full = B) / wrap_elements(full = C) +
  plot_layout(heights = c(1.9, 1.9, 2.1)) + plot_annotation(tag_levels = "A") & TAGS
save_figure(fig5, "Fig5", FIG_W, 6.0)

# ---- Fig. 6 --------------------------------------------------------------------
prof <- read_run("fig6a_profiles.csv")
prof$cond <- factor(prof$ps_nsb, levels = c(0, 0.05, 0.1, 0.2))
prof$lty  <- factor(prof$k, levels = c(1, 0.1, 3))
A <- row_profiles(prof, lab("conc_uM"), setNames(HUE4, levels(prof$cond)), outside = TRUE,
                  colour_labels = c(expression(PS[NSB]*" = 0 mL/h"), expression(PS[NSB]*" = 0.05 mL/h"),
                                    expression(PS[NSB]*" = 0.1 mL/h"), expression(PS[NSB]*" = 0.2 mL/h")),
                  linetypes = c("1" = "solid", "0.1" = "dotted", "3" = "dotdash"),
                  linetype_labels = c(expression(K[media/plastic]*" = 1"), expression(K[media/plastic]*" = 0.1"),
                                      expression(K[media/plastic]*" = 3")))
B <- row_heatmaps(read_run("fig6b_exposure_grid.csv"), "PS_NSB", "K", lab("ps_nsb"), lab("k_media"), x_n_breaks = 4)
C <- pair_bias_narrow(read_run("fig456c_bias_NSB.csv"))
D <- plot_fig6d(read_run("fig6d_panel_predictions.csv"), read_run("fig6d_panel_observed_mean.csv"),
                read_run("fig6d_panel_observed_raw.csv"), base_size = BASE)
fig6 <- wrap_elements(full = A) / wrap_elements(full = B) / (wrap_elements(full = C) | wrap_elements(full = D)) +
  plot_layout(heights = c(1.9, 1.9, 2.2)) + plot_annotation(tag_levels = "A") & TAGS
save_figure(fig6, "Fig6", FIG_W, 6.1)

# ---- Fig. 7 --------------------------------------------------------------------
fold <- read_run("fig7a_fold_errors.csv")
radar <- function(m, title, legend) {
  wide <- radar_table(fold, m)
  plot_radar(wide, RADAR_AXES, 1, radar_max(max(as.matrix(wide[, -1]))), ARCHETYPE_COLOURS,
             ARCHETYPE_LABELS, title, show_legend = legend, base_size = BASE)
}
A <- radar("CL", "Liver clearance", TRUE) | radar("PS", "Gut permeability", FALSE)
drugs <- read_run("fig7b_alldrugs_results.csv")
B <- (plot_drug_scatter(data.frame(name = drugs$Compound, x = drugs$true_CL, y = drugs$est_CL),
                        lab("CLint", "Reference"), lab("CLint", "Estimated"), base_size = BASE) |
        plot_drug_scatter(data.frame(name = drugs$Compound, x = drugs$true_P, y = drugs$est_P),
                          lab("P", "Reference"), lab("P", "Estimated"), base_size = BASE)) +
  plot_layout(guides = "collect") & theme(legend.position = "top")
fig7 <- wrap_elements(full = A) / wrap_elements(full = B) +
  plot_layout(heights = c(2.3, 3.3)) + plot_annotation(tag_levels = "A") & TAGS
save_figure(fig7, "Fig7", FIG_W, 5.8)
pc_log("done: ", pc_run_dir())
