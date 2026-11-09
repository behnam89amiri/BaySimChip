# =============================================================================
# R/plotting.R - figure helpers shared by the figure scripts and
# scripts/figures_composite.R
#
# Labels follow the notation of the paper:
# CL_int, P (10^-6 cm/s), PS_NSB, K_media/plastic, K_ev,apical / K_ev,liver,
# µM. Panel PNGs are written with ragg; composite figures as 600 dpi PNG
# (ragg) and PDF (cairo_pdf).
# =============================================================================

# ---- labels ------------------------------------------------------------------
LAB <- list(
  time        = "Time (h)",
  time_phase  = "Time within phase (h)",
  conc_uM     = "Concentration (µM)",
  conc_au     = "Concentration (a.u.)",
  kev_apical  = quote(K[ev*","*apical]~"(mL/h)"),
  kev_liver   = quote(K[ev*","*liver]~"(mL/h)"),
  ps_nsb      = quote(PS[NSB]~"(mL/h)"),
  k_media     = quote(K[media/plastic]),
  P           = quote(P~"("*10^-6~"cm/s)"),
  CLint       = quote(CL[int]~"(µL/min/"*10^6~"hepatocytes)"),
  CLint2      = quote(atop(CL[int], "(µL/min/"*10^6~"hepatocytes)")),   # two lines (narrow panels)
  ratio_CL    = quote(CL[int*","*est]/CL[int*","*true]),
  ratio_P     = quote(P[est]/P[true]),
  AUC         = "AUC",
  Cmax        = quote(C[max]),
  Tmax        = quote(T[max]),
  T12         = quote(T[1/2]),
  n_apical    = "# samples gut apical",
  n_liver     = "# samples liver"
)
lab <- function(name, prefix = NULL) {
  x <- LAB[[name]]
  if (is.null(prefix)) return(if (is.character(x)) x else as.expression(x))
  if (is.character(x)) return(paste(prefix, x))
  as.expression(bquote(.(prefix)~.(x)))
}

# ggplot2 default hue colours
HUE2 <- scales::hue_pal()(2)
HUE4 <- scales::hue_pal()(4)

theme_paper <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey92", colour = "grey30"))
}

# ---- saving ------------------------------------------------------------------

# Single panel: PNG (ragg)
save_panel <- function(plot, name, width, height, dpi = 300) {
  ggplot2::ggsave(pc_out(paste0(name, ".png")), plot, width = width, height = height,
                  units = "in", dpi = dpi, device = ragg::agg_png, bg = "white")
}

# Composite figure: 600 dpi PNG (ragg) + PDF (cairo_pdf)
save_figure <- function(plot, name, width, height) {
  ggplot2::ggsave(pc_out(paste0(name, ".png")), plot, width = width, height = height,
                  units = "in", dpi = 600, device = ragg::agg_png, bg = "white")
  ggplot2::ggsave(pc_out(paste0(name, ".pdf")), plot, width = width, height = height,
                  units = "in", device = grDevices::cairo_pdf, bg = "white")
}

# ---- concentration-time profiles (Figs. 4A, 5A, 6A) --------------------------

# df: time_h, compartment, conc, cond (factor; colour); optional `lty`
# (factor; line type)
plot_profiles <- function(df, compartment, title, ylab, colours,
                          colour_labels = ggplot2::waiver(), colour_title = NULL,
                          linetypes = NULL, linetype_labels = ggplot2::waiver(), linetype_title = NULL,
                          legend_inside = NULL, base_size = 11, linewidth = 0.8, legend_outside = FALSE) {
  d <- df[df$compartment == compartment, ]
  if (is.null(linetypes)) {
    p <- ggplot2::ggplot(d, ggplot2::aes(time_h, conc, colour = cond)) +
      ggplot2::geom_line(linewidth = linewidth)
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(time_h, conc, colour = cond, linetype = lty, group = run)) +
      ggplot2::geom_line(linewidth = linewidth) +
      ggplot2::scale_linetype_manual(values = linetypes, labels = linetype_labels, name = linetype_title)
  }
  p <- p +
    ggplot2::scale_colour_manual(values = colours, labels = colour_labels, name = colour_title) +
    ggplot2::labs(title = title, x = lab("time"), y = ylab) +
    theme_paper(base_size)
  if (legend_outside) {                      # legend right of the panel (collected by patchwork)
    p + ggplot2::theme(legend.position = "right",
                       legend.key.height = grid::unit(0.9, "lines"),
                       legend.spacing.y = grid::unit(0.1, "lines"),
                       legend.title = ggplot2::element_text(size = base_size * 0.85),
                       legend.text = ggplot2::element_text(size = base_size * 0.8))
  } else if (is.null(legend_inside)) {
    p + ggplot2::theme(legend.position = "none")
  } else {
    p + ggplot2::theme(legend.position = "inside", legend.position.inside = legend_inside,
                       legend.justification = c(1, 1),
                       legend.background = ggplot2::element_blank(),
                       legend.key = ggplot2::element_blank(),
                       legend.key.height = grid::unit(0.9, "lines"),
                       legend.spacing.y = grid::unit(0.1, "lines"),
                       legend.title = ggplot2::element_text(size = base_size * 0.85),
                       legend.text = ggplot2::element_text(size = base_size * 0.8))
  }
}

# ---- exposure heatmaps (Figs. 4B, 5B, 6B) ------------------------------------

# Colour limits = data range, "pretty" breaks; geom_raster(interpolate = TRUE).
# Labels: scales::label_number().
# bar_horizontal = TRUE puts the colour bar below the panel (print-size composites);
# n_breaks: colour-bar ticks; x_n_breaks: number of x-axis ticks (NULL = default)
plot_heatmap <- function(df, xcol, ycol, zcol, title, xlab, ylab, base_size = 11, bar_height_mm = 30,
                         bar_horizontal = FALSE, n_breaks = 6, x_n_breaks = NULL) {
  z <- df[[zcol]]
  lims <- range(z, na.rm = TRUE)
  brks <- pretty(lims, n = n_breaks)
  ggplot2::ggplot(df, ggplot2::aes(.data[[xcol]], .data[[ycol]], fill = .data[[zcol]])) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::scale_fill_viridis_c(option = "viridis", limits = lims, breaks = brks,
                                  labels = scales::label_number(), oob = scales::squish) +
    colour_guide(bar_horizontal, bar_height_mm) +
    (if (is.null(x_n_breaks)) NULL else ggplot2::scale_x_continuous(n.breaks = x_n_breaks)) +
    ggplot2::labs(title = title, x = xlab, y = ylab, fill = NULL) +
    theme_paper(base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    colour_bar_theme(bar_horizontal, base_size)
}

# colour bar to the right of the panel (default) or below it (horizontal)
colour_guide <- function(horizontal, bar_mm) {
  if (horizontal) {
    ggplot2::guides(fill = ggplot2::guide_colorbar(direction = "horizontal", barwidth = grid::unit(bar_mm, "mm"),
                                                   barheight = grid::unit(2.2, "mm")))
  } else {
    ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = grid::unit(bar_mm, "mm"),
                                                   barwidth = grid::unit(4, "mm")))
  }
}
colour_bar_theme <- function(horizontal, base_size) {
  if (!horizontal) return(NULL)
  ggplot2::theme(legend.position = "bottom", legend.margin = ggplot2::margin(0, 0, 0, 0),
                 legend.box.margin = ggplot2::margin(-4, 0, 0, 0),
                 legend.text = ggplot2::element_text(size = base_size * 0.75))
}

# ---- bias maps (Figs. 4C, 5C, 6C) --------------------------------------------

# res: P_1e6_cm_s, CL_uL_min_1e6 and ratio_CL / ratio_PS (log-spaced grid)
plot_bias_map <- function(res, zcol, base_size = 11, bar_height_mm = 30, bar_horizontal = FALSE, n_breaks = 6,
                          ylab = lab("CLint")) {
  title <- if (zcol == "ratio_CL") lab("ratio_CL") else lab("ratio_P")
  z <- res[[zcol]]
  lims <- range(z, na.rm = TRUE)
  ggplot2::ggplot(res, ggplot2::aes(P_1e6_cm_s, CL_uL_min_1e6, fill = .data[[zcol]])) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_viridis_c(option = "viridis", limits = lims, breaks = pretty(lims, n = n_breaks),
                                  labels = scales::label_number(), oob = scales::squish) +
    colour_guide(bar_horizontal, bar_height_mm) +
    ggplot2::labs(title = title, x = lab("P"), y = ylab, fill = NULL) +
    theme_paper(base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    colour_bar_theme(bar_horizontal, base_size)
}

# ---- Fig. 3 ------------------------------------------------------------------

# Arm codes for display: "..._100uM" -> "..._100µM"
arm_label <- function(x) sub("uM$", "µM", x)

# A: measured (points) and fitted (lines) concentrations for selected arms
plot_fig3_fits <- function(obs, sim, arms, base_size = 10) {
  obs <- obs[obs$arm %in% arms, ]
  sim <- sim[sim$arm %in% arms, ]
  obs$arm <- factor(arm_label(obs$arm), levels = arm_label(arms))
  sim$arm <- factor(arm_label(sim$arm), levels = arm_label(arms))
  ggplot2::ggplot() +
    ggplot2::geom_line(data = sim, ggplot2::aes(time_h, conc), colour = HUE2[2], linewidth = 0.8) +
    ggplot2::geom_point(data = obs, ggplot2::aes(time_h, conc), colour = HUE2[1], size = 1.4) +
    ggplot2::facet_wrap(~arm, nrow = 1, scales = "free") +
    ggplot2::labs(x = lab("time"), y = lab("conc_uM")) +
    theme_paper(base_size)
}

# B: arm-wise estimates with +/- RSE whiskers and literature comparators.
# est: long table SheetName, Parameter, Value, RSE, is_ref
plot_fig3_estimates <- function(est, parameter, order, xlab, base_size = 10) {
  d <- est[est$Parameter == parameter, ]
  d$SheetName <- factor(arm_label(d$SheetName), levels = arm_label(order))
  ggplot2::ggplot(d, ggplot2::aes(x = SheetName, y = Value)) +
    ggplot2::geom_point(ggplot2::aes(colour = is_ref, shape = is_ref), size = 2.5) +
    ggplot2::geom_errorbar(data = d[!is.na(d$RSE), ],
                           ggplot2::aes(ymin = Value - Value * RSE / 100, ymax = Value + Value * RSE / 100),
                           width = 0.2) +
    ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "#1f77b4", `TRUE` = "gray40"),
                                 labels = c("Estimate", "Literature"), drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                                labels = c("Estimate", "Literature"), drop = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = xlab, colour = "Data type", shape = "Data type") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = base_size * 0.75),
                   legend.position = "top")
}

# C: systematic comparison boxplots. comp: long table
# Parameter, ComparisonType, Group, Value. The points are jittered with a fixed
# seed (position_jitter(seed = 1)), so the figure is reproducible.
# free_x = TRUE: each facet shows only its own two groups (print-size composite)
plot_fig3_systematic <- function(comp, parameter, panel_label, ylab, show_x, base_size = 10, free_x = FALSE) {
  pal6 <- c("Daily" = "#4C78A8", "Single" = "#F58518", "Dynamic" = "#54A24B",
            "Static" = "#B279A2", "Multi-Organ" = "#72B7B2", "Single-Organ" = "#E45756")
  d <- comp[comp$Parameter == parameter, ]
  d$PanelLabel <- panel_label
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Group, y = Value, fill = Group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.9) +
    ggplot2::geom_jitter(ggplot2::aes(fill = Group), position = ggplot2::position_jitter(width = 0.12, seed = 1),
                         size = if (free_x) 1.3 else 1.7, alpha = 0.85,
                         shape = 21, colour = "black", stroke = 0.4, show.legend = FALSE) +
    ggplot2::facet_grid(PanelLabel ~ ComparisonType, scales = if (free_x) "free" else "free_y") +
    ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
    ggplot2::scale_fill_manual(values = pal6, drop = FALSE) +
    ggplot2::labs(title = NULL, x = NULL, y = ylab) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = "none",
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
                   panel.grid.minor = ggplot2::element_blank())
  if (show_x && free_x) {
    p
  } else if (show_x) {
    p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  } else {
    p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  }
}

# ---- Fig. 6D -------------------------------------------------------------------

# pred: phase, time_fit, c1, included_in_fit; obs_mean / obs_raw: phase,
# time_fit, conc, included_in_fit. Colours: blue = phases used in the fit,
# orange = phase not used.
plot_fig6d <- function(pred, obs_mean, obs_raw, base_size = 10, legend_inside = TRUE) {
  lv <- c("TRUE", "FALSE")
  lbl <- c("TRUE" = "Fitted", "FALSE" = "Not fitted")
  pred$fit <- factor(as.character(pred$included_in_fit), levels = lv)
  obs_mean$fit <- factor(as.character(obs_mean$included_in_fit), levels = lv)
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = obs_raw, ggplot2::aes(time_fit, conc), colour = "grey70", size = 1.6) +
    ggplot2::geom_point(data = obs_mean, ggplot2::aes(time_fit, conc, shape = fit), size = 2) +
    ggplot2::geom_line(data = pred, ggplot2::aes(time_fit, c1, colour = fit), linewidth = 1) +
    ggplot2::facet_wrap(~phase, scales = "free_x", nrow = 1) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"), labels = lbl, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), labels = lbl, drop = FALSE) +
    ggplot2::labs(x = lab("time_phase"), y = lab("conc_uM"), colour = NULL, shape = NULL) +
    theme_paper(base_size)
  if (legend_inside) {
    p + ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.99, 0.95),
                       legend.justification = c(1, 1),
                       legend.background = ggplot2::element_blank(), legend.key = ggplot2::element_blank())
  } else {
    p
  }
}

# ---- Fig. 7A: radar chart --------------------------------------------------------

RADAR_AXES <- c("Evaporation", "Sampling", "NSB")
ARCHETYPE_COLOURS <- stats::setNames(HUE4, c("highCL-highP", "highCL-lowP", "lowCL-highP", "lowCL-lowP"))
ARCHETYPE_LABELS <- c(expression("high "*CL[int]*" – high "*P), expression("high "*CL[int]*" – low "*P),
                      expression("low "*CL[int]*" – high "*P), expression("low "*CL[int]*" – low "*P))

# Radial maximum slightly above the largest fold error (at least 1.1)
radar_max <- function(max_fold) {
  step <- if (!is.finite(max_fold) || max_fold <= 2) 0.1 else if (max_fold <= 5) 0.25 else if (max_fold <= 10) 0.5 else 1
  max(ceiling(max_fold / step) * step, 1.1)
}

# Wide table (name x axis) of fold errors for one metric ("CL" or "PS")
radar_table <- function(fold, metric) {
  wide <- tidyr::pivot_wider(fold[, c("name", "axis", paste0("fold_", metric))],
                             names_from = axis, values_from = paste0("fold_", metric))
  wide <- as.data.frame(wide[order(wide$name), ])
  wide[is.na(wide)] <- 1
  wide
}

# Implementation adapted from fmsb::radarchart(axistype = 1, seg = 4, centerzero =
# FALSE) in ggplot2 (identical geometry: value v -> radius
# 1/(seg+1) + (v - min)/(max - min) * seg/(seg+1); first axis at the top,
# counter-clockwise), so that the radar can be combined with patchwork and
# exported as vector graphics. wide: data frame with column `name` and one
# column per axis (fold errors).
plot_radar <- function(wide, axes, min_val, max_val, colours, series_labels,
                       title, seg = 4, show_legend = FALSE, base_size = 10, fill_alpha = 0.2) {
  n <- length(axes)
  theta <- (seq(90, 450, length.out = n + 1) * pi / 180)[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  ring <- do.call(rbind, lapply(0:seg, function(i) {
    r <- (i + 1) / (seg + 1)
    data.frame(ring = i, x = c(xx, xx[1]) * r, y = c(yy, yy[1]) * r)
  }))
  spokes <- data.frame(x = xx / (seg + 1), y = yy / (seg + 1), xend = xx, yend = yy)
  ax_lab <- data.frame(x = -0.05, y = (0:seg + 1) / (seg + 1),
                       label = round(seq(min_val, max_val, length.out = seg + 1), 2))
  v_lab <- data.frame(x = xx * 1.2, y = yy * 1.2, label = axes)
  series <- do.call(rbind, lapply(seq_len(nrow(wide)), function(k) {
    v <- as.numeric(wide[k, axes])
    r <- 1 / (seg + 1) + (v - min_val) / (max_val - min_val) * seg / (seg + 1)
    data.frame(name = wide$name[k], x = xx * r, y = yy * r, ord = seq_len(n))
  }))
  series$name <- factor(series$name, levels = names(colours))
  ggplot2::ggplot(data.frame(panel = title)) +
    ggplot2::geom_path(data = ring, ggplot2::aes(x, y, group = ring), colour = "grey70", linewidth = 0.35) +
    ggplot2::geom_segment(data = spokes, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          colour = "grey70", linewidth = 0.35) +
    ggplot2::geom_polygon(data = series, ggplot2::aes(x, y, group = name, colour = name, fill = name),
                          linewidth = 0.8, alpha = fill_alpha) +
    ggplot2::geom_point(data = series, ggplot2::aes(x, y, colour = name), size = 1.6) +
    ggplot2::geom_text(data = ax_lab, ggplot2::aes(x, y, label = label), colour = "blue",
                       size = base_size * 0.3, hjust = 0.5) +
    ggplot2::geom_text(data = v_lab, ggplot2::aes(x, y, label = label), size = base_size * 0.36) +
    ggplot2::scale_colour_manual(values = colours, labels = series_labels, name = NULL) +
    ggplot2::scale_fill_manual(values = colours, labels = series_labels, name = NULL) +
    ggplot2::coord_equal(xlim = c(-1.35, 1.35), ylim = c(-0.85, 1.3), clip = "off") +
    ggplot2::facet_wrap(~panel) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey92", colour = "grey70"),
                   strip.text = ggplot2::element_text(face = "bold", margin = ggplot2::margin(3, 3, 3, 3)),
                   panel.border = ggplot2::element_rect(fill = NA, colour = "grey70"),
                   legend.position = if (show_legend) "inside" else "none",
                   legend.position.inside = c(0.02, 0.98),
                   legend.justification = c(0, 1),
                   legend.text = ggplot2::element_text(size = base_size * 0.75),
                   legend.key.size = grid::unit(0.8, "lines"))
}

# ---- Fig. 7B: estimated vs reference values ---------------------------------------

# d: name, x (reference / "observed"), y (estimated / "predicted")
plot_drug_scatter <- function(d, xlab, ylab, base_size = 11) {
  d <- d[d$x > 0 & d$y > 0, ]
  rng <- range(c(d$x, d$y), na.rm = TRUE)
  lims <- c(rng[1] / 1.2, rng[2] * 1.2)
  ratio <- d$y / d$x
  d$fold_class <- ifelse(ratio >= 1 / 2 & ratio <= 2, "Within 2-fold",
                         ifelse(ratio >= 1 / 3 & ratio <= 3, "Within 3-fold", "Outside 3-fold"))
  d$panel <- "Estimated vs reference"
  ggplot2::ggplot(d, ggplot2::aes(x, y, label = name)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.8) +
    ggplot2::geom_abline(slope = 1, intercept = log10(2), colour = "skyblue4", linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = -log10(2), colour = "skyblue4", linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = log10(3), colour = "goldenrod4", linetype = "dotted", linewidth = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = -log10(3), colour = "goldenrod4", linetype = "dotted", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = fold_class), size = 2.2, alpha = 0.9) +
    ggrepel::geom_text_repel(size = base_size * 0.25, max.overlaps = 100, min.segment.length = 1, seed = 1) +
    ggplot2::scale_x_log10(limits = lims) +
    ggplot2::scale_y_log10(limits = lims) +
    ggplot2::scale_colour_manual(values = c("Within 2-fold" = "darkgreen", "Within 3-fold" = "orange2",
                                            "Outside 3-fold" = "firebrick")) +
    ggplot2::labs(x = xlab, y = ylab, colour = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(legend.position = "top", plot.margin = grid::unit(c(5, 15, 5, 5), "pt"))
}
