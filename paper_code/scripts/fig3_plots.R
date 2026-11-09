# =============================================================================
# fig3_plots.R - Fig. 3 panels and composite (diclofenac PK across OoC designs)
#
#   A: measured (points) and fitted (lines) liver/circuit concentrations of
#      four representative arms
#   B: arm-wise estimates of P and CL_int (+/- RSE) with literature comparators
#   C: estimates grouped by regimen, flow and organ coupling
#
# Inputs : fig3_results.csv, fig3_fit_profiles.csv, fig3_observed.csv
#          (written by fig3_fits.R into the same run folder)
# Outputs: fig3B_estimates.png, fig3C_systematic.png, Fig3.png, Fig3.pdf
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
FIG3A_ARMS <- c("D_G-L_SD_PO_100uM", "S_L_SD_IV_100uM", "D_G_SD_PO_100uM", "D_G-L_DD_PO_40uM")
LIT_ROWS   <- c("Invitro-dual-organ", "Invitro-single-organ")
CL_ORDER <- c("D_G-L_SD_PO_1000uM", "D_G-L_SD_PO_100uM", "D_G-L_SD_PO_40uM",
              "D_G-L_DD_PO_1000uM", "D_G-L_DD_PO_100uM", "D_G-L_DD_PO_40uM",
              "D_L_SD_IV_100uM", "D_L_DD_IV_100uM", "S_L_SD_IV_100uM", LIT_ROWS)
PS_ORDER <- c("D_G-L_SD_PO_1000uM", "D_G-L_SD_PO_100uM", "D_G-L_SD_PO_40uM",
              "D_G-L_DD_PO_1000uM", "D_G-L_DD_PO_100uM", "D_G-L_DD_PO_40uM",
              "D_G_SD_PO_100uM", "D_G_DD_PO_100uM", "S_G_SD_PO_100uM", LIT_ROWS)

# ---- data ----------------------------------------------------------------------
res <- read.csv(pc_in_run("fig3_results.csv"))
sim <- read.csv(pc_in_run("fig3_fit_profiles.csv"))
obs <- read.csv(pc_in_run("fig3_observed.csv"))
# display names: the dynamic liver-only chips were dosed apically into a cell-free
# insert (study protocol), so their route code is PO (data sheets keep IV)
DISPLAY <- c(D_L_SD_IV_100uM = "D_L_SD_PO_100uM", D_L_DD_IV_100uM = "D_L_DD_PO_100uM")
rename_arms <- function(x) ifelse(x %in% names(DISPLAY), DISPLAY[x], x)
res$SheetName <- rename_arms(res$SheetName)
sim$arm <- rename_arms(sim$arm)
obs$arm <- rename_arms(obs$arm)
CL_ORDER <- rename_arms(CL_ORDER)
PS_ORDER <- rename_arms(PS_ORDER)

est <- rbind(
  data.frame(SheetName = res$SheetName, Parameter = "CL", Value = res$Clint, RSE = res$RSECL),
  data.frame(SheetName = res$SheetName, Parameter = "PS", Value = res$P, RSE = res$RSEPS)
)
est$is_ref <- est$SheetName %in% LIT_ROWS
est <- est[!is.na(est$Value), ]

# arm codes W_X_Y_Z_Dose -> flow (D/S), system (G-L, G, L), regimen (SD/DD)
meta <- regmatches(est$SheetName, regexec("^([DS])_([GL-]+)_([SD]{2})_", est$SheetName))
arms <- est[!est$is_ref, ]
mm <- do.call(rbind, meta[!est$is_ref])
arms$Flow    <- c(D = "Dynamic", S = "Static")[mm[, 2]]
arms$System  <- ifelse(mm[, 3] == "G-L", "Multi-Organ", "Single-Organ")
arms$Regimen <- c(SD = "Single", DD = "Daily")[mm[, 4]]
comp <- rbind(
  data.frame(Parameter = arms$Parameter, ComparisonType = "Daily vs Single",
             Group = arms$Regimen, Value = arms$Value),
  data.frame(Parameter = arms$Parameter, ComparisonType = "Dynamic vs Static",
             Group = arms$Flow, Value = arms$Value),
  data.frame(Parameter = arms$Parameter, ComparisonType = "Multi-organ vs Single-organ",
             Group = arms$System, Value = arms$Value)
)
comp$Group <- factor(comp$Group, levels = c("Daily", "Single", "Dynamic", "Static", "Multi-Organ", "Single-Organ"))

# ---- panels ------------------------------------------------------------------------
pA <- plot_fig3_fits(obs, sim, FIG3A_ARMS, base_size = 10)
pB_P  <- plot_fig3_estimates(est, "PS", PS_ORDER, lab("P"))
pB_CL <- plot_fig3_estimates(est, "CL", CL_ORDER, lab("CLint"))
pB <- (pB_P | pB_CL) + plot_layout(guides = "collect") & theme(legend.position = "top")
pC_P  <- plot_fig3_systematic(comp, "PS", "Gut permeability", lab("P"), show_x = FALSE)
pC_CL <- plot_fig3_systematic(comp, "CL", "Liver clearance", lab("CLint"), show_x = TRUE)
pC <- pC_P / pC_CL

save_panel(pB, "fig3B_estimates", width = 7.5, height = 5.5)
save_panel(pC, "fig3C_systematic", width = 8.5, height = 5.5)

# composite at print size: RSC double-column width (17.1 cm), base font 8 pt, one row per panel
FIG_W <- 6.73; BASE <- 8
cA <- plot_fig3_fits(obs, sim, FIG3A_ARMS, base_size = BASE)
cB <- (plot_fig3_estimates(est, "PS", PS_ORDER, lab("P"), base_size = BASE) |
         plot_fig3_estimates(est, "CL", CL_ORDER, lab("CLint"), base_size = BASE)) +
  plot_layout(guides = "collect") & theme(legend.position = "top")
cC <- plot_fig3_systematic(comp, "PS", "Gut permeability", lab("P"), show_x = FALSE, base_size = BASE, free_x = TRUE) /
  plot_fig3_systematic(comp, "CL", "Liver clearance", lab("CLint2"), show_x = TRUE, base_size = BASE, free_x = TRUE)
fig3 <- wrap_elements(full = cA) / wrap_elements(full = cB) / wrap_elements(full = cC) +
  plot_layout(heights = c(1.7, 2.5, 2.4)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = BASE + 3, face = "bold"))
save_figure(fig3, "Fig3", width = FIG_W, height = 6.6)
pc_log("done: ", pc_run_dir())
