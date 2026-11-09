# =============================================================================
# si_s6_nsb_washout.R - SI Section S6: NSB characterised by washout
#
# Empty HUMIMIC Chip2 chips (PDMS and non-PDMS), diclofenac and FIAU at 10 and
# 100 µM: medium concentration during the washout in compound-free medium
# (desorption) and the PDMS / non-PDMS ratio of the washout plateau. The tidy
# data table is built from the two workbooks.
#
# Linear time axis; no plot titles (the caption carries them).
# Plateau = mean of all washout samples at t >= 8 h (8 and 72 h).
#
# Inputs : data/empty_chip_nsb/empty_chip_nsb_dcf.xlsx, empty_chip_nsb_fiau.xlsx
# Outputs: si_s6_nsb_tidy_data.csv, si_s6_washout_summary.csv,
#          si_s6_washout.png/.pdf (A: desorption curves, B: plateau ratios)
# Runtime: ~10 s
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ---- parameters ----------------------------------------------------------------
XLSX_FILES <- c(pc_data("empty_chip_nsb", "empty_chip_nsb_dcf.xlsx"),
                pc_data("empty_chip_nsb", "empty_chip_nsb_fiau.xlsx"))
MIN_NOMINAL_UM <- 1        # 0.1 µM is below the LLOQ in the washout
PLATEAU_FROM_H <- 8        # plateau = all washout samples at t >= 8 h
COLS_MATERIAL <- c("PDMS" = "#E74C3C", "Non-PDMS" = "#2C7FB8")
COLS_COMPOUND <- c("Diclofenac" = "#27AE60", "FIAU" = "#9B59B6")

# ---- tidy data (all sheets except README and the controls) -----------------------
tidy <- bind_rows(lapply(XLSX_FILES, function(f) {
  bind_rows(lapply(setdiff(excel_sheets(f), "README"), function(s) {
    d <- read_xlsx(f, sheet = s, col_types = "text")
    data.frame(compound = d$compound, material = d$pdms_status, nominal_um = as.numeric(d$conc_uM),
               phase = d$phase, time_h = as.numeric(d$time_h), conc_raw = d$conc)
  }))
}))
tidy <- tidy %>%
  filter(nominal_um > 0) %>%
  mutate(is_lloq = as.integer(grepl("<", conc_raw)),
         conc = ifelse(is_lloq == 1, 0, suppressWarnings(as.numeric(conc_raw)))) %>%   # "< LLOQ" -> 0
  select(compound, material, nominal_um, phase, time_h, conc, is_lloq)
write.csv(tidy, pc_out("si_s6_nsb_tidy_data.csv"), row.names = FALSE)

dat <- tidy %>%
  mutate(material = ifelse(material == "nonPDMS", "Non-PDMS", "PDMS"),
         compound_label = ifelse(compound == "DCF", "Diclofenac", "FIAU"),
         cond_label = sprintf("%s %g µM", compound_label, nominal_um))

# ---- A: desorption curves ----------------------------------------------------------
wash <- dat %>% filter(phase == "Washout", nominal_um >= MIN_NOMINAL_UM)
wash_mean <- wash %>%
  group_by(compound, compound_label, material, nominal_um, cond_label, time_h) %>%
  summarise(C = mean(conc), SD = sd(conc), n = n(), .groups = "drop")
p_curves <- ggplot(wash_mean, aes(time_h, C, colour = material, fill = material)) +
  geom_ribbon(aes(ymin = C - SD, ymax = C + SD), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.3) +
  facet_wrap(~cond_label, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = COLS_MATERIAL) +
  scale_fill_manual(values = COLS_MATERIAL) +
  labs(x = "Time in washout (h)", y = "Released concentration (µM)", colour = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey95"),
        panel.grid.minor = element_blank())

# ---- B: plateau ratio PDMS / non-PDMS ------------------------------------------------
plateau <- wash %>%
  filter(time_h >= PLATEAU_FROM_H) %>%
  group_by(compound, compound_label, material, nominal_um, cond_label) %>%
  summarise(C_plateau = mean(conc), C_se = sd(conc) / sqrt(n()), .groups = "drop")
fold <- plateau %>%
  select(compound, compound_label, nominal_um, cond_label, material, C_plateau, C_se) %>%
  pivot_wider(names_from = material, values_from = c(C_plateau, C_se)) %>%
  mutate(fold_PDMS_over_nonPDMS = `C_plateau_PDMS` / `C_plateau_Non-PDMS`,
         fold_se = fold_PDMS_over_nonPDMS *
           sqrt((`C_se_PDMS` / `C_plateau_PDMS`)^2 + (`C_se_Non-PDMS` / `C_plateau_Non-PDMS`)^2))
write.csv(fold, pc_out("si_s6_washout_summary.csv"), row.names = FALSE)
print(as.data.frame(fold[, c("cond_label", "C_plateau_PDMS", "C_plateau_Non-PDMS", "fold_PDMS_over_nonPDMS", "fold_se")]), digits = 4)

fold_plot <- fold %>%
  mutate(cond_label = factor(cond_label, levels = cond_label[order(compound_label, nominal_um)]))
p_sum <- ggplot(fold_plot, aes(cond_label, fold_PDMS_over_nonPDMS, fill = compound_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_col(width = 0.6, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = fold_PDMS_over_nonPDMS - fold_se, ymax = fold_PDMS_over_nonPDMS + fold_se), width = 0.2) +
  geom_text(aes(y = fold_PDMS_over_nonPDMS + fold_se, label = sprintf("%.1f×", fold_PDMS_over_nonPDMS)),
            vjust = -0.6, size = 4) +
  scale_fill_manual(values = COLS_COMPOUND, name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Released conc. ratio\n(PDMS / non-PDMS)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "top")

p <- p_curves / p_sum + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "A")
save_figure(p, "si_s6_washout", width = 8, height = 5.8)   # print size (was 11 x 7.4 in)
pc_log("done: ", pc_run_dir())
