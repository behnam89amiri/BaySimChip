# =============================================================================
# exchange_imbalance.R - media-exchange imbalance (main text, section
# "Effect of sampling"; no figure; SI Table S3 "CNBIO-EXCHANGE")
#
# Daily apical dosing (1 µM) for one week with a daily 1.5 mL exchange of the
# liver medium; the added volume equals the removed volume, or is 0.16 mL
# (about 10%) smaller or larger.
#
# Inputs : data/models/Kev2DAnalysis.rds
# Outputs: exchange_profiles.csv, exchange_summary.csv, exchange_<compartment>.png
# Runtime: ~30 s
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages(library(ggplot2))

# ---- parameters ----------------------------------------------------------------
MODEL_FILE <- pc_data("models", "Kev2DAnalysis.rds")
VOLUMES <- c(0.35, 0.1, 2.5, 0.1)      # mL: apical, basolateral, liver, circulation (SI Table S3)
FINAL_TIME_H <- 0.99 * 7 * 24          # one week (166.32 h)
DOSE_INTERVAL_H <- 24                  # daily apical dose, 1 µM
EXCHANGE_VOLUME_ML <- 1.5              # removed from the liver compartment every 24 h
VOLUME_DIFFERENCE_ML <- 0.16           # added volume = 1.5, 1.5 - 0.16, 1.5 + 0.16 mL (~10%)
MOLAR_MASS <- 300                      # placeholder; doses are in µM x mL

# ---- simulation ----------------------------------------------------------------
app <- setup_generic_compound(MODEL_FILE, final_time = FINAL_TIME_H, volumes = VOLUMES)
app$CompData[[1]]$Kev <- 0
app$CompData[[3]]$Kev <- 0
exchange_times <- seq(0, 1.1 * app$FinalTime, by = 24)   # 0 ... 168 h; last one not applied, t = 0 ignored
app$CompData[[3]]$ExMedia <- exchange_table(exchange_times, EXCHANGE_VOLUME_ML, EXCHANGE_VOLUME_ML)
app$MolarMass <- MOLAR_MASS
app$CompData[[1]]$DoseCond <- 2
app$CompData[[1]]$DoseMultData[1, 2] <- DOSE_INTERVAL_H
app$CompData[[1]]$DoseMultData[1, 1] <- app$CompData[[1]]$Vol * 1 * app$MolarMass

conds <- data.frame(cond = c("added = removed", "added 10% less than removed", "added 10% more than removed"),
                    added = c(EXCHANGE_VOLUME_ML, EXCHANGE_VOLUME_ML - VOLUME_DIFFERENCE_ML,
                              EXCHANGE_VOLUME_ML + VOLUME_DIFFERENCE_ML))
prof <- do.call(rbind, lapply(seq_len(nrow(conds)), function(i) {
  app_i <- copy_app(app)
  app_i$CompData[[3]]$ExMedia[, 3] <- conds$added[i]
  app_i <- ConstructDoseMediaEx(app_i)
  cbind(cond = conds$cond[i], simulate_profiles(app_i))
}))
write.csv(prof, pc_out("exchange_profiles.csv"), row.names = FALSE)

# liver concentration on the last day (144 h to end) relative to the balanced case
last_day <- prof[prof$compartment == 3 & prof$time_h >= 144, ]
summ <- do.call(rbind, lapply(conds$cond, function(cn) {
  d <- last_day[last_day$cond == cn, ]
  data.frame(cond = cn, added_mL = conds$added[conds$cond == cn],
             liver_Cmax_last_day = max(d$conc), liver_Cmin_last_day = min(d$conc),
             liver_C_end = d$conc[nrow(d)])
}))
ref <- summ[summ$cond == "added = removed", ]
summ$Cmax_change_pct <- 100 * (summ$liver_Cmax_last_day / ref$liver_Cmax_last_day - 1)
summ$C_end_change_pct <- 100 * (summ$liver_C_end / ref$liver_C_end - 1)
summ$imbalance_pct <- 100 * VOLUME_DIFFERENCE_ML / EXCHANGE_VOLUME_ML * c(0, -1, 1)
write.csv(summ, pc_out("exchange_summary.csv"), row.names = FALSE)
print(summ, digits = 4)

prof$cond <- factor(prof$cond, levels = conds$cond)
for (j in c(1, 3)) {
  p <- plot_profiles(prof, j, title = c("1" = "Gut apical", "3" = "Liver")[[as.character(j)]],
                     ylab = lab("conc_au"), colours = setNames(scales::hue_pal()(3), conds$cond),
                     legend_inside = c(0.98, 0.98), base_size = 12)
  save_panel(p, paste0("exchange_", c("1" = "gut_apical", "3" = "liver")[[as.character(j)]]), width = 8.5, height = 4)
}
pc_log("done: ", pc_run_dir())
