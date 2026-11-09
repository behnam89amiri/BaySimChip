# =============================================================================
# check_app_engine.R - independent check of Fig. 4A with the engine of the
# BaySimChip R Shiny app (BaySimChipR, separate code)
#
# Simulates the app's example "CNBio_evaporation" (generic compound, Kev
# 0.003 / 0.0055 mL/h, 110 h) with and without evaporation and compares the
# liver profile with fig4a_profiles.csv of the same run.
# BaySimChipR is looked up in ../BaySimChipR (next to paper_code) or in
# the environment variable BAYSIMCHIPR_DIR; if it is not found the check is
# skipped (written as "skipped" to the output table).
#
# Inputs : fig4a_profiles.csv (run folder), BaySimChipR/tools/load_engine.R
# Outputs: app_engine_check.csv
# Runtime: ~10 s
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)

# ---- parameters ----------------------------------------------------------------
BAYSIMCHIPR_DIR <- Sys.getenv("BAYSIMCHIPR_DIR", file.path(dirname(pc_root()), "BaySimChipR"))
FINAL_TIME_H <- 110

# ---- check -----------------------------------------------------------------------
loader <- file.path(BAYSIMCHIPR_DIR, "tools", "load_engine.R")
if (!file.exists(loader)) {
  out <- data.frame(condition = NA, quantity = "app engine check", app_engine = NA, paper_code = NA,
                    abs_diff = NA, status = paste("skipped: BaySimChipR not found in", BAYSIMCHIPR_DIR))
} else {
  source(loader)
  ours <- read.csv(pc_in_run("fig4a_profiles.csv"))
  p <- bsc_load_example("CNBio_evaporation")
  p0 <- p
  p0$compartments$kev <- 0
  out <- do.call(rbind, lapply(list(list("w/ evaporation", p), list("w/o evaporation", p0)), function(x) {
    sim <- bsc_simulate(x[[2]], final_time = FINAL_TIME_H)
    liver_app <- sim$conc[, "Liver"]
    o <- ours[ours$cond == x[[1]] & ours$compartment == 3, ]
    common <- intersect(round(sim$time, 6), round(o$time_h, 6))
    a <- liver_app[match(common, round(sim$time, 6))]
    b <- o$conc[match(common, round(o$time_h, 6))]
    data.frame(condition = x[[1]],
               quantity = c("liver Cmax", "liver Tmax (h)", "max |difference| of liver profile (µM)"),
               app_engine = c(max(liver_app), sim$time[which.max(liver_app)], NA),
               paper_code = c(max(o$conc), o$time_h[which.max(o$conc)], NA),
               abs_diff = c(abs(max(liver_app) - max(o$conc)),
                            abs(sim$time[which.max(liver_app)] - o$time_h[which.max(o$conc)]),
                            max(abs(a - b))),
               status = "ok")
  }))
}
write.csv(out, pc_out("app_engine_check.csv"), row.names = FALSE)
print(out, digits = 6)
pc_log("done: ", pc_run_dir())
