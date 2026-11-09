# =============================================================================
# run_all.R - regenerate all figures and result tables (main text Figs. 3-7,
# SI Sections S5-S7) and the cross-check with the engine of the BaySimChip app
#
# Usage, from the paper_code/ folder:
#   Rscript run_all.R                 # everything (about 10-20 min on one core)
#   Rscript run_all.R fig4 fig7       # only scripts whose name starts with these
#
# Each script runs in its own R process (Rscript) with PAPER_CODE_RUN_DIR set,
# so that all outputs of one run go to one new folder: PAPER_CODE_RUN_DIR if
# set, otherwise output/run_YYYY-MM-DD (run_YYYY-MM-DD_2, ... if it exists). Per-script logs go
# to <run>/logs/, run times to <run>/timings.csv, the R session to
# <run>/sessionInfo.txt. A failing script is reported and the others continue.
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- normalizePath(if (length(.f)) dirname(.f) else ".", winslash = "/")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)

SCRIPTS <- c(
  "scripts/fig3_fits.R",
  "scripts/fig3_plots.R",
  "scripts/fig4ab_evaporation.R",
  "scripts/fig5ab_sampling.R",
  "scripts/exchange_imbalance.R",
  "scripts/fig6ab_nsb.R",
  "scripts/fig456c_bias_maps.R",
  "scripts/fig6d_empty_chip.R",
  "scripts/fig7a_archetypes.R",
  "scripts/fig7b_drugs.R",
  "scripts/figures_composite.R",
  "scripts/si_s5_memantine.R",
  "scripts/si_s6_nsb_washout.R",
  "scripts/si_s7_fiau_tktd.R",
  "verification/check_app_engine.R"
)

only <- commandArgs(trailingOnly = TRUE)
if (length(only)) {
  keep <- Reduce(`|`, lapply(only, function(p) startsWith(basename(SCRIPTS), p)))
  SCRIPTS <- SCRIPTS[keep]
}

run_dir <- pc_run_dir()                       # new folder for this run
dir.create(file.path(run_dir, "logs"), showWarnings = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
pc_log("run folder: ", run_dir)

timings <- data.frame(script = character(), seconds = numeric(), status = character())
for (s in SCRIPTS) {
  log_file <- file.path(run_dir, "logs", paste0(sub("\\.R$", "", basename(s)), ".log"))
  pc_log("running ", s)
  t0 <- Sys.time()
  # PAPER_CODE_ROOT / PAPER_CODE_RUN_DIR were set by pc_setup() / pc_run_dir()
  # and are inherited by the child process
  status <- system2(rscript, shQuote(file.path(.root, s)), stdout = log_file, stderr = log_file)
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- identical(as.integer(status), 0L)
  timings <- rbind(timings, data.frame(script = s, seconds = round(sec, 1), status = if (ok) "ok" else "FAILED"))
  pc_log(sprintf("  %s in %.0f s%s", if (ok) "done" else "FAILED", sec, if (ok) "" else paste0(" - see ", log_file)))
}

write.csv(timings, pc_out("timings.csv"), row.names = FALSE)
writeLines(capture.output({
  suppressPackageStartupMessages(lapply(c("deSolve", "minpack.lm", "ggplot2", "patchwork",
                                          "dplyr", "tidyr", "readxl", "ggrepel", "scales", "ragg",
                                          "viridisLite", "pracma"), library, character.only = TRUE))
  print(sessionInfo())
}), pc_out("sessionInfo.txt"))
print(timings)
pc_log(sprintf("total %.1f min; outputs in %s", sum(timings$seconds) / 60, run_dir))
if (any(timings$status != "ok")) quit(status = 1)
