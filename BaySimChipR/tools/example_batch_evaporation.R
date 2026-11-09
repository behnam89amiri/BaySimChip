# =============================================================================
# example_batch_evaporation.R - using the BaySimChip engine without the app
#
# Re-creates the evaporation simulation of the paper (Fig. 4A) with the engine
# functions, as an analysis script would.
#
#   Rscript tools/example_batch_evaporation.R [output_dir]
# =============================================================================

source(file.path("tools", "load_engine.R"))   # run from the BaySimChipR folder

out_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_dir)) out_dir <- tempdir()

# Project made in the app (example "CNBio_evaporation"): CNBio
# gut-liver chip, generic compound (CL 0.1 mL/h, PS 0.01 mL/h), 1 uM apical.
p <- bsc_load_example("CNBio_evaporation")
print(bsc_parameter_table(p)[, c("key", "value", "unit")])

# Scenarios: one column per parameter to change ("element|parameter")
scenarios <- data.frame(
  scenario = c("w/ evaporation", "w/o evaporation"),
  "Apical|kev" = c(0.003, 0),
  "Liver|kev" = c(0.0055, 0),
  check.names = FALSE
)
res <- bsc_run_scenarios(p, scenarios, final_time = 110)

# exposure parameters per scenario and compartment
print(attr(res, "exposure"), digits = 4)

# long data.frame (scenario, time, compartment, conc) -> any plotting tool
utils::write.csv(res, file.path(out_dir, "evaporation_profiles.csv"), row.names = FALSE)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  g <- ggplot(res[res$compartment %in% c("Apical", "Liver"), ],
              aes(time, conc, colour = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~compartment, scales = "free_y") +
    labs(x = "Time (h)", y = "Concentration (µM)", colour = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "evaporation_profiles.png"), g, width = 7, height = 3.5, dpi = 150)
}
cat("Results written to", normalizePath(out_dir), "\n")

# The same engine fits data, e.g.:
#   f <- bsc_fit(bsc_load_example("DCF_gut_liver_SD_1000uM"))
#   f$estimates[, c("label", "estimate", "se", "rse_percent", "ci_lower", "ci_upper")]
