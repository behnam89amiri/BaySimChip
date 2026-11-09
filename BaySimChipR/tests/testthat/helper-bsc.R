# Loaded automatically by testthat before the tests: sources the app code.
suppressPackageStartupMessages({
  library(shiny)
})
bsc_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
for (f in sort(list.files(file.path(bsc_root, "R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f, local = FALSE)
}
options(baysimchip.app_dir = bsc_root)

bsc_ref_path <- file.path(bsc_root, "tests", "reference", "legacy_reference.rds")
bsc_extdata <- function(name) file.path(bsc_root, "inst", "extdata", name)

# CNBio model of the analysis code (inst/extdata/Kev2DAnalysis.rds) with the settings of
# the Fig. 4A / 5A simulations; dose in "ng" = unit of the analysis code
legacy_cnbio <- function(file = "Kev2DAnalysis.rds") {
  p <- bsc_read_project(bsc_extdata(file))
  p$settings$dose_unit <- "ng"
  p <- bsc_set_param(p, "compartment", "Comp2", "cl", 0.1)
  p <- bsc_set_param(p, "perm", "P.S.0", "ps", 0.01)
  vols <- c(Comp0 = 0.35, Comp1 = 0.8, Comp2 = 1.6, Comp3 = 0.2)
  for (nm in names(vols)) p <- bsc_set_param(p, "compartment", nm, "volume", vols[[nm]])
  bsc_set_dosing(p, "Comp0", list(mode = "single", single = 0.35 * 1 * p$settings$molar_mass))
}

# Fitting set-up of the original Fig. 3 analysis for the arms stored in the reference
legacy_fit_project <- function(rr, v_liver = rr$v_liver) {
  p <- bsc_read_project(bsc_extdata("Kev2DAnalysis.rds"))
  p$settings$dose_unit <- "ng"
  mm <- p$settings$molar_mass
  p <- bsc_set_param(p, "compartment", "Comp0", "kev", 0)
  p <- bsc_set_param(p, "compartment", "Comp2", "kev", 0)
  p <- bsc_set_param(p, "compartment", "Comp0", "volume", 0.2)
  p <- bsc_set_param(p, "compartment", "Comp2", "volume", v_liver)
  p <- bsc_set_param(p, "compartment", "Comp3", "volume", 0.2)
  p <- bsc_set_param(p, "compartment", "Comp2", "cl", 0.15, estimate = TRUE)
  p <- bsc_set_param(p, "perm", "P.S.0", "ps", 0.015, estimate = TRUE)
  obs <- rr$observed[rr$observed$location != "Apical", ]
  map <- c(Liver = "Comp2", Basolateral = "Comp1", Apical = "Comp0")
  p$observed <- data.frame(compartment = unname(map[obs$location]), time = obs$time,
                           conc = obs$conc, stringsAsFactors = FALSE)
  if (grepl("_SD_", rr$sheet)) {
    p <- bsc_set_dosing(p, "Comp0", list(mode = "single", single = 0.2 * 1000 * mm))
  } else {
    p <- bsc_set_dosing(p, "Comp0", list(mode = "repeated", repeated_dose = 0.2 * 100 * mm,
                                         repeated_interval = 24))
    p <- bsc_set_dosing(p, "Comp2", list(mode = "single", first_dose_from_data = TRUE))
    p <- bsc_set_exchange(p, "Comp2", data.frame(time = seq(24, 110, by = 24), removed = 0.867,
                                                 added = 0.867))
  }
  p
}

max_rel_diff <- function(a, b, floor = 1e-12) max(abs(a - b) / pmax(abs(b), floor))
