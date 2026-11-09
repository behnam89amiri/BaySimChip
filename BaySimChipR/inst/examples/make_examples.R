# -----------------------------------------------------------------------------
# make_examples.R - builds the example projects in inst/examples/
#
# Usage (from the BaySimChipR root folder):
#   Rscript inst/examples/make_examples.R
#
# The examples are written in the native project format (.rds), starting
# from the CNBio models in inst/extdata/ (.rds).  If the diclofenac data
# (inst/extdata/DCF_*.csv) are missing, they are exported from the data of the
# analysis code (../paper_code/data/fig3/diclofenac_subsets.xlsx) first.
# -----------------------------------------------------------------------------

source(file.path("tools", "load_engine.R"))   # engine
options(baysimchip.app_dir = normalizePath("."))

ext_dir <- file.path("inst", "extdata")
ex_dir <- file.path("inst", "examples")
dir.create(ext_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. diclofenac data (only if missing) ---------------------------------------
data_dir <- file.path("..", "paper_code", "data", "fig3")
dcf_sheets <- c(DCF_gut_liver_SD_1000uM = "D_G-L_SD_PO_1000uM",
                DCF_gut_liver_DD_100uM = "D_G-L_DD_PO_100uM")
xl <- file.path(data_dir, "diclofenac_subsets.xlsx")
for (nm in names(dcf_sheets)) {
  out <- file.path(ext_dir, paste0(nm, ".csv"))
  if (!file.exists(out) && file.exists(xl)) {
    d <- as.data.frame(readxl::read_excel(xl, sheet = dcf_sheets[[nm]], col_types = "text"))
    d <- d[, c("Location", "Time_h", "Conc", "Regimen", "Config", "Route", "Dose_uM")]
    utils::write.csv(d, out, row.names = FALSE)
  }
}

# ---- helpers -----------------------------------------------------------------------
# readable names and a tidy layout for the CNBio TL6 gut-liver models
cnbio_names <- c(Comp0 = "Apical", Comp1 = "Basolateral", Comp2 = "Liver", Comp3 = "Circulation",
                 Comp5 = "NSB_liver", Comp6 = "NSB_apical")
cnbio_xy <- list(Apical = c(-150, -150), Basolateral = c(-150, 20), Liver = c(170, 20),
                 Circulation = c(10, 170), NSB_liver = c(300, 140), NSB_apical = c(-330, -150))

tidy_cnbio <- function(p, name, description) {
  for (old in intersect(names(cnbio_names), p$compartments$name)) {
    p <- bsc_rename(p, "compartment", old, cnbio_names[[old]])
  }
  for (i in seq_len(nrow(p$compartments))) {
    xy <- cnbio_xy[[p$compartments$name[i]]]
    if (!is.null(xy)) {
      p$compartments$x[i] <- xy[1]
      p$compartments$y[i] <- xy[2]
    }
  }
  p$name <- name
  p$description <- description
  p
}

set_vols <- function(p, v) {
  for (nm in names(v)) p <- bsc_set_param(p, "compartment", nm, "volume", v[[nm]])
  p
}

read_dcf <- function(nm) {
  d <- bsc_read_data_file(file.path(ext_dir, paste0(nm, ".csv")))
  # apical data are excluded from the fits (as in the paper)
  bsc_observed_from_table(d, "Time_h", "Conc", comp_col = "Location",
                          comp_map = c(Liver = "Liver", Basolateral = "Basolateral"),
                          blq = "zero")
}

save_example <- function(p, file) {
  v <- bsc_validate(p, "structure")
  stopifnot(!length(v$errors))
  bsc_write_project(p, file.path(ex_dir, paste0(file, ".rds")))
  cat("wrote", file, "\n")
}

paper_defaults <- function(p) {
  p <- set_vols(p, c(Apical = 0.35, Basolateral = 0.8, Liver = 1.6, Circulation = 0.2))
  p <- bsc_set_param(p, "compartment", "Liver", "cl", 0.1)
  p <- bsc_set_param(p, "perm", "P.S.0", "ps", 0.01)
  p$settings$dose_unit <- "nmol"
  # 1 uM in the apical compartment (0.35 mL)
  bsc_set_dosing(p, "Apical", list(mode = "single", single = 0.35))
}

# ---- 2. CNBio gut-liver base model (inst/extdata/Kev2DAnalysis.rds) --------
base <- bsc_read_project(file.path(ext_dir, "Kev2DAnalysis.rds"))
base <- tidy_cnbio(base, "Gut-liver chip (CN Bio TL6)",
                   paste("Gut-liver model of the CN Bio PhysioMimix TL6 chip: apical and basolateral",
                         "gut compartments separated by the epithelium (permeable surface), liver and",
                         "circulation connected by medium flows. Doses in mg, molar mass 300 g/mol."))
base$settings$sim_time <- 72
save_example(base, "CNBio_gut_liver")

# ---- 3. Evaporation (Fig. 4A) -------------------------------------------
ev <- paper_defaults(base)
ev <- bsc_set_param(ev, "compartment", "Apical", "kev", 0.003)
ev <- bsc_set_param(ev, "compartment", "Liver", "kev", 0.0055)
ev$settings$sim_time <- 110
ev$name <- "Gut-liver chip with evaporation"
ev$description <- paste("Generic compound (CL 0.1 mL/h, PS 0.01 mL/h), 1 µM dosed to the apical",
                        "compartment; evaporation 0.003 mL/h (apical) and 0.0055 mL/h (liver).",
                        "Set Kev to 0 and simulate again to compare with 'no evaporation'.")
save_example(ev, "CNBio_evaporation")

# ---- 4. Sampling (Fig. 5A) ------------------------------------------------
sa <- paper_defaults(base)
sa <- bsc_set_param(sa, "compartment", "Apical", "kev", 0)
sa <- bsc_set_param(sa, "compartment", "Liver", "kev", 0)
sa <- bsc_set_exchange(sa, "Liver", data.frame(time = c(0, 1, 2, 4, 8, 24, 48), removed = 0.05, added = 0))
sa <- bsc_set_exchange(sa, "Apical", data.frame(time = c(8, 24), removed = 0.05, added = 0))
sa$settings$sim_time <- 100
sa$name <- "Gut-liver chip with sampling"
sa$description <- paste("Generic compound; 0.05 mL samples from the liver (0-48 h) and from the",
                        "apical compartment (8 and 24 h). A sample at t = 0 has no effect (the",
                        "initial volume is used).")
save_example(sa, "CNBio_sampling")

# ---- 5. Non-specific binding (Fig. 6A) ---------------------------------------
nsb <- bsc_read_project(file.path(ext_dir, "NSBpara2DAnalysis.rds"))
nsb <- bsc_add_compartment(nsb, "Comp6", type = "solid", volume = 1, fu = 1, cl = 0, kev = 0)
nsb <- bsc_add_perm(nsb, "Comp6", "Comp0", ps = 0, name = "P.S.3")
nsb <- tidy_cnbio(nsb, "Gut-liver chip with non-specific binding",
                  paste("Two solid compartments represent the chip material on the liver side and",
                        "on the apical side. NSB exchange rate PS_NSB 0.05 mL/h (liver side) and",
                        "0.01 mL/h (apical side, 1/5 of the area); K_media/plastic = 1 for both."))
nsb <- paper_defaults(nsb)
nsb <- bsc_set_param(nsb, "compartment", "NSB_liver", "fu", 1)
nsb <- bsc_set_param(nsb, "compartment", "NSB_apical", "fu", 1)
nsb <- bsc_set_param(nsb, "perm", "P.S.2", "ps", 0.05)
nsb <- bsc_set_param(nsb, "perm", "P.S.3", "ps", 0.01)
nsb <- bsc_set_exchange(nsb, "Liver", bsc_empty_exchange())
nsb$settings$sim_time <- 100
save_example(nsb, "CNBio_NSB")

# ---- 6. Diclofenac fitting examples (Fig. 3) -----------------------------------
dcf_base <- function(p) {
  p <- set_vols(p, c(Apical = 0.2, Basolateral = 0.8, Liver = 1.8, Circulation = 0.2))
  p <- bsc_set_param(p, "compartment", "Apical", "kev", 0)
  p <- bsc_set_param(p, "compartment", "Liver", "kev", 0)
  p <- bsc_set_param(p, "compartment", "Liver", "cl", 0.15, estimate = TRUE)
  p <- bsc_set_param(p, "perm", "P.S.0", "ps", 0.015, estimate = TRUE)
  p$settings$dose_unit <- "nmol"
  p$settings$molar_mass <- 296.15
  p
}
if (file.exists(file.path(ext_dir, "DCF_gut_liver_SD_1000uM.csv"))) {
  sd <- dcf_base(base)
  sd <- bsc_set_dosing(sd, "Apical", list(mode = "single", single = 0.2 * 1000))
  sd$observed <- read_dcf("DCF_gut_liver_SD_1000uM")
  sd$settings$sim_time <- 72
  sd$name <- "Diclofenac, gut-liver, single dose 1000 µM"
  sd$description <- paste("Diclofenac in the CN Bio TL6 gut-liver chip, arm D_G-L_SD_PO_1000uM",
                          "(liver data; apical data not used, as in the paper). CL (liver) and",
                          "PS (gut epithelium) are estimated. Liver volume 1.8 mL, as used for",
                          "the published fit.")
  save_example(sd, "DCF_gut_liver_SD_1000uM")

  dd <- dcf_base(base)
  dd <- bsc_set_dosing(dd, "Apical", list(mode = "repeated", repeated_dose = 0.2 * 100,
                                         repeated_interval = 24))
  dd <- bsc_set_dosing(dd, "Liver", list(mode = "single", single = 0, first_dose_from_data = TRUE))
  dd <- bsc_set_exchange(dd, "Liver", data.frame(time = seq(24, 96, by = 24), removed = 0.867,
                                                 added = 0.867))
  dd$observed <- read_dcf("DCF_gut_liver_DD_100uM")
  dd$settings$sim_time <- 110
  dd$name <- "Diclofenac, gut-liver, daily dose 100 µM"
  dd$description <- paste("Diclofenac in the CN Bio TL6 gut-liver chip, arm D_G-L_DD_PO_100uM:",
                          "daily apical doses and daily exchange of 0.867 mL liver medium; the",
                          "initial liver concentration is taken from the first measurement.",
                          "CL and PS are estimated.")
  save_example(dd, "DCF_gut_liver_DD_100uM")
}
