# -----------------------------------------------------------------------------
# validation_report.R - prints the comparison of the engine with the analysis
# code of the paper (tests/reference/legacy_reference.rds) and with the Fig. 3
# estimates of the original analysis.
#
#   Rscript tests/validation_report.R        (from the BaySimChipR folder)
# -----------------------------------------------------------------------------

source(file.path("tools", "load_engine.R"))   # engine
options(baysimchip.app_dir = normalizePath("."))
owd <- setwd(file.path("tests", "testthat"))
source("helper-bsc.R")
setwd(owd)
ref <- readRDS(bsc_ref_path)

rel <- function(a, b, floor = 1e-12) max(abs(a - b) / pmax(abs(b), floor))
line <- function(...) cat(sprintf(...), "\n", sep = "")

cat("\n== 1. Simulation: evaporation (Fig. 4; Kev2DAnalysis.rds, CL 0.1, PS 0.01,",
    "V 0.35/0.8/1.6/0.2 mL, 1 uM apical, 110 h) ==\n")
p <- legacy_cnbio()
for (cond in names(ref$evaporation)) {
  r <- ref$evaporation[[cond]]
  pc <- bsc_set_param(p, "compartment", "Comp0", "kev", r$kev_apical)
  pc <- bsc_set_param(pc, "compartment", "Comp2", "kev", r$kev_liver)
  sim <- bsc_simulate(pc, final_time = 110)
  line("  %-16s n = %d time points | max rel. diff liver = %.2e, apical = %.2e, all = %.2e",
       cond, length(sim$time), rel(sim$conc[, "Comp2"], r$yc[, 3]), rel(sim$conc[, "Comp0"], r$yc[, 1]),
       rel(sim$conc, r$yc))
}

cat("\n== 2. Simulation: sampling (Fig. 5) ==\n")
p2 <- bsc_set_param(bsc_set_param(p, "compartment", "Comp0", "kev", 0), "compartment", "Comp2", "kev", 0)
for (cond in names(ref$sampling)) {
  r <- ref$sampling[[cond]]
  pc <- bsc_set_exchange(p2, "Comp2", data.frame(time = c(0, 1, 2, 4, 8, 24, 48), removed = r$samp_vol, added = 0))
  pc <- bsc_set_exchange(pc, "Comp0", data.frame(time = c(8, 24), removed = r$samp_vol, added = 0))
  sim <- bsc_simulate(pc, final_time = 100)
  line("  sample volume %-8s max abs. diff = %.2e", cond, max(abs(sim$conc - r$yc)))
}

cat("\n== 3. Simulation: non-specific binding (Fig. 6) ==\n")
pn <- bsc_read_project(bsc_extdata("NSBpara2DAnalysis.rds"))
pn$settings$dose_unit <- "ng"
pn <- bsc_add_compartment(pn, "Comp6", type = "solid", volume = 1, fu = 0, cl = 0, kev = 0)
pn <- bsc_add_perm(pn, "Comp6", "Comp0", ps = 0, name = "P.S.3")
pn <- bsc_set_params(pn, list("Comp2|cl" = 0.1, "P.S.0|ps" = 0.01, "Comp0|volume" = 0.35,
                              "Comp1|volume" = 0.8, "Comp2|volume" = 1.6, "Comp3|volume" = 0.2))
pn <- bsc_set_dosing(pn, "Comp0", list(mode = "single", single = 0.35 * pn$settings$molar_mass))
for (r in ref$nsb) {
  pc <- bsc_set_params(pn, list("Comp5|fu" = r$K_media_plastic, "Comp6|fu" = r$K_media_plastic,
                                "P.S.2|ps" = r$PS, "P.S.3|ps" = r$PS / 5))
  sim <- bsc_simulate(pc, final_time = 100)
  solid_new <- sweep(sim$conc[, 5:6], 2, rep(r$K_media_plastic, 2), `*`)
  line("  K = %-4g PS = %-5g fluid max rel. diff = %.2e | solid: analysis code / app = fu_solid (max rel. dev. %.1e)",
       r$K_media_plastic, r$PS, rel(sim$conc[, 1:4], r$yc[, 1:4], 1e-9), rel(solid_new, r$yc[, 5:6], 1e-9))
}

cat("\n== 4. Fitting: diclofenac, set-up of the original Fig. 3 analysis (CL liver + PS apical-basolateral) ==\n")
fmt_fit <- function(f) {
  e <- f$estimates
  sprintf("CL = %.8f mL/h (RSE %.4f %%) PS = %.8f mL/h (RSE %.4f %%)",
          e$estimate[1], e$rse_percent[1], e$estimate[2], e$rse_percent[2])
}
for (k in c("fit_sd", "fit_dd", "fit_sd_v18", "fit_dd_v18")) {
  rr <- ref[[k]]
  f <- bsc_fit(legacy_fit_project(rr), final_time = rr$final_time)
  line("  %s (V_liver %.1f mL)", rr$sheet, rr$v_liver)
  line("    new engine : %s", fmt_fit(f))
  line("    analysis   : CL = %.8f mL/h (RSE %.4f %%) PS = %.8f mL/h (RSE %.4f %%)",
       rr$estimates[["CL"]], rr$rse_percent[["CL"]], rr$estimates[["PS"]], rr$rse_percent[["PS"]])
  line("    rel. diff  : CL %.1e, PS %.1e | %d model runs, %.2f s",
       abs(f$estimates$estimate[1] / rr$estimates[["CL"]] - 1),
       abs(f$estimates$estimate[2] / rr$estimates[["PS"]] - 1), f$n_eval, f$elapsed)
  pub <- ref$published[ref$published$SheetName == rr$sheet, ]
  line("    as in Fig. 3: CLint = %.6f uL/min/1e6 cells, P = %.5f 1e-6 cm/s | stored Fig. 3 values: %.6f, %.5f",
       f$estimates$estimate[1] * 1000 / 60 / 0.6, f$estimates$estimate[2] / (3600 * 0.33) * 1e6,
       pub$Clint, pub$P)
}
cat("\nDone.\n")
