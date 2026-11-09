# Equivalence of the engine with the analysis code that produced the figures
# of the paper (reference outputs: tests/reference/legacy_reference.rds,
# produced once with the original analysis scripts).

skip_if_not(file.exists(bsc_ref_path), "legacy reference not available")
ref <- readRDS(bsc_ref_path)

test_that("evaporation scenario (Fig. 4A) is identical to the analysis code", {
  p <- legacy_cnbio()
  for (cond in names(ref$evaporation)) {
    r <- ref$evaporation[[cond]]
    pc <- bsc_set_param(p, "compartment", "Comp0", "kev", r$kev_apical)
    pc <- bsc_set_param(pc, "compartment", "Comp2", "kev", r$kev_liver)
    sim <- bsc_simulate(pc, final_time = 110)
    expect_equal(sim$time, r$t)
    # liver (Comp2) and apical (Comp0) profiles, and all compartments
    expect_lt(max_rel_diff(sim$conc[, "Comp2"], r$yc[, 3]), 1e-12)
    expect_lt(max_rel_diff(sim$conc[, "Comp0"], r$yc[, 1]), 1e-12)
    expect_lt(max(abs(sim$conc - r$yc)), 1e-12)
  }
})

test_that("sampling scenario (Fig. 5A) is identical to the analysis code", {
  p <- legacy_cnbio()
  p <- bsc_set_param(p, "compartment", "Comp0", "kev", 0)
  p <- bsc_set_param(p, "compartment", "Comp2", "kev", 0)
  for (cond in names(ref$sampling)) {
    r <- ref$sampling[[cond]]
    pc <- bsc_set_exchange(p, "Comp2", data.frame(time = c(0, 1, 2, 4, 8, 24, 48),
                                                  removed = r$samp_vol, added = 0))
    pc <- bsc_set_exchange(pc, "Comp0", data.frame(time = c(8, 24), removed = r$samp_vol, added = 0))
    sim <- bsc_simulate(pc, final_time = 100)
    expect_equal(sim$time, r$t)
    expect_lt(max(abs(sim$conc - r$yc)), 1e-12)
  }
})

test_that("NSB scenario (Fig. 6A): fluid profiles identical; solid state of the analysis code = fu_solid x C_solid", {
  p <- bsc_read_project(bsc_extdata("NSBpara2DAnalysis.rds"))
  p$settings$dose_unit <- "ng"
  p <- bsc_add_compartment(p, "Comp6", type = "solid", volume = 1, fu = 0, cl = 0, kev = 0)
  p <- bsc_add_perm(p, "Comp6", "Comp0", ps = 0, name = "P.S.3")
  p <- bsc_set_param(p, "compartment", "Comp2", "cl", 0.1)
  p <- bsc_set_param(p, "perm", "P.S.0", "ps", 0.01)
  vols <- c(Comp0 = 0.35, Comp1 = 0.8, Comp2 = 1.6, Comp3 = 0.2)
  for (nm in names(vols)) p <- bsc_set_param(p, "compartment", nm, "volume", vols[[nm]])
  p <- bsc_set_dosing(p, "Comp0", list(mode = "single", single = 0.35 * p$settings$molar_mass))
  for (r in ref$nsb) {
    pc <- bsc_set_param(p, "compartment", "Comp5", "fu", r$K_media_plastic)
    pc <- bsc_set_param(pc, "compartment", "Comp6", "fu", r$K_media_plastic)
    pc <- bsc_set_param(pc, "perm", "P.S.2", "ps", r$PS)
    pc <- bsc_set_param(pc, "perm", "P.S.3", "ps", r$PS / 5)
    sim <- bsc_simulate(pc, final_time = 100)
    fluid <- 1:4
    expect_lt(max_rel_diff(sim$conc[, fluid], r$yc[, fluid], floor = 1e-9), 1e-6)
    # ConstructODE() of the analysis code builds the fu matrix row-wise (FU[i, j] = fu_i)
    # while the engine and AnalysisODE() use FU[i, j] = fu_j.  For solid
    # compartments this only rescales the solid state: C_analysis = fu_s * C_new.
    solid <- 5:6
    scaled <- sweep(sim$conc[, solid, drop = FALSE], 2, rep(r$K_media_plastic, 2), `*`)
    expect_lt(max_rel_diff(scaled, r$yc[, solid], floor = 1e-9), 1e-6)
  }
})

test_that("fits reproduce the estimates of the original Fig. 3 analysis (V_liver = 1.6 mL)", {
  for (k in c("fit_sd", "fit_dd")) {
    rr <- ref[[k]]
    f <- bsc_fit(legacy_fit_project(rr), final_time = rr$final_time)
    expect_true(f$converged)
    expect_lt(max_rel_diff(f$estimates$estimate, unname(rr$estimates)), 1e-10)
    expect_lt(max_rel_diff(f$estimates$rse_percent, unname(rr$rse_percent)), 1e-8)
  }
})

test_that("fits reproduce the Fig. 3 values of the original analysis (V_liver = 1.8 mL)", {
  pub <- ref$published
  # single dose 1000 uM: exact to all reported digits
  rr <- ref$fit_sd_v18
  f <- bsc_fit(legacy_fit_project(rr))   # default final time (73 h) instead of 75 h
  cl_cells <- f$estimates$estimate[1] * 1000 / 60 / 0.6       # uL/min/1e6 cells (0.6e6 PHH)
  p_cm_s <- f$estimates$estimate[2] / (3600 * 0.33) * 1e6     # 1e-6 cm/s (0.33 cm2)
  row <- pub[pub$SheetName == "D_G-L_SD_PO_1000uM", ]
  expect_lt(abs(cl_cells / row$Clint - 1), 1e-6)
  expect_lt(abs(p_cm_s / row$P - 1), 1e-6)
  expect_lt(abs(f$estimates$rse_percent[1] / row$RSECL - 1), 1e-6)
  expect_lt(abs(f$estimates$rse_percent[2] / row$RSEPS - 1), 1e-6)
  # daily dosing 100 uM: agrees to < 0.03 % (the stored values come from an earlier version of the original script)
  rr <- ref$fit_dd_v18
  f <- bsc_fit(legacy_fit_project(rr), final_time = rr$final_time)
  row <- pub[pub$SheetName == "D_G-L_DD_PO_100uM", ]
  expect_lt(abs(f$estimates$estimate[1] * 1000 / 60 / 0.6 / row$Clint - 1), 3e-4)
  expect_lt(abs(f$estimates$estimate[2] / (3600 * 0.33) * 1e6 / row$P - 1), 3e-4)
})

test_that("exact sampling-time lookup gives the same estimates as the grid lookup of the analysis code", {
  for (k in c("fit_sd", "fit_dd")) {
    rr <- ref[[k]]
    f <- bsc_fit(legacy_fit_project(rr), obs_mode = "exact", final_time = rr$final_time)
    expect_lt(max_rel_diff(f$estimates$estimate, unname(rr$estimates)), 1e-3)
  }
})

test_that("batch API: the evaporation example project reproduces the Fig. 4A simulation of the analysis code", {
  p <- bsc_load_example("CNBio_evaporation")
  sc <- data.frame(scenario = names(ref$evaporation),
                   `Apical|kev` = vapply(ref$evaporation, function(r) r$kev_apical, 0),
                   `Liver|kev` = vapply(ref$evaporation, function(r) r$kev_liver, 0),
                   check.names = FALSE)
  res <- bsc_run_scenarios(p, sc, final_time = 110)
  expect_setequal(unique(res$scenario), names(ref$evaporation))
  for (cond in names(ref$evaporation)) {
    sim <- attr(res, "simulations")[[cond]]
    expect_lt(max(abs(sim$conc - ref$evaporation[[cond]]$yc)), 1e-12)
  }
  expo <- attr(res, "exposure")
  expect_equal(nrow(expo), 8)
  expect_true(all(c("scenario", "compartment", "AUC", "Cmax", "Tmax", "t_half") %in% names(expo)))
})
