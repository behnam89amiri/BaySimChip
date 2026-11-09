# Basic numerical and data-structure tests of the engine (analytical
# solutions, mass balance, events, units, import/export, validation, fitting).

one_comp <- function(volume = 1, cl = 0, fu = 1, kev = 0, dose = 10) {
  p <- bsc_new_project("one compartment")
  p$settings$dose_unit <- "nmol"
  p <- bsc_add_compartment(p, "C", volume = volume, cl = cl, fu = fu, kev = kev)
  bsc_set_dosing(p, "C", list(mode = "single", single = dose))
}

test_that("first-order clearance matches the analytical solution and exposure parameters", {
  p <- one_comp(volume = 2, cl = 0.3, fu = 0.5, dose = 20)
  sim <- bsc_simulate(p, final_time = 48)
  k <- 0.3 * 0.5 / 2
  expect_lt(max_rel_diff(sim$conc[, "C"], 10 * exp(-k * sim$time)), 1e-5)
  e <- bsc_exposure(sim$time, sim$conc)
  tend <- max(sim$time)
  expect_lt(abs(e$AUC / (10 / k * (1 - exp(-k * tend))) - 1), 1e-4)
  expect_equal(e$Cmax, 10)
  expect_equal(e$Tmax, 0)
  expect_lt(abs(e$t_half - log(2) / k), 0.01)
  expect_equal(e$status, "Tmax reached and passed")
})

test_that("evaporation concentrates the compound: C(t) = A0 / (V0 - Kev t)", {
  p <- one_comp(volume = 1, kev = 0.005, dose = 10)
  sim <- bsc_simulate(p, final_time = 100)
  expect_lt(max_rel_diff(sim$conc[, "C"], 10 / (1 - 0.005 * sim$time)), 1e-6)
  expect_lt(max_rel_diff(bsc_total_amount(sim), rep(10, length(sim$time))), 1e-6)
  e <- bsc_exposure(sim$time, sim$conc)
  expect_equal(e$status, "Tmax not reached")
  expect_true(is.na(e$t_half))
})

test_that("closed systems conserve mass, incl. solid compartments with fu != 1", {
  p <- bsc_new_project("closed")
  p$settings$dose_unit <- "nmol"
  p <- bsc_add_compartment(p, "A", volume = 0.5)
  p <- bsc_add_compartment(p, "B", volume = 1.5)
  p <- bsc_add_compartment(p, "S", type = "solid", volume = 1, fu = 0.2)
  p <- bsc_add_flow(p, "A", "B", q = 2)
  p <- bsc_add_flow(p, "B", "A", q = 2)
  p <- bsc_add_perm(p, "A", "S", ps = 0.5)
  p <- bsc_set_dosing(p, "A", list(single = 5))
  sim <- bsc_simulate(p, final_time = 200)
  tot <- bsc_total_amount(sim)
  expect_lt(max(abs(tot / 5 - 1)), 1e-5)
  # equilibrium of the permeable surface: fu_A C_A = fu_S C_S
  n <- length(sim$time)
  expect_lt(abs(sim$conc[n, "S"] * 0.2 / sim$conc[n, "A"] - 1), 1e-4)
})

test_that("media exchange, sampling and doses update the state as in the original code", {
  p <- one_comp(volume = 1, dose = 10)
  p <- bsc_set_exchange(p, "C", data.frame(time = 10, removed = 0.5, added = 0.25))
  p <- bsc_set_dosing(p, "C", list(mode = "custom", custom = data.frame(time = c(0, 10), dose = c(10, 5))))
  sim <- bsc_simulate(p, final_time = 20)
  # constant 10 uM before t = 10; then (10 * 0.5 + 5) / 0.75
  expect_equal(unname(sim$conc[abs(sim$time - 9.99) < 1e-9, "C"]), 10)
  expect_lt(abs(sim$conc[sim$time == 10, "C"] - (10 * 0.5 + 5) / 0.75), 1e-9)
  expect_lt(abs(sim$volume[sim$time == 10, "C"] - 0.75), 1e-12)
  expect_lt(abs(sim$removed_amount[1, 2] - 5), 1e-9)
})

test_that("intervention schedule: merging, repeated doses, final time rules", {
  p <- one_comp(dose = 1)
  p <- bsc_set_dosing(p, "C", list(mode = "repeated", repeated_dose = 2, repeated_interval = 24))
  p <- bsc_set_exchange(p, "C", data.frame(time = c(24, 48, 100), removed = 0.1, added = 0.1))
  ev <- bsc_build_events(p, final_time = 72)
  # doses at 0, 24, 48, 72; exchanges at 24, 48 merged; exchange at 100 dropped;
  # final time equal to the last event -> moved 1 h later (original rule)
  expect_equal(ev$times, c(0, 24, 48, 72, 73))
  expect_equal(unname(ev$doses[1, ]), c(2, 2, 2, 2, 0))
  expect_equal(unname(ev$removed[1, ]), c(0, 0.1, 0.1, 0, 0))
  # custom doses not starting at 0 get an empty event at t = 0
  p2 <- bsc_set_dosing(one_comp(), "C", list(mode = "custom", custom = data.frame(time = 5, dose = 1)))
  ev2 <- bsc_build_events(p2, final_time = 10)
  expect_equal(ev2$times, c(0, 5, 10))
})

test_that("dose units are converted to nmol (mg, ug, ng as in the analysis code, nmol)", {
  expect_equal(bsc_dose_to_nmol(0.01, "mg", 300), 1e6 * 0.01 / 300)
  expect_equal(bsc_dose_to_nmol(105, "ng", 300), 105 / 300)
  expect_equal(bsc_dose_to_nmol(1, "ug", 250), 4)
  expect_equal(bsc_dose_to_nmol(3, "nmol", NA), 3)
  expect_equal(bsc_nmol_to_dose(bsc_dose_to_nmol(0.7, "mg", 296.15), "mg", 296.15), 0.7)
})

test_that("the app reads only .rds / .json", {
  expect_error(bsc_read_project("model.xyz"), "Unsupported file type")
})

test_that("the CNBio models in inst/extdata are read correctly", {
  p <- bsc_read_project(bsc_extdata("Kev2DAnalysis.rds"))
  expect_equal(p$compartments$name, paste0("Comp", 0:3))
  expect_equal(p$compartments$volume, c(0.3, 0.8, 1.6, 0.2))
  expect_equal(p$compartments$cl, c(0, 0, 0.08, 0))
  expect_equal(p$compartments$kev, c(0.003, 0, 0.01, 0))
  expect_equal(p$flows$q, rep(3.6, 3))
  expect_equal(paste(p$flows$from, p$flows$to), c("Comp1 Comp2", "Comp2 Comp3", "Comp3 Comp1"))
  expect_equal(p$perms$ps, 0.03)
  expect_equal(p$settings$molar_mass, 300)
  expect_equal(bsc_get_dosing(p, "Comp0")$single, 0.01)
  expect_false(any(bsc_parameter_table(p)$estimate))
  nsb <- bsc_read_project(bsc_extdata("NSBpara2DAnalysis.rds"))
  expect_equal(nsb$compartments$type, c(rep("fluid", 4), "solid"))
  expect_equal(nrow(bsc_get_exchange(nsb, "Comp2")), 7)
  smp <- bsc_read_project(bsc_extdata("SamplingVol2DAnalysis.rds"))
  expect_equal(bsc_get_exchange(smp, "Comp2")$time, c(0, 1, 2, 4, 8, 24, 48, 72))
  expect_true(bsc_get_dosing(smp, "Comp0")$first_dose_from_data)
})

test_that("projects survive JSON and RDS round trips unchanged", {
  for (ex in bsc_list_examples()) {
    p <- bsc_load_example(ex)
    f1 <- tempfile(fileext = ".json")
    f2 <- tempfile(fileext = ".rds")
    bsc_write_project(p, f1)
    bsc_write_project(p, f2)
    expect_identical(bsc_read_project(f1), p)
    expect_identical(bsc_read_project(f2), p)
  }
})

test_that("structure editing keeps references consistent", {
  p <- bsc_load_example("CNBio_gut_liver")
  p$observed <- data.frame(compartment = "Liver", time = 1, conc = 2)
  p <- bsc_rename(p, "compartment", "Liver", "Hepatic")
  expect_true("Hepatic" %in% c(p$flows$from, p$flows$to))
  expect_equal(p$observed$compartment, "Hepatic")
  expect_true("Hepatic" %in% names(p$dosing))
  p <- bsc_remove_compartments(p, "Hepatic")
  expect_equal(nrow(p$flows), 1)
  expect_equal(nrow(p$observed), 0)
  expect_false("Hepatic" %in% names(p$exchange))
  expect_error(bsc_add_flow(bsc_add_compartment(p, "S", type = "solid"), "S", "Apical"), "fluid")
  expect_error(bsc_add_compartment(p, "Apical"), "already exists")
})

test_that("validation reports missing or invalid values", {
  p <- bsc_load_example("CNBio_gut_liver")
  p <- bsc_set_param(p, "compartment", "Apical", "volume", -1)
  expect_true(any(grepl("Volume in Apical", bsc_validate(p)$errors)))
  expect_error(bsc_simulate(p), "Volume in Apical")
  p <- bsc_load_example("CNBio_gut_liver")
  expect_true(any(grepl("No parameter is marked", bsc_validate(p, "fit")$errors)))
  p <- bsc_set_param(p, "compartment", "Liver", "kev", 0.05)
  expect_error(bsc_simulate(p), "drops to zero")
})

test_that("fitting recovers known parameters from synthetic data (CL, PS, effective fu, NSB PS)", {
  truth <- bsc_load_example("CNBio_NSB")
  truth <- bsc_set_param(truth, "compartment", "NSB_liver", "fu", 0.4)
  truth <- bsc_set_param(truth, "perm", "P.S.2", "ps", 0.08)
  times <- c(0.5, 1, 2, 4, 8, 24, 48, 72, 96)
  sim <- bsc_simulate(truth, final_time = 100)
  idx <- match(times, round(sim$time, 6))
  truth$observed <- rbind(
    data.frame(compartment = "Liver", time = times, conc = sim$conc[idx, "Liver"]),
    data.frame(compartment = "Apical", time = times, conc = sim$conc[idx, "Apical"]))
  p <- truth
  p <- bsc_set_param(p, "compartment", "Liver", "cl", 0.05, estimate = TRUE)
  p <- bsc_set_param(p, "perm", "P.S.0", "ps", 0.02, estimate = TRUE)
  p <- bsc_set_param(p, "compartment", "NSB_liver", "fu", 1, estimate = TRUE)
  p <- bsc_set_param(p, "perm", "P.S.2", "ps", 0.05, estimate = TRUE)
  f <- bsc_fit(p, obs_mode = "exact", final_time = 100, maxiter = 200)
  expect_true(f$converged)
  expect_equal(f$estimates$key, c("compartment|Liver|cl", "compartment|NSB_liver|fu",
                                  "perm|P.S.0|ps", "perm|P.S.2|ps"))
  expect_lt(max_rel_diff(f$estimates$estimate, c(0.1, 0.4, 0.01, 0.08)), 1e-3)
  expect_true(all(f$estimates$ci_lower <= f$estimates$estimate & f$estimates$estimate <= f$estimates$ci_upper))
  p2 <- bsc_apply_estimates(p, f)
  expect_lt(abs(bsc_get_param(p2, "compartment", "Liver", "cl")$value / 0.1 - 1), 1e-3)
  expect_true(bsc_get_param(p2, "compartment", "Liver", "cl")$estimate)
})

test_that("observed data import handles compartment columns and BLQ entries", {
  df <- data.frame(Location = c("Liver", "Liver", "Apical", "Basolateral"),
                   Time_h = c("0", "1", "1", "2"), Conc = c("< 0.05", "2.5", "3", "x"),
                   stringsAsFactors = FALSE)
  g <- bsc_guess_columns(df)
  expect_equal(unlist(g), c(time = "Time_h", conc = "Conc", compartment = "Location"))
  o <- bsc_observed_from_table(df, "Time_h", "Conc", comp_col = "Location",
                               comp_map = c(Liver = "Liver", Apical = ""))
  expect_equal(o$compartment, c("Liver", "Liver"))
  expect_equal(o$conc, c(0, 2.5))
  o2 <- bsc_observed_from_table(df, "Time_h", "Conc", compartment = "Liver", blq = "half")
  expect_equal(o2$conc, c(0.025, 2.5, 3))
})

test_that("editor requests (graphChange) add and delete model elements", {
  p <- bsc_load_example("CNBio_gut_liver")
  r <- bsc_apply_graph_change(p, list(cmd = "addNode", x = 12.4, y = -30, type = "solid"))
  new <- r$selected$name
  expect_equal(r$project$compartments$type[r$project$compartments$name == new], "solid")
  expect_equal(unlist(r$project$compartments[r$project$compartments$name == new, c("x", "y")]),
               c(x = 12, y = -30))
  r2 <- bsc_apply_graph_change(r$project, list(cmd = "addEdge", from = paste0("compartment|", new),
                                               to = "compartment|Liver", type = "perm"))
  expect_equal(r2$selected$kind, "perm")
  expect_equal(nrow(r2$project$perms), 2)
  expect_error(bsc_apply_graph_change(r2$project, list(cmd = "addEdge", type = "flow",
                                                       from = "compartment|Liver",
                                                       to = paste0("compartment|", new))), "fluid")
  r3 <- bsc_apply_graph_change(r2$project, list(cmd = "addEdge", type = "flow",
                                                from = "compartment|Liver", to = "compartment|Apical"))
  expect_equal(r3$selected, list(kind = "flow", name = "Flow3"))
  r4 <- bsc_apply_graph_change(r3$project, list(cmd = "deleteElements",
                                                nodes = list(paste0("compartment|", new)),
                                                edges = list("flow|Flow3")))
  expect_false(new %in% r4$project$compartments$name)
  expect_equal(nrow(r4$project$perms), 1)
  expect_equal(nrow(r4$project$flows), 3)
})

test_that("bsc_set_params accepts full and short parameter keys", {
  p <- bsc_load_example("CNBio_gut_liver")
  p <- bsc_set_params(p, list("Liver|cl" = 0.2, "perm|P.S.0|ps" = 0.05, "Flow0|q" = 2))
  expect_equal(bsc_get_param(p, "compartment", "Liver", "cl")$value, 0.2)
  expect_equal(bsc_get_param(p, "perm", "P.S.0", "ps")$value, 0.05)
  expect_equal(bsc_get_param(p, "flow", "Flow0", "q")$value, 2)
  p <- bsc_set_params(p, c("Liver|cl" = 0.3), estimate = TRUE)
  expect_true(bsc_get_param(p, "compartment", "Liver", "cl")$estimate)
  expect_error(bsc_set_params(p, list("Nowhere|cl" = 1)), "Unknown parameter")
})
