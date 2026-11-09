# Server-side tests of the Shiny modules with shiny::testServer().

new_test_state <- function(example) {
  p <- bsc_load_example(example)
  state <- bsc_new_state(p)
  shiny::isolate(bsc_set_project(state, p))
  state
}
file_input <- function(path, name = basename(path)) {
  data.frame(name = name, size = file.size(path), type = "", datapath = path,
             stringsAsFactors = FALSE)
}

test_that("diagram module: selection, dragging and editor requests (graphChange)", {
  state <- new_test_state("CNBio_gut_liver")
  testServer(mod_diagram_server, args = list(state = state, editable = TRUE), {
    session$flushReact()
    expect_equal(nrow(net_data()$nodes), 4)
    expect_equal(nrow(net_data()$edges), 4)
    expect_error(output$net, NA)
    session$setInputs(net_select = list(id = "compartment|Liver", nonce = 1))
    expect_equal(state$selected, list(kind = "compartment", name = "Liver"))
    session$setInputs(net_select = list(id = "perm|P.S.0", nonce = 2))
    expect_equal(state$selected, list(kind = "perm", name = "P.S.0"))
    session$setInputs(net_drag = list(pos = list(list(id = "compartment|Liver", x = 11, y = 22)), nonce = 3))
    i <- match("Liver", state$project$compartments$name)
    expect_equal(c(state$project$compartments$x[i], state$project$compartments$y[i]), c(11, 22))
    # editor: add a fluid compartment by clicking at (300, -40)
    session$setInputs(net_graphChange = list(cmd = "addNode", x = 300, y = -40, type = "fluid", nonce = 4))
    expect_equal(nrow(state$project$compartments), 5)
    new <- state$selected$name
    expect_equal(state$project$compartments$x[state$project$compartments$name == new], 300)
    expect_equal(nrow(net_data()$nodes), 5)
    # draw a flow and a permeable surface
    session$setInputs(net_graphChange = list(cmd = "addEdge", from = paste0("compartment|", new),
                                             to = "compartment|Liver", type = "flow", nonce = 5))
    expect_equal(nrow(state$project$flows), 4)
    expect_equal(state$selected$kind, "flow")
    session$setInputs(net_graphChange = list(cmd = "addEdge", from = "compartment|Apical",
                                             to = paste0("compartment|", new), type = "perm", nonce = 6))
    expect_equal(nrow(state$project$perms), 2)
    # a refused request (flow into a solid compartment) leaves the project unchanged
    session$setInputs(net_graphChange = list(cmd = "addNode", x = 0, y = 300, type = "solid", nonce = 7))
    solid <- state$selected$name
    before <- state$project
    session$setInputs(net_graphChange = list(cmd = "addEdge", from = "compartment|Liver",
                                             to = paste0("compartment|", solid), type = "flow", nonce = 8))
    expect_identical(state$project, before)
    # delete the new compartment (its flow and permeable surface go with it)
    session$setInputs(net_graphChange = list(cmd = "deleteElements",
                                             nodes = list(paste0("compartment|", new)),
                                             edges = list()))
    expect_false(new %in% state$project$compartments$name)
    expect_equal(nrow(state$project$flows), 3)
    expect_equal(nrow(state$project$perms), 1)
    expect_equal(nrow(net_data()$edges), 4)
  })
})

test_that("model structure module: editor, properties panel, tables, files", {
  state <- new_test_state("CNBio_gut_liver")
  testServer(mod_model_structure_server, args = list(state = state), {
    session$flushReact()
    expect_error(output$comp_table, NA)
    expect_error(output$flow_table, NA)
    expect_error(output$perm_table, NA)
    # drag-and-drop editor (nested diagram module)
    session$setInputs(`diagram-net_graphChange` = list(cmd = "addNode", x = 50, y = 250,
                                                       type = "solid", nonce = 1))
    expect_equal(sum(state$project$compartments$type == "solid"), 1)
    solid <- state$selected$name
    # properties panel (nested element editor): rename, volume, effective fu
    session$setInputs(`editor-element` = paste0("compartment|", solid))
    session$setInputs(`editor-name` = "Plastic", `editor-rename` = 1)
    expect_true("Plastic" %in% state$project$compartments$name)
    expect_equal(state$selected$name, "Plastic")
    session$setInputs(`editor-element` = "compartment|Plastic")
    session$setInputs(`editor-vol` = 0.5, `editor-fu` = 0.25, `editor-fu_unknown` = TRUE)
    r <- state$project$compartments[state$project$compartments$name == "Plastic", ]
    expect_equal(c(r$volume, r$fu), c(0.5, 0.25))
    expect_true(r$fu_est)
    session$setInputs(`diagram-net_graphChange` = list(cmd = "addEdge", from = "compartment|Plastic",
                                                       to = "compartment|Liver", type = "perm",
                                                       nonce = 2))
    expect_equal(nrow(state$project$perms), 2)
    session$setInputs(`editor-element` = paste0("perm|", state$project$perms$name[2]))
    session$setInputs(`editor-ps` = 0.07)
    expect_equal(state$project$perms$ps[2], 0.07)
    # type change: a solid compartment with flows is refused
    session$setInputs(`editor-element` = "compartment|Liver", `editor-type` = "solid")
    expect_equal(state$project$compartments$type[state$project$compartments$name == "Liver"], "fluid")
    # tables (secondary view): add, rename, edit Q, remove with confirmation
    session$setInputs(new_comp_name = "Gut2", new_comp_type = "fluid", add_comp = 1)
    expect_true("Gut2" %in% state$project$compartments$name)
    session$setInputs(new_flow_from = "Gut2", new_flow_to = "Liver", new_flow_q = 2, add_flow = 1)
    expect_equal(nrow(state$project$flows), 4)
    i <- match("Gut2", state$project$compartments$name)
    session$setInputs(comp_table_cell_edit = data.frame(row = i, col = 0, value = "Gut_2"))
    expect_true("Gut_2" %in% state$project$flows$from)
    session$setInputs(comp_table_cell_edit = data.frame(row = i, col = 2, value = "0.4"))
    expect_equal(state$project$compartments$volume[i], 0.4)
    session$setInputs(flow_table_cell_edit = data.frame(row = 4, col = 3, value = "2.5"))
    expect_equal(state$project$flows$q[4], 2.5)
    j <- match("Plastic", state$project$compartments$name)
    session$setInputs(comp_table_rows_selected = j, remove_comp = 1)
    session$setInputs(confirm_remove = 1)
    expect_false("Plastic" %in% state$project$compartments$name)
    expect_equal(nrow(state$project$perms), 1)
    # save (.rds, native) and JSON export, then re-load both
    session$setInputs(name = "My chip")
    f <- output$save
    expect_match(f, "\\.rds$")
    expect_equal(bsc_read_project(f)$name, "My chip")
    expect_true("Gut_2" %in% bsc_read_project(f)$compartments$name)
    fj <- output$export_json
    expect_match(fj, "\\.json$")
    expect_equal(bsc_read_project(fj)$compartments, bsc_read_project(f)$compartments)
    session$setInputs(example = "CNBio_NSB", load_example = 1)
    expect_equal(sum(state$project$compartments$type == "solid"), 2)
    session$setInputs(load = file_input(f))
    expect_equal(state$project$name, "My chip")
    session$setInputs(load = file_input(fj))
    expect_true("Gut_2" %in% state$project$compartments$name)
    # other file types are refused
    before <- state$project
    other <- tempfile(fileext = ".xyz")
    writeBin(as.raw(0:255), other)
    session$setInputs(load = file_input(other))
    expect_identical(state$project, before)
    session$setInputs(new = 1)
    expect_equal(nrow(state$project$compartments), 0)
  })
})

test_that("a late browser echo of the previous element's values is not written into the next element", {
  # Regression test: the server pushes Liver's values to the form; before the browser's
  # echo of those values arrives, the user selects Apical. The late echo (Liver's volume,
  # Kev and 'unknown' flag) must not be copied into Apical.
  state <- new_test_state("CNBio_evaporation")
  testServer(mod_conditions_server, args = list(state = state), {
    session$flushReact()
    p0 <- state$project
    session$setInputs(`editor-element` = "compartment|Liver")   # server pushes Liver's values
    session$setInputs(`editor-element` = "compartment|Apical")  # ... and then Apical's
    # late echoes of the values pushed for Liver arrive now
    session$setInputs(`editor-vol` = 1.6, kev = 0.0055)
    ap <- state$project$compartments[state$project$compartments$name == "Apical", ]
    expect_equal(c(ap$volume, ap$kev), c(0.35, 0.003))
    # a genuine user edit afterwards still goes through
    session$setInputs(`editor-vol` = 0.3)
    ap <- state$project$compartments[state$project$compartments$name == "Apical", ]
    expect_equal(ap$volume, 0.3)
    li <- state$project$compartments[state$project$compartments$name == "Liver", ]
    expect_equal(c(li$volume, li$kev), c(p0$compartments$volume[3], p0$compartments$kev[3]))
  })
})

test_that("experimental conditions module: parameters, unknown flags, dosing, exchange", {
  state <- new_test_state("CNBio_gut_liver")
  testServer(mod_conditions_server, args = list(state = state), {
    session$flushReact()
    session$setInputs(`editor-element` = "compartment|Liver")
    expect_equal(state$selected, list(kind = "compartment", name = "Liver"))
    expect_equal(output$element_kind, "fluid")
    session$setInputs(`editor-vol` = 1.8, `editor-cl` = 0.2, `editor-fu` = 0.9, kev = 0.0055)
    r <- state$project$compartments[state$project$compartments$name == "Liver", ]
    expect_equal(c(r$volume, r$cl, r$fu, r$kev), c(1.8, 0.2, 0.9, 0.0055))
    session$setInputs(`editor-cl_unknown` = TRUE, kev_unknown = TRUE)
    expect_true(bsc_get_param(state$project, "compartment", "Liver", "cl")$estimate)
    expect_true(bsc_get_param(state$project, "compartment", "Liver", "kev")$estimate)
    # dosing
    session$setInputs(dose_mode = "repeated", dose_rep = 0.02, dose_interval = 12)
    d <- bsc_get_dosing(state$project, "Liver")
    expect_equal(c(d$mode, d$repeated_dose, d$repeated_interval), c("repeated", "0.02", "12"))
    session$setInputs(dose_mode = "custom", dose_add = 1)
    session$setInputs(dose_add = 2)
    expect_equal(nrow(bsc_get_dosing(state$project, "Liver")$custom), 3)
    session$setInputs(dose_table_cell_edit = data.frame(row = 3, col = 1, value = "0.5"))
    expect_equal(bsc_get_dosing(state$project, "Liver")$custom$dose[3], 0.5)
    session$setInputs(first_dose = TRUE)
    expect_true(bsc_get_dosing(state$project, "Liver")$first_dose_from_data)
    session$setInputs(dose_unit = "nmol", molar_mass = 296.15)
    expect_equal(state$project$settings$dose_unit, "nmol")
    expect_equal(state$project$settings$molar_mass, 296.15)
    # media exchange table
    n0 <- nrow(bsc_get_exchange(state$project, "Liver"))
    session$setInputs(ex_add = 1)
    expect_equal(nrow(bsc_get_exchange(state$project, "Liver")), n0 + 1)
    session$setInputs(exchange_table_cell_edit = data.frame(row = n0 + 1, col = 1, value = "0.867"))
    expect_equal(bsc_get_exchange(state$project, "Liver")$removed[n0 + 1], 0.867)
    session$setInputs(exchange_table_rows_selected = n0 + 1, ex_remove = 1)
    expect_equal(nrow(bsc_get_exchange(state$project, "Liver")), n0)
    expect_error(output$exchange_table, NA)
    expect_error(output$dose_table, NA)
    # flows and permeable surfaces
    session$setInputs(`editor-element` = "flow|Flow0")
    expect_equal(output$element_kind, "flow")
    session$setInputs(`editor-q` = 5)
    expect_equal(state$project$flows$q[1], 5)
    session$setInputs(`editor-element` = "perm|P.S.0", `editor-ps_unknown` = TRUE)
    expect_true(state$project$perms$ps_est[1])
    # overview table edits
    pt <- bsc_parameter_table(state$project)
    i <- match("compartment|Apical|kev", pt$key)
    session$setInputs(param_table_cell_edit = data.frame(row = i, col = 2, value = "0.004"))
    expect_equal(bsc_get_param(state$project, "compartment", "Apical", "kev")$value, 0.004)
    session$setInputs(param_table_cell_edit = data.frame(row = i, col = 4, value = "yes"))
    expect_true(bsc_get_param(state$project, "compartment", "Apical", "kev")$estimate)
    expect_error(output$param_table, NA)
  })
})

test_that("simulation module: simulate, overlay runs, exposure, downloads", {
  state <- new_test_state("CNBio_evaporation")
  testServer(mod_simulation_server, args = list(state = state), {
    session$flushReact()
    session$setInputs(sim_time = 110, comps = c("Apical", "Liver"), yscale = "log", run_label = "",
                      overlay = TRUE, method = "lsoda", rtol = 1e-6, atol = 1e-6, dt = 0.01,
                      win_from = NA, win_to = NA)
    session$setInputs(simulate = 1)
    expect_length(state$sim_runs, 1)
    expect_equal(state$sim_runs[[1]]$label, "Run 1")
    expect_equal(max(state$sim_runs[[1]]$sim$time), 110.09)
    state$project <- bsc_set_param(state$project, "compartment", "Liver", "kev", 0)
    session$setInputs(run_label = "no liver evaporation", simulate = 2)
    expect_length(state$sim_runs, 2)
    e <- exposure()
    expect_equal(nrow(e), 4)
    expect_true(all(e$AUC > 0))
    expect_error(output$plot, NA)
    expect_error(output$exposure, NA)
    expect_error(output$profile_table, NA)
    d <- utils::read.csv(output$dl_profiles)
    expect_equal(sort(unique(d$run)), sort(c("Run 1", "no liver evaporation")))
    expect_true(nrow(utils::read.csv(output$dl_exposure)) == 4)
    # clicking a compartment in this tab's diagram adds it to the plot selection
    session$setInputs(`diagram-net_select` = list(id = "compartment|Basolateral", nonce = 1))
    expect_equal(state$selected$name, "Basolateral")
    session$setInputs(overlay = FALSE, simulate = 3)
    expect_length(state$sim_runs, 1)
    session$setInputs(clear = 1)
    expect_length(state$sim_runs, 0)
    # invalid model -> no run, no crash
    state$project <- bsc_set_param(state$project, "compartment", "Liver", "volume", NA)
    session$setInputs(simulate = 4)
    expect_length(state$sim_runs, 0)
  })
})

test_that("fitting module: fit, apply estimates, data import, unknown selection", {
  state <- new_test_state("DCF_gut_liver_SD_1000uM")
  testServer(mod_fitting_server, args = list(state = state), {
    session$flushReact()
    expect_error(output$fit_plot, NA)
    session$setInputs(obs_mode = "grid", maxiter = 50, method = "lsoda", yscale = "linear")
    session$setInputs(fit = 1)
    expect_false(is.null(state$fit))
    e <- state$fit$result$estimates
    # Fig. 3 values for this arm
    expect_lt(abs(e$estimate[1] * 1000 / 60 / 0.6 / 2.093318 - 1), 1e-6)
    expect_lt(abs(e$estimate[2] / (3600 * 0.33) * 1e6 / 16.04920 - 1), 1e-6)
    expect_error(output$estimates, NA)
    expect_error(output$fit_summary, NA)
    expect_error(output$resid_plot, NA)
    expect_error(output$exposure, NA)
    expect_true(nrow(utils::read.csv(output$dl_estimates)) == 2)
    # the fit is stored in the native project file and restored on load
    saved <- bsc_project_for_save(state)
    expect_equal(saved$results$fit$result$estimates$estimate, e$estimate)
    st2 <- bsc_new_state()
    shiny::isolate(bsc_set_project(st2, saved))
    expect_equal(shiny::isolate(st2$fit$result$estimates$estimate), e$estimate)
    session$setInputs(apply = 1)
    expect_equal(bsc_get_param(state$project, "compartment", "Liver", "cl")$value, e$estimate[1])
    # import the raw sheet (Location column: Liver / Apical) and map it
    csv <- file.path(bsc_root, "inst", "extdata", "DCF_gut_liver_SD_1000uM.csv")
    session$setInputs(data_file = file_input(csv))
    expect_false(is.null(raw()))
    session$setInputs(col_time = "Time_h", col_conc = "Conc", col_comp = "Location", blq = "zero",
                      import_mode = "replace")
    expect_equal(map_values(), c("Liver", "Apical"))
    session$setInputs(map_1 = "Liver", map_2 = "Apical", do_import = 1)
    obs <- state$project$observed
    expect_equal(sum(obs$compartment == "Liver"), 18)
    expect_equal(sum(obs$compartment == "Apical"), 2)
    # select unknowns through the selectize input; edit bounds
    session$setInputs(unknowns = "compartment|Liver|cl")
    expect_equal(bsc_parameter_table(state$project)$key[bsc_parameter_table(state$project)$estimate],
                 "compartment|Liver|cl")
    session$setInputs(unknown_table_cell_edit = data.frame(row = 1, col = 3, value = "5"))
    expect_equal(state$project$fit$bounds$upper, 5)
    # applying an outdated fit is refused
    session$setInputs(apply = 2)
    expect_equal(bsc_get_param(state$project, "compartment", "Liver", "cl")$value, e$estimate[1])
    # observed table editing
    session$setInputs(obs_add = 1)
    expect_equal(nrow(state$project$observed), 21)
    session$setInputs(obs_table_rows_selected = 21, obs_remove = 1)
    expect_equal(nrow(state$project$observed), 20)
    session$setInputs(obs_clear = 1)
    expect_equal(nrow(state$project$observed), 0)
    # fitting without data is refused gracefully
    session$setInputs(fit = 2)
    expect_equal(state$fit$keys, c("compartment|Liver|cl", "perm|P.S.0|ps"))
  })
})

test_that("the whole app server runs with testServer", {
  app <- shiny::shinyAppDir(bsc_root)
  testServer(app, {
    session$flushReact()
    expect_equal(state$project$name, "Gut-liver chip (CN Bio TL6)")
    session$setInputs(`simulation-simulate` = 1, `simulation-sim_time` = 72,
                      `simulation-comps` = "Liver")
    expect_length(state$sim_runs, 1)
    session$setInputs(`structure-diagram-net_graphChange` = list(cmd = "addNode", x = 0, y = 0,
                                                                 type = "fluid", nonce = 1))
    expect_equal(nrow(state$project$compartments), 5)
    expect_error(output$project_badge, NA)
  })
})
