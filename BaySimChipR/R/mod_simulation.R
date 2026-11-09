# =============================================================================
# mod_simulation.R - tab 3 "Simulation"
#
# Simulation time, "Simulate" / "Clear plot" buttons, profiles of one or
# several compartments, exposure parameters (AUC, Cmax, half-life, Tmax).
# Runs can be overlaid (e.g. to compare conditions) and the results
# downloaded.
# =============================================================================

mod_simulation_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 310,
      title = "Simulation",
      shiny::numericInput(ns("sim_time"), "Simulation time (h)", value = 72, min = 0.1, step = 1),
      shiny::selectizeInput(ns("comps"), "Compartments to plot", choices = NULL, multiple = TRUE,
                            options = list(plugins = list("remove_button"))),
      shiny::radioButtons(ns("yscale"), "Concentration axis",
                          c("Linear" = "linear", "Logarithmic" = "log"), inline = TRUE),
      shiny::textInput(ns("run_label"), "Run label (optional)", placeholder = "e.g. with evaporation"),
      shiny::checkboxInput(ns("overlay"), "Keep previous runs (overlay)", value = TRUE),
      shiny::div(
        class = "d-flex gap-2",
        shiny::actionButton(ns("simulate"), "Simulate", class = "btn-primary", icon = shiny::icon("play")),
        shiny::actionButton(ns("clear"), "Clear plot", icon = shiny::icon("eraser"))
      ),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Exposure window",
          shiny::numericInput(ns("win_from"), "From (h)", value = NA, min = 0),
          shiny::numericInput(ns("win_to"), "To (h)", value = NA, min = 0),
          shiny::p(class = "small text-muted", "Empty = whole simulation.")
        ),
        bslib::accordion_panel(
          "ODE solver",
          shiny::selectInput(ns("method"), "Solver",
                             c("lsoda (default, as the paper's analysis)" = "lsoda",
                               "ode45 (Dormand-Prince)" = "ode45",
                               "bdf (stiff)" = "bdf")),
          shiny::numericInput(ns("rtol"), "Relative tolerance", value = 1e-6, min = 1e-12),
          shiny::numericInput(ns("atol"), "Absolute tolerance", value = 1e-6, min = 1e-14),
          shiny::numericInput(ns("dt"), "Output step (h)", value = 0.01, min = 0.001)
        )
      ),
      shiny::div(
        class = "d-flex flex-wrap gap-2",
        shiny::downloadButton(ns("dl_profiles"), "Profiles (CSV)", class = "btn-sm"),
        shiny::downloadButton(ns("dl_exposure"), "Exposure (CSV)", class = "btn-sm")
      )
    ),
    bslib::layout_columns(
      col_widths = c(8, 4),
      bslib::navset_card_tab(
        full_screen = TRUE,
        title = "Concentration-time profiles",
        bslib::nav_panel("Plot", plotly::plotlyOutput(ns("plot"), height = "430px")),
        bslib::nav_panel("Table", DT::DTOutput(ns("profile_table"), fill = FALSE))
      ),
      bslib::card(
        bslib::card_header("Model"),
        mod_diagram_ui(ns("diagram"), height = "290px"),
        shiny::p(class = "small text-muted mb-0", "Click a compartment to add it to the plot."),
        shiny::uiOutput(ns("notes"))
      )
    ),
    bslib::card(
      bslib::card_header("Exposure parameters"),
      DT::DTOutput(ns("exposure"), fill = FALSE)
    )
  )
}

mod_simulation_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    # clicking a compartment in this tab's diagram adds it to the plot
    mod_diagram_server("diagram", state, on_select = function(sel) {
      if (sel$kind == "compartment" && !sel$name %in% input$comps) {
        shiny::updateSelectizeInput(session, "comps", selected = c(input$comps, sel$name))
      }
    })

    # ---- simulation time is stored in the project ----------------------------------
    shiny::observeEvent(state$project_id, {
      shiny::updateNumericInput(session, "sim_time", value = state$project$settings$sim_time)
      run_counter <<- 0
    })
    shiny::observeEvent(input$sim_time, {
      v <- input$sim_time
      if (is.numeric(v) && length(v) && is.finite(v) && v > 0 &&
          !isTRUE(all.equal(v, state$project$settings$sim_time))) {
        p <- state$project
        p$settings$sim_time <- v
        state$project <- p
      }
    }, ignoreInit = TRUE)

    # ---- compartment choices ---------------------------------------------------------
    comp_info <- shiny::reactiveVal(NULL)
    shiny::observe({
      cp <- state$project$compartments[, c("name", "type")]
      if (!identical(cp, shiny::isolate(comp_info()))) comp_info(cp)
    })
    shiny::observeEvent(comp_info(), {
      cp <- comp_info()
      sel <- intersect(shiny::isolate(input$comps), cp$name)
      if (!length(sel)) sel <- cp$name[cp$type == "fluid"]
      shiny::updateSelectizeInput(session, "comps", choices = cp$name, selected = sel)
    })

    # ---- run / clear ---------------------------------------------------------------------
    run_counter <- 0
    shiny::observeEvent(input$simulate, {
      p <- state$project
      if (!(is.numeric(input$sim_time) && isTRUE(input$sim_time > 0))) {
        bsc_notify("The simulation time must be > 0.", "error")
        return()
      }
      sim <- tryCatch(
        bsc_simulate(p, final_time = input$sim_time + 0.1,
                     method = input$method %||% "lsoda", rtol = input$rtol %||% 1e-6,
                     atol = input$atol %||% 1e-6, dt = input$dt %||% 0.01),
        error = function(e) {
          bsc_notify(paste("Simulation not possible:", conditionMessage(e)), "error")
          NULL
        })
      shiny::req(sim)
      run_counter <<- run_counter + 1
      label <- trimws(input$run_label %||% "")
      if (!nzchar(label)) label <- paste("Run", run_counter)
      run <- list(label = label, sim = sim, colors = bsc_compartment_colors(p),
                  unknown = sum(bsc_parameter_table(p)$estimate), created = Sys.time())
      state$sim_runs <- if (isTRUE(input$overlay)) c(state$sim_runs, list(run)) else list(run)
    })
    shiny::observeEvent(input$clear, {
      state$sim_runs <- list()
      run_counter <<- 0
    })

    # ---- plot -------------------------------------------------------------------------------
    output$plot <- plotly::renderPlotly({
      runs <- state$sim_runs
      comps <- input$comps
      if (!length(runs)) return(bsc_plotly_message("Press 'Simulate' to compute the profiles."))
      if (!length(comps)) return(bsc_plotly_message("Select the compartments to plot."))
      ylog <- identical(input$yscale, "log")
      fig <- plotly::plot_ly()
      for (k in seq_along(runs)) {
        r <- runs[[k]]
        sim <- r$sim
        idx <- bsc_thin_index(sim$time, sim$segment)
        for (cc in intersect(comps, colnames(sim$conc))) {
          y <- sim$conc[idx, cc]
          if (ylog) y[y <= 0] <- NA
          nm <- if (length(runs) > 1) paste0(cc, " · ", r$label) else cc
          fig <- plotly::add_trace(
            fig, x = sim$time[idx], y = y, type = "scatter", mode = "lines", name = nm,
            line = list(color = r$colors[[cc]] %||% BSC_SERIES_COLORS[1], width = 2,
                        dash = BSC_RUN_DASHES[(k - 1) %% length(BSC_RUN_DASHES) + 1]),
            hovertemplate = paste0(nm, ": %{y:.4g} µM<extra></extra>"))
        }
      }
      bsc_plotly_layout(fig, ylog = ylog)
    })

    output$profile_table <- DT::renderDT({
      runs <- state$sim_runs
      shiny::req(length(runs), length(input$comps))
      tabs <- lapply(runs, function(r) {
        t <- bsc_profile_table(r$sim, input$comps, step = 1)
        if (length(runs) > 1) names(t)[-1] <- paste0(names(t)[-1], " · ", r$label)
        t
      })
      df <- Reduce(function(a, b) merge(a, b, by = "Time (h)", all = TRUE), tabs)
      bsc_dt(df, selection = "none", scroll_y = "360px", digits = 5)
    })

    # ---- exposure parameters ---------------------------------------------------------------
    exposure <- shiny::reactive({
      runs <- state$sim_runs
      comps <- input$comps
      shiny::req(length(runs), length(comps))
      win <- c(input$win_from, input$win_to)
      if (length(win) != 2 || any(!is.finite(win))) win <- NULL
      do.call(rbind, lapply(runs, function(r) {
        cc <- intersect(comps, colnames(r$sim$conc))
        if (!length(cc)) return(NULL)
        e <- bsc_exposure(r$sim$time, r$sim$conc[, cc, drop = FALSE], window = win)
        cbind(run = r$label, e, stringsAsFactors = FALSE)
      }))
    })

    output$exposure <- DT::renderDT({
      e <- exposure()
      shiny::req(e)
      note <- e$status
      note[is.na(e$t_half)] <- paste0(note[is.na(e$t_half)], "; Cmax/2 not reached")
      df <- data.frame(Run = e$run, Compartment = e$compartment,
                       `AUC (µM·h)` = e$AUC, `Cmax (µM)` = e$Cmax, `Tmax (h)` = e$Tmax,
                       `Half-life (h)` = e$t_half, Note = note, check.names = FALSE)
      bsc_dt(df, selection = "none", scroll_y = "220px", digits = 4)
    })

    output$notes <- shiny::renderUI({
      runs <- state$sim_runs
      if (!length(runs)) return(NULL)
      n_unknown <- runs[[length(runs)]]$unknown
      if (n_unknown > 0) {
        shiny::p(class = "small text-warning mb-0 mt-2",
                 shiny::icon("triangle-exclamation"),
                 sprintf("%d parameter(s) are marked unknown; the last run used their current values.",
                         n_unknown))
      }
    })

    # ---- downloads ----------------------------------------------------------------------------
    output$dl_profiles <- shiny::downloadHandler(
      filename = function() paste0("BaySimChip_profiles_", format(Sys.Date()), ".csv"),
      content = function(file) {
        runs <- state$sim_runs
        df <- do.call(rbind, lapply(runs, function(r) {
          cbind(run = r$label, bsc_simulation_long(r$sim), stringsAsFactors = FALSE)
        }))
        if (is.null(df)) df <- data.frame(run = character(0), time = numeric(0),
                                          compartment = character(0), conc = numeric(0))
        names(df)[names(df) == "time"] <- "time_h"
        names(df)[names(df) == "conc"] <- "conc_uM"
        utils::write.csv(df, file, row.names = FALSE)
      })
    output$dl_exposure <- shiny::downloadHandler(
      filename = function() paste0("BaySimChip_exposure_", format(Sys.Date()), ".csv"),
      content = function(file) {
        e <- tryCatch(exposure(), error = function(err) NULL)
        if (is.null(e)) e <- data.frame()
        utils::write.csv(e, file, row.names = FALSE)
      })
  })
}
