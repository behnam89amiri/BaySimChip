# =============================================================================
# mod_fitting.R - tab 4 "Fitting and analysis"
#
# Experimental c-t data are loaded for one or several compartments (CSV or
# Excel, optionally with a compartment column), the unknown parameters are
# listed with initial values and bounds, "Fit model" estimates them
# (Levenberg-Marquardt), and estimates are reported with SE, RSE and 95 % CI,
# with plots of data and fitted profiles and the exposure parameters.
# =============================================================================

mod_fitting_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(5, 7),
    shiny::div(
      bslib::card(
        bslib::card_header("Experimental data"),
        shiny::fileInput(ns("data_file"), "Load concentration-time data (.csv, .txt, .xlsx)",
                         accept = c(".csv", ".txt", ".tsv", ".xlsx", ".xls"), width = "100%"),
        shiny::uiOutput(ns("import_ui")),
        DT::DTOutput(ns("obs_table"), fill = FALSE),
        shiny::div(
          class = "d-flex flex-wrap gap-2 mt-2",
          shiny::actionButton(ns("obs_add"), "Add row", icon = shiny::icon("plus"), class = "btn-sm"),
          shiny::actionButton(ns("obs_remove"), "Remove selected", icon = shiny::icon("trash"),
                              class = "btn-sm"),
          shiny::actionButton(ns("obs_clear"), "Clear all", icon = shiny::icon("eraser"),
                              class = "btn-sm")
        ),
        shiny::p(class = "small text-muted mt-1 mb-0",
                 "Time in h, concentration in µM. Double-click a cell to edit; new rows are added",
                 "to the compartment selected in the diagram.")
      ),
      bslib::card(
        bslib::card_header("Parameters to estimate"),
        shiny::selectizeInput(ns("unknowns"), NULL, choices = NULL, multiple = TRUE, width = "100%",
                              options = list(plugins = list("remove_button"),
                                             placeholder = "Select unknown parameters, e.g. CL, PS")),
        DT::DTOutput(ns("unknown_table"), fill = FALSE),
        shiny::p(class = "small text-muted mt-1",
                 "Initial value = current model value (edit it here or in 'Experimental conditions').",
                 "Default bounds: 1e-4 × and 1000 × the initial value."),
        bslib::accordion(
          open = FALSE,
          bslib::accordion_panel(
            "Fit settings",
            shiny::radioButtons(ns("obs_mode"), "Model value at the sampling times",
                                c("Nearest point of a 0.1 h grid (as the paper's analysis)" = "grid",
                                  "Exact sampling times" = "exact")),
            shiny::numericInput(ns("maxiter"), "Maximum iterations (Levenberg-Marquardt)",
                                value = 50, min = 1, step = 10),
            shiny::selectInput(ns("method"), "ODE solver",
                               c("lsoda (default)" = "lsoda", "ode45 (Dormand-Prince)" = "ode45",
                                 "bdf (stiff)" = "bdf"))
          )
        ),
        shiny::div(class = "mt-2",
                   shiny::actionButton(ns("fit"), "Fit model", class = "btn-primary",
                                       icon = shiny::icon("chart-line")))
      )
    ),
    shiny::div(
      bslib::layout_columns(
        col_widths = c(8, 4),
        bslib::navset_card_tab(
          full_screen = TRUE,
          title = "Model fit",
          bslib::nav_panel(
            "Plot",
            shiny::radioButtons(ns("yscale"), NULL, c("Linear" = "linear", "Logarithmic" = "log"),
                                inline = TRUE),
            plotly::plotlyOutput(ns("fit_plot"), height = "380px")),
          bslib::nav_panel("Residuals", plotly::plotlyOutput(ns("resid_plot"), height = "400px")),
          bslib::nav_panel("Table", DT::DTOutput(ns("resid_table"), fill = FALSE))
        ),
        bslib::card(
          bslib::card_header("Model"),
          mod_diagram_ui(ns("diagram"), height = "300px")
        )
      ),
      bslib::card(
        bslib::card_header("Estimated parameters"),
        shiny::uiOutput(ns("fit_summary")),
        DT::DTOutput(ns("estimates"), fill = FALSE),
        shiny::div(
          class = "d-flex flex-wrap gap-2 mt-2",
          shiny::actionButton(ns("apply"), "Apply estimates to model", icon = shiny::icon("check")),
          shiny::downloadButton(ns("dl_estimates"), "Estimates (CSV)", class = "btn-sm"),
          shiny::downloadButton(ns("dl_fitted"), "Fitted profiles (CSV)", class = "btn-sm")
        )
      ),
      bslib::card(
        bslib::card_header("Exposure parameters (fitted model)"),
        DT::DTOutput(ns("exposure"), fill = FALSE)
      )
    )
  )
}

mod_fitting_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    mod_diagram_server("diagram", state)

    modify <- function(fun) {
      res <- tryCatch(fun(state$project), error = function(e) {
        bsc_notify(conditionMessage(e), "error")
        NULL
      })
      if (!is.null(res)) state$project <- res
      invisible(!is.null(res))
    }
    selected_comp <- function() {
      sel <- state$selected
      cp <- state$project$compartments
      if (!is.null(sel) && sel$kind == "compartment" && sel$name %in% cp$name) return(sel$name)
      if (nrow(cp)) cp$name[cp$type == "fluid"][1] else NA_character_
    }

    # ---- data import --------------------------------------------------------------------
    raw <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$data_file, {
      f <- input$data_file
      r <- tryCatch({
        sheets <- bsc_data_sheets(f$datapath, f$name)
        list(df = bsc_read_data_file(f$datapath, sheet = sheets[1], name = f$name),
             name = f$name, path = f$datapath, sheets = sheets, sheet = sheets[1])
      }, error = function(e) {
        bsc_notify(paste("Could not read the file:", conditionMessage(e)), "error")
        NULL
      })
      raw(r)
    })
    shiny::observeEvent(input$sheet, {
      r <- raw()
      shiny::req(r, input$sheet)
      if (!identical(input$sheet, r$sheet)) {
        r$df <- bsc_read_data_file(r$path, sheet = input$sheet, name = r$name)
        r$sheet <- input$sheet
        raw(r)
      }
    })

    output$import_ui <- shiny::renderUI({
      r <- raw()
      shiny::req(r)
      g <- bsc_guess_columns(r$df)
      nm <- names(r$df)
      shiny::div(
        class = "bsc-import border rounded p-2 mb-2",
        shiny::strong(r$name),
        if (length(r$sheets) > 1) {
          shiny::selectInput(ns("sheet"), "Sheet", r$sheets, selected = r$sheet, width = "100%")
        },
        shiny::div(
          class = "d-flex gap-2 flex-wrap",
          shiny::selectInput(ns("col_time"), "Time (h)", nm, selected = g$time, width = "150px"),
          shiny::selectInput(ns("col_conc"), "Concentration (µM)", nm, selected = g$conc,
                             width = "170px"),
          shiny::selectInput(ns("col_comp"), "Compartment column",
                             c("(none - one compartment)" = "", nm),
                             selected = if (is.na(g$compartment)) "" else g$compartment,
                             width = "200px")
        ),
        shiny::uiOutput(ns("map_ui")),
        shiny::radioButtons(ns("blq"), "Entries below the LLOQ ('< x')",
                            c("set to 0" = "zero", "drop" = "na", "LLOQ / 2" = "half"), inline = TRUE),
        shiny::radioButtons(ns("import_mode"), NULL,
                            c("Replace data of these compartments" = "replace", "Append" = "append"),
                            inline = TRUE),
        shiny::actionButton(ns("do_import"), "Import data", icon = shiny::icon("file-import"),
                            class = "btn-primary btn-sm")
      )
    })

    map_values <- shiny::reactive({
      r <- raw()
      cc <- input$col_comp
      if (is.null(r) || is.null(cc) || !nzchar(cc) || !cc %in% names(r$df)) return(character(0))
      v <- unique(stats::na.omit(trimws(r$df[[cc]])))
      v[nzchar(v)]
    })
    guess_comp <- function(value, comps) {
      i <- match(tolower(value), tolower(comps))
      if (!is.na(i)) return(comps[i])
      hit <- comps[grepl(tolower(value), tolower(comps), fixed = TRUE) |
                     vapply(tolower(comps), function(cn) grepl(cn, tolower(value), fixed = TRUE), TRUE)]
      if (length(hit) == 1) hit else ""
    }
    output$map_ui <- shiny::renderUI({
      r <- raw()
      shiny::req(r)
      comps <- state$project$compartments$name
      vals <- map_values()
      if (!length(vals)) {
        return(shiny::selectInput(ns("assign_comp"), "Compartment of all rows", comps,
                                  selected = shiny::isolate(selected_comp()), width = "220px"))
      }
      if (length(vals) > 25) {
        return(shiny::p(class = "text-danger small", "The compartment column has more than 25",
                        "different values - choose another column."))
      }
      shiny::tagList(
        shiny::p(class = "small mb-1", "Model compartment of each value of the compartment column:"),
        shiny::div(class = "d-flex gap-2 flex-wrap", lapply(seq_along(vals), function(k) {
          shiny::selectInput(ns(paste0("map_", k)), vals[k], c("(ignore)" = "", comps),
                             selected = guess_comp(vals[k], comps), width = "170px")
        }))
      )
    })

    shiny::observeEvent(input$do_import, {
      r <- raw()
      shiny::req(r)
      vals <- map_values()
      map <- NULL
      if (length(vals)) {
        map <- stats::setNames(vapply(seq_along(vals), function(k) input[[paste0("map_", k)]] %||% "", ""),
                               vals)
      }
      obs <- tryCatch(
        bsc_observed_from_table(r$df, input$col_time, input$col_conc,
                                comp_col = if (length(vals)) input$col_comp else NULL,
                                comp_map = map, compartment = input$assign_comp,
                                blq = input$blq %||% "zero"),
        error = function(e) {
          bsc_notify(conditionMessage(e), "error")
          NULL
        })
      shiny::req(obs)
      if (!nrow(obs)) {
        bsc_notify("No data points imported - check the column choices and the mapping.", "warning")
        return()
      }
      modify(function(p) {
        if (identical(input$import_mode, "replace")) {
          p$observed <- p$observed[!p$observed$compartment %in% unique(obs$compartment), , drop = FALSE]
        }
        p$observed <- rbind(p$observed, obs)
        # make sure the simulation covers the data
        if (max(obs$time) > p$settings$sim_time) p$settings$sim_time <- ceiling(max(obs$time))
        p
      })
      bsc_notify(sprintf("Imported %d data points (%s).", nrow(obs),
                         paste(unique(obs$compartment), collapse = ", ")))
      raw(NULL)
    })

    # ---- observed data table ------------------------------------------------------------
    obs_tab <- shiny::reactiveVal(NULL)
    shiny::observe({
      v <- state$project$observed
      if (!identical(v, shiny::isolate(obs_tab()))) obs_tab(v)
    })
    output$obs_table <- DT::renderDT({
      o <- obs_tab()
      shiny::req(!is.null(o))
      df <- data.frame(Compartment = o$compartment, `Time (h)` = o$time, `Conc (µM)` = o$conc,
                       check.names = FALSE, stringsAsFactors = FALSE)
      bsc_dt(df, editable_cols = names(df), scroll_y = "220px")
    })
    shiny::observeEvent(input$obs_table_cell_edit, {
      info <- input$obs_table_cell_edit
      modify(function(p) {
        i <- info$row
        j <- info$col + 1
        if (j == 1) {
          v <- trimws(as.character(info$value))
          if (!v %in% p$compartments$name) stop("Unknown compartment '", v, "'.", call. = FALSE)
          p$observed$compartment[i] <- v
        } else {
          p$observed[[c("", "time", "conc")[j]]][i] <- bsc_cell_number(info$value)
        }
        p
      })
    })
    shiny::observeEvent(input$obs_add, modify(function(p) {
      cn <- selected_comp()
      if (is.na(cn)) stop("Add a compartment to the model first.", call. = FALSE)
      prev <- p$observed$time[p$observed$compartment == cn]
      p$observed <- rbind(p$observed, data.frame(compartment = cn,
                                                 time = if (length(prev)) max(prev) + 1 else 0,
                                                 conc = 0, stringsAsFactors = FALSE))
      p
    }))
    shiny::observeEvent(input$obs_remove, modify(function(p) {
      r <- input$obs_table_rows_selected
      if (length(r)) p$observed <- p$observed[-r, , drop = FALSE]
      rownames(p$observed) <- NULL
      p
    }))
    shiny::observeEvent(input$obs_clear, modify(function(p) {
      p$observed <- bsc_empty_observed()
      p
    }))

    # ---- unknown parameters ----------------------------------------------------------------
    ptab <- shiny::reactiveVal(NULL)
    guard <- bsc_echo_guard()   # ignore the browser's echo of the selection pushed below
    shiny::observe({
      v <- bsc_parameter_table(state$project)[, c("key", "label", "estimate")]
      if (!identical(v, shiny::isolate(ptab()))) ptab(v)
    })
    shiny::observeEvent(ptab(), {
      v <- ptab()
      ch <- stats::setNames(v$key, v$label)
      guard$note("unknowns", v$key[v$estimate])
      shiny::updateSelectizeInput(session, "unknowns", choices = ch, selected = v$key[v$estimate])
    })
    shiny::observeEvent(input$unknowns, {
      if (guard$is_echo("unknowns", input$unknowns %||% character(0))) return()
      v <- ptab()
      shiny::req(v)
      want <- v$key %in% (input$unknowns %||% character(0))
      if (identical(want, v$estimate)) return()
      modify(function(p) {
        pt <- bsc_parameter_table(p)
        for (i in which(pt$estimate != (pt$key %in% input$unknowns))) {
          p <- bsc_set_param(p, pt$kind[i], pt$element[i], pt$parameter[i],
                             estimate = pt$key[i] %in% input$unknowns)
        }
        p
      })
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    fitpar <- shiny::reactive(bsc_fit_parameters(state$project))
    output$unknown_table <- DT::renderDT({
      fp <- fitpar()
      df <- data.frame(Parameter = fp$label, Initial = fp$initial, Lower = fp$lower,
                       Upper = fp$upper, Unit = fp$unit, stringsAsFactors = FALSE)
      bsc_dt(df, editable_cols = c("Initial", "Lower", "Upper"), selection = "none",
             scroll_y = "160px", digits = 4)
    })
    shiny::observeEvent(input$unknown_table_cell_edit, {
      info <- input$unknown_table_cell_edit
      fp <- fitpar()
      i <- info$row
      shiny::req(i >= 1, i <= nrow(fp))
      col <- c("Parameter", "Initial", "Lower", "Upper", "Unit")[info$col + 1]
      modify(function(p) {
        v <- bsc_cell_number(info$value)
        if (col == "Initial") {
          return(bsc_set_param(p, fp$kind[i], fp$element[i], fp$parameter[i], value = v))
        }
        if (!col %in% c("Lower", "Upper")) return(p)
        b <- p$fit$bounds
        k <- which(b$kind == fp$kind[i] & b$element == fp$element[i] & b$parameter == fp$parameter[i])
        if (!length(k)) {
          b <- rbind(b, data.frame(kind = fp$kind[i], element = fp$element[i],
                                   parameter = fp$parameter[i], lower = NA_real_, upper = NA_real_,
                                   stringsAsFactors = FALSE))
          k <- nrow(b)
        }
        b[[tolower(col)]][k] <- v
        p$fit$bounds <- b
        p
      })
    })

    # ---- fit settings stored in the project ------------------------------------------------------
    shiny::observeEvent(state$project_id, {
      shiny::updateRadioButtons(session, "obs_mode", selected = state$project$fit$obs_mode %||% "grid")
      shiny::updateNumericInput(session, "maxiter", value = state$project$fit$maxiter %||% 50)
    })
    shiny::observeEvent(list(input$obs_mode, input$maxiter), {
      p <- state$project
      om <- input$obs_mode %||% "grid"
      mi <- as.integer(input$maxiter %||% 50)
      if (!identical(p$fit$obs_mode, om) || !identical(as.integer(p$fit$maxiter), mi)) {
        p$fit$obs_mode <- om
        p$fit$maxiter <- mi
        state$project <- p
      }
    }, ignoreInit = TRUE)

    # ---- fitting ---------------------------------------------------------------------------------
    shiny::observeEvent(input$fit, {
      p <- state$project
      v <- bsc_validate(p, "fit")
      if (length(v$errors)) {
        bsc_notify(paste(v$errors, collapse = "\n"), "error")
        return()
      }
      prog <- tryCatch(shiny::Progress$new(session, min = 0, max = 1), error = function(e) NULL)
      if (!is.null(prog)) {
        on.exit(prog$close(), add = TRUE)
        prog$set(message = "Fitting the model", value = 0.05)
      }
      progress <- function(n) {
        if (!is.null(prog)) prog$set(value = min(0.95, 0.05 + n / 150),
                                     detail = sprintf("%d model evaluations", n))
      }
      res <- tryCatch(
        bsc_fit(p, obs_mode = input$obs_mode %||% "grid",
                maxiter = as.integer(input$maxiter %||% 50),
                method = input$method %||% "lsoda", progress = progress),
        error = function(e) {
          bsc_notify(paste("Fit failed:", conditionMessage(e)), "error")
          NULL
        })
      shiny::req(res)
      state$fit <- list(result = res, keys = res$estimates$key, colors = bsc_compartment_colors(p),
                        created = Sys.time())
      if (!res$converged) {
        bsc_notify(paste("The optimiser stopped without convergence:", res$message), "warning")
      }
    })

    shiny::observeEvent(input$apply, {
      f <- state$fit
      if (is.null(f)) {
        bsc_notify("Fit the model first.", "warning")
        return()
      }
      cur <- bsc_parameter_table(state$project)
      if (!identical(cur$key[cur$estimate], f$keys)) {
        bsc_notify("The unknown parameters changed since the fit - please fit again.", "error")
        return()
      }
      modify(function(p) bsc_apply_estimates(p, f$result))
      bsc_notify("Estimated values written into the model (they stay marked as unknown).")
    })

    # ---- results -------------------------------------------------------------------------------------
    preview <- shiny::reactive({
      p <- state$project
      shiny::req(nrow(p$observed) > 0)
      fp <- bsc_fit_parameters(p)
      if (nrow(fp)) p <- bsc_apply_estimates(p, fp$initial)
      tryCatch(bsc_simulate(p, final_time = max(p$settings$sim_time + 0.1, max(p$observed$time) + 1)),
               error = function(e) NULL)
    }) |> shiny::debounce(500)

    output$fit_plot <- plotly::renderPlotly({
      p <- state$project
      obs <- p$observed
      f <- state$fit
      if (!nrow(obs) && is.null(f)) {
        return(bsc_plotly_message("Load experimental data to start."))
      }
      comps <- unique(c(obs$compartment, if (!is.null(f)) f$result$residuals$compartment))
      comps <- intersect(p$compartments$name, comps)
      cols <- bsc_compartment_colors(p)
      ylog <- identical(input$yscale, "log")
      fig <- plotly::plot_ly()
      sim <- if (!is.null(f)) f$result$simulation else preview()
      line_name <- if (!is.null(f)) "fit" else "current model"
      if (!is.null(sim)) {
        idx <- bsc_thin_index(sim$time, sim$segment)
        for (cc in intersect(comps, colnames(sim$conc))) {
          y <- sim$conc[idx, cc]
          if (ylog) y[y <= 0] <- NA
          fig <- plotly::add_trace(
            fig, x = sim$time[idx], y = y, type = "scatter", mode = "lines",
            name = paste0(cc, " · ", line_name), legendgroup = cc,
            line = list(color = cols[[cc]], width = 2, dash = if (is.null(f)) "dash" else "solid"),
            hovertemplate = paste0(cc, " ", line_name, ": %{y:.4g} µM<extra></extra>"))
        }
      }
      for (cc in comps) {
        o <- obs[obs$compartment == cc, , drop = FALSE]
        if (!nrow(o)) next
        y <- o$conc
        if (ylog) y[y <= 0] <- NA
        fig <- plotly::add_trace(
          fig, x = o$time, y = y, type = "scatter", mode = "markers",
          name = paste0(cc, " · observed"), legendgroup = cc,
          marker = list(color = cols[[cc]], size = 9, line = list(color = "#ffffff", width = 1.5)),
          hovertemplate = paste0(cc, " observed: %{y:.4g} µM at %{x} h<extra></extra>"))
      }
      bsc_plotly_layout(fig, ylog = ylog) |> plotly::layout(hovermode = "closest")
    })

    output$resid_plot <- plotly::renderPlotly({
      f <- state$fit
      if (is.null(f)) return(bsc_plotly_message("Fit the model to see the residuals."))
      rd <- f$result$residuals
      cols <- f$colors
      fig <- plotly::plot_ly()
      for (cc in unique(rd$compartment)) {
        o <- rd[rd$compartment == cc, ]
        fig <- plotly::add_trace(
          fig, x = o$time, y = o$residual, type = "scatter", mode = "markers", name = cc,
          marker = list(color = cols[[cc]] %||% BSC_SERIES_COLORS[1], size = 9,
                        line = list(color = "#ffffff", width = 1.5)),
          hovertemplate = paste0(cc, ": %{y:.3g} µM at %{x} h<extra></extra>"))
      }
      bsc_plotly_layout(fig, ylab = "Residual: predicted - observed (µM)") |>
        plotly::layout(hovermode = "closest",
                       shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = 0,
                                          y1 = 0, line = list(color = BSC_INK$axis, width = 1))))
    })

    output$resid_table <- DT::renderDT({
      f <- state$fit
      shiny::req(f)
      rd <- f$result$residuals
      df <- data.frame(Compartment = rd$compartment, `Time (h)` = rd$time,
                       `Observed (µM)` = rd$observed, `Predicted (µM)` = rd$predicted,
                       `Residual (µM)` = rd$residual, check.names = FALSE)
      bsc_dt(df, selection = "none", scroll_y = "360px", digits = 5)
    })

    output$fit_summary <- shiny::renderUI({
      f <- state$fit
      if (is.null(f)) {
        return(shiny::p(class = "text-muted", "No fit yet. Select the unknown parameters and press",
                        "'Fit model'."))
      }
      r <- f$result
      status <- if (r$converged) {
        shiny::span(class = "text-success", shiny::icon("circle-check"), "Converged")
      } else {
        shiny::span(class = "text-danger", shiny::icon("circle-exclamation"), "Not converged")
      }
      bound <- if (any(r$estimates$at_bound)) {
        shiny::div(class = "text-warning small", shiny::icon("triangle-exclamation"),
                   "At least one estimate is at a bound - widen the bounds or check identifiability.")
      }
      shiny::div(
        class = "small mb-2",
        status, sprintf(" (%s) after %d iterations, %d model runs, %.1f s.", r$message, r$iterations,
                        r$n_eval, r$elapsed),
        shiny::br(),
        sprintf("n = %d observations, p = %d parameters, SSR = %s, residual SD = %s µM, AIC = %s.",
                r$n_obs, r$n_par, bsc_fmt(r$ssr), bsc_fmt(r$sigma), bsc_fmt(r$aic)),
        shiny::br(),
        sprintf("SE from the Jacobian at the solution; %g %% CI = estimate ± t(%d) × SE.",
                100 * r$conf_level, r$df),
        " RSE > 50 % is highlighted (poorly identifiable).",
        bound
      )
    })

    output$estimates <- DT::renderDT({
      f <- state$fit
      shiny::req(f)
      e <- f$result$estimates
      df <- data.frame(Parameter = e$label, Estimate = e$estimate, SE = e$se,
                       `RSE (%)` = e$rse_percent, `CI 2.5%` = e$ci_lower, `CI 97.5%` = e$ci_upper,
                       Unit = e$unit, Initial = e$initial,
                       Bound = ifelse(e$at_bound, "at bound", ""),
                       check.names = FALSE, stringsAsFactors = FALSE)
      bsc_dt(df, selection = "none", scroll_y = "200px", digits = 4) |>
        DT::formatStyle("RSE (%)", color = DT::styleInterval(50, c(BSC_INK$primary, "#d03b3b")),
                        fontWeight = DT::styleInterval(50, c("normal", "bold")))
    })

    output$exposure <- DT::renderDT({
      f <- state$fit
      shiny::req(f, f$result$simulation)
      sim <- f$result$simulation
      comps <- intersect(unique(f$result$residuals$compartment), colnames(sim$conc))
      e <- bsc_exposure(sim$time, sim$conc[, comps, drop = FALSE])
      note <- e$status
      note[is.na(e$t_half)] <- paste0(note[is.na(e$t_half)], "; Cmax/2 not reached")
      df <- data.frame(Compartment = e$compartment, `AUC (µM·h)` = e$AUC,
                       `Cmax (µM)` = e$Cmax, `Tmax (h)` = e$Tmax, `Half-life (h)` = e$t_half,
                       Note = note, check.names = FALSE)
      bsc_dt(df, selection = "none", scroll_y = "160px", digits = 4)
    })

    output$dl_estimates <- shiny::downloadHandler(
      filename = function() paste0("BaySimChip_estimates_", format(Sys.Date()), ".csv"),
      content = function(file) {
        f <- state$fit
        e <- if (is.null(f)) data.frame() else f$result$estimates
        utils::write.csv(e, file, row.names = FALSE)
      })
    output$dl_fitted <- shiny::downloadHandler(
      filename = function() paste0("BaySimChip_fitted_profiles_", format(Sys.Date()), ".csv"),
      content = function(file) {
        f <- state$fit
        df <- if (is.null(f) || is.null(f$result$simulation)) data.frame() else
          bsc_simulation_long(f$result$simulation)
        utils::write.csv(df, file, row.names = FALSE)
      })
  })
}
