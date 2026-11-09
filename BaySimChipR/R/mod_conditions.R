# =============================================================================
# mod_conditions.R - tab 2 "Experimental conditions"
#
# Panel "Chip and biological parameters" for the selected
# compartment (volume, clearance, fu), solid compartment (partition coefficient
# K_media/plastic), flow (flow rate) or permeable surface (permeability x area,
# PS; PS_NSB for exchange with a solid compartment); panel "Dosing"
# (single / repetitive / custom administration, dose in mg, compound molar
# mass, "use first experimental concentration"); panel "Sampling and media
# exchange" (evaporation rate, table of time / extracted / added volume).
# Each parameter has an explicit "unknown (estimate)" check box and its value
# is the initial guess.
# =============================================================================

.bsc_table_buttons <- function(ns, prefix, load_label = "Load (.xlsx/.csv)") {
  shiny::div(
    class = "d-flex flex-wrap gap-2 align-items-start mt-2",
    shiny::actionButton(ns(paste0(prefix, "_add")), "Add row", icon = shiny::icon("plus"),
                        class = "btn-sm"),
    shiny::actionButton(ns(paste0(prefix, "_remove")), "Remove selected", icon = shiny::icon("trash"),
                        class = "btn-sm"),
    shiny::actionButton(ns(paste0(prefix, "_clear")), "Clear", icon = shiny::icon("eraser"),
                        class = "btn-sm"),
    shiny::div(class = "bsc-inline-file",
               shiny::fileInput(ns(paste0(prefix, "_file")), NULL, buttonLabel = load_label,
                                accept = c(".xlsx", ".xls", ".csv", ".txt"), width = "230px"))
  )
}

mod_conditions_ui <- function(id) {
  ns <- shiny::NS(id)
  unit_choices <- stats::setNames(names(BSC_DOSE_UNITS), BSC_DOSE_UNIT_LABELS)
  # wide screens: parameters | dosing and sampling | model and parameter table;
  # medium screens: two columns, dosing and sampling below the parameters
  bslib::layout_columns(
    col_widths = bslib::breakpoints(sm = 12, lg = c(6, 6), xxl = c(8, 4)),
    bslib::layout_columns(
      col_widths = bslib::breakpoints(sm = 12, xxl = c(6, 6)),
      shiny::div(
        bslib::card(
          bslib::card_header("Chip and biological parameters"),
          mod_element_editor_ui(ns("editor"), show_type = FALSE, show_kev = FALSE)
        )
      ),
      shiny::div(
        shiny::conditionalPanel(
          "output.element_kind == 'fluid'", ns = ns,
          bslib::card(
            bslib::card_header("Dosing"),
            shiny::radioButtons(ns("dose_mode"), NULL,
                                c("Single dose" = "single", "Repetitive administration" = "repeated",
                                  "Custom administration" = "custom"), inline = TRUE),
            shiny::conditionalPanel(
              "input.dose_mode == 'single'", ns = ns,
              shiny::numericInput(ns("dose_single"), "Dose at t = 0", value = 0, min = 0, width = "200px")),
            shiny::conditionalPanel(
              "input.dose_mode == 'repeated'", ns = ns,
              shiny::div(class = "d-flex gap-3 flex-wrap",
                         shiny::numericInput(ns("dose_rep"), "Dose", value = 0, min = 0, width = "160px"),
                         shiny::numericInput(ns("dose_interval"), "Every (h)", value = 24, min = 0.01,
                                             width = "160px"))),
            shiny::conditionalPanel(
              "input.dose_mode == 'custom'", ns = ns,
              DT::DTOutput(ns("dose_table"), fill = FALSE),
              .bsc_table_buttons(ns, "dose")),
            shiny::checkboxInput(ns("first_dose"),
                                 "Use first experimental concentration (if available) for dose calculation",
                                 width = "100%"),
            shiny::div(
              class = "d-flex gap-3 flex-wrap",
              shiny::selectInput(ns("dose_unit"), "Dose unit (all compartments)", unit_choices,
                                 width = "200px"),
              shiny::numericInput(ns("molar_mass"), "Compound molar mass (g/mol)", value = 300,
                                  min = 0, width = "220px"))
          ),
          bslib::card(
            bslib::card_header("Sampling and media exchange"),
            .bsc_param_row(ns("kev"), ns("kev_unknown"), "Evaporation rate, Kev (mL/h)"),
            DT::DTOutput(ns("exchange_table"), fill = FALSE),
            shiny::p(class = "small text-muted mt-1 mb-0",
                     "Time (h), volume extracted (mL) and volume added (mL). Sampling = extraction",
                     "without addition; a media exchange removes and adds medium (drug-free)."),
            .bsc_table_buttons(ns, "ex")
          )
        ),
        shiny::conditionalPanel(
          "output.element_kind != 'fluid'", ns = ns,
          shiny::p(class = "small text-muted mt-2",
                   "Dosing, sampling and media exchange are set per fluid compartment.",
                   "Select a compartment in the diagram or in the list.")
        )
      )
    ),
    shiny::div(
      bslib::card(
        bslib::card_header("Model"),
        mod_diagram_ui(ns("diagram"), height = "330px")
      ),
      bslib::card(
        bslib::card_header("All parameters"),
        DT::DTOutput(ns("param_table"), fill = FALSE),
        shiny::p(class = "small text-muted mt-1 mb-0",
                 "Double-click a value or the 'Unknown' cell (yes / no) to edit. Click a row to select",
                 "the element.")
      )
    )
  )
}

mod_conditions_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    mod_diagram_server("diagram", state)
    current <- mod_element_editor_server("editor", state)

    modify <- function(fun) {
      res <- tryCatch(fun(state$project), error = function(e) {
        bsc_notify(conditionMessage(e), "error")
        NULL
      })
      if (!is.null(res)) state$project <- res
      invisible(!is.null(res))
    }
    same <- function(a, b) {
      a <- if (length(a)) a[[1]] else NA
      b <- if (length(b)) b[[1]] else NA
      (is.na(a) && is.na(b)) || (!is.na(a) && !is.na(b) && isTRUE(all.equal(a, b)))
    }

    guard <- bsc_echo_guard()   # ignore the browser's echo of values pushed to the form

    output$element_kind <- shiny::renderText({
      el <- current()
      if (is.null(el)) "" else if (el$kind == "compartment") el$type else el$kind
    })
    shiny::outputOptions(output, "element_kind", suspendWhenHidden = FALSE)

    # ---- project -> dosing / evaporation fields ----------------------------------------
    shiny::observe({
      el <- current()
      shiny::req(el, el$kind == "compartment")
      p <- state$project
      upd_num <- function(id, value) {
        if (!same(value, shiny::isolate(input[[id]]))) {
          guard$note(id, value)
          shiny::updateNumericInput(session, id, value = value)
        }
      }
      upd_chk <- function(id, value) {
        if (!identical(isTRUE(value), isTRUE(shiny::isolate(input[[id]])))) {
          guard$note(id, isTRUE(value))
          shiny::updateCheckboxInput(session, id, value = isTRUE(value))
        }
      }
      r <- p$compartments[p$compartments$name == el$name, ]
      upd_num("kev", r$kev)
      upd_chk("kev_unknown", r$kev_est)
      d <- bsc_get_dosing(p, el$name)
      if (!identical(d$mode, shiny::isolate(input$dose_mode))) {
        guard$note("dose_mode", d$mode)
        shiny::updateRadioButtons(session, "dose_mode", selected = d$mode)
      }
      upd_num("dose_single", d$single)
      upd_num("dose_rep", d$repeated_dose)
      upd_num("dose_interval", d$repeated_interval)
      upd_chk("first_dose", d$first_dose_from_data)
    })

    shiny::observeEvent(state$project_id, {
      s <- state$project$settings
      guard$note("dose_unit", s$dose_unit)
      guard$note("molar_mass", s$molar_mass)
      shiny::updateSelectInput(session, "dose_unit", selected = s$dose_unit)
      shiny::updateNumericInput(session, "molar_mass", value = s$molar_mass)
    })

    # ---- evaporation rate (shown in the media exchange panel) ---------
    set_kev <- function(value = NULL, estimate = NULL) {
      el <- shiny::isolate(current())
      if (is.null(el) || el$kind != "compartment") return(invisible())
      cur <- bsc_get_param(state$project, "compartment", el$name, "kev")
      if (!is.null(value)) {
        value <- if (length(value)) as.numeric(value) else NA_real_
        if (same(cur$value, value)) value <- NULL
      }
      if (!is.null(estimate) && identical(cur$estimate, isTRUE(estimate))) estimate <- NULL
      if (is.null(value) && is.null(estimate)) return(invisible())
      modify(function(p) bsc_set_param(p, "compartment", el$name, "kev", value = value, estimate = estimate))
    }
    shiny::observeEvent(input$kev, {
      if (guard$is_echo("kev", input$kev)) return()
      set_kev(value = input$kev)
    }, ignoreInit = TRUE, ignoreNULL = FALSE)
    shiny::observeEvent(input$kev_unknown, {
      if (guard$is_echo("kev_unknown", isTRUE(input$kev_unknown))) return()
      set_kev(estimate = input$kev_unknown)
    }, ignoreInit = TRUE)

    # ---- dosing ---------------------------------------------------------------------------
    set_dosing <- function(field, value) {
      el <- shiny::isolate(current())
      if (is.null(el) || el$kind != "compartment") return(invisible())
      d <- bsc_get_dosing(state$project, el$name)
      if (is.numeric(d[[field]]) && same(d[[field]], value)) return(invisible())
      if (identical(d[[field]], value)) return(invisible())
      modify(function(p) bsc_set_dosing(p, el$name, stats::setNames(list(value), field)))
    }
    shiny::observeEvent(input$dose_mode, {
      if (guard$is_echo("dose_mode", input$dose_mode)) return()
      set_dosing("mode", input$dose_mode)
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$dose_single, {
      if (guard$is_echo("dose_single", input$dose_single)) return()
      set_dosing("single", as.numeric(input$dose_single))
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$dose_rep, {
      if (guard$is_echo("dose_rep", input$dose_rep)) return()
      set_dosing("repeated_dose", as.numeric(input$dose_rep))
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$dose_interval, {
      if (guard$is_echo("dose_interval", input$dose_interval)) return()
      set_dosing("repeated_interval", as.numeric(input$dose_interval))
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$first_dose, {
      if (guard$is_echo("first_dose", isTRUE(input$first_dose))) return()
      set_dosing("first_dose_from_data", isTRUE(input$first_dose))
      el <- current()
      if (isTRUE(input$first_dose) && !is.null(el) &&
          !is.finite(bsc_first_observed_conc(state$project, el$name))) {
        bsc_notify(paste0("No observed data for ", el$name, " yet - the dose is computed once data",
                          " are loaded in the 'Fitting and analysis' tab."), "warning")
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$dose_unit, {
      if (guard$is_echo("dose_unit", input$dose_unit)) return()
      if (!identical(input$dose_unit, state$project$settings$dose_unit)) {
        p <- state$project
        p$settings$dose_unit <- input$dose_unit
        state$project <- p
      }
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$molar_mass, {
      if (guard$is_echo("molar_mass", input$molar_mass)) return()
      if (!same(input$molar_mass, state$project$settings$molar_mass)) {
        p <- state$project
        p$settings$molar_mass <- as.numeric(input$molar_mass)
        state$project <- p
      }
    }, ignoreInit = TRUE)

    # ---- custom dosing table -----------------------------------------------------------------
    custom_tab <- shiny::reactiveVal(NULL)
    shiny::observe({
      el <- current()
      v <- if (!is.null(el) && el$kind == "compartment") bsc_get_dosing(state$project, el$name)$custom else NULL
      if (!identical(v, shiny::isolate(custom_tab()))) custom_tab(v)
    })
    output$dose_table <- DT::renderDT({
      cu <- custom_tab()
      shiny::req(!is.null(cu))
      unit <- BSC_DOSE_UNIT_LABELS[[state$project$settings$dose_unit]]
      df <- data.frame(cu$time, cu$dose)
      names(df) <- c("Time (h)", paste0("Dose (", unit, ")"))
      bsc_dt(df, editable_cols = names(df), scroll_y = "160px")
    })
    edit_custom <- function(fun) {
      el <- current()
      shiny::req(el, el$kind == "compartment")
      modify(function(p) {
        cu <- fun(bsc_get_dosing(p, el$name)$custom)
        bsc_set_dosing(p, el$name, list(custom = cu))
      })
    }
    shiny::observeEvent(input$dose_table_cell_edit, {
      info <- input$dose_table_cell_edit
      edit_custom(function(cu) {
        cu[info$row, info$col + 1] <- bsc_cell_number(info$value)
        cu
      })
    })
    shiny::observeEvent(input$dose_add, edit_custom(function(cu) {
      t_new <- if (nrow(cu)) max(cu$time, na.rm = TRUE) + 24 else 0
      d_new <- if (nrow(cu)) cu$dose[nrow(cu)] else 0
      rbind(cu, data.frame(time = t_new, dose = d_new))
    }))
    shiny::observeEvent(input$dose_remove, edit_custom(function(cu) {
      r <- input$dose_table_rows_selected
      if (length(r)) cu[-r, , drop = FALSE] else cu
    }))
    shiny::observeEvent(input$dose_clear, edit_custom(function(cu) cu[0, , drop = FALSE]))
    shiny::observeEvent(input$dose_file, {
      f <- input$dose_file
      edit_custom(function(cu) {
        m <- bsc_read_numeric_table(f$datapath, 2, name = f$name)
        data.frame(time = m[, 1], dose = m[, 2])
      })
    })

    # ---- media exchange table ------------------------------------------------------------------
    exchange_tab <- shiny::reactiveVal(NULL)
    shiny::observe({
      el <- current()
      v <- if (!is.null(el) && el$kind == "compartment") bsc_get_exchange(state$project, el$name) else NULL
      if (!identical(v, shiny::isolate(exchange_tab()))) exchange_tab(v)
    })
    output$exchange_table <- DT::renderDT({
      ex <- exchange_tab()
      shiny::req(!is.null(ex))
      df <- data.frame(`Time (h)` = ex$time, `Extracted volume (mL)` = ex$removed,
                       `Added volume (mL)` = ex$added, check.names = FALSE)
      bsc_dt(df, editable_cols = names(df), scroll_y = "180px")
    })
    edit_exchange <- function(fun) {
      el <- current()
      shiny::req(el, el$kind == "compartment")
      modify(function(p) bsc_set_exchange(p, el$name, fun(bsc_get_exchange(p, el$name))))
    }
    shiny::observeEvent(input$exchange_table_cell_edit, {
      info <- input$exchange_table_cell_edit
      edit_exchange(function(ex) {
        ex[info$row, info$col + 1] <- bsc_cell_number(info$value)
        ex
      })
    })
    shiny::observeEvent(input$ex_add, edit_exchange(function(ex) {
      t_new <- if (nrow(ex)) max(ex$time, na.rm = TRUE) + 24 else 24
      last <- if (nrow(ex)) ex[nrow(ex), ] else data.frame(time = 0, removed = 0, added = 0)
      rbind(ex, data.frame(time = t_new, removed = last$removed, added = last$added))
    }))
    shiny::observeEvent(input$ex_remove, edit_exchange(function(ex) {
      r <- input$exchange_table_rows_selected
      if (length(r)) ex[-r, , drop = FALSE] else ex
    }))
    shiny::observeEvent(input$ex_clear, edit_exchange(function(ex) ex[0, , drop = FALSE]))
    shiny::observeEvent(input$ex_file, {
      f <- input$ex_file
      edit_exchange(function(ex) {
        m <- bsc_read_numeric_table(f$datapath, 3, name = f$name)
        data.frame(time = m[, 1], removed = m[, 2], added = m[, 3])
      })
    })

    # ---- overview of all parameters ----------------------------------------------------------
    param_tab <- shiny::reactiveVal(NULL)
    shiny::observe({
      v <- bsc_parameter_table(state$project)
      if (!identical(v, shiny::isolate(param_tab()))) param_tab(v)
    })
    output$param_table <- DT::renderDT({
      pt <- param_tab()
      shiny::req(pt)
      short <- c(volume = "Volume", fu = "fu", cl = "CL", kev = "Kev", q = "Q", ps = "PS")
      par <- unname(short[pt$parameter])
      solid <- pt$kind == "compartment" & pt$parameter == "fu" &
        pt$element %in% state$project$compartments$name[state$project$compartments$type == "solid"]
      par[solid] <- "K_media/plastic"
      nsb <- pt$kind == "perm" &
        vapply(pt$element, function(n) bsc_is_nsb_surface(state$project, n), logical(1))
      par[nsb] <- "PS_NSB"
      df <- data.frame(Element = pt$element, Parameter = par, Value = pt$value, Unit = pt$unit,
                       Unknown = ifelse(pt$estimate, "yes", ""), stringsAsFactors = FALSE)
      bsc_dt(df, editable_cols = c("Value", "Unknown"), selection = "single", scroll_y = "330px") |>
        DT::formatStyle("Unknown", fontWeight = "bold", color = "#a8431a")
    })
    shiny::observeEvent(input$param_table_cell_edit, {
      info <- input$param_table_cell_edit
      pt <- param_tab()
      i <- info$row
      shiny::req(i >= 1, i <= nrow(pt))
      col <- c("Element", "Parameter", "Value", "Unit", "Unknown")[info$col + 1]
      modify(function(p) {
        if (col == "Value") {
          bsc_set_param(p, pt$kind[i], pt$element[i], pt$parameter[i],
                        value = bsc_cell_number(info$value))
        } else if (col == "Unknown") {
          bsc_set_param(p, pt$kind[i], pt$element[i], pt$parameter[i],
                        estimate = bsc_cell_logical(info$value))
        } else p
      })
    })
    shiny::observeEvent(input$param_table_rows_selected, {
      pt <- param_tab()
      r <- input$param_table_rows_selected
      if (length(r) == 1 && r <= nrow(pt)) state$selected <- list(kind = pt$kind[r], name = pt$element[r])
    })
  })
}
