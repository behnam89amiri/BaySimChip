# =============================================================================
# mod_element_editor.R - properties of the selected model element
#
# Side panel of the drag-and-drop editor (Model structure tab) and the
# "Chip and biological parameters" panel of the Experimental conditions tab
# (it shows the fields of the element selected in the diagram):
#   compartment:        name, type, volume, CL, fu, Kev
#   solid compartment:  name, type, volume, partition coefficient K_media/plastic
#                       (stored in the same field as fu: C_medium = K * C_solid)
#   fluid flow:         name, Q
#   permeable surface:  name, PS (PS_NSB if it connects to a solid compartment)
# Every parameter has an "unknown (estimate)" check box; the value of an
# unknown parameter is its initial guess for fitting.
# =============================================================================

.bsc_param_row <- function(id_value, id_unknown, label) {
  shiny::div(
    class = "bsc-param-row d-flex align-items-end gap-3",
    shiny::div(class = "bsc-param-value",
               shiny::numericInput(id_value, label, value = NA, min = 0, width = "100%")),
    shiny::div(class = "bsc-param-unknown mb-3",
               shiny::tags$span(title = "Estimate this parameter in 'Fitting and analysis'",
                                shiny::checkboxInput(id_unknown, "unknown", value = FALSE)))
  )
}

# Labels (HTML) - the same notation as in the manuscript
.bsc_label_fu    <- "Fraction unbound, fu (&ndash;)"
.bsc_label_kmp   <- "Partition coefficient, K<sub>media/plastic</sub> (&ndash;)"
.bsc_label_ps    <- "Permeability &times; area, PS (mL/h)"
.bsc_label_psnsb <- "NSB exchange rate, PS<sub>NSB</sub> (mL/h)"

mod_element_editor_ui <- function(id, show_type = TRUE, show_kev = TRUE) {
  ns <- shiny::NS(id)
  is_comp <- "output.kind == 'fluid' || output.kind == 'solid'"
  shiny::tagList(
    shiny::selectInput(ns("element"), "Selected element (or click it in the diagram)",
                       choices = NULL, width = "100%"),
    shiny::uiOutput(ns("info")),
    shiny::conditionalPanel(
      "output.kind != ''", ns = ns,
      shiny::div(
        class = "d-flex gap-2 align-items-end",
        shiny::div(class = "flex-grow-1", shiny::textInput(ns("name"), "Name", width = "100%")),
        shiny::div(class = "mb-3", shiny::actionButton(ns("rename"), "Rename", class = "btn-sm"))
      )
    ),
    if (show_type) shiny::conditionalPanel(
      is_comp, ns = ns,
      shiny::selectInput(ns("type"), "Type",
                         c("Compartment (fluid)" = "fluid",
                           "Solid compartment (non-specific binding)" = "solid"),
                         width = "100%")),
    shiny::conditionalPanel(
      is_comp, ns = ns,
      .bsc_param_row(ns("vol"), ns("vol_unknown"), "Volume (mL)"),
      shiny::conditionalPanel(
        "output.kind == 'fluid'", ns = ns,
        .bsc_param_row(ns("cl"), ns("cl_unknown"), "Clearance, CL (mL/h)")),
      .bsc_param_row(ns("fu"), ns("fu_unknown"), shiny::HTML(.bsc_label_fu)),
      if (show_kev) shiny::conditionalPanel(
        "output.kind == 'fluid'", ns = ns,
        .bsc_param_row(ns("kev"), ns("kev_unknown"), "Evaporation rate, Kev (mL/h)"))
    ),
    shiny::conditionalPanel(
      "output.kind == 'flow'", ns = ns,
      .bsc_param_row(ns("q"), ns("q_unknown"), "Flow rate, Q (mL/h)")),
    shiny::conditionalPanel(
      "output.kind == 'perm'", ns = ns,
      .bsc_param_row(ns("ps"), ns("ps_unknown"), shiny::HTML(.bsc_label_ps))),
    shiny::p(class = "small text-muted mb-0",
             "Parameters marked as unknown are estimated in 'Fitting and analysis';",
             "their current value is the initial guess.")
  )
}

#' @return reactive with the selected element: list(kind, name, type) or NULL
mod_element_editor_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
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

    # ---- element selection (synced with the diagrams / other tabs) ----------------
    choices <- shiny::reactiveVal(NULL)
    shiny::observe({
      ch <- bsc_element_choices(state$project)
      if (!identical(ch, shiny::isolate(choices()))) choices(ch)
    })
    shiny::observeEvent(choices(), {
      all <- unlist(choices(), use.names = FALSE)
      sel <- bsc_sel_key(shiny::isolate(state$selected))
      if (!sel %in% all) sel <- if (length(all)) all[1] else character(0)
      shiny::updateSelectInput(session, "element", choices = choices(), selected = sel)
    }, ignoreNULL = FALSE)
    shiny::observeEvent(state$selected, {
      key <- bsc_sel_key(state$selected)
      if (!identical(key, input$element) &&
          key %in% unlist(shiny::isolate(choices()), use.names = FALSE)) {
        shiny::updateSelectInput(session, "element", selected = key)
      }
    })
    shiny::observeEvent(input$element, {
      sel <- bsc_sel_parse(input$element)
      if (!is.null(sel) && !identical(sel, state$selected)) state$selected <- sel
    })

    current <- shiny::reactive({
      sel <- bsc_sel_parse(input$element)
      p <- state$project
      if (is.null(sel) || !bsc_element_exists(p, sel)) return(NULL)
      if (sel$kind == "compartment") sel$type <- p$compartments$type[p$compartments$name == sel$name]
      sel
    })

    output$kind <- shiny::renderText({
      el <- current()
      if (is.null(el)) "" else if (el$kind == "compartment") el$type else el$kind
    })
    shiny::outputOptions(output, "kind", suspendWhenHidden = FALSE)

    output$info <- shiny::renderUI({
      el <- current()
      if (is.null(el)) {
        return(shiny::p(class = "text-muted", "No element selected. Build the model in the",
                        "diagram of the 'Model structure' tab."))
      }
      p <- state$project
      txt <- switch(el$kind,
                    compartment = if (el$type == "solid") {
                      shiny::HTML(paste(
                        "Solid compartment: chip material that binds the compound non-specifically (NSB).",
                        "At equilibrium C<sub>medium</sub> = K<sub>media/plastic</sub> &times; C<sub>solid</sub>,",
                        "so a smaller K<sub>media/plastic</sub> means stronger binding. The volume is a",
                        "reference volume of the material; changing it only rescales K<sub>media/plastic</sub>."))
                    } else "Fluid compartment (medium / tissue).",
                    flow = {
                      f <- p$flows[p$flows$name == el$name, ]
                      paste0("Fluid flow from ", f$from, " to ", f$to, ".")
                    },
                    perm = {
                      s <- p$perms[p$perms$name == el$name, ]
                      if (bsc_is_nsb_surface(p, el$name)) {
                        paste0("Exchange surface between ", s$between1, " & ", s$between2,
                               " (binding to the chip material, NSB).")
                      } else {
                        paste0("Permeable surface between ", s$between1, " & ", s$between2, ".")
                      }
                    })
      shiny::p(class = "small text-muted", txt)
    })

    # ---- project -> form (only fields whose stored value differs) -------------------
    last_name_key <- ""
    last_fu_label <- ""
    last_ps_label <- ""
    set_label <- function(id, html) {
      session$sendCustomMessage("bsc-set-label", list(id = session$ns(id), html = html))
    }
    shiny::observe({
      el <- current()
      shiny::req(el)
      p <- state$project
      key <- bsc_sel_key(el)
      if (!identical(key, last_name_key)) {
        shiny::updateTextInput(session, "name", value = el$name)
        last_name_key <<- key
      }
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
      if (el$kind == "compartment") {
        r <- p$compartments[p$compartments$name == el$name, ]
        if (!identical(r$type, shiny::isolate(input$type))) {
          guard$note("type", r$type)
          shiny::updateSelectInput(session, "type", selected = r$type)
        }
        upd_num("vol", r$volume)
        upd_chk("vol_unknown", r$volume_est)
        upd_num("cl", r$cl)
        upd_chk("cl_unknown", r$cl_est)
        upd_num("fu", r$fu)
        upd_chk("fu_unknown", r$fu_est)
        upd_num("kev", r$kev)
        upd_chk("kev_unknown", r$kev_est)
        lbl <- if (r$type == "solid") .bsc_label_kmp else .bsc_label_fu
        if (!identical(lbl, last_fu_label)) {
          set_label("fu", lbl)
          last_fu_label <<- lbl
        }
      } else if (el$kind == "flow") {
        r <- p$flows[p$flows$name == el$name, ]
        upd_num("q", r$q)
        upd_chk("q_unknown", r$q_est)
      } else if (el$kind == "perm") {
        r <- p$perms[p$perms$name == el$name, ]
        upd_num("ps", r$ps)
        upd_chk("ps_unknown", r$ps_est)
        lbl <- if (bsc_is_nsb_surface(p, el$name)) .bsc_label_psnsb else .bsc_label_ps
        if (!identical(lbl, last_ps_label)) {
          set_label("ps", lbl)
          last_ps_label <<- lbl
        }
      }
    })

    # ---- form -> project ----------------------------------------------------------------
    shiny::observeEvent(input$rename, {
      el <- current()
      shiny::req(el)
      new <- trimws(input$name %||% "")
      if (identical(new, el$name)) return()
      ok <- modify(function(p) bsc_rename(p, el$kind, el$name, new))
      if (ok) {
        state$selected <- list(kind = el$kind, name = new)
      } else {
        shiny::updateTextInput(session, "name", value = el$name)
      }
    })

    shiny::observeEvent(input$type, {
      if (guard$is_echo("type", input$type)) return()
      el <- current()
      if (is.null(el) || el$kind != "compartment" || identical(input$type, el$type)) return()
      ok <- modify(function(p) bsc_set_compartment_type(p, el$name, input$type))
      if (!ok) shiny::updateSelectInput(session, "type", selected = el$type)
    }, ignoreInit = TRUE)

    set_param <- function(kind, parameter, value = NULL, estimate = NULL) {
      el <- shiny::isolate(current())
      if (is.null(el) || el$kind != kind) return(invisible())
      cur <- bsc_get_param(state$project, kind, el$name, parameter)
      if (!is.null(value)) {
        value <- if (length(value)) as.numeric(value) else NA_real_
        if (same(cur$value, value)) value <- NULL
      }
      if (!is.null(estimate) && identical(cur$estimate, isTRUE(estimate))) estimate <- NULL
      if (is.null(value) && is.null(estimate)) return(invisible())
      modify(function(p) bsc_set_param(p, kind, el$name, parameter, value = value, estimate = estimate))
    }
    bind_param <- function(id, kind, parameter) {
      id_unknown <- paste0(id, "_unknown")
      shiny::observeEvent(input[[id]], {
        if (guard$is_echo(id, input[[id]])) return()
        set_param(kind, parameter, value = input[[id]])
      }, ignoreInit = TRUE, ignoreNULL = FALSE)
      shiny::observeEvent(input[[id_unknown]], {
        if (guard$is_echo(id_unknown, isTRUE(input[[id_unknown]]))) return()
        set_param(kind, parameter, estimate = input[[id_unknown]])
      }, ignoreInit = TRUE)
    }
    bind_param("vol", "compartment", "volume")
    bind_param("cl", "compartment", "cl")
    bind_param("fu", "compartment", "fu")
    bind_param("kev", "compartment", "kev")
    bind_param("q", "flow", "q")
    bind_param("ps", "perm", "ps")

    current
  })
}
