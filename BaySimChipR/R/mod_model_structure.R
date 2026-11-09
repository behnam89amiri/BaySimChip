# =============================================================================
# mod_model_structure.R - tab 1 "Model structure"
#
# Block-diagram editor with the buttons New model / Load model / Save model,
# Compartment, Solid compartment, Permeable surface, Fluid flow and "Remove
# selected":
#   * centre: drag-and-drop editor (visNetwork in manipulation mode) with
#     these buttons; compartments are placed by clicking in the canvas, flows
#     and permeable surfaces are drawn by dragging between compartments,
#     compartments can be moved, the selection can be removed;
#   * right: properties of the selected element (name, type, volume, CL, fu,
#     Kev; Q; PS);
#   * left: project files - native format .rds (model, conditions, layout,
#     data and the last fit), JSON export, bundled examples;
#   * bottom: editable tables as a secondary view, synced with the diagram.
# =============================================================================

mod_model_structure_ui <- function(id) {
  ns <- shiny::NS(id)
  examples <- bsc_example_choices()
  shiny::tagList(
    bslib::layout_columns(
      col_widths = c(3, 6, 3),
      bslib::card(
        bslib::card_header("Project"),
        shiny::div(
          class = "d-flex flex-wrap gap-2 mb-2",
          shiny::actionButton(ns("new"), "New model", icon = shiny::icon("file"), class = "btn-sm"),
          shiny::downloadButton(ns("save"), "Save model", icon = shiny::icon("floppy-disk"),
                                class = "btn-sm"),
          shiny::downloadButton(ns("export_json"), "Export JSON", icon = shiny::icon("file-code"),
                                class = "btn-sm btn-outline-secondary")
        ),
        shiny::fileInput(ns("load"), "Load model (.rds, .json)", accept = c(".rds", ".json"),
                         width = "100%"),
        shiny::selectInput(ns("example"), "Example projects", choices = examples, width = "100%"),
        shiny::actionButton(ns("load_example"), "Load example", icon = shiny::icon("folder-open"),
                            class = "btn-sm mb-3"),
        shiny::textInput(ns("name"), "Model name", width = "100%"),
        shiny::textAreaInput(ns("description"), "Description", rows = 3, width = "100%")
      ),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("Model"),
        mod_diagram_ui(ns("diagram"), height = "520px", editable = TRUE),
        shiny::p(class = "small text-muted mb-0",
                 "Click a button above, then click in the diagram (compartments) or drag from one",
                 "compartment to another (flows, permeable surfaces). Drag compartments to arrange",
                 "the layout; select an element to edit its properties.")
      ),
      bslib::card(
        bslib::card_header("Properties"),
        mod_element_editor_ui(ns("editor"), show_type = TRUE, show_kev = TRUE)
      )
    ),
    bslib::accordion(
      open = FALSE,
      class = "mb-3",
      bslib::accordion_panel(
        "Model tables (alternative to the diagram)",
        icon = shiny::icon("table"),
        bslib::navset_underline(
          id = ns("edit_tabs"),
          bslib::nav_panel(
            "Compartments",
            DT::DTOutput(ns("comp_table"), fill = FALSE),
            shiny::p(class = "small text-muted mt-1",
                     "Double-click a cell to rename a compartment or change its type (fluid / solid)."),
            shiny::div(
              class = "d-flex flex-wrap gap-2 align-items-end",
              shiny::textInput(ns("new_comp_name"), "Name", placeholder = "e.g. Liver", width = "140px"),
              shiny::selectInput(ns("new_comp_type"), "Type",
                                 c("Compartment (fluid)" = "fluid", "Solid compartment (NSB)" = "solid"),
                                 width = "200px"),
              shiny::div(class = "mb-3",
                         shiny::actionButton(ns("add_comp"), "Add", icon = shiny::icon("plus")),
                         shiny::actionButton(ns("remove_comp"), "Remove selected",
                                             icon = shiny::icon("trash")))
            )
          ),
          bslib::nav_panel(
            "Fluid flows",
            DT::DTOutput(ns("flow_table"), fill = FALSE),
            shiny::p(class = "small text-muted mt-1",
                     "Flows connect fluid compartments (direction from → to). Q in mL/h."),
            shiny::div(
              class = "d-flex flex-wrap gap-2 align-items-end",
              shiny::selectInput(ns("new_flow_from"), "From", choices = NULL, width = "130px"),
              shiny::selectInput(ns("new_flow_to"), "To", choices = NULL, width = "130px"),
              shiny::numericInput(ns("new_flow_q"), "Q (mL/h)", value = 1, min = 0, width = "100px"),
              shiny::div(class = "mb-3",
                         shiny::actionButton(ns("add_flow"), "Add", icon = shiny::icon("plus")),
                         shiny::actionButton(ns("remove_flow"), "Remove selected",
                                             icon = shiny::icon("trash")))
            )
          ),
          bslib::nav_panel(
            "Permeable surfaces",
            DT::DTOutput(ns("perm_table"), fill = FALSE),
            shiny::p(class = "small text-muted mt-1",
                     "Permeable surfaces exchange compound between two compartments (fluid or solid).",
                     "PS = permeability × surface area in mL/h."),
            shiny::div(
              class = "d-flex flex-wrap gap-2 align-items-end",
              shiny::selectInput(ns("new_perm_1"), "Between", choices = NULL, width = "130px"),
              shiny::selectInput(ns("new_perm_2"), "and", choices = NULL, width = "130px"),
              shiny::numericInput(ns("new_perm_ps"), "PS (mL/h)", value = 0.01, min = 0,
                                  width = "100px"),
              shiny::div(class = "mb-3",
                         shiny::actionButton(ns("add_perm"), "Add", icon = shiny::icon("plus")),
                         shiny::actionButton(ns("remove_perm"), "Remove selected",
                                             icon = shiny::icon("trash")))
            )
          )
        )
      )
    )
  )
}

mod_model_structure_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    mod_diagram_server("diagram", state, editable = TRUE)
    mod_element_editor_server("editor", state)

    # apply a change to the project; errors are reported, not raised
    modify <- function(fun) {
      res <- tryCatch(fun(state$project), error = function(e) {
        bsc_notify(conditionMessage(e), "error")
        NULL
      })
      if (!is.null(res)) state$project <- res
      invisible(!is.null(res))
    }

    # ---- project name / description ---------------------------------------------
    shiny::observeEvent(state$project_id, {
      shiny::updateTextInput(session, "name", value = state$project$name)
      shiny::updateTextAreaInput(session, "description", value = state$project$description)
    })
    shiny::observeEvent(input$name, {
      if (!identical(input$name, state$project$name)) {
        p <- state$project
        p$name <- input$name
        state$project <- p
      }
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$description, {
      if (!identical(input$description, state$project$description)) {
        p <- state$project
        p$description <- input$description
        state$project <- p
      }
    }, ignoreInit = TRUE)

    # ---- new / load / example / save --------------------------------------------------
    shiny::observeEvent(input$new, {
      bsc_set_project(state, bsc_new_project())
      bsc_notify("New empty model: add compartments with the buttons above the diagram.")
    })

    shiny::observeEvent(input$load, {
      f <- input$load
      p <- tryCatch(bsc_read_project(f$datapath, name = f$name), error = function(e) {
        bsc_notify(paste("Could not load the model:", conditionMessage(e)), "error")
        NULL
      })
      shiny::req(p)
      bsc_set_project(state, p)
      bsc_notify(paste0("Loaded '", f$name, "'."))
    })

    shiny::observeEvent(input$load_example, {
      shiny::req(input$example)
      p <- tryCatch(bsc_load_example(input$example), error = function(e) {
        bsc_notify(conditionMessage(e), "error")
        NULL
      })
      shiny::req(p)
      bsc_set_project(state, p)
      bsc_notify(paste0("Loaded example '", input$example, "'."))
    })

    file_stem <- function() {
      nm <- gsub("[^A-Za-z0-9_.-]+", "_", state$project$name)
      if (nzchar(nm)) nm else "BaySimChip_model"
    }
    output$save <- shiny::downloadHandler(
      filename = function() paste0(file_stem(), ".rds"),
      content = function(file) bsc_write_project(bsc_project_for_save(state), file, format = "rds"),
      contentType = "application/octet-stream"
    )
    output$export_json <- shiny::downloadHandler(
      filename = function() paste0(file_stem(), ".json"),
      content = function(file) bsc_write_project(state$project, file, format = "json"),
      contentType = "application/json"
    )

    # ---- compartments table ---------------------------------------------------------------
    output$comp_table <- DT::renderDT({
      cp <- state$project$compartments
      fl <- state$project$flows
      pm <- state$project$perms
      df <- data.frame(Name = cp$name, Type = cp$type,
                       `Volume (mL)` = cp$volume,
                       Flows = vapply(cp$name, function(n) sum(fl$from == n | fl$to == n), 0),
                       `Permeable surfaces` = vapply(cp$name, function(n) {
                         sum(pm$between1 == n | pm$between2 == n)
                       }, 0),
                       check.names = FALSE, stringsAsFactors = FALSE)
      bsc_dt(df, editable_cols = c("Name", "Type", "Volume (mL)"))
    })

    shiny::observeEvent(input$comp_table_cell_edit, {
      info <- input$comp_table_cell_edit
      cp <- state$project$compartments
      i <- info$row
      j <- info$col + 1
      shiny::req(i >= 1, i <= nrow(cp))
      val <- trimws(as.character(info$value))
      if (j == 1) {
        if (modify(function(p) bsc_rename(p, "compartment", cp$name[i], val))) {
          state$selected <- list(kind = "compartment", name = val)
        }
      } else if (j == 2) {
        modify(function(p) bsc_set_compartment_type(p, cp$name[i], tolower(val)))
      } else if (j == 3) {
        modify(function(p) bsc_set_param(p, "compartment", cp$name[i], "volume", bsc_cell_number(val)))
      }
    })

    shiny::observeEvent(input$add_comp, {
      ok <- modify(function(p) {
        bsc_auto_layout(bsc_add_compartment(p, name = input$new_comp_name, type = input$new_comp_type))
      })
      if (ok) {
        cp <- state$project$compartments
        state$selected <- list(kind = "compartment", name = cp$name[nrow(cp)])
        shiny::updateTextInput(session, "new_comp_name", value = "")
      }
    })

    # rows selected in a table, or the element selected in the diagram
    selected_names <- function(kind, rows) {
      tbl <- switch(kind, compartment = state$project$compartments, flow = state$project$flows,
                    perm = state$project$perms)
      out <- tbl$name[rows[rows <= nrow(tbl)]]
      sel <- state$selected
      if (!length(out) && !is.null(sel) && sel$kind == kind) out <- sel$name
      out
    }

    pending <- shiny::reactiveVal(NULL)
    confirm_remove <- function(kind, names) {
      if (!length(names)) {
        bsc_notify("Select the items to remove first (table rows or diagram).", "warning")
        return(invisible())
      }
      pending(list(kind = kind, names = names))
      what <- switch(kind, compartment = "compartment(s)", flow = "fluid flow(s)",
                     perm = "permeable surface(s)")
      extra <- if (kind == "compartment") {
        " Connected flows, permeable surfaces, dosing and data are removed as well."
      } else ""
      shiny::showModal(shiny::modalDialog(
        title = "Remove selected items",
        paste0("Would you like to delete the ", what, " ", paste(names, collapse = ", "), "?", extra),
        footer = shiny::tagList(shiny::modalButton("Cancel"),
                                shiny::actionButton(ns("confirm_remove"), "Delete",
                                                    class = "btn-danger")),
        easyClose = TRUE))
    }

    shiny::observeEvent(input$remove_comp, {
      confirm_remove("compartment", selected_names("compartment", input$comp_table_rows_selected))
    })
    shiny::observeEvent(input$remove_flow, {
      confirm_remove("flow", selected_names("flow", input$flow_table_rows_selected))
    })
    shiny::observeEvent(input$remove_perm, {
      confirm_remove("perm", selected_names("perm", input$perm_table_rows_selected))
    })
    shiny::observeEvent(input$confirm_remove, {
      x <- pending()
      shiny::req(x)
      modify(function(p) {
        switch(x$kind,
               compartment = bsc_remove_compartments(p, x$names),
               flow = bsc_remove_flows(p, x$names),
               perm = bsc_remove_perms(p, x$names))
      })
      if (!is.null(state$selected) && !bsc_element_exists(state$project, state$selected)) {
        state$selected <- NULL
      }
      pending(NULL)
      shiny::removeModal()
    })

    # ---- flows table -------------------------------------------------------------------------
    output$flow_table <- DT::renderDT({
      fl <- state$project$flows
      df <- data.frame(Name = fl$name, From = fl$from, To = fl$to, `Q (mL/h)` = fl$q,
                       check.names = FALSE, stringsAsFactors = FALSE)
      bsc_dt(df, editable_cols = names(df))
    })

    shiny::observeEvent(input$flow_table_cell_edit, {
      info <- input$flow_table_cell_edit
      fl <- state$project$flows
      i <- info$row
      j <- info$col + 1
      shiny::req(i >= 1, i <= nrow(fl))
      modify(function(p) {
        if (j == 1) return(bsc_rename(p, "flow", fl$name[i], info$value))
        if (j == 4) return(bsc_set_param(p, "flow", fl$name[i], "q", bsc_cell_number(info$value)))
        col <- c("", "from", "to")[j]
        p$flows[[col]][i] <- trimws(as.character(info$value))
        v <- bsc_validate(p, "structure")
        if (length(v$errors)) stop(paste(v$errors, collapse = " "), call. = FALSE)
        p
      })
    })

    shiny::observeEvent(input$add_flow, {
      ok <- modify(function(p) bsc_add_flow(p, input$new_flow_from, input$new_flow_to,
                                            q = input$new_flow_q))
      if (ok) {
        fl <- state$project$flows
        state$selected <- list(kind = "flow", name = fl$name[nrow(fl)])
      }
    })

    # ---- permeable surfaces table ----------------------------------------------------------------
    output$perm_table <- DT::renderDT({
      pm <- state$project$perms
      df <- data.frame(Name = pm$name, `Between 1` = pm$between1, `Between 2` = pm$between2,
                       `PS (mL/h)` = pm$ps, check.names = FALSE, stringsAsFactors = FALSE)
      bsc_dt(df, editable_cols = names(df))
    })

    shiny::observeEvent(input$perm_table_cell_edit, {
      info <- input$perm_table_cell_edit
      pm <- state$project$perms
      i <- info$row
      j <- info$col + 1
      shiny::req(i >= 1, i <= nrow(pm))
      modify(function(p) {
        if (j == 1) return(bsc_rename(p, "perm", pm$name[i], info$value))
        if (j == 4) return(bsc_set_param(p, "perm", pm$name[i], "ps", bsc_cell_number(info$value)))
        col <- c("", "between1", "between2")[j]
        p$perms[[col]][i] <- trimws(as.character(info$value))
        v <- bsc_validate(p, "structure")
        if (length(v$errors)) stop(paste(v$errors, collapse = " "), call. = FALSE)
        p
      })
    })

    shiny::observeEvent(input$add_perm, {
      ok <- modify(function(p) bsc_add_perm(p, input$new_perm_1, input$new_perm_2,
                                            ps = input$new_perm_ps))
      if (ok) {
        pm <- state$project$perms
        state$selected <- list(kind = "perm", name = pm$name[nrow(pm)])
      }
    })

    # ---- choices of the add forms ------------------------------------------------------------------
    comp_names <- shiny::reactiveVal(NULL)
    shiny::observe({
      cp <- state$project$compartments
      v <- list(all = cp$name, fluid = cp$name[cp$type == "fluid"])
      if (!identical(v, shiny::isolate(comp_names()))) comp_names(v)
    })
    shiny::observeEvent(comp_names(), {
      v <- comp_names()
      shiny::updateSelectInput(session, "new_flow_from", choices = v$fluid)
      shiny::updateSelectInput(session, "new_flow_to", choices = v$fluid,
                               selected = if (length(v$fluid) > 1) v$fluid[2] else NULL)
      shiny::updateSelectInput(session, "new_perm_1", choices = v$all)
      shiny::updateSelectInput(session, "new_perm_2", choices = v$all,
                               selected = if (length(v$all) > 1) v$all[2] else NULL)
    })

    # ---- selection sync between tables and diagram ------------------------------------------------
    shiny::observeEvent(input$comp_table_rows_selected, {
      r <- input$comp_table_rows_selected
      cp <- state$project$compartments
      if (length(r) == 1 && r <= nrow(cp)) state$selected <- list(kind = "compartment", name = cp$name[r])
    })
    shiny::observeEvent(input$flow_table_rows_selected, {
      r <- input$flow_table_rows_selected
      fl <- state$project$flows
      if (length(r) == 1 && r <= nrow(fl)) state$selected <- list(kind = "flow", name = fl$name[r])
    })
    shiny::observeEvent(input$perm_table_rows_selected, {
      r <- input$perm_table_rows_selected
      pm <- state$project$perms
      if (length(r) == 1 && r <= nrow(pm)) state$selected <- list(kind = "perm", name = pm$name[r])
    })
    shiny::observeEvent(state$selected, {
      sel <- state$selected
      shiny::req(sel)
      tbl <- switch(sel$kind, compartment = "comp_table", flow = "flow_table", perm = "perm_table")
      rows <- switch(sel$kind, compartment = state$project$compartments$name,
                     flow = state$project$flows$name, perm = state$project$perms$name)
      i <- match(sel$name, rows)
      cur <- input[[paste0(tbl, "_rows_selected")]]
      if (!is.na(i) && !identical(as.integer(cur), as.integer(i))) {
        DT::selectRows(DT::dataTableProxy(tbl), i)
      }
      tab <- switch(sel$kind, compartment = "Compartments", flow = "Fluid flows",
                    perm = "Permeable surfaces")
      bslib::nav_select("edit_tabs", tab)
    })
  })
}
