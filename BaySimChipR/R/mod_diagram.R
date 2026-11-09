# =============================================================================
# mod_diagram.R - model diagram shared by all tabs (and the drag-and-drop editor)
#
# The model diagram is visible on every tab and clicking a compartment
# selects it: each tab embeds this sub-module.
#
#  * selecting a compartment / flow / permeable surface sets the shared
#    selection (state$selected);
#  * dragging compartments stores their x/y positions in the project, so the
#    layout is saved and restored with the model;
#  * editable = TRUE (Model structure tab) adds the drag-and-drop editor:
#    toolbar buttons put the vis.js network into manipulation mode (add
#    compartment / solid compartment by clicking in the canvas, draw a fluid
#    flow or a permeable surface by dragging from one compartment to another,
#    delete the selection).  The requests arrive as input$net_graphChange and
#    are applied to the project, which is the single source of truth.
#
# The widget is drawn once per project; later changes (from the editor, the
# tables or other tabs) are applied incrementally through a visNetwork proxy
# so that zoom, pan and selection are preserved.
# =============================================================================

mod_diagram_ui <- function(id, height = "400px", editable = FALSE) {
  ns <- shiny::NS(id)
  out_id <- ns("net")
  shiny::div(
    class = if (editable) "bsc-diagram bsc-editor" else "bsc-diagram",
    if (editable) bsc_editor_toolbar(out_id),
    if (editable) shiny::div(id = paste0(out_id, "_status"), class = "bsc-status small"),
    visNetwork::visNetworkOutput(out_id, height = height),
    bsc_network_legend()
  )
}

#' @param state shared reactive state (see bsc_new_state())
#' @param editable handle editor requests (input$net_graphChange)
#' @param on_select optional function(selection) called after a user selection
#' @param defaults default values of new flows / permeable surfaces
mod_diagram_server <- function(id, state, editable = FALSE, on_select = NULL,
                               defaults = list(q = 1, ps = 0.01)) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    net_data <- shiny::reactiveVal(NULL)
    redraw <- shiny::reactiveVal(0L)   # forces a full re-draw when incremented
    client <- NULL            # diagram data currently shown in the browser
    rendered_pid <- -1L       # project_id of the last full render

    shiny::observe({
      d <- bsc_network_data(state$project)
      if (!identical(d, shiny::isolate(net_data()))) net_data(d)
    })

    # full render only for a new project; everything else goes via the proxy
    output$net <- visNetwork::renderVisNetwork({
      d <- net_data()
      pid <- state$project_id
      redraw()
      shiny::req(d)
      if (!is.null(client) && identical(pid, rendered_pid)) shiny::req(FALSE, cancelOutput = TRUE)
      client <<- d
      rendered_pid <<- pid
      session$onFlushed(function() select_in_diagram(shiny::isolate(state$selected)), once = TRUE)
      bsc_network_widget(d, output_id = ns("net"), editable = editable)
    })

    shiny::observeEvent(net_data(), {
      d <- net_data()
      if (is.null(client) || identical(d, client)) return()
      bsc_network_sync(visNetwork::visNetworkProxy(ns("net")), client, d)
      client <<- d
      select_in_diagram(shiny::isolate(state$selected))
    })

    select_in_diagram <- function(sel) {
      if (is.null(client)) return(invisible())
      proxy <- visNetwork::visNetworkProxy(ns("net"))
      if (is.null(sel)) {
        visNetwork::visUnselectAll(proxy)
      } else if (sel$kind == "compartment") {
        visNetwork::visSelectNodes(proxy, id = paste0("compartment|", sel$name), highlightEdges = FALSE)
      } else {
        visNetwork::visSelectEdges(proxy, id = paste0(sel$kind, "|", sel$name))
      }
    }
    shiny::observeEvent(state$selected, select_in_diagram(state$selected),
                        ignoreNULL = FALSE, ignoreInit = TRUE)

    # ---- user selection --------------------------------------------------------------
    shiny::observeEvent(input$net_select, {
      sel <- bsc_sel_parse(input$net_select$id)
      shiny::req(sel)
      state$selected <- sel
      if (is.function(on_select)) on_select(sel)
    })

    # ---- node positions after dragging ------------------------------------------------
    shiny::observeEvent(input$net_drag, {
      pos <- input$net_drag$pos
      shiny::req(length(pos))
      p <- state$project
      changed <- FALSE
      for (q in pos) {
        sel <- bsc_sel_parse(q$id)
        if (is.null(sel) || sel$kind != "compartment") next
        i <- match(sel$name, p$compartments$name)
        if (is.na(i)) next
        if (!identical(p$compartments$x[i], as.numeric(q$x)) ||
            !identical(p$compartments$y[i], as.numeric(q$y))) {
          p$compartments$x[i] <- as.numeric(q$x)
          p$compartments$y[i] <- as.numeric(q$y)
          changed <- TRUE
        }
      }
      if (changed) state$project <- p
    })

    # ---- editor requests (vis.js manipulation mode) ----------------------------------
    if (editable) {
      shiny::observeEvent(input$net_graphChange, {
        ev <- input$net_graphChange
        res <- tryCatch(bsc_apply_graph_change(state$project, ev, defaults), error = function(e) {
          bsc_notify(conditionMessage(e), "error")
          NULL
        })
        if (is.null(res)) {
          # the browser may have changed already (e.g. a refused deletion):
          # redraw it from the project
          client <<- NULL
          redraw(redraw() + 1L)
          return()
        }
        state$project <- res$project
        if (!is.null(res$selected)) {
          state$selected <- res$selected
        } else if (!is.null(state$selected) && !bsc_element_exists(res$project, state$selected)) {
          state$selected <- NULL
        }
      })
    }

    invisible(net_data)
  })
}

#' Does a selected element (list(kind, name)) exist in the project?
bsc_element_exists <- function(p, sel) {
  if (is.null(sel)) return(FALSE)
  tbl <- switch(sel$kind, compartment = p$compartments, flow = p$flows, perm = p$perms, NULL)
  !is.null(tbl) && sel$name %in% tbl$name
}

#' Apply an editor request (input$<id>_graphChange) to a project
#'
#' @param ev list(cmd = "addNode" (x, y, type) | "addEdge" (from, to, type) |
#'   "deleteElements" (nodes, edges)); element ids are "kind|name"
#' @return list(project, selected) - the new element is selected
bsc_apply_graph_change <- function(p, ev, defaults = list(q = 1, ps = 0.01)) {
  cmd <- ev$cmd %||% ""
  if (cmd == "addNode") {
    type <- if (identical(ev$type, "solid")) "solid" else "fluid"
    p <- bsc_add_compartment(p, type = type, x = round(as.numeric(ev$x)), y = round(as.numeric(ev$y)))
    name <- p$compartments$name[nrow(p$compartments)]
    return(list(project = p, selected = list(kind = "compartment", name = name)))
  }
  if (cmd == "addEdge") {
    from <- bsc_sel_parse(ev$from)
    to <- bsc_sel_parse(ev$to)
    if (is.null(from) || is.null(to) || from$kind != "compartment" || to$kind != "compartment") {
      stop("Connections must be drawn between two compartments.", call. = FALSE)
    }
    if (identical(ev$type, "perm")) {
      p <- bsc_add_perm(p, from$name, to$name, ps = defaults$ps)
      return(list(project = p, selected = list(kind = "perm", name = p$perms$name[nrow(p$perms)])))
    }
    p <- bsc_add_flow(p, from$name, to$name, q = defaults$q)
    return(list(project = p, selected = list(kind = "flow", name = p$flows$name[nrow(p$flows)])))
  }
  if (cmd == "deleteElements") {
    ids <- c(unlist(ev$nodes), unlist(ev$edges))
    sel <- lapply(ids, bsc_sel_parse)
    kinds <- vapply(sel, function(s) s$kind, "")
    names <- vapply(sel, function(s) s$name, "")
    p <- bsc_remove_compartments(p, names[kinds == "compartment"])
    p <- bsc_remove_flows(p, names[kinds == "flow"])
    p <- bsc_remove_perms(p, names[kinds == "perm"])
    return(list(project = p, selected = NULL))
  }
  stop("Unknown editor command '", cmd, "'.", call. = FALSE)
}
