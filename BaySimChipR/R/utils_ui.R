# =============================================================================
# utils_ui.R - helpers shared by the Shiny modules (no numerics here)
# =============================================================================

# ---- colours ---------------------------------------------------------------
# Categorical palette (validated for colour-vision deficiency on the light
# surface, adjacent pairs); slots are assigned to compartments in table order
# so that a compartment keeps its colour in every plot.
BSC_SERIES_COLORS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                       "#e87ba4", "#008300", "#4a3aa7", "#e34948")
BSC_INK <- list(primary = "#0b0b0b", secondary = "#52514e", muted = "#898781",
                grid = "#e1e0d9", axis = "#c3c2b7", surface = "#fcfcfb")
BSC_NODE_FLUID <- list(background = "#2a78d6", border = "#1c5cab", font = "#ffffff")
BSC_NODE_SOLID <- list(background = "#c3c2b7", border = "#52514e", font = "#0b0b0b")
BSC_HIGHLIGHT <- list(background = "#eb6834", border = "#a8431a")
BSC_FLOW_COLOR <- "#1c5cab"
BSC_PERM_COLOR <- "#52514e"
BSC_RUN_DASHES <- c("solid", "dash", "dot", "dashdot", "longdash", "longdashdot")

#' Colour of each compartment (fixed by position in the compartment table)
bsc_compartment_colors <- function(p) {
  nm <- p$compartments$name
  cols <- rep(BSC_SERIES_COLORS, length.out = length(nm))
  stats::setNames(cols, nm)
}

# ---- notifications ---------------------------------------------------------
bsc_notify <- function(msg, type = c("message", "warning", "error"), duration = 6) {
  type <- match.arg(type)
  try(shiny::showNotification(msg, type = type, duration = if (type == "error") 10 else duration),
      silent = TRUE)
  invisible(msg)
}

# ---- number formatting -----------------------------------------------------
bsc_fmt <- function(x, digits = 4) {
  out <- formatC(signif(x, digits), digits = digits, format = "g", flag = "#")
  out <- sub("\\.$", "", sub("0+(e|$)", "\\1", out))
  out[is.na(x)] <- ""
  trimws(out)
}

# ---- shared reactive state -------------------------------------------------

#' Reactive state shared by all modules
#'
#' project    the current project (see engine_project.R)
#' project_id incremented whenever a different project is loaded (UI reset)
#' selected   list(kind = "compartment"|"flow"|"perm", name) or NULL
#' sim_runs   list of simulation runs shown in the Simulation tab
#' fit        last fit result (bsc_fit) or NULL
bsc_new_state <- function(project = bsc_new_project()) {
  shiny::reactiveValues(project = project, project_id = 1L, selected = NULL,
                        sim_runs = list(), fit = NULL)
}

#' Replace the project (New / Load / Example); a fit stored in the project
#' (p$results$fit, written by "Save model") is restored if it still matches
#' the unknown parameters
bsc_set_project <- function(state, p) {
  p <- bsc_auto_layout(bsc_normalize_project(p))
  fit <- p$results$fit
  pt <- bsc_parameter_table(p)
  if (!is.null(fit) && !identical(fit$keys, pt$key[pt$estimate])) fit <- NULL
  state$project <- p
  state$sim_runs <- list()
  state$fit <- fit
  first <- if (nrow(p$compartments)) list(kind = "compartment", name = p$compartments$name[1]) else NULL
  state$selected <- first
  state$project_id <- shiny::isolate(state$project_id) + 1L
  invisible(p)
}

#' Project with the last fit attached (saved in the native .rds format)
bsc_project_for_save <- function(state) {
  p <- shiny::isolate(state$project)
  f <- shiny::isolate(state$fit)
  if (is.null(f)) {
    p$results <- NULL
  } else {
    f$result$nls <- NULL   # optimiser internals (closures) are not stored
    p$results <- list(fit = f)
  }
  p
}

# key <-> selection
bsc_sel_key <- function(sel) if (is.null(sel)) "" else paste(sel$kind, sel$name, sep = "|")
bsc_sel_parse <- function(key) {
  if (is.null(key) || !nzchar(key)) return(NULL)
  parts <- strsplit(key, "|", fixed = TRUE)[[1]]
  list(kind = parts[1], name = paste(parts[-1], collapse = "|"))
}

#' Grouped choices for a selectInput listing all model elements
bsc_element_choices <- function(p) {
  ch <- list()
  if (nrow(p$compartments)) {
    lab <- ifelse(p$compartments$type == "solid", paste0(p$compartments$name, " (solid)"),
                  p$compartments$name)
    ch[["Compartments"]] <- stats::setNames(paste0("compartment|", p$compartments$name), lab)
  }
  if (nrow(p$flows)) {
    ch[["Fluid flows"]] <- stats::setNames(paste0("flow|", p$flows$name),
                                           paste0(p$flows$name, " (", p$flows$from, " → ", p$flows$to, ")"))
  }
  if (nrow(p$perms)) {
    ch[["Permeable surfaces"]] <- stats::setNames(paste0("perm|", p$perms$name),
                                                  paste0(p$perms$name, " (", p$perms$between1, " ↔ ",
                                                         p$perms$between2, ")"))
  }
  ch
}

# ---- DT helpers ------------------------------------------------------------

#' Compact DT table; `editable_cols` are the (display) columns users may edit
bsc_dt <- function(df, editable_cols = character(0), selection = "multiple",
                   scroll_y = "240px", digits = NULL, dom = "t") {
  disable <- which(!names(df) %in% editable_cols) - 1L
  editable <- if (length(editable_cols)) {
    list(target = "cell", disable = list(columns = disable))
  } else FALSE
  dt <- DT::datatable(
    df, rownames = FALSE, selection = selection, editable = editable,
    class = "compact hover cell-border",
    options = list(dom = dom, paging = FALSE, ordering = FALSE, scrollY = scroll_y,
                   scrollCollapse = TRUE, autoWidth = FALSE,
                   language = list(emptyTable = "No rows - use the buttons below to add some.")))
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  if (!is.null(digits) && length(num_cols)) dt <- DT::formatSignif(dt, num_cols, digits = digits)
  dt
}

#' Value of an edited DT cell as number (empty -> NA)
bsc_cell_number <- function(v) {
  v <- trimws(as.character(v))
  if (!nzchar(v)) return(NA_real_)
  out <- suppressWarnings(as.numeric(sub(",", ".", v, fixed = TRUE)))
  if (is.na(out)) stop("'", v, "' is not a number.", call. = FALSE)
  out
}

bsc_cell_logical <- function(v) {
  v <- tolower(trimws(as.character(v)))
  if (v %in% c("yes", "y", "true", "t", "1", "x", "unknown", "estimate")) return(TRUE)
  if (v %in% c("no", "n", "false", "f", "0", "", "known")) return(FALSE)
  stop("Use 'yes' or 'no'.", call. = FALSE)
}

# ---- model diagram (visNetwork) ------------------------------------------------

#' Nodes and edges of the model diagram
bsc_network_data <- function(p) {
  p <- bsc_auto_layout(p)
  cp <- p$compartments
  if (!nrow(cp)) {
    return(list(nodes = data.frame(id = character(0), label = character(0)),
                edges = data.frame(id = character(0), from = character(0), to = character(0))))
  }
  nodes <- data.frame(
    id = paste0("compartment|", cp$name),
    label = cp$name,
    group = cp$type,
    x = cp$x, y = cp$y,
    title = ifelse(cp$type == "solid",
                   paste0("<b>", cp$name, "</b><br>Solid compartment (non-specific binding)"),
                   paste0("<b>", cp$name, "</b><br>Fluid compartment")),
    stringsAsFactors = FALSE)
  fl <- p$flows
  pm <- p$perms
  edges <- rbind(
    if (nrow(fl)) data.frame(
      id = paste0("flow|", fl$name), from = paste0("compartment|", fl$from),
      to = paste0("compartment|", fl$to), label = fl$name, arrows = "to", dashes = FALSE,
      color.color = BSC_FLOW_COLOR, color.highlight = BSC_HIGHLIGHT$background,
      width = 2, smooth.type = "curvedCW",
      title = paste0("<b>", fl$name, "</b><br>Fluid flow ", fl$from, " → ", fl$to),
      stringsAsFactors = FALSE),
    if (nrow(pm)) data.frame(
      id = paste0("perm|", pm$name), from = paste0("compartment|", pm$between1),
      to = paste0("compartment|", pm$between2), label = pm$name, arrows = "", dashes = TRUE,
      color.color = BSC_PERM_COLOR, color.highlight = BSC_HIGHLIGHT$background,
      width = 3, smooth.type = "curvedCCW",
      title = paste0("<b>", pm$name, "</b><br>Permeable surface ", pm$between1, " ↔ ",
                     pm$between2),
      stringsAsFactors = FALSE))
  if (is.null(edges)) {
    edges <- data.frame(id = character(0), from = character(0), to = character(0))
  }
  list(nodes = nodes, edges = edges)
}

#' visNetwork widget of the model diagram
#'
#' @param output_id namespaced id of the visNetwork output; the widget sends
#'   input$<output_id>_select (clicked element), input$<output_id>_drag (node
#'   positions after dragging) and, when editable, input$<output_id>_graphChange
#'   (requests to add / delete elements, as visNetwork's manipulation mode)
#' @param editable enable the drag-and-drop editor (vis.js manipulation mode)
bsc_network_widget <- function(data, output_id, editable = FALSE, height = "400px") {
  js_select <- sprintf(
    "function(p) {
       var id = null;
       if (p.nodes.length > 0) { id = p.nodes[0]; } else if (p.edges.length > 0) { id = p.edges[0]; }
       if (id !== null) { Shiny.setInputValue('%s_select', {id: id, nonce: Math.random()}, {priority: 'event'}); }
     }", output_id)
  js_drag <- sprintf(
    "function(p) {
       if (p.nodes.length > 0) {
         var pos = this.getPositions(p.nodes); var out = [];
         for (var k in pos) { out.push({id: k, x: Math.round(pos[k].x), y: Math.round(pos[k].y)}); }
         Shiny.setInputValue('%s_drag', {pos: out, nonce: Math.random()}, {priority: 'event'});
       }
     }", output_id)
  font <- "system-ui, -apple-system, 'Segoe UI', sans-serif"
  w <- visNetwork::visNetwork(data$nodes, data$edges, height = height, width = "100%") |>
    visNetwork::visGroups(groupname = "fluid", shape = "box", borderWidth = 1.5,
                          margin = 12, shapeProperties = list(borderRadius = 8),
                          font = list(color = BSC_NODE_FLUID$font, size = 16, face = font),
                          color = list(background = BSC_NODE_FLUID$background,
                                       border = BSC_NODE_FLUID$border,
                                       highlight = BSC_HIGHLIGHT, hover = BSC_NODE_FLUID)) |>
    visNetwork::visGroups(groupname = "solid", shape = "box", borderWidth = 1.5,
                          margin = 8, shapeProperties = list(borderRadius = 2, borderDashes = c(5, 4)),
                          font = list(color = BSC_NODE_SOLID$font, size = 14, face = font),
                          color = list(background = BSC_NODE_SOLID$background,
                                       border = BSC_NODE_SOLID$border,
                                       highlight = BSC_HIGHLIGHT, hover = BSC_NODE_SOLID)) |>
    visNetwork::visEdges(font = list(size = 12, color = BSC_INK$secondary, face = font,
                                     strokeWidth = 4, strokeColor = "#ffffff", align = "middle"),
                         arrows = list(to = list(scaleFactor = 0.8)),
                         smooth = list(enabled = TRUE, roundness = 0.15),
                         selectionWidth = 2, hoverWidth = 0.5) |>
    visNetwork::visPhysics(enabled = FALSE) |>
    visNetwork::visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE,
                               hover = TRUE, tooltipDelay = 300, multiselect = FALSE) |>
    visNetwork::visEvents(select = js_select, dragEnd = js_drag)
  if (editable) {
    # vis.js manipulation mode.  The add callbacks only forward the request
    # (position / end points and the element type chosen in the toolbar) to
    # the server as input$<id>_graphChange; the server creates the element in
    # the project and pushes it to the diagram.  Deletion uses visNetwork's
    # default handler (confirmation + graphChange "deleteElements").
    js_add_node <- htmlwidgets::JS(sprintf(
      "function(data, callback) {
         var el = document.getElementById('graph%1$s');
         var type = (el && el.bscNodeType) ? el.bscNodeType : 'fluid';
         Shiny.setInputValue('%1$s_graphChange',
           {cmd: 'addNode', x: data.x, y: data.y, type: type, nonce: Math.random()}, {priority: 'event'});
         callback(null);
         if (window.bscEditorDone) { bscEditorDone('%1$s'); }
       }", output_id))
    js_add_edge <- htmlwidgets::JS(sprintf(
      "function(data, callback) {
         var el = document.getElementById('graph%1$s');
         var type = (el && el.bscEdgeType) ? el.bscEdgeType : 'flow';
         if (data.from !== data.to) {
           Shiny.setInputValue('%1$s_graphChange',
             {cmd: 'addEdge', from: data.from, to: data.to, type: type, nonce: Math.random()}, {priority: 'event'});
         }
         callback(null);
         if (window.bscEditorDone) { bscEditorDone('%1$s'); }
       }", output_id))
    w <- visNetwork::visOptions(
      w, manipulation = list(enabled = TRUE, initiallyActive = FALSE,
                             addNode = js_add_node, addEdge = js_add_edge, editEdge = FALSE,
                             deleteNode = TRUE, deleteEdge = TRUE))
  }
  w
}

#' Apply the difference between two diagram data sets through a visNetwork
#' proxy (keeps zoom, pan and selection instead of re-drawing the widget)
bsc_network_sync <- function(proxy, old, new) {
  row_keys <- function(df) {
    if (!nrow(df)) return(character(0))
    df <- df[, sort(names(df)), drop = FALSE]
    do.call(paste, c(lapply(df, as.character), sep = "\r"))
  }
  diff_rows <- function(o, n) {
    if (!nrow(n)) return(n)
    if (!nrow(o) || !identical(sort(names(o)), sort(names(n)))) return(n)
    n[!row_keys(n) %in% row_keys(o), , drop = FALSE]
  }
  rm_nodes <- setdiff(old$nodes$id, new$nodes$id)
  rm_edges <- setdiff(old$edges$id, new$edges$id)
  if (length(rm_edges)) visNetwork::visRemoveEdges(proxy, id = rm_edges)
  if (length(rm_nodes)) visNetwork::visRemoveNodes(proxy, id = rm_nodes)
  up_nodes <- diff_rows(old$nodes, new$nodes)
  if (nrow(up_nodes)) visNetwork::visUpdateNodes(proxy, nodes = up_nodes)
  up_edges <- diff_rows(old$edges, new$edges)
  if (nrow(up_edges)) visNetwork::visUpdateEdges(proxy, edges = up_edges)
  invisible(list(removed_nodes = rm_nodes, removed_edges = rm_edges,
                 updated_nodes = up_nodes$id, updated_edges = up_edges$id))
}

#' Toolbar of the drag-and-drop editor
bsc_editor_toolbar <- function(output_id) {
  btn <- function(label, icon, js, cls = "btn-outline-primary", title = NULL) {
    shiny::tags$button(type = "button", class = paste("btn btn-sm", cls), onclick = js,
                       title = title, shiny::icon(icon), label)
  }
  shiny::div(
    class = "bsc-toolbar d-flex flex-wrap gap-2 align-items-center",
    btn("Compartment", "square", sprintf("bscEditorMode('%s','node','fluid')", output_id),
        title = "Add a fluid compartment: click, then click in the diagram"),
    btn("Solid compartment", "square-minus", sprintf("bscEditorMode('%s','node','solid')", output_id),
        title = "Add a solid compartment (chip material, non-specific binding)"),
    btn("Fluid flow", "arrow-right-long", sprintf("bscEditorMode('%s','edge','flow')", output_id),
        title = "Add a flow: click, then drag from the source to the target compartment"),
    btn("Permeable surface", "grip-lines-vertical", sprintf("bscEditorMode('%s','edge','perm')", output_id),
        title = "Add a permeable surface: click, then drag between two compartments"),
    btn("Remove selected", "trash", sprintf("bscEditorMode('%s','delete')", output_id),
        cls = "btn-outline-danger"),
    btn("", "expand", sprintf("bscEditorMode('%s','fit')", output_id), cls = "btn-outline-secondary",
        title = "Fit the diagram to the view")
  )
}

#' Legend shown under the diagram
bsc_network_legend <- function() {
  item <- function(cls, label) {
    shiny::span(class = "bsc-legend-item", shiny::span(class = cls), label)
  }
  shiny::div(
    class = "bsc-legend",
    item("bsc-swatch bsc-swatch-fluid", "Fluid compartment"),
    item("bsc-swatch bsc-swatch-solid", "Solid compartment (NSB)"),
    item("bsc-line bsc-line-flow", "Fluid flow"),
    item("bsc-line bsc-line-perm", "Permeable surface")
  )
}

# ---- plots -----------------------------------------------------------------

#' Indices to keep when drawing long profiles (keeps every interval start)
bsc_thin_index <- function(time, segment = NULL, max_points = 2500) {
  n <- length(time)
  if (n <= max_points) return(seq_len(n))
  k <- ceiling(n / max_points)
  idx <- seq(1, n, by = k)
  if (!is.null(segment)) {
    starts <- which(c(TRUE, diff(segment) != 0))
    idx <- c(idx, starts, pmax(starts - 1L, 1L))
  }
  sort(unique(c(idx, n)))
}

#' Common plotly layout of concentration-time plots
bsc_plotly_layout <- function(fig, ylog = FALSE, xlab = "Time (h)", ylab = "Concentration (µM)") {
  axis <- list(gridcolor = BSC_INK$grid, gridwidth = 1, zeroline = FALSE,
               linecolor = BSC_INK$axis, linewidth = 1, showline = TRUE,
               tickfont = list(color = BSC_INK$secondary, size = 12),
               titlefont = list(color = BSC_INK$primary, size = 13))
  plotly::layout(
    fig,
    xaxis = c(axis, list(title = xlab, rangemode = "tozero")),
    yaxis = c(axis, list(title = ylab, type = if (ylog) "log" else "linear",
                         exponentformat = "power")),
    hovermode = "x unified",
    legend = list(orientation = "h", y = -0.2, font = list(color = BSC_INK$secondary, size = 12)),
    font = list(family = "system-ui, -apple-system, 'Segoe UI', sans-serif"),
    plot_bgcolor = "#ffffff", paper_bgcolor = "#ffffff",
    margin = list(l = 60, r = 20, t = 20, b = 50)
  ) |>
    plotly::config(displaylogo = FALSE,
                   modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"),
                   toImageButtonOptions = list(format = "png", scale = 3,
                                               filename = "BaySimChip_profiles"))
}

#' Empty plot with a message
bsc_plotly_message <- function(msg) {
  plotly::plot_ly(type = "scatter", mode = "markers", x = numeric(0), y = numeric(0)) |>
    plotly::layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
                   annotations = list(list(text = msg, showarrow = FALSE, xref = "paper",
                                           yref = "paper", x = 0.5, y = 0.5,
                                           font = list(color = BSC_INK$muted, size = 14))),
                   plot_bgcolor = "#ffffff", paper_bgcolor = "#ffffff") |>
    plotly::config(displayModeBar = FALSE)
}

#' Profiles on a coarse grid for the table view of a plot
bsc_profile_table <- function(sim, compartments, step = 1) {
  if (is.null(sim)) return(data.frame())
  grid <- seq(0, max(sim$time), by = step)
  idx <- vapply(grid, function(g) which.min(abs(sim$time - g)), integer(1))
  out <- data.frame(`Time (h)` = round(sim$time[idx], 4), check.names = FALSE)
  for (cc in intersect(compartments, colnames(sim$conc))) {
    out[[paste0(cc, " (µM)")]] <- signif(sim$conc[idx, cc], 5)
  }
  out
}

# ---- echo guard for forms bound to the selected element / loaded project --------------
# The server pushes values into a form with update*Input(); the browser then reports the
# new value back as an ordinary input event (an "echo"). If the selection - or the whole
# project - changes before that echo arrives, a binding that writes input values into the
# currently selected element would copy the previous element's values into the new one.
# bsc_echo_guard() remembers the values pushed recently, per input id:
#   guard$note(id, value)     call whenever the server pushes a value to the form;
#   guard$is_echo(id, value)  TRUE if an incoming value only echoes such a push (the event
#                             is then ignored). When the echo of the latest push arrives,
#                             browser and server are in sync again and the memory for that
#                             input is cleared, so later edits by the user always go through.
bsc_echo_guard <- function(max_age = 5) {
  pushed <- new.env(parent = emptyenv())
  prune <- function(id) {
    lst <- pushed[[id]]
    if (is.null(lst)) return(list())
    now <- Sys.time()
    lst <- Filter(function(z) as.numeric(difftime(now, z$t, units = "secs")) < max_age, lst)
    pushed[[id]] <- lst
    lst
  }
  list(
    note = function(id, value) {
      lst <- prune(id)
      pushed[[id]] <- c(lst, list(list(value = value, t = Sys.time())))
      invisible()
    },
    is_echo = function(id, value) {
      lst <- prune(id)
      if (!length(lst)) return(FALSE)
      hit <- which(vapply(lst, function(z) .bsc_same_value(z$value, value), logical(1)))
      if (!length(hit)) return(FALSE)
      if (max(hit) == length(lst)) pushed[[id]] <- list()
      TRUE
    }
  )
}

.bsc_same_value <- function(a, b) {
  if (is.null(a)) a <- character(0)
  if (is.null(b)) b <- character(0)
  if (length(a) != 1 || length(b) != 1) {            # sets, e.g. a multiple selectize input
    return(identical(sort(as.character(a)), sort(as.character(b))))
  }
  if (is.logical(a) || is.logical(b)) return(identical(isTRUE(a), isTRUE(b)))
  na <- suppressWarnings(as.numeric(a))
  nb <- suppressWarnings(as.numeric(b))
  if (is.na(na) || is.na(nb)) return(identical(as.character(a), as.character(b)))
  isTRUE(all.equal(na, nb))
}
