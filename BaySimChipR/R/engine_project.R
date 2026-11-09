# =============================================================================
# engine_project.R - BaySimChip project data structure
#
# A *project* is a plain R list (no environments, no globals) that holds the
# model structure, the experimental conditions, observed data and fit
# settings.  All engine functions take a project as input and return new
# values; nothing is modified in place.
#
# Units (same as the analysis code of the paper):
#   time h, volume mL, concentration uM (= nmol/mL), CL / Q / PS mL/h,
#   evaporation rate Kev mL/h, dose in `settings$dose_unit` (default "mg",
#   converted to nmol with the compound molar mass in g/mol).
#
# Structure
#   name, description          character
#   settings                   list(molar_mass, dose_unit, sim_time)
#   compartments               data.frame(name, type ("fluid"/"solid"),
#                                volume, volume_est, fu, fu_est, cl, cl_est,
#                                kev, kev_est, x, y)
#   flows                      data.frame(name, from, to, q, q_est)
#   perms                      data.frame(name, between1, between2, ps, ps_est)
#   dosing                     named list (by compartment) of
#                                list(mode = "single"|"repeated"|"custom",
#                                     single, repeated_dose, repeated_interval,
#                                     custom = data.frame(time, dose),
#                                     first_dose_from_data)
#   exchange                   named list (by compartment) of
#                                data.frame(time, removed, added)
#   observed                   data.frame(compartment, time, conc)
#   fit                        list(bounds = data.frame(kind, element,
#                                parameter, lower, upper), obs_mode, maxiter)
#   results (optional)         stored results, e.g. the last fit (only kept
#                                in the native .rds format)
#
# A parameter whose `*_est` flag is TRUE is "unknown": it is estimated by the
# fitting module and its current value is used as the initial guess.
#
# Files: the native format is .rds (bsc_write_project / bsc_read_project);
# JSON is offered as an interoperable export.
#
# The engine (R/engine_*.R) has no Shiny dependency and can be used from
# analysis scripts: source("tools/load_engine.R") - see README.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Shallow list update: elements of `val` replace those of `x` (unlike
# utils::modifyList this never merges nested lists/data.frames; NULLs ignored)
.bsc_merge <- function(x, val) {
  for (nm in names(val)) if (!is.null(val[[nm]])) x[[nm]] <- val[[nm]]
  x
}

#' Dose units understood by the engine and their factor to "ng"
#' (amount in nmol = dose * factor / molar mass; "nmol" needs no molar mass).
BSC_DOSE_UNITS <- c(mg = 1e6, ug = 1e3, ng = 1, nmol = NA_real_)

BSC_DOSE_UNIT_LABELS <- c(mg = "mg", ug = "µg", ng = "ng", nmol = "nmol")

BSC_COMPARTMENT_TYPES <- c("fluid", "solid")

BSC_DOSE_MODES <- c("single", "repeated", "custom")

# ---------------------------------------------------------------------------
# Empty tables (templates used for construction and JSON parsing)
# ---------------------------------------------------------------------------
bsc_empty_compartments <- function() {
  data.frame(name = character(0), type = character(0),
             volume = numeric(0), volume_est = logical(0),
             fu = numeric(0), fu_est = logical(0),
             cl = numeric(0), cl_est = logical(0),
             kev = numeric(0), kev_est = logical(0),
             x = numeric(0), y = numeric(0),
             stringsAsFactors = FALSE)
}

bsc_empty_flows <- function() {
  data.frame(name = character(0), from = character(0), to = character(0),
             q = numeric(0), q_est = logical(0), stringsAsFactors = FALSE)
}

bsc_empty_perms <- function() {
  data.frame(name = character(0), between1 = character(0),
             between2 = character(0), ps = numeric(0), ps_est = logical(0),
             stringsAsFactors = FALSE)
}

bsc_empty_observed <- function() {
  data.frame(compartment = character(0), time = numeric(0), conc = numeric(0),
             stringsAsFactors = FALSE)
}

bsc_empty_exchange <- function() {
  data.frame(time = numeric(0), removed = numeric(0), added = numeric(0))
}

bsc_empty_custom_dose <- function() {
  data.frame(time = numeric(0), dose = numeric(0))
}

bsc_empty_bounds <- function() {
  data.frame(kind = character(0), element = character(0),
             parameter = character(0), lower = numeric(0), upper = numeric(0),
             stringsAsFactors = FALSE)
}

#' Default dosing specification of a compartment
bsc_default_dosing <- function() {
  list(mode = "single", single = 0, repeated_dose = 0, repeated_interval = 24,
       custom = bsc_empty_custom_dose(), first_dose_from_data = FALSE)
}

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

#' Create an empty project
bsc_new_project <- function(name = "New model", description = "") {
  list(
    format = "BaySimChipR project",
    format_version = 1L,
    name = name,
    description = description,
    settings = list(molar_mass = 300, dose_unit = "mg", sim_time = 72),
    compartments = bsc_empty_compartments(),
    flows = bsc_empty_flows(),
    perms = bsc_empty_perms(),
    dosing = list(),
    exchange = list(),
    observed = bsc_empty_observed(),
    fit = list(bounds = bsc_empty_bounds(), obs_mode = "grid", maxiter = 50L)
  )
}

# Next free default name "<prefix><k>"
.bsc_next_name <- function(existing, prefix) {
  k <- 0L
  while (paste0(prefix, k) %in% existing) k <- k + 1L
  paste0(prefix, k)
}

#' Add a compartment
#'
#' @param type "fluid" (medium, organ) or "solid" (chip material for
#'   non-specific binding; its `fu` holds the partition coefficient K_media/plastic).
#' @return the modified project
bsc_add_compartment <- function(p, name = NULL, type = c("fluid", "solid"),
                                volume = 1, fu = 1, cl = 0, kev = 0,
                                x = NA_real_, y = NA_real_) {
  type <- match.arg(type)
  if (is.null(name) || !nzchar(trimws(name))) {
    name <- .bsc_next_name(c(p$compartments$name, p$flows$name, p$perms$name), "Comp")
  }
  name <- trimws(name)
  if (name %in% p$compartments$name) {
    stop("A compartment called '", name, "' already exists.", call. = FALSE)
  }
  row <- data.frame(name = name, type = type,
                    volume = as.numeric(volume), volume_est = FALSE,
                    fu = as.numeric(fu), fu_est = FALSE,
                    cl = as.numeric(cl), cl_est = FALSE,
                    kev = as.numeric(kev), kev_est = FALSE,
                    x = as.numeric(x), y = as.numeric(y),
                    stringsAsFactors = FALSE)
  p$compartments <- rbind(p$compartments, row)
  p$dosing[[name]] <- bsc_default_dosing()
  p$exchange[[name]] <- bsc_empty_exchange()
  p
}

#' Add a fluid flow (from -> to). Flows can only connect fluid compartments.
bsc_add_flow <- function(p, from, to, q = 1, name = NULL) {
  .bsc_check_comp(p, c(from, to))
  if (identical(from, to)) stop("A flow needs two different compartments.", call. = FALSE)
  types <- p$compartments$type[match(c(from, to), p$compartments$name)]
  if (any(types != "fluid")) {
    stop("Flows can only connect fluid (non-solid) compartments.", call. = FALSE)
  }
  if (is.null(name) || !nzchar(trimws(name))) {
    name <- .bsc_next_name(p$flows$name, "Flow")
  }
  if (name %in% p$flows$name) stop("A flow called '", name, "' already exists.", call. = FALSE)
  p$flows <- rbind(p$flows, data.frame(name = name, from = from, to = to,
                                       q = as.numeric(q), q_est = FALSE,
                                       stringsAsFactors = FALSE))
  p
}

#' Add a permeable surface between two compartments (fluid or solid)
bsc_add_perm <- function(p, between1, between2, ps = 0.01, name = NULL) {
  .bsc_check_comp(p, c(between1, between2))
  if (identical(between1, between2)) {
    stop("A permeable surface needs two different compartments.", call. = FALSE)
  }
  if (is.null(name) || !nzchar(trimws(name))) {
    name <- .bsc_next_name(p$perms$name, "P.S.")
  }
  if (name %in% p$perms$name) {
    stop("A permeable surface called '", name, "' already exists.", call. = FALSE)
  }
  p$perms <- rbind(p$perms, data.frame(name = name, between1 = between1,
                                       between2 = between2, ps = as.numeric(ps),
                                       ps_est = FALSE, stringsAsFactors = FALSE))
  p
}

.bsc_check_comp <- function(p, names) {
  bad <- setdiff(names, p$compartments$name)
  if (length(bad)) stop("Unknown compartment(s): ", paste(bad, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' Remove compartments together with their flows, permeable surfaces,
#' dosing / exchange settings and observed data
bsc_remove_compartments <- function(p, names) {
  names <- intersect(names, p$compartments$name)
  if (!length(names)) return(p)
  p$compartments <- p$compartments[!p$compartments$name %in% names, , drop = FALSE]
  p$flows <- p$flows[!(p$flows$from %in% names | p$flows$to %in% names), , drop = FALSE]
  p$perms <- p$perms[!(p$perms$between1 %in% names | p$perms$between2 %in% names), , drop = FALSE]
  p$dosing <- p$dosing[setdiff(names(p$dosing), names)]
  p$exchange <- p$exchange[setdiff(names(p$exchange), names)]
  p$observed <- p$observed[!p$observed$compartment %in% names, , drop = FALSE]
  p <- .bsc_prune_bounds(p)
  rownames(p$compartments) <- NULL
  p
}

bsc_remove_flows <- function(p, names) {
  p$flows <- p$flows[!p$flows$name %in% names, , drop = FALSE]
  rownames(p$flows) <- NULL
  .bsc_prune_bounds(p)
}

bsc_remove_perms <- function(p, names) {
  p$perms <- p$perms[!p$perms$name %in% names, , drop = FALSE]
  rownames(p$perms) <- NULL
  .bsc_prune_bounds(p)
}

# drop bounds that refer to elements that no longer exist
.bsc_prune_bounds <- function(p) {
  b <- p$fit$bounds
  if (is.null(b) || !nrow(b)) return(p)
  keep <- (b$kind == "compartment" & b$element %in% p$compartments$name) |
    (b$kind == "flow" & b$element %in% p$flows$name) |
    (b$kind == "perm" & b$element %in% p$perms$name)
  p$fit$bounds <- b[keep, , drop = FALSE]
  p
}

#' Rename a compartment, flow or permeable surface (references are updated)
bsc_rename <- function(p, kind = c("compartment", "flow", "perm"), old, new) {
  kind <- match.arg(kind)
  new <- trimws(new)
  if (!nzchar(new)) stop("Names cannot be empty.", call. = FALSE)
  if (identical(old, new)) return(p)
  tbl <- switch(kind, compartment = "compartments", flow = "flows", perm = "perms")
  if (!old %in% p[[tbl]]$name) stop("'", old, "' does not exist.", call. = FALSE)
  if (new %in% p[[tbl]]$name) stop("The name '", new, "' is already used.", call. = FALSE)
  p[[tbl]]$name[p[[tbl]]$name == old] <- new
  if (kind == "compartment") {
    p$flows$from[p$flows$from == old] <- new
    p$flows$to[p$flows$to == old] <- new
    p$perms$between1[p$perms$between1 == old] <- new
    p$perms$between2[p$perms$between2 == old] <- new
    names(p$dosing)[names(p$dosing) == old] <- new
    names(p$exchange)[names(p$exchange) == old] <- new
    p$observed$compartment[p$observed$compartment == old] <- new
  }
  b <- p$fit$bounds
  if (!is.null(b) && nrow(b)) {
    b$element[b$kind == kind & b$element == old] <- new
    p$fit$bounds <- b
  }
  p
}

#' Change the type of a compartment ("fluid" or "solid")
#'
#' Solid compartments represent chip material (non-specific binding); they
#' cannot be connected by fluid flows, have no clearance and are not dosed.
bsc_set_compartment_type <- function(p, name, type) {
  type <- tolower(trimws(type))
  if (!type %in% BSC_COMPARTMENT_TYPES) stop("Type must be 'fluid' or 'solid'.", call. = FALSE)
  i <- match(name, p$compartments$name)
  if (is.na(i)) stop("Unknown compartment '", name, "'.", call. = FALSE)
  if (type == "solid" && any(p$flows$from == name | p$flows$to == name)) {
    stop("'", name, "' is connected by fluid flows; remove them before making it solid.",
         call. = FALSE)
  }
  p$compartments$type[i] <- type
  if (type == "solid") {
    p$compartments$cl[i] <- 0
    p$compartments$cl_est[i] <- FALSE
  }
  p
}

# ---------------------------------------------------------------------------
# Parameter access
# ---------------------------------------------------------------------------

# parameter name -> (table, value column, flag column)
.bsc_param_columns <- function(kind, parameter) {
  ok <- switch(kind,
               compartment = parameter %in% c("volume", "fu", "cl", "kev"),
               flow = parameter == "q",
               perm = parameter == "ps",
               FALSE)
  if (!isTRUE(ok)) stop("Unknown parameter '", parameter, "' for a ", kind, ".", call. = FALSE)
  list(table = switch(kind, compartment = "compartments", flow = "flows", perm = "perms"),
       value = parameter, flag = paste0(parameter, "_est"))
}

#' Set the value and/or the "estimate" flag of a model parameter
#'
#' @param kind "compartment", "flow" or "perm"
#' @param element element name
#' @param parameter "volume", "fu", "cl", "kev" (compartments), "q" (flows),
#'   "ps" (permeable surfaces)
bsc_set_param <- function(p, kind, element, parameter, value = NULL, estimate = NULL) {
  cols <- .bsc_param_columns(kind, parameter)
  tbl <- p[[cols$table]]
  i <- match(element, tbl$name)
  if (is.na(i)) stop("Unknown ", kind, " '", element, "'.", call. = FALSE)
  if (!is.null(value)) tbl[[cols$value]][i] <- as.numeric(value)
  if (!is.null(estimate)) tbl[[cols$flag]][i] <- isTRUE(as.logical(estimate))
  p[[cols$table]] <- tbl
  p
}

bsc_get_param <- function(p, kind, element, parameter) {
  cols <- .bsc_param_columns(kind, parameter)
  tbl <- p[[cols$table]]
  i <- match(element, tbl$name)
  if (is.na(i)) stop("Unknown ", kind, " '", element, "'.", call. = FALSE)
  list(value = tbl[[cols$value]][i], estimate = isTRUE(tbl[[cols$flag]][i]))
}

#' Dosing specification of one compartment
bsc_get_dosing <- function(p, comp) {
  d <- p$dosing[[comp]]
  if (is.null(d)) d <- bsc_default_dosing()
  d
}

bsc_set_dosing <- function(p, comp, dosing) {
  .bsc_check_comp(p, comp)
  d <- .bsc_merge(bsc_get_dosing(p, comp), dosing)
  if (!d$mode %in% BSC_DOSE_MODES) stop("Unknown dosing mode '", d$mode, "'.", call. = FALSE)
  d$custom <- .bsc_as_df(d$custom, bsc_empty_custom_dose())
  p$dosing[[comp]] <- d
  p
}

#' Media exchange / sampling table of one compartment
#' (time h, removed volume mL, added volume mL)
bsc_get_exchange <- function(p, comp) {
  e <- p$exchange[[comp]]
  if (is.null(e)) e <- bsc_empty_exchange()
  e
}

bsc_set_exchange <- function(p, comp, exchange) {
  .bsc_check_comp(p, comp)
  p$exchange[[comp]] <- .bsc_as_df(exchange, bsc_empty_exchange())
  p
}

# Coerce a data.frame / matrix / list of rows to the column layout of `template`
.bsc_as_df <- function(x, template) {
  cols <- names(template)
  if (is.null(x) || length(x) == 0) return(template)
  if (is.matrix(x)) {
    x <- as.data.frame(x)
    names(x)[seq_len(min(ncol(x), length(cols)))] <- cols[seq_len(min(ncol(x), length(cols)))]
  }
  if (is.list(x) && !is.data.frame(x)) {
    # list of records (as read from JSON)
    x <- do.call(rbind, lapply(x, function(r) {
      as.data.frame(lapply(stats::setNames(cols, cols), function(cn) {
        v <- r[[cn]]
        if (is.null(v)) NA else v
      }), stringsAsFactors = FALSE)
    }))
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (cn in cols) {
    if (!cn %in% names(x)) x[[cn]] <- rep(NA, nrow(x))
    proto <- template[[cn]]
    if (is.numeric(proto)) x[[cn]] <- suppressWarnings(as.numeric(x[[cn]]))
    if (is.logical(proto)) x[[cn]] <- as.logical(x[[cn]])
    if (is.character(proto)) x[[cn]] <- as.character(x[[cn]])
  }
  x <- x[, cols, drop = FALSE]
  rownames(x) <- NULL
  x
}

# ---------------------------------------------------------------------------
# Doses
# ---------------------------------------------------------------------------

#' Convert doses to amounts in nmol (amount / volume in mL = concentration uM)
#'
#' mg: 1e6 * dose / molar mass; ng: dose / molar mass (the unit of the
#' analysis code).
bsc_dose_to_nmol <- function(dose, unit = "mg", molar_mass = 300) {
  unit <- match.arg(unit, names(BSC_DOSE_UNITS))
  if (unit == "nmol") return(dose)
  f <- BSC_DOSE_UNITS[[unit]]
  if (f == 1) dose / molar_mass else f * dose / molar_mass
}

#' Inverse of [bsc_dose_to_nmol()]
bsc_nmol_to_dose <- function(amount, unit = "mg", molar_mass = 300) {
  unit <- match.arg(unit, names(BSC_DOSE_UNITS))
  if (unit == "nmol") return(amount)
  amount * molar_mass / BSC_DOSE_UNITS[[unit]]
}

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

#' Check a project for problems
#'
#' @param purpose "structure" (only model structure), "simulate" (all values
#'   needed) or "fit" (values or initial guesses needed + observed data)
#' @return list(errors = character, warnings = character)
bsc_validate <- function(p, purpose = c("simulate", "structure", "fit")) {
  purpose <- match.arg(purpose)
  err <- character(0)
  wrn <- character(0)
  cp <- p$compartments
  if (!nrow(cp)) err <- c(err, "The model has no compartments.")
  if (anyDuplicated(cp$name)) err <- c(err, "Compartment names must be unique.")
  if (any(!nzchar(cp$name))) err <- c(err, "Compartment names cannot be empty.")
  if (any(!cp$type %in% BSC_COMPARTMENT_TYPES)) {
    err <- c(err, "Compartment type must be 'fluid' or 'solid'.")
  }
  fl <- p$flows
  if (nrow(fl)) {
    if (anyDuplicated(fl$name)) err <- c(err, "Flow names must be unique.")
    bad <- !(fl$from %in% cp$name & fl$to %in% cp$name)
    if (any(bad)) err <- c(err, paste0("Flow '", fl$name[bad], "' refers to an unknown compartment."))
    same <- fl$from == fl$to
    if (any(same)) err <- c(err, paste0("Flow '", fl$name[same], "' starts and ends in the same compartment."))
    tp <- cp$type[match(c(fl$from, fl$to), cp$name)]
    if (any(tp %in% "solid")) err <- c(err, "Flows cannot involve solid compartments.")
  }
  pm <- p$perms
  if (nrow(pm)) {
    if (anyDuplicated(pm$name)) err <- c(err, "Permeable surface names must be unique.")
    bad <- !(pm$between1 %in% cp$name & pm$between2 %in% cp$name)
    if (any(bad)) err <- c(err, paste0("Permeable surface '", pm$name[bad], "' refers to an unknown compartment."))
  }
  if (purpose %in% c("simulate", "fit") && nrow(cp)) {
    pt <- bsc_parameter_table(p)
    # for fitting, unknown parameters may be empty (initial guess 1 is used)
    known <- if (purpose == "fit") !pt$estimate else rep(TRUE, nrow(pt))
    miss <- pt[known & (is.na(pt$value) | !is.finite(pt$value)), , drop = FALSE]
    if (nrow(miss)) {
      err <- c(err, paste0(miss$label, " needs a value",
                           ifelse(miss$estimate, " (it is marked as unknown; its value is the initial guess).", ".")))
    }
    neg <- pt[!is.na(pt$value) & pt$value < 0, , drop = FALSE]
    if (nrow(neg)) err <- c(err, paste0(neg$label, " cannot be negative."))
    v0 <- pt[known & pt$parameter == "volume" & !is.na(pt$value) & pt$value <= 0, , drop = FALSE]
    if (nrow(v0)) err <- c(err, paste0(v0$label, " must be > 0."))
    for (cn in cp$name) {
      d <- bsc_get_dosing(p, cn)
      if (d$mode == "repeated" && !(is.finite(d$repeated_interval) && d$repeated_interval > 0)) {
        err <- c(err, paste0("Dosing interval in ", cn, " must be > 0."))
      }
      if (d$mode == "custom" && nrow(d$custom) &&
          any(!is.finite(d$custom$time) | !is.finite(d$custom$dose))) {
        err <- c(err, paste0("Custom dosing table of ", cn, " contains empty cells."))
      }
      ex <- bsc_get_exchange(p, cn)
      if (nrow(ex) && any(!is.finite(as.matrix(ex)))) {
        err <- c(err, paste0("Media exchange table of ", cn, " contains empty cells."))
      }
      if (nrow(ex) && any(ex$time < 0, na.rm = TRUE)) {
        err <- c(err, paste0("Media exchange times in ", cn, " must be >= 0."))
      }
    }
    if (!p$settings$dose_unit %in% names(BSC_DOSE_UNITS)) err <- c(err, "Unknown dose unit.")
    if (p$settings$dose_unit != "nmol" &&
        !(is.finite(p$settings$molar_mass) && p$settings$molar_mass > 0)) {
      err <- c(err, "The compound molar mass must be > 0.")
    }
  }
  if (purpose == "fit") {
    pt <- bsc_parameter_table(p)
    if (!any(pt$estimate)) {
      err <- c(err, "No parameter is marked as unknown. Mark at least one parameter for estimation.")
    }
    obs <- p$observed
    obs <- obs[obs$compartment %in% cp$name & is.finite(obs$time) & is.finite(obs$conc), , drop = FALSE]
    if (!nrow(obs)) err <- c(err, "No observed data. Load concentration-time data first.")
    if (nrow(obs) && nrow(obs) <= sum(pt$estimate)) {
      err <- c(err, "Fewer observations than unknown parameters.")
    }
  }
  list(errors = unique(err), warnings = unique(wrn))
}

# ---------------------------------------------------------------------------
# Saving and loading (native .rds, JSON export)
# ---------------------------------------------------------------------------

#' Serialise a project to JSON text
bsc_project_to_json <- function(p, pretty = TRUE) {
  p <- bsc_normalize_project(p)
  comps <- lapply(seq_len(nrow(p$compartments)), function(i) {
    r <- as.list(p$compartments[i, , drop = FALSE])
    d <- bsc_get_dosing(p, r$name)
    r$dosing <- d
    r$exchange <- bsc_get_exchange(p, r$name)
    r
  })
  out <- list(
    format = "BaySimChipR project",
    format_version = 1L,
    name = p$name,
    description = p$description,
    settings = p$settings,
    compartments = comps,
    flows = p$flows,
    perms = p$perms,
    observed = p$observed,
    fit = p$fit
  )
  jsonlite::toJSON(out, auto_unbox = TRUE, pretty = pretty, digits = NA,
                   na = "null", null = "null", dataframe = "rows")
}

#' Parse a project from JSON text
bsc_project_from_json <- function(txt) {
  j <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
  if (!identical(j$format, "BaySimChipR project")) {
    stop("This file is not a BaySimChipR project.", call. = FALSE)
  }
  p <- bsc_new_project(name = j$name %||% "Imported model", description = j$description %||% "")
  s <- j$settings %||% list()
  p$settings$molar_mass <- as.numeric(s$molar_mass %||% 300)
  p$settings$dose_unit <- as.character(s$dose_unit %||% "mg")
  p$settings$sim_time <- as.numeric(s$sim_time %||% 72)
  comps <- j$compartments %||% list()
  p$compartments <- .bsc_as_df(comps, bsc_empty_compartments())
  for (r in comps) {
    d <- r$dosing %||% list()
    dd <- bsc_default_dosing()
    if (!is.null(d$mode)) dd$mode <- as.character(d$mode)
    for (f in c("single", "repeated_dose", "repeated_interval")) {
      if (!is.null(d[[f]])) dd[[f]] <- as.numeric(d[[f]])
    }
    dd$custom <- .bsc_as_df(d$custom, bsc_empty_custom_dose())
    dd$first_dose_from_data <- isTRUE(d$first_dose_from_data)
    p$dosing[[r$name]] <- dd
    p$exchange[[r$name]] <- .bsc_as_df(r$exchange, bsc_empty_exchange())
  }
  p$flows <- .bsc_as_df(j$flows, bsc_empty_flows())
  p$perms <- .bsc_as_df(j$perms, bsc_empty_perms())
  p$observed <- .bsc_as_df(j$observed, bsc_empty_observed())
  f <- j$fit %||% list()
  p$fit$bounds <- .bsc_as_df(f$bounds, bsc_empty_bounds())
  p$fit$obs_mode <- as.character(f$obs_mode %||% "grid")
  p$fit$maxiter <- as.integer(f$maxiter %||% 50L)
  bsc_normalize_project(p)
}

#' Make sure all fields exist and have the expected types
bsc_normalize_project <- function(p) {
  tmpl <- bsc_new_project()
  for (nm in names(tmpl)) if (is.null(p[[nm]])) p[[nm]] <- tmpl[[nm]]
  p$settings <- .bsc_merge(tmpl$settings, p$settings)
  p$compartments <- .bsc_as_df(p$compartments, bsc_empty_compartments())
  for (f in c("volume_est", "fu_est", "cl_est", "kev_est")) {
    p$compartments[[f]][is.na(p$compartments[[f]])] <- FALSE
  }
  p$compartments$type[is.na(p$compartments$type)] <- "fluid"
  p$flows <- .bsc_as_df(p$flows, bsc_empty_flows())
  p$flows$q_est[is.na(p$flows$q_est)] <- FALSE
  p$perms <- .bsc_as_df(p$perms, bsc_empty_perms())
  p$perms$ps_est[is.na(p$perms$ps_est)] <- FALSE
  p$observed <- .bsc_as_df(p$observed, bsc_empty_observed())
  for (cn in p$compartments$name) {
    p$dosing[[cn]] <- .bsc_merge(bsc_default_dosing(), p$dosing[[cn]] %||% list())
    p$dosing[[cn]]$custom <- .bsc_as_df(p$dosing[[cn]]$custom, bsc_empty_custom_dose())
    p$exchange[[cn]] <- .bsc_as_df(p$exchange[[cn]], bsc_empty_exchange())
  }
  p$dosing <- p$dosing[p$compartments$name]
  p$exchange <- p$exchange[p$compartments$name]
  if (!length(p$compartments$name)) {
    p$dosing <- list()
    p$exchange <- list()
  }
  p$fit <- .bsc_merge(tmpl$fit, p$fit %||% list())
  p$fit$bounds <- .bsc_as_df(p$fit$bounds, bsc_empty_bounds())
  p
}

#' Write a project to a file
#'
#' The native format is ".rds" (R data; keeps everything incl. stored
#' results).  ".json" is an interoperable export of the model, conditions,
#' observed data and fit settings (stored results are not exported).
#'
#' @param format "rds" or "json"; default taken from the file extension
#'   ("rds" if there is none)
bsc_write_project <- function(p, path, format = NULL) {
  if (is.null(format)) {
    ext <- tolower(tools::file_ext(path))
    format <- if (ext == "json") "json" else "rds"
  }
  format <- match.arg(format, c("rds", "json"))
  if (format == "rds") {
    saveRDS(bsc_normalize_project(p), path)
  } else {
    writeLines(bsc_project_to_json(p), path, useBytes = TRUE)
  }
  invisible(path)
}

#' Read a project from ".rds" (native) or ".json"
#'
#' @param name original file name (e.g. of a Shiny upload) used to detect
#'   the format; defaults to `path`
bsc_read_project <- function(path, name = NULL) {
  ext <- tolower(tools::file_ext(if (is.null(name)) path else name))
  switch(ext,
         rds = {
           p <- readRDS(path)
           if (!is.list(p) || !identical(p$format, "BaySimChipR project")) {
             stop("This .rds file does not contain a BaySimChipR project.", call. = FALSE)
           }
           bsc_normalize_project(p)
         },
         json = bsc_project_from_json(paste(readLines(path, warn = FALSE, encoding = "UTF-8"),
                                            collapse = "\n")),
         stop("Unsupported file type '.", ext, "' (use .rds or .json).", call. = FALSE))
}
