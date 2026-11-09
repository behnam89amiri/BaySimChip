# =============================================================================
# engine_batch.R - helpers for batch use of the engine (analysis scripts)
#
# The engine functions work on plain project lists and need no Shiny, so the
# analyses of the paper (parameter sweeps, repeated fits, ...) can call them
# directly:
#
#   source("BaySimChipR/tools/load_engine.R")
#   p   <- bsc_read_project("my_chip.rds")
#   sc  <- data.frame(scenario = c("with evaporation", "without"),
#                     "Apical|kev" = c(0.003, 0), check.names = FALSE)
#   res <- bsc_run_scenarios(p, sc, final_time = 110)   # long data.frame
#   attr(res, "exposure")                                # AUC, Cmax, ...
# =============================================================================

#' Row of the parameter table addressed by a key
#'
#' Accepts the full key "kind|element|parameter" (e.g. "compartment|Liver|cl")
#' or the short form "element|parameter" (e.g. "Liver|cl", "P.S.0|ps").
.bsc_resolve_key <- function(pt, key) {
  i <- match(key, pt$key)
  if (!is.na(i)) return(i)
  i <- which(paste(pt$element, pt$parameter, sep = "|") == key)
  if (length(i) == 1) return(i)
  if (!length(i)) stop("Unknown parameter '", key, "'. Use keys like 'Liver|cl' or ",
                       "'compartment|Liver|cl' (see bsc_parameter_table()).", call. = FALSE)
  stop("Parameter '", key, "' is ambiguous; use the full key 'kind|element|parameter'.",
       call. = FALSE)
}

#' Set several parameter values (and optionally the "unknown" flags) at once
#'
#' @param values named list or vector, names = parameter keys (see
#'   `.bsc_resolve_key`)
#' @param estimate optional logical: mark these parameters as unknown (TRUE)
#'   or known (FALSE)
bsc_set_params <- function(p, values, estimate = NULL) {
  pt <- bsc_parameter_table(p)
  for (nm in names(values)) {
    i <- .bsc_resolve_key(pt, nm)
    p <- bsc_set_param(p, pt$kind[i], pt$element[i], pt$parameter[i],
                       value = values[[nm]], estimate = estimate)
  }
  p
}

#' Simulate several scenarios (variants) of a project
#'
#' @param p base project
#' @param scenarios either a data.frame with one row per scenario, an optional
#'   column `scenario` (labels) and one column per parameter key to change; or
#'   a named list of functions function(project) -> project
#' @param compartments compartments to return (default all)
#' @param ... passed to [bsc_simulate()] (final_time, method, rtol, ...)
#' @return long data.frame(scenario, time, compartment, conc) with the
#'   attributes "exposure" (exposure parameters per scenario and compartment)
#'   and "simulations" (the bsc_simulation objects)
bsc_run_scenarios <- function(p, scenarios, compartments = NULL, ...) {
  if (is.data.frame(scenarios)) {
    lab <- if ("scenario" %in% names(scenarios)) as.character(scenarios$scenario) else
      paste("Scenario", seq_len(nrow(scenarios)))
    pars <- setdiff(names(scenarios), "scenario")
    mods <- lapply(seq_len(nrow(scenarios)), function(i) {
      vals <- as.list(scenarios[i, pars, drop = FALSE])
      function(q) bsc_set_params(q, vals)
    })
    names(mods) <- lab
  } else if (is.list(scenarios) && length(scenarios) &&
             all(vapply(scenarios, is.function, logical(1)))) {
    mods <- scenarios
    if (is.null(names(mods))) names(mods) <- paste("Scenario", seq_along(mods))
  } else {
    stop("'scenarios' must be a data.frame or a list of functions.", call. = FALSE)
  }
  if (anyDuplicated(names(mods))) stop("Scenario labels must be unique.", call. = FALSE)
  sims <- lapply(mods, function(f) bsc_simulate(f(p), ...))
  long <- do.call(rbind, lapply(names(sims), function(nm) {
    cbind(scenario = nm, bsc_simulation_long(sims[[nm]], compartments), stringsAsFactors = FALSE)
  }))
  expo <- do.call(rbind, lapply(names(sims), function(nm) {
    s <- sims[[nm]]
    cc <- if (is.null(compartments)) colnames(s$conc) else intersect(compartments, colnames(s$conc))
    cbind(scenario = nm, bsc_exposure(s$time, s$conc[, cc, drop = FALSE]), stringsAsFactors = FALSE)
  }))
  rownames(long) <- NULL
  attr(long, "exposure") <- expo
  attr(long, "simulations") <- sims
  long
}
