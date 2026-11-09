# =============================================================================
# engine_model.R - parameters, ODE coefficient matrix and event schedule
#
# Ports of ConstructODE() / AnalysisODE() (matrix assembly) and
# ConstructDoseMediaEx() (intervention schedule) of the paper's analysis code
# (paper_code/R/engine.R),
# rewritten as pure functions.
#
# Model (per compartment i, concentration C_i in uM, volume V_i(t) in mL):
#
#   d(V_i C_i)/dt = sum_j Q_ji C_j - sum_j Q_ij C_i                  (flows)
#                 + sum_j PS_ij (fu_j C_j - fu_i C_i)                  (perm. surfaces)
#                 - CL_i fu_i C_i                                      (clearance)
#   dV_i/dt       = -Kev_i                                             (evaporation)
#
# which gives  dC/dt = (A / V(t)) C  with
#   A = (-diag(CL) + PERM) .* FU + FLOW + diag(Kev),   FU[i, j] = fu_j
# and V(t) = V_start - Kev (t - t_start) within each interval between two
# interventions (doses, media exchange, sampling).
# =============================================================================

#' Table of all model parameters in the canonical order
#'
#' The order is the one used by the analysis code for the
#' vector of unknown parameters: for every compartment (in table order) volume,
#' fu, CL, Kev; then every flow Q; then every permeable surface PS.
#'
#' @return data.frame(key, kind, element, parameter, label, unit, value, estimate)
bsc_parameter_table <- function(p) {
  cp <- p$compartments
  rows <- list()
  for (i in seq_len(nrow(cp))) {
    solid <- identical(cp$type[i], "solid")
    for (par in c("volume", "fu", "cl", "kev")) {
      rows[[length(rows) + 1]] <- data.frame(
        kind = "compartment", element = cp$name[i], parameter = par,
        label = .bsc_param_label(par, cp$name[i], solid),
        unit = .bsc_param_unit(par),
        value = cp[[par]][i], estimate = isTRUE(cp[[paste0(par, "_est")]][i]),
        stringsAsFactors = FALSE)
    }
  }
  fl <- p$flows
  for (i in seq_len(nrow(fl))) {
    rows[[length(rows) + 1]] <- data.frame(
      kind = "flow", element = fl$name[i], parameter = "q",
      label = paste0("Flow rate from ", fl$from[i], " to ", fl$to[i]),
      unit = "mL/h", value = fl$q[i], estimate = isTRUE(fl$q_est[i]),
      stringsAsFactors = FALSE)
  }
  pm <- p$perms
  for (i in seq_len(nrow(pm))) {
    rows[[length(rows) + 1]] <- data.frame(
      kind = "perm", element = pm$name[i], parameter = "ps",
      label = if (bsc_is_nsb_surface(p, pm$name[i])) {
        paste0("NSB exchange rate PS_NSB between ", pm$between1[i], " & ", pm$between2[i])
      } else {
        paste0("Permeability * surface area (PS) between ", pm$between1[i], " & ", pm$between2[i])
      },
      unit = "mL/h", value = pm$ps[i], estimate = isTRUE(pm$ps_est[i]),
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) {
    return(data.frame(key = character(0), kind = character(0), element = character(0),
                      parameter = character(0), label = character(0), unit = character(0),
                      value = numeric(0), estimate = logical(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  out$key <- bsc_param_key(out$kind, out$element, out$parameter)
  out <- out[, c("key", "kind", "element", "parameter", "label", "unit", "value", "estimate")]
  rownames(out) <- NULL
  out
}

#' Unique key of a parameter, e.g. "compartment|Liver|cl"
bsc_param_key <- function(kind, element, parameter) paste(kind, element, parameter, sep = "|")

.bsc_param_label <- function(par, comp, solid = FALSE) {
  switch(par,
         volume = paste("Volume in", comp),
         fu = if (solid) paste("Partition coefficient K_media/plastic of", comp) else paste("fu in", comp),
         cl = paste("Clearance in", comp),
         kev = paste("Evaporation rate in", comp))
}

#' TRUE if a permeable surface connects a compartment with a solid compartment,
#' i.e. it describes exchange with the chip material (non-specific binding, NSB).
#' Its PS is then the NSB exchange rate PS_NSB; the partition coefficient
#' K_media/plastic is stored as the "fu" of the solid compartment
#' (equilibrium: C_medium = K_media/plastic * C_solid).
bsc_is_nsb_surface <- function(p, name) {
  pm <- p$perms[p$perms$name == name, , drop = FALSE]
  if (!nrow(pm)) return(FALSE)
  solids <- p$compartments$name[p$compartments$type == "solid"]
  any(c(pm$between1[1], pm$between2[1]) %in% solids)
}

.bsc_param_unit <- function(par) {
  switch(par, volume = "mL", fu = "-", cl = "mL/h", kev = "mL/h", q = "mL/h", ps = "mL/h")
}

#' Assemble the linear ODE system for given parameter values
#'
#' @param p project
#' @param theta optional numeric vector with values for the parameters flagged
#'   as unknown (in the order of [bsc_parameter_table()]); if NULL the stored
#'   values (initial guesses for unknowns) are used.
#' @return list(A = coefficient matrix, v0, kev, fu, cl, names, n)
bsc_assemble_system <- function(p, theta = NULL) {
  pt <- bsc_parameter_table(p)
  if (!is.null(theta)) {
    idx <- which(pt$estimate)
    if (length(theta) != length(idx)) {
      stop("theta has ", length(theta), " values but ", length(idx),
           " parameters are marked as unknown.", call. = FALSE)
    }
    pt$value[idx] <- unname(theta)
  }
  cp <- p$compartments
  n <- nrow(cp)
  get_comp <- function(par) pt$value[pt$kind == "compartment" & pt$parameter == par]
  v0 <- get_comp("volume")
  fu <- get_comp("fu")
  cl <- get_comp("cl")
  kev <- get_comp("kev")

  # Clearance and fu matrices: FU[, j] = fu_j (AnalysisODE() of the analysis code:
  # Compfu[, i] <- fu_i)
  CompClearance <- matrix(0, n, n)
  diag(CompClearance) <- cl
  Compfu <- matrix(0, n, n)
  for (i in seq_len(n)) Compfu[, i] <- fu[i]

  # Flows: outflux on the diagonal of the source, influx in [to, from]
  ConnectionFlow <- matrix(0, n, n)
  fl <- p$flows
  qv <- pt$value[pt$kind == "flow"]
  for (i in seq_len(nrow(fl))) {
    a <- match(fl$from[i], cp$name)
    b <- match(fl$to[i], cp$name)
    ConnectionFlow[a, a] <- ConnectionFlow[a, a] - qv[i]
    ConnectionFlow[b, a] <- ConnectionFlow[b, a] + qv[i]
  }

  # Permeable surfaces: symmetric exchange between two compartments
  ConnectionPerm <- matrix(0, n, n)
  pm <- p$perms
  psv <- pt$value[pt$kind == "perm"]
  for (i in seq_len(nrow(pm))) {
    a <- match(pm$between1[i], cp$name)
    b <- match(pm$between2[i], cp$name)
    ConnectionPerm[a, a] <- ConnectionPerm[a, a] - psv[i]
    ConnectionPerm[a, b] <- ConnectionPerm[a, b] + psv[i]
    ConnectionPerm[b, b] <- ConnectionPerm[b, b] - psv[i]
    ConnectionPerm[b, a] <- ConnectionPerm[b, a] + psv[i]
  }

  A <- (-CompClearance + ConnectionPerm) * Compfu + ConnectionFlow + diag(kev, nrow = n)
  dimnames(A) <- list(cp$name, cp$name)
  list(A = A, v0 = v0, kev = kev, fu = fu, cl = cl, names = cp$name, n = n)
}

#' First observed concentration of a compartment (mean of the replicates at
#' the earliest sampling time), as used by the dosing option "Use first
#' experimental concentration for dose calculation".
bsc_first_observed_conc <- function(p, comp) {
  obs <- p$observed
  obs <- obs[obs$compartment == comp & is.finite(obs$time) & is.finite(obs$conc), , drop = FALSE]
  if (!nrow(obs)) return(NA_real_)
  t1 <- min(obs$time)
  mean(obs$conc[obs$time == t1])
}

#' Doses of one compartment in the project dose unit
#'
#' @return list(times, doses)
.bsc_dose_events <- function(p, comp, final_time, v0_comp) {
  d <- bsc_get_dosing(p, comp)
  s <- p$settings
  if (isTRUE(d$first_dose_from_data) && d$mode %in% c("single", "repeated")) {
    c1 <- bsc_first_observed_conc(p, comp)
    if (is.finite(c1)) {
      # dose that gives the first observed concentration in the initial volume
      # (analysis code: mean(y[t == 0]) * MolarMass * Vol, dose in ng)
      dose <- bsc_nmol_to_dose(c1 * v0_comp, s$dose_unit, s$molar_mass)
      if (s$dose_unit == "ng") dose <- c1 * s$molar_mass * v0_comp
      if (d$mode == "single") d$single <- dose else d$repeated_dose <- dose
    }
  }
  switch(d$mode,
         single = list(times = 0, doses = d$single),
         repeated = {
           tt <- seq(0, final_time, by = d$repeated_interval)
           list(times = tt, doses = rep(d$repeated_dose, length(tt)))
         },
         custom = {
           cu <- d$custom
           cu <- cu[is.finite(cu$time) & is.finite(cu$dose), , drop = FALSE]
           list(times = cu$time, doses = cu$dose)
         })
}

#' Build the intervention schedule (doses, media exchange, sampling)
#'
#' Port of ConstructDoseMediaEx(): all dose and exchange events of all
#' compartments are collected, sorted (stable) and events at identical times
#' are merged by summing.  The final time is appended (if it coincides with
#' the last event it is moved 1 h later, as in the original code).
#'
#' Deviations (only for inputs the original code could not handle):
#' events after `final_time` are dropped and an empty event at t = 0 is added
#' if no event starts at 0.
#'
#' @param v0 initial volumes (needed for "first dose from data")
#' @return list(times = intervention times incl. final time,
#'   doses = n x K matrix of amounts in nmol, removed, added = n x K volumes)
bsc_build_events <- function(p, final_time, v0 = NULL) {
  cp <- p$compartments
  n <- nrow(cp)
  if (is.null(v0)) v0 <- cp$volume
  times <- numeric(0)
  doses <- matrix(0, n, 0)
  removed <- matrix(0, n, 0)
  added <- matrix(0, n, 0)
  for (i in seq_len(n)) {
    de <- .bsc_dose_events(p, cp$name[i], final_time, v0[i])
    m <- matrix(0, n, length(de$times))
    m[i, ] <- de$doses
    times <- c(times, de$times)
    doses <- cbind(doses, m)
    removed <- cbind(removed, matrix(0, n, length(de$times)))
    added <- cbind(added, matrix(0, n, length(de$times)))

    ex <- bsc_get_exchange(p, cp$name[i])
    ex <- ex[is.finite(ex$time), , drop = FALSE]
    mr <- matrix(0, n, nrow(ex))
    ma <- matrix(0, n, nrow(ex))
    mr[i, ] <- ex$removed
    ma[i, ] <- ex$added
    times <- c(times, ex$time)
    removed <- cbind(removed, mr)
    added <- cbind(added, ma)
    doses <- cbind(doses, matrix(0, n, nrow(ex)))
  }
  keep <- times <= final_time
  times <- times[keep]
  doses <- doses[, keep, drop = FALSE]
  removed <- removed[, keep, drop = FALSE]
  added <- added[, keep, drop = FALSE]
  if (!length(times) || min(times) > 0) {
    times <- c(0, times)
    doses <- cbind(0, doses)
    removed <- cbind(0, removed)
    added <- cbind(0, added)
  }

  # stable sort and merge of simultaneous events
  indx <- order(times)
  st <- times[indx]
  sd <- doses[, indx, drop = FALSE]
  sr <- removed[, indx, drop = FALSE]
  sa <- added[, indx, drop = FALSE]
  mt <- st[1]
  md <- sd[, 1, drop = FALSE]
  mr <- sr[, 1, drop = FALSE]
  ma <- sa[, 1, drop = FALSE]
  if (length(st) > 1) {
    for (jj in 2:length(st)) {
      if (st[jj] - st[jj - 1] != 0) {
        mt <- c(mt, st[jj])
        md <- cbind(md, sd[, jj, drop = FALSE])
        mr <- cbind(mr, sr[, jj, drop = FALSE])
        ma <- cbind(ma, sa[, jj, drop = FALSE])
      } else {
        k <- ncol(md)
        md[, k] <- md[, k] + sd[, jj]
        mr[, k] <- mr[, k] + sr[, jj]
        ma[, k] <- ma[, k] + sa[, jj]
      }
    }
  }
  it <- c(mt, final_time)
  if (it[length(it)] == it[length(it) - 1]) it[length(it)] <- it[length(it)] + 1

  s <- p$settings
  dose_nmol <- cbind(md, matrix(0, nrow = n, ncol = 1))
  dose_nmol <- bsc_dose_to_nmol(dose_nmol, s$dose_unit, s$molar_mass)
  list(times = it,
       doses = dose_nmol,
       removed = cbind(mr, matrix(0, nrow = n, ncol = 1)),
       added = cbind(ma, matrix(0, nrow = n, ncol = 1)))
}
