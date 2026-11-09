# =============================================================================
# engine_solver.R - piecewise ODE integration, simulation, exposure parameters
#
# Port of CalCulate() / ODEmodel() (simulation) and RtoODE() in AnalysisODE()
# (fitting) of the analysis code (paper_code/R/engine.R).  The system is integrated with deSolve
# piecewise between intervention times.  Within an interval the volumes
# decrease linearly by evaporation, V(t) = V_start - Kev (t - t_start), and
# dC/dt = (A / V(t)) C.  At each intervention time the state is updated:
#
#   V_new = V_end - removed + added
#   C_new = (C_end (V_end - removed) + dose) / (V_end - removed + added)
#
# Doses/removals/additions at the first intervention time (t = 0) only enter
# through the initial concentration dose / V0 (as in the original code).
# =============================================================================

#' ODE right-hand side within one interval (identical arithmetic to ODEmodel() of the analysis code)
.bsc_rhs <- function(t, C, parms) {
  vVector <- parms$v_start - parms$kev * (t - parms$t_start)
  VMatrix <- matrix(vVector, nrow = parms$n, ncol = parms$n)
  list((parms$A / VMatrix) %*% C)
}

#' Output times of one interval
#'
#' grid = "sim": every `dt` h from the interval start to end - dt, plus the end
#'   point (the end point is not kept in the output) - as CalCulate().
#' grid = "fit": every `dt` h from start to end, both ends kept - RtoODE().
#' grid = "exact": start, end and the sampling times inside the interval.
.bsc_interval_times <- function(a, b, grid, dt, sample_times = NULL) {
  switch(grid,
         sim = {
           to <- b - dt
           if (to < a) unique(c(a, b)) else c(seq(from = a, to = to, by = dt), b)
         },
         fit = {
           tt <- seq(a, b, by = dt)
           if (length(tt) < 2) c(a, b) else tt
         },
         exact = sort(unique(c(a, sample_times[sample_times > a & sample_times < b], b))))
}

#' Integrate the linear ODE system piecewise between intervention times
#'
#' @param A coefficient matrix from [bsc_assemble_system()]
#' @param v0,kev initial volumes and evaporation rates
#' @param events schedule from [bsc_build_events()]
#' @param grid "sim", "fit" or "exact" (see `.bsc_interval_times`)
#' @param dt output step (default 0.01 h for "sim", 0.1 h for "fit")
#' @param sample_times sampling times for grid = "exact"
#' @param method,rtol,atol passed to [deSolve::ode()]; the defaults are the
#'   deSolve defaults used by the analysis code (lsoda, 1e-6, 1e-6).
#' @return list(time, conc (matrix), volume (matrix), segment, keep_end,
#'   removed_amount (matrix n x K, nmol removed at each intervention))
bsc_integrate <- function(A, v0, kev, events, grid = c("sim", "fit", "exact"), dt = NULL,
                          sample_times = NULL, method = "lsoda", rtol = 1e-6, atol = 1e-6,
                          comp_names = NULL) {
  grid <- match.arg(grid)
  if (is.null(dt)) dt <- if (grid == "sim") 0.01 else 0.1
  n <- length(v0)
  tk <- events$times
  K <- length(tk)
  if (is.null(comp_names)) comp_names <- paste0("C", seq_len(n))
  if (any(!is.finite(v0)) || any(v0 <= 0)) {
    stop("All initial volumes must be > 0.", call. = FALSE)
  }
  InitialConcentration <- events$doses[, 1] / v0
  Volumes <- matrix(v0, nrow = n, ncol = K)
  removed_amount <- matrix(0, n, K)
  out_t <- vector("list", K - 1)
  out_c <- vector("list", K - 1)
  out_v <- vector("list", K - 1)
  out_s <- vector("list", K - 1)

  for (iDose in seq_len(K - 1)) {
    a <- tk[iDose]
    b <- tk[iDose + 1]
    InitVolumesofInterval <- Volumes[, iDose]
    EndVolumesofInterval <- InitVolumesofInterval - kev * (b - a)
    if (any(EndVolumesofInterval <= 0)) {
      j <- which(EndVolumesofInterval <= 0)[1]
      stop(sprintf("The volume of %s drops to zero before t = %g h (evaporation too large?).",
                   comp_names[j], b), call. = FALSE)
    }
    times <- .bsc_interval_times(a, b, grid, dt, sample_times)
    parms <- list(A = A, v_start = InitVolumesofInterval, kev = kev, t_start = a, n = n)
    if (grid == "exact") {
      # deSolve limits the step size to the largest gap between output times;
      # with sparse sampling times this would allow very large steps and make
      # finite-difference Jacobians noisy, so cap it at the fit grid step.
      sol <- deSolve::ode(y = InitialConcentration, times = times, func = .bsc_rhs,
                          parms = parms, method = method, rtol = rtol, atol = atol,
                          hmax = dt)
    } else {
      sol <- deSolve::ode(y = InitialConcentration, times = times, func = .bsc_rhs,
                          parms = parms, method = method, rtol = rtol, atol = atol)
    }
    if (nrow(sol) < length(times) || any(!is.finite(sol[, -1]))) {
      stop(sprintf("The ODE solver failed in the interval %g-%g h.", a, b), call. = FALSE)
    }
    EndConcentration <- sol[nrow(sol), -1]

    # state update at the next intervention time
    rem <- events$removed[, iDose + 1]
    add <- events$added[, iDose + 1]
    remaining <- EndVolumesofInterval - rem
    if (any(remaining + add <= 0)) {
      j <- which(remaining + add <= 0)[1]
      stop(sprintf("The volume of %s is <= 0 after the media exchange at t = %g h.",
                   comp_names[j], b), call. = FALSE)
    }
    removed_amount[, iDose + 1] <- EndConcentration * rem
    Volumes[, iDose + 1] <- EndVolumesofInterval - rem + add
    InitialConcentration <- (EndConcentration * remaining + events$doses[, iDose + 1]) /
      (EndVolumesofInterval - rem + add)

    keep <- if (grid == "sim") seq_len(nrow(sol) - 1) else seq_len(nrow(sol))
    tt <- sol[keep, 1]
    out_t[[iDose]] <- tt
    out_c[[iDose]] <- sol[keep, -1, drop = FALSE]
    # volumes at the output times: V_start - Kev (t - t_start)
    out_v[[iDose]] <- matrix(InitVolumesofInterval, nrow = length(tt), ncol = n, byrow = TRUE) -
      outer(tt - a, kev)
    out_s[[iDose]] <- rep(iDose, length(tt))
  }
  conc <- do.call(rbind, out_c)
  vol <- do.call(rbind, out_v)
  colnames(conc) <- comp_names
  colnames(vol) <- comp_names
  rownames(conc) <- NULL
  list(time = unlist(out_t, use.names = FALSE), conc = conc, volume = vol,
       segment = unlist(out_s, use.names = FALSE), volumes_at_events = Volumes,
       removed_amount = removed_amount)
}

#' Simulate concentration-time profiles of a project
#'
#' @param p project
#' @param final_time end of the simulation (h).  Default: settings$sim_time +
#'   0.1, so that the requested end time is part of the 0.01 h output grid.
#' @param theta optional values for the unknown parameters (see
#'   [bsc_assemble_system()])
#' @param dt output step (h)
#' @param method,rtol,atol ODE solver settings (deSolve)
#' @return object of class "bsc_simulation": list(time, conc, volume, segment,
#'   events, system, final_time)
bsc_simulate <- function(p, final_time = NULL, theta = NULL, dt = 0.01,
                         method = "lsoda", rtol = 1e-6, atol = 1e-6, validate = TRUE) {
  if (is.null(final_time)) final_time <- p$settings$sim_time + 0.1
  if (!(is.finite(final_time) && final_time > 0)) stop("The simulation time must be > 0.", call. = FALSE)
  if (validate) {
    v <- bsc_validate(p, if (is.null(theta)) "simulate" else "fit")
    errs <- v$errors[!grepl("^No (parameter|observed)|Fewer observations", v$errors)]
    if (length(errs)) stop(paste(errs, collapse = "\n"), call. = FALSE)
  }
  sys <- bsc_assemble_system(p, theta)
  ev <- bsc_build_events(p, final_time, v0 = sys$v0)
  out <- bsc_integrate(sys$A, sys$v0, sys$kev, ev, grid = "sim", dt = dt,
                       method = method, rtol = rtol, atol = atol, comp_names = sys$names)
  structure(c(out, list(events = ev, system = sys, final_time = final_time,
                        solver = list(method = method, rtol = rtol, atol = atol, dt = dt))),
            class = "bsc_simulation")
}

#' Simulation result as a long data.frame (time, compartment, conc)
bsc_simulation_long <- function(sim, compartments = NULL) {
  cn <- colnames(sim$conc)
  if (is.null(compartments)) compartments <- cn
  compartments <- intersect(compartments, cn)
  do.call(rbind, lapply(compartments, function(cc) {
    data.frame(time = sim$time, compartment = cc, conc = sim$conc[, cc],
               stringsAsFactors = FALSE)
  }))
}

#' Total amount of compound (nmol) in all compartments over time
bsc_total_amount <- function(sim) rowSums(sim$conc * sim$volume)

#' Exposure parameters (AUC, Cmax, Tmax, half-life) of simulated profiles
#'
#' Definitions: AUC by the
#' trapezoidal rule over the simulated grid, Cmax/Tmax at the global maximum,
#' half-life = time from Tmax until the concentration first falls to Cmax/2
#' (linear interpolation; NA if this does not happen within the simulation).
#'
#' @param time numeric vector
#' @param conc matrix (time x compartments) with column names
#' @param window optional c(from, to) restricting the analysis to a time window
#' @return data.frame(compartment, AUC, Cmax, Tmax, t_half, status)
bsc_exposure <- function(time, conc, window = NULL) {
  conc <- as.matrix(conc)
  if (is.null(colnames(conc))) colnames(conc) <- paste0("C", seq_len(ncol(conc)))
  sel <- rep(TRUE, length(time))
  if (!is.null(window) && length(window) == 2 && all(is.finite(window))) {
    sel <- time >= window[1] & time <= window[2]
  }
  tt <- time[sel]
  res <- lapply(colnames(conc), function(cc) {
    y <- conc[sel, cc]
    if (length(y) < 2) {
      return(data.frame(compartment = cc, AUC = NA_real_, Cmax = NA_real_, Tmax = NA_real_,
                        t_half = NA_real_, status = "too few points", stringsAsFactors = FALSE))
    }
    auc <- sum(diff(tt) * (y[-1] + y[-length(y)]) / 2)
    imax <- which.max(y)
    cmax <- y[imax]
    tmax <- tt[imax]
    mono <- all(diff(y) >= 0)
    status <- if (imax == length(y) && mono) {
      "Tmax not reached"
    } else if (imax == length(y)) {
      "Tmax is a local maximum"
    } else {
      "Tmax reached and passed"
    }
    th <- NA_real_
    if (imax < length(y) && cmax > 0) {
      half <- cmax / 2
      post <- (imax + 1):length(y)
      k <- which(y[post] <= half)[1]
      if (!is.na(k)) {
        j2 <- post[k]
        j1 <- j2 - 1
        tc <- if (y[j1] == y[j2]) tt[j2] else tt[j1] + (half - y[j1]) * (tt[j2] - tt[j1]) / (y[j2] - y[j1])
        th <- tc - tmax
      }
    }
    data.frame(compartment = cc, AUC = auc, Cmax = cmax, Tmax = tmax, t_half = th,
               status = status, stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}
