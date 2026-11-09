# =============================================================================
# engine_fit.R - parameter estimation
#
# Port of AnalysisODE() of the paper's analysis code (paper_code/R/engine.R): the parameters flagged as
# unknown are estimated by nonlinear least squares with the
# Levenberg-Marquardt algorithm (minpack.lm::nls.lm).  Residuals are
#   model prediction - observed concentration
# (unweighted).  Standard errors come from the Jacobian at the solution
# (summary.nls.lm, as in the analysis code); 95 % confidence intervals:
# estimate +/- t(0.975, n - p) * SE.
#
# Observation lookup ("obs_mode"):
#   "grid"  (default) - identical to the analysis code: the ODE is solved on a
#           0.1 h grid per interval and each observation is compared with the
#           nearest grid point (an observation at an intervention time is
#           compared with the value just *before* the intervention).
#   "exact" - the ODE solution is evaluated exactly at the sampling times
#           (same before-intervention convention).
# =============================================================================

#' Unknown parameters with initial values and bounds
#'
#' Initial value = current parameter value (values <= 0 or missing are
#' replaced by 1, the default).  Default bounds are those of the
#' analysis code: 1e-4 x and 1000 x the initial value; user bounds are taken
#' from p$fit$bounds.
#'
#' @return data.frame(key, kind, element, parameter, label, unit, initial,
#'   lower, upper)
bsc_fit_parameters <- function(p) {
  pt <- bsc_parameter_table(p)
  u <- pt[pt$estimate, , drop = FALSE]
  if (!nrow(u)) {
    return(data.frame(key = character(0), kind = character(0), element = character(0),
                      parameter = character(0), label = character(0), unit = character(0),
                      initial = numeric(0), lower = numeric(0), upper = numeric(0),
                      stringsAsFactors = FALSE))
  }
  initial <- u$value
  initial[!is.finite(initial) | initial <= 0] <- 1
  lower <- 1e-4 * initial
  upper <- 1000 * initial
  b <- p$fit$bounds
  if (!is.null(b) && nrow(b)) {
    bk <- bsc_param_key(b$kind, b$element, b$parameter)
    m <- match(u$key, bk)
    has <- !is.na(m)
    lb <- b$lower[m[has]]
    ub <- b$upper[m[has]]
    lower[has][is.finite(lb)] <- lb[is.finite(lb)]
    upper[has][is.finite(ub)] <- ub[is.finite(ub)]
  }
  data.frame(key = u$key, kind = u$kind, element = u$element, parameter = u$parameter,
             label = u$label, unit = u$unit, initial = initial, lower = lower, upper = upper,
             stringsAsFactors = FALSE)
}

#' Observed data prepared for fitting (ordered by compartment, then file order)
.bsc_fit_observations <- function(p) {
  obs <- p$observed
  cp <- p$compartments$name
  obs <- obs[obs$compartment %in% cp & is.finite(obs$time) & is.finite(obs$conc), , drop = FALSE]
  ci <- match(obs$compartment, cp)
  obs <- obs[order(ci, seq_len(nrow(obs))), , drop = FALSE]
  obs$comp_index <- match(obs$compartment, cp)
  rownames(obs) <- NULL
  obs
}

#' Model predictions at the observation times
.bsc_predict_obs <- function(sol, obs, mode) {
  if (mode == "grid") {
    # nearest output time; which.min returns the first index, i.e. the end of
    # the previous interval when the observation is at an intervention time
    idx <- vapply(obs$time, function(tp) which.min(abs(sol$time - tp)), integer(1))
    sol$conc[cbind(idx, obs$comp_index)]
  } else {
    idx <- match(obs$time, sol$time)
    sol$conc[cbind(idx, obs$comp_index)]
  }
}

#' Integrate for fitting purposes (grid or exact observation lookup)
.bsc_integrate_for_fit <- function(sys, events, obs, mode, dt, method, rtol, atol) {
  if (mode == "grid") {
    return(bsc_integrate(sys$A, sys$v0, sys$kev, events, grid = "fit", dt = dt,
                         method = method, rtol = rtol, atol = atol, comp_names = sys$names))
  }
  # exact: integrate interval by interval, keeping start, sampling times and
  # end of every interval; the first match of an observation time is its
  # pre-intervention value
  bsc_integrate(sys$A, sys$v0, sys$kev, events, grid = "exact", sample_times = unique(obs$time),
                method = method, rtol = rtol, atol = atol, comp_names = sys$names)
}

#' Estimate the unknown parameters of a project from its observed data
#'
#' @param p project with parameters flagged as unknown (`*_est` columns) and
#'   observed data (p$observed)
#' @param obs_mode "grid" (as the analysis code) or "exact" (see file header)
#' @param final_time end time of the model run (default: the larger of the
#'   simulation time + 0.1 and the last observation + 1 h)
#' @param maxiter maximum Levenberg-Marquardt iterations (nls.lm default 50)
#' @param conf_level confidence level of the intervals
#' @param progress optional function(n_eval) called after every model run
#' @return object of class "bsc_fit" (see README)
bsc_fit <- function(p, obs_mode = NULL, final_time = NULL, maxiter = NULL,
                    method = "lsoda", rtol = 1e-6, atol = 1e-6, dt_fit = 0.1,
                    conf_level = 0.95, progress = NULL) {
  obs_mode <- match.arg(obs_mode %||% p$fit$obs_mode %||% "grid", c("grid", "exact"))
  maxiter <- as.integer(maxiter %||% p$fit$maxiter %||% 50L)
  v <- bsc_validate(p, "fit")
  if (length(v$errors)) stop(paste(v$errors, collapse = "\n"), call. = FALSE)

  par_tab <- bsc_fit_parameters(p)
  obs <- .bsc_fit_observations(p)
  if (is.null(final_time)) {
    final_time <- max(p$settings$sim_time + 0.1, max(obs$time) + 1)
  }
  # initial guesses are written into the project so that the events (e.g.
  # first dose from data) use the same volumes as the analysis code
  p0 <- bsc_apply_estimates(p, par_tab$initial)
  sys0 <- bsc_assemble_system(p0)
  events <- bsc_build_events(p0, final_time, v0 = sys0$v0)

  n_eval <- 0L
  big <- 1e6 * (1 + abs(obs$conc))
  residuals_function <- function(r) {
    n_eval <<- n_eval + 1L
    if (is.function(progress)) progress(n_eval)
    res <- tryCatch({
      sys <- bsc_assemble_system(p0, r)
      sol <- .bsc_integrate_for_fit(sys, events, obs, obs_mode, dt_fit, method, rtol, atol)
      .bsc_predict_obs(sol, obs, obs_mode) - obs$conc
    }, error = function(e) big)
    if (any(!is.finite(res))) res[!is.finite(res)] <- big[!is.finite(res)]
    res
  }

  r0 <- par_tab$initial
  lower <- par_tab$lower
  upper <- par_tab$upper
  r0 <- pmin(pmax(r0, lower), upper)
  t_start <- Sys.time()
  fit <- minpack.lm::nls.lm(par = r0, fn = residuals_function, lower = lower, upper = upper,
                            control = minpack.lm::nls.lm.control(maxiter = maxiter))
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  est <- fit$par
  np <- length(est)
  nobs <- length(fit$fvec)
  df <- nobs - np
  ssr <- sum(fit$fvec^2)
  # standard errors exactly as the analysis code: coef(summary(nls.lm))
  sm <- tryCatch(summary(fit), error = function(e) NULL)
  se <- if (is.null(sm)) rep(NA_real_, np) else unname(stats::coef(sm)[, 2])
  covm <- tryCatch(chol2inv(chol(fit$hessian)) * ssr / df, error = function(e) NULL)
  corm <- if (is.null(covm)) NULL else stats::cov2cor(covm)
  if (!is.null(corm)) dimnames(corm) <- list(par_tab$label, par_tab$label)
  tq <- if (df > 0) stats::qt(1 - (1 - conf_level) / 2, df) else NA_real_
  at_bound <- abs(est - lower) <= 1e-8 * pmax(abs(lower), 1e-12) |
    abs(est - upper) <= 1e-8 * pmax(abs(upper), 1e-12)

  estimates <- data.frame(
    key = par_tab$key, kind = par_tab$kind, element = par_tab$element,
    parameter = par_tab$parameter, label = par_tab$label, unit = par_tab$unit,
    initial = r0, estimate = est, se = se, rse_percent = 100 * se / est,
    ci_lower = est - tq * se, ci_upper = est + tq * se,
    lower = lower, upper = upper, at_bound = at_bound,
    stringsAsFactors = FALSE)

  # predictions at the observations with the final estimates
  sys <- bsc_assemble_system(p0, est)
  sol <- .bsc_integrate_for_fit(sys, events, obs, obs_mode, dt_fit, method, rtol, atol)
  pred <- .bsc_predict_obs(sol, obs, obs_mode)
  resid_df <- data.frame(compartment = obs$compartment, time = obs$time, observed = obs$conc,
                         predicted = pred, residual = pred - obs$conc, stringsAsFactors = FALSE)

  # fitted profiles on the simulation grid
  p_fit <- bsc_apply_estimates(p, est)
  sim <- tryCatch(bsc_simulate(p_fit, final_time = final_time, method = method,
                               rtol = rtol, atol = atol, validate = FALSE),
                  error = function(e) NULL)

  structure(list(
    estimates = estimates,
    n_obs = nobs, n_par = np, df = df, ssr = ssr,
    sigma = if (df > 0) sqrt(ssr / df) else NA_real_,
    aic = nobs * log(ssr / nobs) + 2 * np,
    correlation = corm,
    iterations = fit$niter, n_eval = n_eval,
    info = fit$info, message = fit$message,
    converged = fit$info %in% 1:4,
    residuals = resid_df,
    simulation = sim,
    obs_mode = obs_mode, final_time = final_time, conf_level = conf_level,
    elapsed = elapsed,
    nls = fit
  ), class = "bsc_fit")
}

#' Write estimated values into the project (the "unknown" flags are kept)
#'
#' @param est a "bsc_fit" object or a numeric vector in the order of the
#'   unknown parameters
bsc_apply_estimates <- function(p, est) {
  if (inherits(est, "bsc_fit")) est <- est$estimates$estimate
  pt <- bsc_parameter_table(p)
  u <- pt[pt$estimate, , drop = FALSE]
  if (length(est) != nrow(u)) stop("Number of estimates does not match the unknown parameters.",
                                   call. = FALSE)
  for (i in seq_len(nrow(u))) {
    p <- bsc_set_param(p, u$kind[i], u$element[i], u$parameter[i], value = est[i])
  }
  p
}
