# =============================================================================
# R/engine.R - BaySimChip R engine used for the analyses of the paper
#
# Model import (load_model), schedule of doses and media exchanges
# (ConstructDoseMediaEx), ODE system (ConstructODE, ODEmodel), simulation
# (CalCulate), parameter estimation (AnalysisODE) and exposure metrics
# (exposure_metrics, ExposureParas).
#
# Units: time h, volume mL, concentration µM (amounts = µM x mL), flows, PS,
# CL and Kev in mL/h.
#
# The model ("app") is a list. Its elements CompData (compartments),
# FlowData (flows) and PermData (permeable surfaces) are lists of
# environments, exactly as created by the model import. Environments are
# reference objects: `b <- a` does not copy them. Use copy_app() whenever an
# independent copy is needed.
#
# Compartment order of the CNBio MPS-TL6 models used in the paper:
#   1 gut apical ("Comp0"), 2 gut basolateral ("Comp1"), 3 liver ("Comp2"),
#   4 circulation ("Comp3"); NSB models add 5 = solid (device material)
#   coupled to the liver/basolateral circuit ("Comp5") and, added by
#   add_apical_nsb(), 6 = solid coupled to the apical compartment ("Comp6").
# PermData[[1]] is the gut barrier (PS); PermData[[2]] / [[3]] are the NSB
# surfaces (PS_NSB). For a solid compartment the field `fu` holds
# K_media/plastic.
#
# Behaviours of the original implementation that are kept (see README.md):
#   * ConstructODE() builds the fu matrix transposed with respect to
#     AnalysisODE() (Compfu[i, j] = fu_i instead of fu_j). Fluid
#     profiles are identical; the simulated concentration of a solid
#     compartment is K_media/plastic x the concentration of the AnalysisODE
#     formulation. With K = 1 (Figs. 6C, 7A) there is no difference at all.
#   * Media exchange / sampling events at t = 0 are ignored (only doses at
#     t = 0 are applied); the initial volumes are the nominal volumes.
#   * DoseCond = 3 (custom dosing) reads the custom dosing table by rows; it
#     is not used in the paper.
# =============================================================================

# ---- model import ----------------------------------------------------------

# Read a BaySimChip model (data/models/*.rds, see data/README.md) into the
# `app` list.
load_model <- function(model_file) {
  readRDS(model_file)
}

# Independent (deep) copy of an app, including its environments
copy_app <- function(app) unserialize(serialize(app, NULL))

# Copy of one environment
clone_env <- function(e) {
  stopifnot(is.environment(e))
  new <- new.env(parent = parent.env(e))
  for (nm in ls(e, all.names = TRUE)) {
    val <- get(nm, envir = e, inherits = FALSE)
    if (is.environment(val)) val <- clone_env(val)
    assign(nm, val, envir = new)
  }
  new
}

# ---- dosing and media-exchange schedule ------------------------------------

# Merge all dose and media-exchange events of all compartments into one sorted
# list of intervention times. Adds
# app$InterventionTimes, app$DosesAll (amounts / MolarMass), app$ExMediaAll and
# app$AddMediaAll (one column per intervention time; FinalTime appended).
ConstructDoseMediaEx <- function(app) {
  TmpInterventionTimes <- c()
  TmpDosesAll <- NULL
  TmpExMediaAll <- NULL
  TmpAddMediaAll <- NULL

  for (i in seq_along(app$CompData)) {
    # dose events of compartment i
    if (app$CompData[[i]]$DoseCond == 1) {
      tmpInterventionTimes <- 0
      tmpDosesAll <- app$CompData[[i]]$DoseSingData
    } else if (app$CompData[[i]]$DoseCond == 2) {
      tmpInterventionTimes <- seq(0, app$FinalTime, by = app$CompData[[i]]$DoseMultData[1, 2])
      tmpDosesAll <- rep(app$CompData[[i]]$DoseMultData[1, 1], length(tmpInterventionTimes))
    } else if (app$CompData[[i]]$DoseCond == 3) {
      tmpInterventionTimes <- app$CompData[[i]]$DoseCustData[1, ]
      tmpDosesAll <- app$CompData[[i]]$DoseCustData[2, ]
    }

    TmpInterventionTimes <- c(TmpInterventionTimes, tmpInterventionTimes)
    ttmpDosesAll <- matrix(0, nrow = app$NumofComp, ncol = length(tmpInterventionTimes))
    ttmpDosesAll[i, ] <- tmpDosesAll
    if (is.null(TmpDosesAll)) {
      TmpDosesAll <- ttmpDosesAll
      TmpExMediaAll <- matrix(0, nrow = app$NumofComp, ncol = length(tmpInterventionTimes))
      TmpAddMediaAll <- matrix(0, nrow = app$NumofComp, ncol = length(tmpInterventionTimes))
    } else {
      TmpDosesAll <- cbind(TmpDosesAll, ttmpDosesAll)
      TmpExMediaAll <- cbind(TmpExMediaAll, matrix(0, nrow = app$NumofComp, ncol = length(tmpInterventionTimes)))
      TmpAddMediaAll <- cbind(TmpAddMediaAll, matrix(0, nrow = app$NumofComp, ncol = length(tmpInterventionTimes)))
    }

    # media exchange / sampling events of compartment i
    tmpInterventionTimes <- app$CompData[[i]]$ExMedia[, 1]
    tmpExMediaAll  <- app$CompData[[i]]$ExMedia[, 2]
    tmpAddMediaAll <- app$CompData[[i]]$ExMedia[, 3]
    ttmpExMediaAll <- matrix(0, nrow = app$NumofComp, ncol = length(tmpInterventionTimes))
    ttmpExMediaAll[i, ] <- tmpExMediaAll
    ttmpAddMediaAll <- matrix(0, nrow = app$NumofComp, ncol = length(tmpInterventionTimes))
    ttmpAddMediaAll[i, ] <- tmpAddMediaAll
    TmpInterventionTimes <- c(TmpInterventionTimes, tmpInterventionTimes)
    if (is.null(TmpExMediaAll)) {
      TmpExMediaAll <- ttmpExMediaAll
      TmpAddMediaAll <- ttmpAddMediaAll
      TmpDosesAll <- cbind(TmpDosesAll, matrix(0, nrow = app$NumofComp, ncol = length(tmpInterventionTimes)))
    } else {
      TmpExMediaAll <- cbind(TmpExMediaAll, ttmpExMediaAll)
      TmpAddMediaAll <- cbind(TmpAddMediaAll, ttmpAddMediaAll)
      TmpDosesAll <- cbind(TmpDosesAll, matrix(0, nrow = app$NumofComp, ncol = length(tmpInterventionTimes)))
    }
  }

  # sort by time and merge events that happen at the same time
  indx <- order(TmpInterventionTimes)
  TTmpInterventionTimes <- TmpInterventionTimes[indx]
  TTmpExMediaAll  <- TmpExMediaAll[, indx, drop = FALSE]
  TTmpAddMediaAll <- TmpAddMediaAll[, indx, drop = FALSE]
  TTmpDosesAll    <- TmpDosesAll[, indx, drop = FALSE]
  TTTmpExMediaAll  <- TTmpExMediaAll[, 1, drop = FALSE]
  TTTmpAddMediaAll <- TTmpAddMediaAll[, 1, drop = FALSE]
  TTTmpDosesAll    <- TTmpDosesAll[, 1, drop = FALSE]
  TTTmpInterventionTimes <- TTmpInterventionTimes[1]

  for (jj in 2:length(TTmpInterventionTimes)) {
    if (TTmpInterventionTimes[jj] - TTmpInterventionTimes[jj - 1] != 0) {
      TTTmpInterventionTimes <- c(TTTmpInterventionTimes, TTmpInterventionTimes[jj])
      TTTmpExMediaAll  <- cbind(TTTmpExMediaAll, TTmpExMediaAll[, jj, drop = FALSE])
      TTTmpAddMediaAll <- cbind(TTTmpAddMediaAll, TTmpAddMediaAll[, jj, drop = FALSE])
      TTTmpDosesAll    <- cbind(TTTmpDosesAll, TTmpDosesAll[, jj, drop = FALSE])
    } else {
      TTTmpExMediaAll[, ncol(TTTmpExMediaAll)]   <- TTTmpExMediaAll[, ncol(TTTmpExMediaAll)] + TTmpExMediaAll[, jj]
      TTTmpAddMediaAll[, ncol(TTTmpAddMediaAll)] <- TTTmpAddMediaAll[, ncol(TTTmpAddMediaAll)] + TTmpAddMediaAll[, jj]
      TTTmpDosesAll[, ncol(TTTmpDosesAll)]       <- TTTmpDosesAll[, ncol(TTTmpDosesAll)] + TTmpDosesAll[, jj]
    }
  }

  app$InterventionTimes <- c(TTTmpInterventionTimes, app$FinalTime)
  if (tail(app$InterventionTimes, n = 1) == app$InterventionTimes[length(app$InterventionTimes) - 1]) {
    app$InterventionTimes[length(app$InterventionTimes)] <- app$InterventionTimes[length(app$InterventionTimes)] + 1
  }
  app$DosesAll    <- cbind(TTTmpDosesAll, matrix(0, nrow = nrow(TTTmpDosesAll), ncol = 1))
  app$DosesAll    <- app$DosesAll / app$MolarMass
  app$ExMediaAll  <- cbind(TTTmpExMediaAll, matrix(0, nrow = nrow(TTTmpExMediaAll), ncol = 1))
  app$AddMediaAll <- cbind(TTTmpAddMediaAll, matrix(0, nrow = nrow(TTTmpAddMediaAll), ncol = 1))
  app
}

# ---- ODE system with known parameters --------------------------------------

# Build the coefficient matrix of the linear ODE system from compartments,
# flows and permeable surfaces. Parameters that are NULL
# are counted as unknown (app$NumofUnknownPara) for AnalysisODE().
ConstructODE <- function(app) {
  app$ConnectionFlow <- matrix(0, nrow = app$NumofComp, ncol = app$NumofComp)
  app$ConnectionPerm <- matrix(0, nrow = app$NumofComp, ncol = app$NumofComp)
  app$ConnectionFlowVector <- matrix(0, nrow = length(app$FlowData), ncol = 3)
  app$ConnectionPermVector <- matrix(0, nrow = length(app$PermData), ncol = 3)
  app$v0Vector  <- numeric(app$NumofComp)
  app$fuVector  <- numeric(app$NumofComp)
  app$KevVector <- numeric(app$NumofComp)
  app$ClVector  <- numeric(app$NumofComp)
  app$CompNameID <- vector("list", app$NumofComp)
  app$NumofUnknownPara <- 0

  for (i in seq_along(app$CompData)) {
    app$CompNameID[[i]] <- app$CompData[[i]]$Name

    if (is.null(app$CompData[[i]]$Vol)) {
      app$NumofUnknownPara <- app$NumofUnknownPara + 1
      if (app$Estimated[[1]]$Done == 1) {
        app$v0Vector[i] <- app$Estimated[[app$NumofUnknownPara]]$Values
        app$Estimated[[app$NumofUnknownPara]]$Names <- paste("Volume in", app$CompData[[i]]$Name)
      }
    } else {
      app$v0Vector[i] <- as.numeric(app$CompData[[i]]$Vol)
    }

    if (is.null(app$CompData[[i]]$fu)) {
      app$NumofUnknownPara <- app$NumofUnknownPara + 1
      if (app$Estimated[[1]]$Done == 1) {
        app$fuVector[i] <- app$Estimated[[app$NumofUnknownPara]]$Values
        app$Estimated[[app$NumofUnknownPara]]$Names <- paste("fu in", app$CompData[[i]]$Name)
      }
    } else {
      app$fuVector[i] <- as.numeric(app$CompData[[i]]$fu)
    }

    if (is.null(app$CompData[[i]]$CL)) {
      app$NumofUnknownPara <- app$NumofUnknownPara + 1
      if (app$Estimated[[1]]$Done == 1) {
        app$ClVector[i] <- app$Estimated[[app$NumofUnknownPara]]$Values
        app$Estimated[[app$NumofUnknownPara]]$Names <- paste("Clearance in", app$CompData[[i]]$Name)
      }
    } else {
      app$ClVector[i] <- as.numeric(app$CompData[[i]]$CL)
    }

    if (is.null(app$CompData[[i]]$Kev)) {
      app$NumofUnknownPara <- app$NumofUnknownPara + 1
      if (app$Estimated[[1]]$Done == 1) {
        app$KevVector[i] <- app$Estimated[[app$NumofUnknownPara]]$Values
        app$Estimated[[app$NumofUnknownPara]]$Names <- paste("Evaporation rate in", app$CompData[[i]]$Name)
      }
    } else {
      app$KevVector[i] <- as.numeric(app$CompData[[i]]$Kev)
    }
  }

  app$CompClearance <- diag(app$ClVector)
  # orientation of the original implementation: Compfu[i, j] = fu_i (see file header)
  app$Compfu <- matrix(rep(app$fuVector, each = app$NumofComp), nrow = app$NumofComp, byrow = TRUE)

  # flows (compartments are matched by name with grep(), i.e. substring matching)
  for (i in seq_along(app$FlowData)) {
    from_value <- as.character(app$FlowData[[i]]$from[[1]][1, 1])
    to_value   <- as.character(app$FlowData[[i]]$to[[1]][1, 1])
    from_idx <- grep(from_value, sapply(app$CompNameID, function(x) x))
    to_idx   <- grep(to_value, sapply(app$CompNameID, function(x) x))
    app$ConnectionFlowVector[i, 1] <- from_idx
    app$ConnectionFlowVector[i, 2] <- to_idx
    if (is.null(app$FlowData[[i]]$Q)) {
      app$NumofUnknownPara <- app$NumofUnknownPara + 1
      if (app$Estimated[[1]]$Done == 1) {
        app$ConnectionFlowVector[i, 3] <- app$Estimated[[app$NumofUnknownPara]]$Values
        app$Estimated[[app$NumofUnknownPara]]$Names <- paste("Flow rate from", from_value, "to", to_value)
      }
    } else {
      app$ConnectionFlowVector[i, 3] <- as.numeric(app$FlowData[[i]]$Q)
    }
  }
  for (i in seq_len(nrow(app$ConnectionFlowVector))) {
    from_idx <- app$ConnectionFlowVector[i, 1]
    to_idx   <- app$ConnectionFlowVector[i, 2]
    flow_value <- app$ConnectionFlowVector[i, 3]
    app$ConnectionFlow[from_idx, from_idx] <- app$ConnectionFlow[from_idx, from_idx] - flow_value  # outflux
    app$ConnectionFlow[to_idx, from_idx]   <- app$ConnectionFlow[to_idx, from_idx] + flow_value    # influx
  }

  # permeable surfaces (gut barrier PS, NSB surfaces PS_NSB)
  for (i in seq_along(app$PermData)) {
    between1_value <- as.character(app$PermData[[i]]$between1[[1]][1, 1])
    between2_value <- as.character(app$PermData[[i]]$between2[[1]][1, 1])
    comp1_idx <- grep(between1_value, sapply(app$CompNameID, function(x) x))
    comp2_idx <- grep(between2_value, sapply(app$CompNameID, function(x) x))
    app$ConnectionPermVector[i, 1] <- comp1_idx
    app$ConnectionPermVector[i, 2] <- comp2_idx
    if (is.null(app$PermData[[i]]$PS)) {
      app$NumofUnknownPara <- app$NumofUnknownPara + 1
      if (app$Estimated[[1]]$Done == 1) {
        app$ConnectionPermVector[i, 3] <- app$Estimated[[app$NumofUnknownPara]]$Values
        app$Estimated[[app$NumofUnknownPara]]$Names <- paste("Permeability * surface area between", between1_value, "&", between2_value)
      }
    } else {
      app$ConnectionPermVector[i, 3] <- as.numeric(app$PermData[[i]]$PS)
    }
  }
  for (i in seq_len(nrow(app$ConnectionPermVector))) {
    idx1 <- app$ConnectionPermVector[i, 1]
    idx2 <- app$ConnectionPermVector[i, 2]
    ps_value <- app$ConnectionPermVector[i, 3]
    app$ConnectionPerm[idx1, idx1] <- app$ConnectionPerm[idx1, idx1] - ps_value
    app$ConnectionPerm[idx1, idx2] <- app$ConnectionPerm[idx1, idx2] + ps_value
    app$ConnectionPerm[idx2, idx2] <- app$ConnectionPerm[idx2, idx2] - ps_value
    app$ConnectionPerm[idx2, idx1] <- app$ConnectionPerm[idx2, idx1] + ps_value
  }

  # dC/dt = (coefMatrixKnown / V(t)) %*% C ; V(t) = V0 - Kev * t within an interval
  app$coefMatrixKnown <- (-app$CompClearance + app$ConnectionPerm) * app$Compfu + app$ConnectionFlow + diag(app$KevVector)
  app
}

# Right-hand side for deSolve::ode()
ODEmodel <- function(t, C, app) {
  app$vVector <- app$InitVolumesofInterval - app$KevVector * (t - app$InitTimeofInterval)
  VMatrix <- matrix(app$vVector, nrow = app$NumofComp, ncol = app$NumofComp)  # VMatrix[i, j] = V_i
  app$CoefMatrix <- app$coefMatrixKnown / VMatrix
  SystemODE <- app$CoefMatrix %*% C
  list(SystemODE)
}

# Simulate concentration-time profiles. The ODE is solved
# piecewise between intervention times on a 0.01 h grid; at each intervention
# the volumes are updated (evaporation, removal, addition) and doses added.
# Adds app$t (time, h) and app$yc (concentration matrix, one column per
# compartment). The last grid point of each interval is not stored.
CalCulate <- function(app) {
  InitialConcentration <- app$DosesAll[, 1] / app$v0Vector
  Volumes <- matrix(rep(app$v0Vector, times = length(app$InterventionTimes)),
                    nrow = length(app$v0Vector),
                    ncol = length(app$InterventionTimes))
  Time <- c()
  Conc <- matrix(ncol = length(InitialConcentration), nrow = 0)

  for (iDose in 1:(length(app$InterventionTimes) - 1)) {
    tt <- c(app$InterventionTimes[iDose], app$InterventionTimes[iDose + 1])
    delta_t <- tt[2] - tt[1]
    app$EndVolumesofInterval  <- Volumes[, iDose] - app$KevVector * delta_t
    app$InitVolumesofInterval <- Volumes[, iDose]
    app$InitTimeofInterval    <- app$InterventionTimes[iDose]

    # 0.01 h output grid; very short intervals use their end points only
    from <- tt[1]
    to   <- tt[2] - 0.01
    if (!is.finite(from) || !is.finite(tt[2])) {
      stop("Non-finite intervention times: tt = ", paste(tt, collapse = ", "))
    }
    if (to < from) {
      ttt <- c(from, tt[2])
    } else {
      ttt <- seq(from = from, to = to, by = 0.01)
    }
    times <- c(ttt, tt[2])

    sol <- deSolve::ode(y = InitialConcentration, times = times, func = ODEmodel, parms = app)
    EndConcentration <- sol[nrow(sol), -1]

    Volumes[, iDose + 1] <- app$EndVolumesofInterval - app$ExMediaAll[, iDose + 1] + app$AddMediaAll[, iDose + 1]
    InitialConcentration <- (EndConcentration * (app$EndVolumesofInterval - app$ExMediaAll[, iDose + 1]) + app$DosesAll[, iDose + 1]) /
      (app$EndVolumesofInterval - app$ExMediaAll[, iDose + 1] + app$AddMediaAll[, iDose + 1])

    Time <- c(Time, sol[1:(nrow(sol) - 1), 1])
    Conc <- rbind(Conc, sol[1:(nrow(sol) - 1), -1])
  }

  app$t  <- Time
  app$yc <- Conc
  app
}

# ---- parameter estimation --------------------------------------------------

# Estimate the unknown (NULL) parameters by least squares with
# minpack.lm::nls.lm (Levenberg-Marquardt), bounds 1e-4 x r0 to 1000 x r0.
# Unknowns are ordered: compartment parameters (Vol, fu,
# CL, Kev of compartment 1, 2, ...), then flows, then permeable surfaces.
# Needs app$expData[[i]]$t / $y (observed time / concentration of compartment
# i) and app$r0 (start values). The model is solved on a 0.1 h grid and the
# prediction at an observation time is the grid point closest to it (the first
# one at an intervention time, i.e. the value before the event).
# Adds app$ParameterFit (estimate, SE, t, p) and app$Estimated[[k]]$Values.
AnalysisODE <- function(app) {
  app$SampTimes <- c()
  app$compartmentSample <- c()
  app$SampData <- c()
  for (ii in seq_along(app$expData)) {
    app$SampTimes <- c(app$SampTimes, app$expData[[ii]]$t)
    app$compartmentSample <- c(app$compartmentSample, rep(ii, length(app$expData[[ii]]$t)))
    app$SampData <- c(app$SampData, app$expData[[ii]]$y)
  }

  rdown <- 0.0001 * app$r0
  rup   <- 1000 * app$r0
  r0    <- app$r0

  residuals_function <- function(r, app) {
    app$CompNameID <- vector("list", app$NumofComp)
    thisUnknownPara <- 0
    CompClearance <- matrix(0, nrow = app$NumofComp, ncol = app$NumofComp)
    Compfu        <- matrix(0, nrow = app$NumofComp, ncol = app$NumofComp)
    app$v0VectorUnknown  <- numeric(app$NumofComp)
    app$KevVectorUnknown <- numeric(app$NumofComp)

    for (i in seq_along(app$CompData)) {
      app$CompNameID[[i]] <- app$CompData[[i]]$Name
      if (is.null(app$CompData[[i]]$Vol) || app$CompData[[i]]$Vol == "") {
        thisUnknownPara <- thisUnknownPara + 1
        app$v0VectorUnknown[i] <- r[thisUnknownPara]
      } else {
        app$v0VectorUnknown[i] <- as.numeric(app$CompData[[i]]$Vol)
      }
      if (is.null(app$CompData[[i]]$fu) || app$CompData[[i]]$fu == "") {
        thisUnknownPara <- thisUnknownPara + 1
        Compfu[, i] <- r[thisUnknownPara]
      } else {
        Compfu[, i] <- as.numeric(app$CompData[[i]]$fu)
      }
      if (is.null(app$CompData[[i]]$CL) || app$CompData[[i]]$CL == "") {
        thisUnknownPara <- thisUnknownPara + 1
        CompClearance[i, i] <- r[thisUnknownPara]
      } else {
        CompClearance[i, i] <- as.numeric(app$CompData[[i]]$CL)
      }
      if (is.null(app$CompData[[i]]$Kev) || app$CompData[[i]]$Kev == "") {
        thisUnknownPara <- thisUnknownPara + 1
        app$KevVectorUnknown[i] <- r[thisUnknownPara]
      } else {
        app$KevVectorUnknown[i] <- as.numeric(app$CompData[[i]]$Kev)
      }
    }

    ConnectionFlow <- matrix(0, nrow = app$NumofComp, ncol = app$NumofComp)
    ConnectionPerm <- matrix(0, nrow = app$NumofComp, ncol = app$NumofComp)
    ConnectionFlowVector <- matrix(0, nrow = length(app$FlowData), ncol = 2)
    ConnectionFlowVectorValue <- numeric(length(app$FlowData))
    ConnectionPermVector <- matrix(0, nrow = length(app$PermData), ncol = 2)
    ConnectionPermVectorValue <- numeric(length(app$PermData))

    for (i in seq_along(app$FlowData)) {
      from_value <- as.character(app$FlowData[[i]]$from[[1]][1, 1])
      to_value   <- as.character(app$FlowData[[i]]$to[[1]][1, 1])
      from_idx <- grep(from_value, sapply(app$CompNameID, function(x) x))
      to_idx   <- grep(to_value, sapply(app$CompNameID, function(x) x))
      ConnectionFlowVector[i, 1] <- from_idx
      ConnectionFlowVector[i, 2] <- to_idx
      if (is.null(app$FlowData[[i]]$Q) || app$FlowData[[i]]$Q == "") {
        thisUnknownPara <- thisUnknownPara + 1
        ConnectionFlowVectorValue[i] <- r[thisUnknownPara]
      } else {
        ConnectionFlowVectorValue[i] <- as.numeric(app$FlowData[[i]]$Q)
      }
    }

    for (i in seq_along(app$PermData)) {
      between1_value <- as.character(app$PermData[[i]]$between1[[1]][1, 1])
      between2_value <- as.character(app$PermData[[i]]$between2[[1]][1, 1])
      comp1_idx <- grep(between1_value, sapply(app$CompNameID, function(x) x))
      comp2_idx <- grep(between2_value, sapply(app$CompNameID, function(x) x))
      ConnectionPermVector[i, 1] <- comp1_idx
      ConnectionPermVector[i, 2] <- comp2_idx
      if (is.null(app$PermData[[i]]$PS) || app$PermData[[i]]$PS == "") {
        thisUnknownPara <- thisUnknownPara + 1
        ConnectionPermVectorValue[i] <- r[thisUnknownPara]
      } else {
        ConnectionPermVectorValue[i] <- as.numeric(app$PermData[[i]]$PS)
      }
    }

    for (i in seq_len(nrow(ConnectionFlowVector))) {
      from_idx <- ConnectionFlowVector[i, 1]
      to_idx   <- ConnectionFlowVector[i, 2]
      flow_value <- ConnectionFlowVectorValue[i]
      ConnectionFlow[from_idx, from_idx] <- ConnectionFlow[from_idx, from_idx] - flow_value
      ConnectionFlow[to_idx, from_idx]   <- ConnectionFlow[to_idx, from_idx] + flow_value
    }
    for (i in seq_len(nrow(ConnectionPermVector))) {
      idx1 <- ConnectionPermVector[i, 1]
      idx2 <- ConnectionPermVector[i, 2]
      perm_value <- ConnectionPermVectorValue[i]
      ConnectionPerm[idx1, idx1] <- ConnectionPerm[idx1, idx1] - perm_value
      ConnectionPerm[idx1, idx2] <- ConnectionPerm[idx1, idx2] + perm_value
      ConnectionPerm[idx2, idx2] <- ConnectionPerm[idx2, idx2] - perm_value
      ConnectionPerm[idx2, idx1] <- ConnectionPerm[idx2, idx1] + perm_value
    }
    app$coefMatrixUnknown <- (-CompClearance + ConnectionPerm) * Compfu + ConnectionFlow + diag(app$KevVectorUnknown)

    app <- RtoODE(app)

    model_predictions <- numeric(length(app$SampTimes))
    for (i in seq_along(app$SampTimes)) {
      idx <- which.min(abs(app$SolTime - app$SampTimes[i]))
      model_predictions[i] <- app$SolConc[idx, app$compartmentSample[i]]
    }
    model_predictions - app$SampData
  }

  # piecewise solution on a 0.1 h grid (same volume / dose bookkeeping as CalCulate)
  RtoODE <- function(app) {
    initialConcentration <- app$DosesAll[, 1] / app$v0VectorUnknown
    Volumes <- matrix(rep(app$v0VectorUnknown, length(app$InterventionTimes)), nrow = app$NumofComp)
    app$SolTime <- c()
    app$SolConc <- matrix(nrow = 0, ncol = app$NumofComp)

    for (iDose in 1:(length(app$InterventionTimes) - 1)) {
      tt <- c(app$InterventionTimes[iDose], app$InterventionTimes[iDose + 1])
      delta_t <- tt[2] - tt[1]
      app$EndVolumesofInterval  <- Volumes[, iDose] - app$KevVectorUnknown * delta_t
      app$InitVolumesofInterval <- Volumes[, iDose]
      app$InitTimeofInterval    <- app$InterventionTimes[iDose]

      ode_func <- function(t, C, parms) {
        vVector <- app$InitVolumesofInterval - app$KevVectorUnknown * (t - app$InitTimeofInterval)
        VMatrix <- matrix(rep(vVector, each = app$NumofComp), nrow = app$NumofComp)
        VMatrix <- t(VMatrix)
        SystemODE <- (app$coefMatrixUnknown / VMatrix) %*% C
        list(as.vector(SystemODE))
      }

      times <- seq(tt[1], tt[2], by = 0.1)
      sol <- deSolve::ode(y = initialConcentration, times = times, func = ode_func, parms = NULL)
      app$SolTime <- c(app$SolTime, sol[, "time"])
      app$SolConc <- rbind(app$SolConc, sol[, -1, drop = FALSE])
      EndConcentration <- sol[nrow(sol), -1]

      Volumes[, iDose + 1] <- app$EndVolumesofInterval - app$ExMediaAll[, iDose + 1] + app$AddMediaAll[, iDose + 1]
      initialConcentration <- (EndConcentration * (app$EndVolumesofInterval - app$ExMediaAll[, iDose + 1]) + app$DosesAll[, iDose + 1]) /
        (app$EndVolumesofInterval - app$ExMediaAll[, iDose + 1] + app$AddMediaAll[, iDose + 1])
    }
    app
  }

  fit <- minpack.lm::nls.lm(par = r0, fn = residuals_function, lower = rdown, upper = rup, app = app)
  app$ParameterFit <- coef(summary(fit))
  x <- fit$par
  for (jjj in seq_along(x)) {
    app$Estimated[[jjj]]$Values <- x[jjj]
  }
  app
}

# ---- exposure metrics ------------------------------------------------------

# Exposure metrics of one profile (used for the
# heatmaps of Figs. 4B, 5B and 6B). No smoothing.
#   Cmax, Tmax : maximum and its time
#   AUC        : trapezoidal rule over the simulated time range
#   T12        : profile-based half time = time from Tmax until the profile
#                falls to Cmax/2 (linear interpolation); if not reached, the
#                terminal log-linear slope of the last `tail_n` points is used
exposure_metrics <- function(time, conc, tail_n = 5) {
  stopifnot(length(time) == length(conc))
  stopifnot(length(time) >= 2)
  ord <- order(time)
  time <- as.numeric(time[ord])
  conc <- as.numeric(conc[ord])

  Cmax <- max(conc, na.rm = TRUE)
  i_max <- which.max(conc)
  Tmax <- time[i_max]
  dt <- diff(time)
  auc <- sum(dt * (conc[-length(conc)] + conc[-1]) / 2, na.rm = TRUE)

  half <- 0.5 * Cmax
  t_half_abs <- NA_real_
  T12 <- NA_real_
  if (i_max < length(conc)) {
    t2 <- time[(i_max + 1):length(time)]
    c2 <- conc[(i_max + 1):length(conc)]
    idx <- which(c2 <= half)[1]
    if (!is.na(idx)) {
      t_curr <- t2[idx]
      c_curr <- c2[idx]
      if (idx == 1) {
        t_prev <- Tmax
        c_prev <- Cmax
      } else {
        t_prev <- t2[idx - 1]
        c_prev <- c2[idx - 1]
      }
      if (!isTRUE(all.equal(c_prev, c_curr))) {
        frac <- (half - c_prev) / (c_curr - c_prev)
        t_half_abs <- t_prev + frac * (t_curr - t_prev)
      } else {
        t_half_abs <- t_curr
      }
      T12 <- t_half_abs - Tmax
    }
  }

  if (is.na(T12)) {
    t_tail <- time[time > Tmax]
    c_tail <- conc[time > Tmax]
    keep <- which(is.finite(c_tail) & c_tail > 0)
    t_tail <- t_tail[keep]
    c_tail <- c_tail[keep]
    if (length(c_tail) >= 3) {
      n_use <- min(tail_n, length(c_tail))
      t_fit <- tail(t_tail, n_use)
      c_fit <- tail(c_tail, n_use)
      fit <- stats::lm(log(c_fit) ~ t_fit)
      k <- -as.numeric(stats::coef(fit)[2])
      if (is.finite(k) && k > 0) {
        T12 <- log(2) / k
        t_half_abs <- Tmax + T12
      } else {
        T12 <- NA_real_
        t_half_abs <- NA_real_
      }
    }
  }
  list(Cmax = Cmax, Tmax = Tmax, AUC = auc, T12 = T12, t_half_abs = t_half_abs)
}

# Exposure metrics of all compartments after a moving-average smoothing
# (window = 1% of the points; edges set to 0). NOT used by any figure of the
# paper: the heatmap scripts of Figs. 4B/5B/6B call exposure_metrics(). Kept
# for completeness (needs the package pracma).
ExposureParas <- function(app) {
  app$AUC   <- rep(0, length(app$Comp))
  app$Cmax  <- rep(0, length(app$Comp))
  app$Tcmax <- rep(0, length(app$Comp))
  app$Thalf <- rep(0, length(app$Comp))
  app$TmaxTime <- vector("list", length(app$Comp))

  for (ii in 1:length(app$Comp)) {
    window_size <- round(length(app$yc[, ii]) / 100)
    app$yc[, ii] <- stats::filter(app$yc[, ii], rep(1 / window_size, window_size), sides = 2)
    app$yc[, ii][is.na(app$yc[, ii])] <- 0
    app$AUC[ii] <- pracma::trapz(app$t, app$yc[, ii])
    app$Cmax[ii] <- max(app$yc[, ii])
    indmax <- which.max(app$yc[, ii])
    app$Tcmax[ii] <- app$t[indmax]
    if (length(app$TmaxTime) < ii) app$TmaxTime[[ii]] <- list()
    if (indmax == length(app$yc[, ii]) && all(diff(app$yc[, ii]) >= 0)) {
      app$TmaxTime[[ii]]$reached <- "TmaxNotReached"
    } else if (indmax == length(app$yc[, ii]) && !all(diff(app$yc[, ii]) >= 0)) {
      app$TmaxTime[[ii]]$reached <- "TmaxIsLocalMax"
    } else if (indmax != length(app$yc[, ii]) && !all(diff(app$yc[, ii]) >= 0)) {
      app$TmaxTime[[ii]]$reached <- "TmaxReached&Passed"
    }
    tryCatch({
      app$Thalf[ii] <- pracma::interp1(app$yc[(indmax + 1):length(app$yc[, ii]), ii],
                                       app$t[(indmax + 1):length(app$t)],
                                       app$Cmax[ii] / 2) - app$Tcmax[ii]
    }, error = function(e) NULL)
  }
  app
}
