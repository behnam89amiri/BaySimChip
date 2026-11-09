# =============================================================================
# R/scenarios.R - set-up of the virtual experiments used in the figures
#
# Each function sets the parameters of one virtual experiment (values, and the
# order of the assignments where it matters, as used for the figures of the
# paper).
# =============================================================================

# Generic compound in the CNBio MPS-TL6 gut-liver model (Figs. 4-6 A/B and the
# media-exchange scenario):
#   CL_int (liver, compartment 3) = 0.1 mL/h, PS (gut barrier) = 0.01 mL/h,
#   single apical dose giving 1 µM (amount = V_apical x 1 µM x MolarMass,
#   the engine divides by MolarMass again).
# The model files store CL 0.08 / PS 0.03 mL/h (earlier defaults); they are
# always overwritten here and never used.
setup_generic_compound <- function(model_file,
                                   final_time,
                                   CL = 0.1,
                                   PS = 0.01,
                                   volumes = c(0.35, 0.8, 1.6, 0.2),
                                   dose_uM = 1) {
  app <- load_model(model_file)
  app$FinalTime <- final_time
  app$CompData[[3]]$CL <- CL
  app$PermData[[1]]$PS <- PS
  for (i in seq_along(volumes)) app$CompData[[i]]$Vol <- volumes[i]
  app$CompData[[1]]$DoseSingData <- app$CompData[[1]]$Vol * dose_uM * app$MolarMass
  app
}

# Add the apical solid (NSB) compartment to the NSB model NSBpara2DAnalysis.rds:
# compartment 6 ("Comp6", a copy of the basolateral solid compartment 5) and
# permeable surface 3 between Comp6 and the gut apical compartment ("Comp0").
# PS_NSB of surface 3 is set by the caller (apical PS_NSB = basolateral
# PS_NSB / ratio).
add_apical_nsb <- function(app) {
  app$CompData[[6]] <- clone_env(app$CompData[[5]])
  app$CompData[[6]]$Name <- list(Name = matrix("Comp6", nrow = 1, ncol = 1))
  surf <- new.env(parent = emptyenv())
  surf$Name     <- list(Name = matrix("P.S.3", nrow = 1, ncol = 1))
  surf$between1 <- list(Name = matrix("Comp6", nrow = 1, ncol = 1))
  surf$between2 <- list(Name = matrix("Comp0", nrow = 1, ncol = 1))
  surf$PS <- 0
  app$PermData[[3]] <- surf
  app$NumofComp <- 6
  app
}

# Set NSB parameters of both solid compartments: K_media/plastic (stored in
# `fu` of compartments 5 and 6) and PS_NSB (basolateral surface 2; apical
# surface 3 = PS_NSB / nsb_ratio)
set_nsb <- function(app, k_media_plastic, ps_nsb, nsb_ratio) {
  app$CompData[[5]]$fu <- k_media_plastic
  app$PermData[[2]]$PS <- ps_nsb
  app$CompData[[6]]$fu <- k_media_plastic
  app$PermData[[3]]$PS <- ps_nsb / nsb_ratio
  app
}

# Media-exchange / sampling table of one compartment from a list of event
# times: rows = times[1:(n-1)] (the last listed time is not applied, as in the
# original analysis), removed volume `removed`, added volume `added` (mL)
exchange_table <- function(times, removed, added = 0) {
  n <- length(times) - 1
  ex <- matrix(0, nrow = n, ncol = 3)
  ex[, 1] <- times[seq_len(n)]
  ex[, 2] <- removed
  ex[, 3] <- added
  ex
}

# Simulate an app and return a long data frame (time, compartment, conc)
simulate_profiles <- function(app, compartments = seq_len(app$NumofComp)) {
  app <- ConstructODE(app)
  app <- CalCulate(app)
  do.call(rbind, lapply(compartments, function(j) {
    data.frame(time_h = app$t, compartment = j, conc = app$yc[, j])
  }))
}

# Liver (compartment 3) exposure metrics of a simulated app
liver_exposure <- function(app) {
  app <- ConstructODE(app)
  app <- CalCulate(app)
  m <- exposure_metrics(app$t, app$yc[, 3])
  data.frame(Cmax = m$Cmax, Tmax = m$Tmax, AUC = m$AUC, T12 = m$T12)
}

# Unit conversions used for reporting (SI S1 "Units and reported quantities")
#   P [cm/s] = PS [mL/h] / (3600 s/h x 0.33 cm^2 gut insert area)
#   CL_int [µL/min per 10^6 hepatocytes] = CL [mL/h] x 1000/60 / 0.6
#   (0.6 x 10^6 hepatocytes per liver compartment)
GUT_AREA_CM2 <- 0.33
HEPATOCYTES_1E6 <- 0.6
ps_to_P_1e6cms <- function(PS) PS * (1000000 / (3600 * GUT_AREA_CM2))
cl_to_uLmin1e6 <- function(CL) CL * (10 / 6) * (1000 / 60)
