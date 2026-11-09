# =============================================================================
# fig3_fits.R - Fig. 3: diclofenac PK fits per study arm (CN Bio TL6 study)
#
# Fits hepatic intrinsic clearance CL_int (liver, compartment 3) and/or gut
# permeability PS (gut barrier, PermData[[1]]) to the measured medium
# concentrations of each experimental arm.
#
# Inputs : data/models/Kev2DAnalysis.rds (CNBio MPS-TL6 gut-liver model)
#          data/fig3/diclofenac_subsets.xlsx (one sheet per arm)
#          data/fig3/literature_comparators.csv (grey triangles in Fig. 3B)
# Outputs: fig3_estimates.csv  - per arm: CL (µL/min), PS (cm/s), RSE (%)
#          fig3_results.csv    - CL_int (µL/min/10^6 hepatocytes),
#                                P (10^-6 cm/s), RSE, + literature rows
#                                (plotted in Fig. 3B/C)
#          fig3_fit_profiles.csv, fig3_observed.csv - data for Fig. 3A
#          fig3_fit_<arm>.png  - one fit plot per arm
# Runtime: ~1 min
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})

# ---- parameters ----------------------------------------------------------------
MODEL_FILE  <- pc_data("models", "Kev2DAnalysis.rds")
DATA_FILE <- pc_data("fig3", "diclofenac_subsets.xlsx")
LIT_FILE  <- pc_data("fig3", "literature_comparators.csv")

# Liver compartment volume (mL) used for the fits of Fig. 3
LIVER_VOLUME_ML  <- 1.8
APICAL_VOLUME_ML <- 0.2      # gut apical (protocol: 0.2 mL)
CIRC_VOLUME_ML   <- 0.2      # circulation / channels; basolateral 0.8 mL from the model file
R0_CL_PS   <- c(0.15, 0.015) # start values CL_int, PS (mL/h)
R0_PS_ONLY <- 0.015          # start value PS (mL/h), gut-only arms
DAILY_FINAL_TIME_H    <- 110 # daily-dosing arms
MEDIA_EXCHANGE_TIMES  <- seq(24, 110, by = 24)  # h, daily 1/3-volume exchange (liver circuit)
EXCHANGE_VOLUME_ML    <- 0.867                  # removed = added (protocol: 867 µL)
STATIC_FLOW_ML_H      <- 20  # flow used for the static arms

# Arms (order of the result table). S_L_DD_IV_100uM is not fitted: data from one
# well only, an inconsistent first daily-phase sample and an undocumented daily dose
# in the static wells.
# Daily-dose (DD) arms: the study protocol has a medium exchange on day 7 after the
# pre-dose sample, so the daily phase starts with drug-free basolateral/liver
# medium; the day-7 (t = 0) sample belongs to the preceding dose and is not fitted.
ARMS <- c("D_G-L_SD_PO_1000uM", "D_G-L_SD_PO_100uM", "D_G-L_SD_PO_40uM",
          "D_G_SD_PO_100uM", "D_L_SD_IV_100uM", "S_G_SD_PO_100uM", "S_L_SD_IV_100uM",
          "D_G-L_DD_PO_1000uM", "D_G-L_DD_PO_100uM", "D_G-L_DD_PO_40uM",
          "D_G_DD_PO_100uM", "D_L_DD_IV_100uM")

# Location in the data -> compartment (apical data are not used in the fits)
LOCATION_COMPARTMENT <- c(Apical = 1, Basolateral = 2, Liver = 3)

# ---- one arm -------------------------------------------------------------------

# Read and clean the data of one arm: "< 0.05" -> 0, no NA, no apical samples
read_arm_data <- function(arm) {
  dat <- read_xlsx(DATA_FILE, sheet = arm)
  dat$Conc <- ifelse(grepl("<\\s*0\\.05", dat$Conc), 0, suppressWarnings(as.numeric(dat$Conc)))
  dat <- dat[!is.na(dat$Time_h) & !is.na(dat$Conc), ]
  if (grepl("_DD_", arm)) dat <- dat[as.numeric(dat$Time_h) > 0, ]   # pre-exchange day-7 sample (see above)
  dat[dat$Location != "Apical", ]
}

# Build the arm-specific model (dose amounts are V x concentration x MolarMass,
# the engine divides by MolarMass)
build_arm_app <- function(arm, dat) {
  app <- load_model(MODEL_FILE)
  for (loc in unique(dat$Location)) {
    idx <- LOCATION_COMPARTMENT[[loc]]
    sub <- dat[dat$Location == loc, ]
    app$expData[[idx]] <- list(t = as.numeric(sub$Time_h), y = as.numeric(sub$Conc))
  }

  # default: estimate CL_int (liver) and PS (gut barrier); no evaporation
  app$Estimated <- list(list(Done = 0), list(Done = 0))
  app$r0 <- R0_CL_PS
  app$CompData[[3]]$CL <- NULL
  app$PermData[[1]]$PS <- NULL
  app$CompData[[1]]$Kev <- 0
  app$CompData[[3]]$Kev <- 0
  app$CompData[[1]]$Vol <- APICAL_VOLUME_ML
  app$CompData[[3]]$Vol <- LIVER_VOLUME_ML
  app$CompData[[4]]$Vol <- CIRC_VOLUME_ML

  dose_uM <- as.numeric(sub(".*_([0-9]+)uM$", "\\1", arm))
  single_apical <- function(app) {
    app$CompData[[1]]$DoseCond <- 1
    app$CompData[[1]]$DoseSingData <- app$CompData[[1]]$Vol * dose_uM * app$MolarMass
    app
  }
  ps_only <- function(app) {             # gut-only arms: no liver, PS only
    app$Estimated <- list(list(Done = 0))
    app$r0 <- R0_PS_ONLY
    app$CompData[[3]]$CL <- 0
    app$PermData[[1]]$PS <- NULL
    app
  }
  daily <- function(app) {               # daily apical dosing after the day-7 exchange (C(0) = 0), daily 1/3 exchange
    app$CompData[[1]]$DoseCond <- 2
    app$CompData[[1]]$DoseMultData[1, 1] <- app$CompData[[1]]$Vol * dose_uM * app$MolarMass
    app$CompData[[1]]$DoseMultData[1, 2] <- 24
    app$FinalTime <- DAILY_FINAL_TIME_H
    app$CompData[[3]]$ExMedia <- cbind(MEDIA_EXCHANGE_TIMES,
                                       rep(EXCHANGE_VOLUME_ML, length(MEDIA_EXCHANGE_TIMES)),
                                       rep(EXCHANGE_VOLUME_ML, length(MEDIA_EXCHANGE_TIMES)))
    app
  }
  static_flows <- function(app) {
    for (k in 1:3) app$FlowData[[k]]$Q <- STATIC_FLOW_ML_H
    app
  }

  if (arm %in% c("D_G-L_SD_PO_1000uM", "D_G-L_SD_PO_100uM", "D_G-L_SD_PO_40uM")) {
    app <- single_apical(app)
  } else if (arm == "D_G_SD_PO_100uM") {
    app <- ps_only(single_apical(app))
  } else if (arm == "D_L_SD_IV_100uM") {
    # the dose is placed in the apical compartment and PS is co-estimated
    app <- single_apical(app)
  } else if (arm == "S_G_SD_PO_100uM") {
    app$CompData[[2]]$Vol <- 0.8
    app$CompData[[3]]$Vol <- 0.1
    app$CompData[[4]]$Vol <- 0.1
    app <- ps_only(single_apical(static_flows(app)))
  } else if (arm == "S_L_SD_IV_100uM") {
    app$CompData[[1]]$Vol <- 0.01
    app$CompData[[2]]$Vol <- 0.01
    app$CompData[[3]]$Vol <- 0.75
    app$CompData[[4]]$Vol <- 0.01
    app <- static_flows(app)
    app$CompData[[1]]$DoseSingData <- 0
    app$CompData[[3]]$DoseCond <- 1
    app$CompData[[3]]$DoseSingData <- app$CompData[[3]]$Vol * dose_uM * app$MolarMass
  } else if (arm %in% c("D_G-L_DD_PO_1000uM", "D_G-L_DD_PO_100uM", "D_G-L_DD_PO_40uM", "D_L_DD_IV_100uM")) {
    app <- daily(app)                    # D_L_DD_IV_100uM: apical dosing, as for D_L_SD_IV_100uM
  } else if (arm == "D_G_DD_PO_100uM") {
    app <- daily(ps_only(app))
  } else {
    stop("Unknown arm: ", arm)
  }
  app
}

fit_arm <- function(arm) {
  dat <- read_arm_data(arm)
  app <- build_arm_app(arm, dat)
  app <- ConstructDoseMediaEx(app)
  app <- ConstructODE(app)
  app <- AnalysisODE(app)

  RSECL <- NaN
  RSEPS <- NaN
  pf <- app$ParameterFit
  if (grepl("G-L", arm) || grepl("_L_", arm)) {          # CL_int and PS estimated
    app$CompData[[3]]$CL <- app$Estimated[[1]]$Values
    RSECL <- 100 * pf[1, 2] / pf[1, 1]
    app$PermData[[1]]$PS <- app$Estimated[[2]]$Values
    RSEPS <- 100 * pf[2, 2] / pf[2, 1]
  } else {                                                 # gut-only: PS estimated
    app$PermData[[1]]$PS <- app$Estimated[[1]]$Values
    RSEPS <- 100 * pf[1, 2] / pf[1, 1]
  }

  est <- data.frame(SheetName = arm,
                    CL = app$CompData[[3]]$CL * 1000 / 60,          # mL/h -> µL/min
                    RSECL = RSECL,
                    PS = app$PermData[[1]]$PS / (3600 * GUT_AREA_CM2),  # mL/h -> cm/s
                    RSEPS = RSEPS)

  # simulation with the estimates (Fig. 3A lines)
  app <- ConstructODE(app)
  app <- CalCulate(app)
  sim <- data.frame(arm = arm, time_h = app$t, compartment = 3, conc = app$yc[, 3])
  obs <- data.frame(arm = arm, location = dat$Location, time_h = as.numeric(dat$Time_h), conc = dat$Conc)
  list(est = est, sim = sim, obs = obs)
}

# ---- run all arms ----------------------------------------------------------------
pc_log("Fig. 3 fits, liver volume ", LIVER_VOLUME_ML, " mL")
res <- lapply(ARMS, function(arm) {
  pc_log("  ", arm)
  fit_arm(arm)
})
est <- do.call(rbind, lapply(res, `[[`, "est"))
sim <- do.call(rbind, lapply(res, `[[`, "sim"))
obs <- do.call(rbind, lapply(res, `[[`, "obs"))

# reported units: CL_int per 10^6 hepatocytes, P in 10^-6 cm/s;
# only the parameters that are identifiable for the arm type are reported
is_gut_only   <- grepl("_G_", est$SheetName)
is_liver_only <- grepl("_L_", est$SheetName)
results <- data.frame(SheetName = est$SheetName,
                      Clint = ifelse(is_gut_only, NA, est$CL / HEPATOCYTES_1E6),
                      RSECL = ifelse(is_gut_only, NA, est$RSECL),
                      P     = ifelse(is_liver_only, NA, est$PS * 1e6),
                      RSEPS = ifelse(is_liver_only, NA, est$RSEPS))
lit <- read.csv(LIT_FILE)
results <- rbind(results, data.frame(SheetName = lit$SheetName, Clint = lit$Clint, RSECL = NA,
                                     P = lit$P, RSEPS = NA))

write.csv(est, pc_out("fig3_estimates.csv"), row.names = FALSE)
write.csv(results, pc_out("fig3_results.csv"), row.names = FALSE)
write.csv(sim, pc_out("fig3_fit_profiles.csv"), row.names = FALSE)
write.csv(obs, pc_out("fig3_observed.csv"), row.names = FALSE)

for (arm in ARMS) {
  p <- plot_fig3_fits(obs, sim, arm, base_size = 12)
  save_panel(p, paste0("fig3_fit_", arm), width = 6, height = 4)
}
print(results, digits = 6)
pc_log("done: ", pc_run_dir())
