# =============================================================================
# fig7b_drugs.R - Fig. 7B: CL_int and P estimated with evaporation, sampling
# and NSB omitted from the analysis model, for 29 drugs with literature Caco-2
# P and hepatocyte CL_int (Li et al. 2007), versus the input ("observed") values
#
# For each drug: liver data (0-72 h) are simulated with evaporation, sampling
# and NSB at the nominal values of Fig. 7A plus 35% proportional noise, then
# CL_int and PS are re-estimated with all three processes removed from the
# analysis model (as for the bias maps of Figs. 4C-6C).
#
# Inputs : data/models/NSBpara2DAnalysis.rds, data/fig7/Li2007_Papp_CLint.csv
# Outputs: fig7b_alldrugs_results.csv, fig7b_CL.png, fig7b_P.png
# Runtime: ~2 min
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages(library(ggplot2))

# ---- parameters ----------------------------------------------------------------
MODEL_FILE  <- pc_data("models", "NSBpara2DAnalysis.rds")
DRUG_FILE <- pc_data("fig7", "Li2007_Papp_CLint.csv")  # CL, PS in mL/h (model units)
# nominal platform processes of Fig. 7A (scripts/fig7a_archetypes.R, Table S3)
KEV_APICAL <- 0.0026        # mL/h
KEV_LIVER  <- 0.005         # mL/h
K_MEDIA_PLASTIC <- 1
PS_NSB <- 0.05              # mL/h, basolateral
NSB_RATIO <- 5              # apical PS_NSB = PS_NSB / 5
SAMPLE_VOLUME <- 0.05       # mL
TIMES_LIVER  <- c(0, 1, 2, 4, 8, 24, 48, 72)   # observations; withdrawals 1-48 h (0 h ignored, 72 h not applied)
TIMES_APICAL <- c(24, 48, 72)                  # withdrawals 24, 48 h
FINAL_TIME_H <- 72
PROP_NOISE_SD <- 0.35       # proportional noise of the virtual data
ANALYSIS <- "all_neglected" # evaporation, sampling and NSB removed (run_alldrugs())
SEED <- 17                  # representative realisation: over seeds 1-60, CL_int outside 2-fold
                            # median 9 (range 5-14), outside 3-fold median 4 (1-7); seed 17: 9 and 4;
                            # P outside 2-fold 3, outside 3-fold 1
NAME_FIXES <- c(Gancyclovir = "Ganciclovir")   # spelling in the input file

# ---- run -----------------------------------------------------------------------
drugs <- read.csv(DRUG_FILE)
app <- build_misspec_app(MODEL_FILE, kev_apical = KEV_APICAL, kev_liver = KEV_LIVER,
                         ps_nsb = PS_NSB, k_media_plastic = K_MEDIA_PLASTIC, nsb_ratio = NSB_RATIO,
                         samp_vol = SAMPLE_VOLUME, times_liver = TIMES_LIVER, times_apical = TIMES_APICAL,
                         final_time = FINAL_TIME_H)
set.seed(SEED)
final <- run_alldrugs(app, drugs, PROP_NOISE_SD, analysis = ANALYSIS)
fix <- final$Compound %in% names(NAME_FIXES)
final$Compound[fix] <- NAME_FIXES[final$Compound[fix]]
write.csv(final, pc_out("fig7b_alldrugs_results.csv"), row.names = FALSE)

# ---- plots -----------------------------------------------------------------------
pCL <- plot_drug_scatter(data.frame(name = final$Compound, x = final$true_CL, y = final$est_CL),
                         lab("CLint", "Reference"), lab("CLint", "Estimated"))
pP  <- plot_drug_scatter(data.frame(name = final$Compound, x = final$true_P, y = final$est_P),
                         lab("P", "Reference"), lab("P", "Estimated"))
save_panel(pCL, "fig7b_CL", width = 6, height = 5)
save_panel(pP, "fig7b_P", width = 6, height = 5)
lr_CL <- abs(log(final$est_CL / final$true_CL))
lr_P  <- abs(log(final$est_P / final$true_P))
pc_log(sprintf("within 2-fold: CL %d/%d, P %d/%d; outside 3-fold: CL %d, P %d",
               sum(lr_CL <= log(2), na.rm = TRUE), nrow(final), sum(lr_P <= log(2), na.rm = TRUE), nrow(final),
               sum(lr_CL > log(3), na.rm = TRUE), sum(lr_P > log(3), na.rm = TRUE)))
pc_log("done: ", pc_run_dir())
