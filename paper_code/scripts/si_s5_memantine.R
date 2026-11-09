# =============================================================================
# si_s5_memantine.R - SI Section S5: memantine PK in the gut-liver MPS
#
# Joint fit of gut permeability PS and liver clearance CL_int to the liver
# concentrations of the two oral arms (1 and 10 µM); the intravenous arm
# (0.09 µM) is predicted as a check.
#
# Model (self-contained, no model file): apical compartment and one central
# compartment (basolateral + liver, well mixed; the flow itself is not modelled):
#   dC_ap/dt = -PS (C_ap - C_c) / V_ap
#   dC_c/dt  = (PS (C_ap - C_c) - CL C_c) / V_c
# solved with a fixed-step RK4 (dt = 0.02 h). Fit: minpack.lm::nls.lm on
# log(PS), log(CL); residuals scaled by the maximum of each arm; RSE (%) =
# 100 x SE of the log parameter (vcov of the fit).
#
# Inputs : data/si_s5_memantine/MNT-{PO-1,PO-10,IV-009}-{Ga,Lv}.xlsx
#          (Ga = gut apical, Lv = liver; time (h), concentration (µM))
# Outputs: si_s5_memantine_parameters.csv, si_s5_memantine_fit.png/.pdf
# Runtime: ~10 s
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
})

# ---- parameters ----------------------------------------------------------------
DATA_DIR <- pc_data("si_s5_memantine")
# volumes of the study protocol (gut apical 300 µL, gut basolateral 800 µL,
# liver 1800 µL; liver + gut basolateral = 2.6 mL)
V_AP <- 0.3                 # mL, gut apical
V_BL <- 0.8                 # mL, gut basolateral
V_LV <- 1.8                 # mL, liver
V_C  <- V_BL + V_LV         # central pool (basolateral + liver)
DT_H <- 0.02                # RK4 step (h)
START <- c(logPS = log(0.5), logCL = log(0.5))
LOWER <- c(log(1e-3), log(1e-3))
UPPER <- c(log(50), log(50))
MAXITER <- 200
ARMS <- data.frame(id = c("PO_1", "PO_10", "IV"),
                   prefix = c("MNT-PO-1", "MNT-PO-10", "MNT-IV-009"),
                   route = c("PO", "PO", "IV"),
                   label = c("PO 1 µM (fit)", "PO 10 µM (fit)", "IV 0.09 µM (validation)"))

# ---- data ----------------------------------------------------------------------
read_loc <- function(prefix, loc) {
  d <- suppressMessages(read_excel(file.path(DATA_DIR, sprintf("%s-%s.xlsx", prefix, loc)), col_names = FALSE))
  d <- d[, 1:2]
  colnames(d) <- c("time_h", "conc")
  d$time_h <- suppressWarnings(as.numeric(d$time_h))
  d$conc   <- suppressWarnings(as.numeric(d$conc))
  d[!is.na(d$time_h) & !is.na(d$conc), ]
}
liver_data <- list()
apical0 <- list()
for (k in seq_len(nrow(ARMS))) {
  lv <- read_loc(ARMS$prefix[k], "Lv")
  liver_data[[ARMS$id[k]]] <- lv %>% group_by(time_h) %>% summarise(conc = mean(conc), .groups = "drop")
  if (ARMS$route[k] == "PO") {
    ga <- read_loc(ARMS$prefix[k], "Ga")
    apical0[[ARMS$id[k]]] <- mean(ga$conc[ga$time_h == min(ga$time_h)])   # initial apical concentration
  }
}
iv_lv0 <- liver_data[["IV"]]$conc[liver_data[["IV"]]$time_h == 0]          # IV: initial central concentration

# ---- model ---------------------------------------------------------------------
deriv <- function(y, PS, CL) c(-PS * (y[1] - y[2]) / V_AP, (PS * (y[1] - y[2]) - CL * y[2]) / V_C)
simulate <- function(PS, CL, route, c0, t_out) {
  y <- if (route == "PO") c(c0, 0) else c(0, c0)
  t_req <- sort(unique(c(0, t_out)))
  grid <- seq(0, max(t_req), by = DT_H)
  c_ap <- numeric(length(grid)); c_c <- numeric(length(grid))
  c_ap[1] <- y[1]; c_c[1] <- y[2]
  for (i in 2:length(grid)) {
    k1 <- deriv(y, PS, CL)
    k2 <- deriv(y + DT_H / 2 * k1, PS, CL)
    k3 <- deriv(y + DT_H / 2 * k2, PS, CL)
    k4 <- deriv(y + DT_H * k3, PS, CL)
    y <- y + DT_H / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    c_ap[i] <- y[1]; c_c[i] <- y[2]
  }
  data.frame(time = t_req, C_ap = approx(grid, c_ap, xout = t_req)$y, C_lv = approx(grid, c_c, xout = t_req)$y)
}

# ---- fit -------------------------------------------------------------------------
resid_fun <- function(logp) {
  PS <- exp(logp[1]); CL <- exp(logp[2])
  res <- c()
  for (id in c("PO_1", "PO_10")) {
    ld <- liver_data[[id]]
    sim <- simulate(PS, CL, "PO", apical0[[id]], ld$time_h)
    pred <- approx(sim$time, sim$C_lv, xout = ld$time_h, rule = 2)$y
    res <- c(res, (pred - ld$conc) / max(ld$conc))
  }
  res
}
fit <- minpack.lm::nls.lm(par = START, fn = resid_fun, lower = LOWER, upper = UPPER,
                          control = minpack.lm::nls.lm.control(maxiter = MAXITER))
PS_hat <- exp(coef(fit)[1]); CL_hat <- exp(coef(fit)[2])
se <- tryCatch(sqrt(diag(vcov(fit))), error = function(e) c(NA, NA))
params <- data.frame(parameter = c("PS (gut permeability)", "P (gut permeability)", "CL (liver clearance)"),
                     value = c(PS_hat, PS_hat / (3600 * GUT_AREA_CM2), CL_hat),
                     unit = c("mL/h", "cm/s", "mL/h"),
                     RSE_pct = c(100 * se[1], 100 * se[1], 100 * se[2]))
write.csv(params, pc_out("si_s5_memantine_parameters.csv"), row.names = FALSE)
print(params, digits = 5)

# ---- figure --------------------------------------------------------------------
preds <- bind_rows(lapply(seq_len(nrow(ARMS)), function(k) {
  c0 <- if (ARMS$route[k] == "PO") apical0[[ARMS$id[k]]] else iv_lv0
  sim <- simulate(PS_hat, CL_hat, ARMS$route[k], c0, seq(0, 72, length.out = 200))
  data.frame(time_h = sim$time, conc = sim$C_lv, arm_label = ARMS$label[k])
}))
obs <- bind_rows(lapply(seq_len(nrow(ARMS)), function(k) mutate(liver_data[[ARMS$id[k]]], arm_label = ARMS$label[k])))
cols <- setNames(c("#2C7FB8", "#E67E22", "#7B4FA0"), ARMS$label)
p <- ggplot() +
  geom_line(data = preds, aes(time_h, conc, colour = arm_label), linewidth = 1) +
  geom_point(data = obs, aes(time_h, conc, colour = arm_label), size = 2.6) +
  facet_wrap(~arm_label, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = cols) +
  labs(x = lab("time"), y = "Liver concentration (µM)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey95"))
save_figure(p, "si_s5_memantine_fit", width = 8, height = 2.7)   # print size (was 11 x 3.4 in)
pc_log("done: ", pc_run_dir())
