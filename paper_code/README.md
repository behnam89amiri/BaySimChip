# paper_code - analyses and figures of the paper

R code and processed data that reproduce the analyses and figures of the paper
"Quantitative modeling of platform-specific processes to improve
pharmacokinetic studies in organ-on-a-chip systems" (submitted): main-text
Figs. 3-7 and the SI figures S3-S5 (SI Sections S5-S7: memantine PK in the
gut-liver MPS, non-specific binding characterised by washout, FIAU
toxicokinetics-toxicodynamics). The simulations and fits use the BaySimChip R
engine in `R/engine.R`. The BaySimChip app (`../BaySimChipR/`) has its own
implementation of the engine, which is cross-checked against this one
(`verification/check_app_engine.R`).

## How to run

R >= 4.2 with the packages listed under "Software". From the `paper_code/`
folder:

```
Rscript run_all.R                   # all scripts (about 10-20 min on one core)
Rscript run_all.R fig4 fig7         # only the scripts whose name starts with fig4 / fig7
```

Each run writes to a new folder `output/run_YYYY-MM-DD` (`_2`, `_3`, ... if the
folder exists), or to the folder given in the environment variable
`PAPER_CODE_RUN_DIR`: result tables (CSV), panel PNGs, composite figures
(600 dpi PNG + PDF), `logs/`, `timings.csv` and `sessionInfo.txt`. Nothing is
overwritten: writing a file that already exists stops the script.

Single scripts can be run from `paper_code/` as well. Set `PAPER_CODE_RUN_DIR`
so that scripts that build on each other use the same folder:

```
# bash                                         # PowerShell
export PAPER_CODE_RUN_DIR=output/my_run        $env:PAPER_CODE_RUN_DIR = "output/my_run"
Rscript scripts/fig3_fits.R
Rscript scripts/fig3_plots.R
```

Scripts never call `setwd()`; paths are relative to the `paper_code` root
(helpers `pc_data()`, `pc_out()` in `R/setup.R`).

Dependencies between scripts (same run folder): `fig3_plots.R` needs
`fig3_fits.R`; `figures_composite.R` needs `fig4ab_*`, `fig5ab_*`, `fig6ab_*`,
`fig456c_bias_maps.R`, `fig6d_*`, `fig7a_*` and `fig7b_*`;
`verification/check_app_engine.R` needs `fig4ab_*` (and `../BaySimChipR/`).
All other scripts are independent.

## Folder layout

| Path | Content |
|---|---|
| `run_all.R` | runs all scripts in order, one R process each |
| `R/setup.R` | root/path helpers, output folder, never-overwrite rule |
| `R/engine.R` | BaySimChip R engine (model import, dosing/exchange schedule, ODE system, simulation, least-squares fit, exposure metrics) |
| `R/scenarios.R` | set-up of the virtual experiments (generic compound, apical NSB compartment, sampling tables), unit conversions |
| `R/misspecification.R` | virtual "truth" data and refits with an incomplete analysis model (Figs. 4C-6C, 7) |
| `R/plotting.R` | labels, panel and composite plotting functions |
| `R/fiau_tktd.R` | functions of SI Section S7: data readers, dosing schedule and exposure model of the three plates, damage model (exact integration), multistart fit with profile likelihood, static Hill fit |
| `scripts/` | one script per figure or panel group (table below) |
| `verification/check_app_engine.R` | Fig. 4A recomputed with the engine of the BaySimChip app |
| `data/` | processed data and the BaySimChip models, see `data/README.md` |
| `output/v1.0.0/` | outputs of the run behind the paper (see "Reproducibility") |

## Figures, scripts, outputs

| Figure | Script | Main outputs |
|---|---|---|
| Fig. 3A-C | `scripts/fig3_fits.R` | `fig3_estimates.csv`, `fig3_results.csv`, `fig3_fit_profiles.csv`, `fig3_observed.csv`, `fig3_fit_<arm>.png` |
|  | `scripts/fig3_plots.R` | `fig3B_estimates.png`, `fig3C_systematic.png`, **`Fig3.png/.pdf`** |
| Fig. 4A/B | `scripts/fig4ab_evaporation.R` | `fig4a_profiles.csv`, `fig4a_liver_cmax.csv`, `fig4b_exposure_grid.csv`, panels |
| Fig. 5A/B | `scripts/fig5ab_sampling.R` | `fig5a_profiles.csv`, `fig5b_exposure_grid.csv`, `fig5b_T12_change.csv`, panels |
| text (media exchange) | `scripts/exchange_imbalance.R` | `exchange_profiles.csv`, `exchange_summary.csv` |
| Fig. 6A/B | `scripts/fig6ab_nsb.R` | `fig6a_profiles.csv`, `fig6b_exposure_grid.csv`, panels |
| Figs. 4C, 5C, 6C | `scripts/fig456c_bias_maps.R` (cases Evap, Samp, NSB) | `fig456c_bias_<case>.csv`, panels |
| Fig. 6D | `scripts/fig6d_empty_chip.R` | `fig6d_fit_results_<fit>.csv`, `fig6d_panel_*.csv`, `fig6d_panel.png` |
| Fig. 7A | `scripts/fig7a_archetypes.R` | `fig7a_misspecification_results.csv`, `fig7a_fold_errors.csv`, radar panels |
| Fig. 7B | `scripts/fig7b_drugs.R` | `fig7b_alldrugs_results.csv`, panels |
| Figs. 4-7 composites | `scripts/figures_composite.R` | **`Fig4` ... `Fig7` .png/.pdf** |
| SI S5 (Fig. S3) | `scripts/si_s5_memantine.R` | `si_s5_memantine_parameters.csv`, **`si_s5_memantine_fit.png/.pdf`** |
| SI S6 (Fig. S4) | `scripts/si_s6_nsb_washout.R` | `si_s6_nsb_tidy_data.csv`, `si_s6_washout_summary.csv`, **`si_s6_washout.png/.pdf`** |
| SI S7 (Fig. S5) | `scripts/si_s7_fiau_tktd.R` (+ `R/fiau_tktd.R`) | `si_s7_pk_tidy.csv`, `si_s7_pk_parameters.csv`, `si_s7_albumin_dosemeans.csv`, `si_s7_tktd_parameters.csv`, `si_s7_tktd_fit_summary.csv`, `si_s7_static_ec50.csv`, `si_s7_ec50_drift.csv`, **`si_s7_tktd.png/.pdf`** |
| cross-check | `verification/check_app_engine.R` | `app_engine_check.csv` |

Fig. 2 (screenshots of the app) is made with `../BaySimChipR/tools/make_screenshots.mjs`.

## Model and numerical settings

Units: h, mL, µM; flows, PS, PS_NSB, CL_int and K_ev in mL/h. Compartments of
the CN Bio TL6 model: 1 gut apical, 2 gut basolateral, 3 liver, 4 circulation
(Q = 3.6 mL/h between 2 -> 3 -> 4 -> 2); NSB model: + 5 solid compartment
exchanging with the liver compartment of the circuit and 6 solid compartment
exchanging with the apical compartment (apical PS_NSB = basolateral / 5 in
Figs. 6 and 7). For a solid compartment the field `fu` holds K_media/plastic.
Reported units: P = PS / (3600 x 0.33 cm²); CL_int per 10^6 hepatocytes =
CL x 1000/60 / 0.6.

- Simulation (`CalCulate`): `deSolve::ode()` default method `lsoda`, default
  `rtol = atol = 1e-6`, solved piecewise between dosing/exchange events,
  output every 0.01 h.
- Fit (`AnalysisODE`): `minpack.lm::nls.lm` (Levenberg-Marquardt), bounds
  1e-4 x start to 1000 x start, default control (maxiter 50); model solved
  with lsoda (default tolerances) on a 0.1 h grid, prediction = nearest grid
  point; SE from `summary()` of the nls.lm fit (Jacobian at the optimum),
  RSE = 100 x SE / estimate.
- Fig. 6D (`fig6d_empty_chip.R`): lsoda (default tolerances) + nls.lm
  (maxiter 200).
- SI S5 memantine: fixed-step RK4 (dt 0.02 h) + nls.lm on log parameters; RSE
  from `vcov()` on the log scale. SI S7 TKTD: the damage ODE is solved exactly
  between exposure changes (piecewise-constant medium concentration), `nls.lm`
  with 72 start values, parameters on log / logit scale with bounds, RSE from
  the Jacobian and 95% profile-likelihood intervals; the exposure is a
  one-compartment medium model fitted to the LC-MS data on the log scale.
- Exposure metrics of the heatmaps (Figs. 4B, 5B, 6B): `exposure_metrics()` on
  the unsmoothed 0.01 h profile (AUC by trapezoidal rule, Cmax/Tmax at the
  maximum, T1/2 = time from Tmax until Cmax/2 is reached).
- Behaviours of the original implementation that are kept: `ConstructODE()`
  builds the fu matrix transposed with respect to `AnalysisODE()` (fluid
  profiles are identical; the concentration of a solid compartment is scaled by
  K_media/plastic); media exchange / sampling events at t = 0 are ignored (doses
  at t = 0 are applied); in the sampling schedules the last listed time is not
  applied as a withdrawal; custom dosing tables are read by rows (not used).
- Fig. 3: liver compartment volume 1.8 mL (parameter at the top of
  `fig3_fits.R`); the arms with daily dosing start the daily phase with
  drug-free basolateral/liver medium after the day-7 medium exchange.
- Random numbers: virtual data with noise use fixed seeds (Figs. 4C-6C seed 1,
  Fig. 7A seed 128, Fig. 7B seed 17); the jitter of the points in Fig. 3C uses
  seed 1.

## Reproducibility

`output/v1.0.0/` holds the outputs of this version of the code, produced with

```
PAPER_CODE_RUN_DIR=output/v1.0.0 Rscript run_all.R
```

on Windows 11 with R 4.2.3 in 11 min (packages and versions:
`output/v1.0.0/sessionInfo.txt`; run times: `output/v1.0.0/timings.csv`; the
per-script logs are not included).
The figures of the paper are the composite files `Fig3` ... `Fig7` and
`si_s5_memantine_fit`, `si_s6_washout`, `si_s7_tktd` (.png/.pdf) of that
folder; the numbers reported in the paper are in the CSV files. All scripts use
fixed seeds, so repeated runs on the same platform give the same CSV files.
Other versions of R or of the packages, or another platform, can change the
last digits of numerical results (ODE solver, optimiser).

The engine of the BaySimChip app gives the same Fig. 4A liver profiles as
`R/engine.R` (`app_engine_check.csv`).

## Software

R 4.2.3. Packages: deSolve 1.40, minpack.lm 1.2-4, ggplot2 3.5.2, patchwork
1.2.0, scales 1.4.0, ragg 1.3.0, dplyr 1.1.4, tidyr 1.3.1, readxl 1.4.5,
ggrepel 0.9.5, viridisLite 0.4.2, pracma 2.4.4 (only for `ExposureParas()`,
which no figure uses). The scripts do not install packages. Licences of the
packages: see the README in the repository root.
