# BaySimChip (R Shiny)

Pharmacokinetic (PK) modelling of organ-on-chip (OoC) systems: build a
compartmental model of the chip, define the experimental conditions (dosing,
sampling, media exchange, evaporation, non-specific binding), simulate
concentration-time profiles and estimate unknown parameters (e.g. intrinsic
clearance and permeability) from measured data.

BaySimChip is an R Shiny app for compartmental PK models of OoC systems. It has
four tabs (see below) and uses the engine that was used for the analyses of the
paper (refactored into pure functions, see [Validation](#validation)).

## Quick start

```r
# R >= 4.1 (developed and tested with R 4.2.3)
install.packages(c("shiny", "bslib", "DT", "visNetwork", "plotly",
                   "deSolve", "minpack.lm", "jsonlite", "readxl"))
shiny::runApp("path/to/BaySimChipR", launch.browser = TRUE)
```

From a terminal, in the folder that contains `BaySimChipR/`:

```
Rscript -e "shiny::runApp('BaySimChipR', launch.browser = TRUE)"
```

The app opens with the example "Gut-liver chip (CN Bio TL6)". Further examples (evaporation,
sampling, non-specific binding, two diclofenac fitting examples) are loaded
from *Model structure > Example projects*.

Optional packages: `ggplot2` (batch example), `testthat`, `callr`, `curl`
(tests).

## The four tabs

| Tab | Content |
|---|---|
| **Model structure** | *New* / *Load* / *Save model*. Drag-and-drop editor (visNetwork in manipulation mode) with the buttons *Compartment*, *Solid compartment*, *Fluid flow*, *Permeable surface* and *Remove selected*: click *Compartment* / *Solid compartment*, then click in the diagram; click *Fluid flow* / *Permeable surface*, then drag from one compartment to the other; drag compartments to arrange the layout; *Remove selected* (with confirmation). A properties panel edits the selected element (name, type, volume, CL, fu, Kev; Q; PS). Editable tables (compartments, flows, permeable surfaces) are a second, synced view. Projects are saved as `.rds` (native) and can be exported as JSON. |
| **Experimental conditions** | Panels *Chip and biological parameters* (volume, clearance, fu; K_media/plastic of solid compartments; flow rate; permeability x area), *Dosing* (single / repetitive / custom administration, compound molar mass, *use first experimental concentration*) and *Sampling and media exchange* (evaporation rate, table of time / extracted / added volume) for the element selected in the diagram or list; every parameter has an *unknown (estimate)* check box (its value is the initial guess). Dose unit selectable (mg, µg, ng, nmol). Custom-dose and exchange tables can be loaded from Excel/CSV. Overview table of all parameters. |
| **Simulation** | Simulation time, *Simulate*, *Clear plot*; several compartments per plot, overlay of runs (e.g. with / without evaporation), linear / log axis, table view, exposure parameters (AUC, Cmax, half-life, Tmax) per run and compartment (optionally in a time window), CSV downloads, ODE solver settings. |
| **Fitting and analysis** | CSV/Excel import for one or several compartments (column mapping, values "< LLOQ" set to 0 / dropped / LLOQ/2), editable data table, selection of the unknown parameters with initial values and bounds, Levenberg-Marquardt fit (`minpack.lm::nls.lm`), estimates with SE, RSE % and 95 % CI, fit / residual plots, exposure parameters of the fitted model, *Apply estimates to model*. |

The model diagram is shown on every tab; clicking an element selects it
everywhere.

## Files

* **`.rds` (native)** - *Save model* writes the whole project: structure,
  diagram layout, experimental conditions, observed data, fit settings and the
  last fit. *Load model* reads `.rds` and `.json`.
* **`.json`** - *Export JSON*: model, conditions and data in a readable,
  interoperable format (without stored results).
* The example projects in `inst/examples/` are built from the CNBio models in
  `inst/extdata/` (`.rds`) with `inst/examples/make_examples.R`.

## Using the engine in analysis scripts (no app)

The engine (`R/engine_*.R`) is a set of pure functions on a plain project list
(no Shiny, no global state), so the paper's batch analyses can call it directly:

```r
source("path/to/BaySimChipR/tools/load_engine.R")

p   <- bsc_read_project("my_chip.rds")              # or bsc_load_example("CNBio_evaporation")
p   <- bsc_set_params(p, list("Liver|cl" = 0.1, "P.S.0|ps" = 0.01))
sim <- bsc_simulate(p, final_time = 110)              # sim$time, sim$conc, sim$volume
bsc_exposure(sim$time, sim$conc)                      # AUC, Cmax, Tmax, half-life

sc  <- data.frame(scenario = c("with evaporation", "without"),
                  "Apical|kev" = c(0.003, 0), "Liver|kev" = c(0.0055, 0), check.names = FALSE)
res <- bsc_run_scenarios(p, sc, final_time = 110)     # long data.frame + attr(res, "exposure")

p   <- bsc_set_params(p, list("Liver|cl" = 0.15, "P.S.0|ps" = 0.015), estimate = TRUE)
p$observed <- data.frame(compartment = "Liver", time = c(1, 4, 24), conc = c(8, 19, 40))
fit <- bsc_fit(p)                                     # fit$estimates: estimate, se, rse_percent, ci_*
```

`tools/example_batch_evaporation.R` re-creates the evaporation simulation of
Fig. 4 this way. Main functions: `bsc_new_project`, `bsc_add_compartment`,
`bsc_add_flow`, `bsc_add_perm`, `bsc_set_param(s)`, `bsc_set_dosing`,
`bsc_set_exchange`, `bsc_parameter_table`, `bsc_validate`, `bsc_simulate`,
`bsc_run_scenarios`, `bsc_exposure`, `bsc_fit`, `bsc_apply_estimates`,
`bsc_read_project`, `bsc_write_project`.

## Model and numerics

Units: h, mL, µM (= nmol/mL); CL, Q, PS and Kev in mL/h; doses in the chosen
unit (mg are converted with the compound molar mass).

For compartment *i* with concentration *C_i* and volume *V_i(t)*:

```
d(V_i C_i)/dt = sum_j Q_ji C_j - sum_j Q_ij C_i          (fluid flows)
              + sum_j PS_ij (fu_j C_j - fu_i C_i)         (permeable surfaces)
              - CL_i fu_i C_i                             (clearance)
dV_i/dt       = -Kev_i                                    (evaporation)
```

i.e. `dC/dt = (A / V(t)) C` with `A = (-diag(CL) + PERM) .* FU + FLOW + diag(Kev)`,
`FU[i, j] = fu_j`. Solid compartments (chip material) exchange with fluid
compartments through permeable surfaces (non-specific binding, NSB). The
"fu" field of a solid compartment holds the partition coefficient
`K_media/plastic`, and the PS of such a surface is the NSB exchange rate
`PS_NSB`: at equilibrium `C_medium = K_media/plastic * C_solid` (fu = 1 in the
medium), so a smaller `K_media/plastic` means stronger binding. The app labels
these parameters accordingly, using the same notation as the manuscript.

* The system is integrated with `deSolve::ode` (default `lsoda`, rtol = atol =
  1e-6; `ode45` available) piecewise between intervention times (doses, media
  exchange, sampling). Within an interval `V(t) = V_start - Kev (t - t_start)`.
  At an intervention: `V = V_end - removed + added`,
  `C = (C_end (V_end - removed) + dose) / (V_end - removed + added)`.
  Events at the same time are merged; doses at t = 0 give the initial
  concentration.
* Profiles are reported on a 0.01 h grid up to the simulation time + 0.1 h.
* Fitting: residuals `prediction - observation` (unweighted) minimised with
  Levenberg-Marquardt (`minpack.lm::nls.lm`, default 50 iterations). By default
  the prediction is the nearest point of a 0.1 h solution grid (identical to
  the R code used for the paper); *exact sampling times* is an option.
  Standard errors from the Jacobian at the solution (`summary.nls.lm`), 95 % CI
  = estimate ± t(0.975, n - p) SE (n observations, p parameters), RSE = SE /
  estimate.
  Default bounds 1e-4 x and 1000 x the initial value.

## Validation

The engine was checked against the analysis code of the paper.
`tests/reference/legacy_reference.rds` holds reference outputs of the original
analysis scripts (`paper_code/R/engine.R` is a cleaned version of them with
identical results); `tests/validation_report.R` prints the comparison. Results
(R 4.2.3, deSolve 1.40, minpack.lm 1.2-4):

| Check | Result |
|---|---|
| Evaporation simulation (Fig. 4A: CL 0.1, PS 0.01 mL/h, V 0.35/0.8/1.6/0.2 mL, Kev apical 0.003 / liver 0.0055 mL/h, 110 h), liver and apical profiles | max. relative difference **0** (bit-identical), with and without evaporation |
| Sampling simulation (Fig. 5A, 0 - 0.15 mL samples) | max. difference 0 |
| NSB simulation (Fig. 6A), fluid compartments | max. relative difference 0 (K = 1), 4e-15 (K = 0.1), 3e-9 (K = 3) |
| Fit D_G-L_SD_PO_1000uM (V_liver 1.6 mL) | CL 0.07720042, PS 0.01729556 mL/h, RSE 7.560 / 8.936 % - identical to the original analysis |
| Same arm with V_liver 1.8 mL (as in Fig. 3) | CL 0.07535944 mL/h = 2.093318 µL/min/10^6 cells, P = 16.04919e-6 cm/s, RSE 7.784 / 9.092 % = Fig. 3 to all digits |
| Fit D_G-L_DD_PO_100uM (daily dosing and media exchange; the daily phase starts from the first measured liver concentration, as in the original analysis) | identical to the original analysis (V_liver 1.6 and 1.8 mL). Fig. 3 uses a revised set-up for the daily-dose arms (drug-free liver medium at the start of the daily phase, see `paper_code/scripts/fig3_fits.R`). |

Run the tests (from `tests/testthat`):

```
Rscript -e "testthat::test_dir('.')"
```

They cover the engine (analytical solutions, mass balance, events, units,
files, synthetic fits incl. K_media/plastic and PS_NSB), the
equivalence with the analysis code of the paper, every Shiny module (`shiny::testServer`,
incl. the editor requests), the whole app server, and an HTTP smoke test that
starts the app in a background process.

## Relation to the analysis code of the paper (`paper_code/`)

* The fit objective is the one of the analysis code (`AnalysisODE()` in
  `paper_code/R/engine.R`): nearest 0.1 h grid point, observations at an
  intervention time compared with the value just before it.
* `FU[i, j] = fu_j` as in `AnalysisODE()`. `ConstructODE()` (simulations)
  uses `fu_i`; for the paper's models (fu = 1 in all fluid compartments) this
  only rescales the state of solid compartments (value of the analysis code =
  fu_solid x C_solid); fluid profiles are identical.
* Custom dosing tables have column 1 time, column 2 dose;
  `ConstructDoseMediaEx()` of the analysis code reads rows.
* Exposure parameters: no smoothing; half-life = time from Tmax to Cmax/2, NA
  if not reached.
* Events after the end of the simulation are ignored; an empty event at t = 0
  is added if needed. A removal at t = 0 has no effect (the initial volume is
  used), as in the analysis code.
* Compartment names are matched exactly (the analysis code uses substring
  matching).

## Folder layout

```
app.R                     Shiny app (page with the four tabs)
R/engine_project.R        project structure, editing, validation, .rds/.json
R/engine_model.R          parameter table, ODE matrix, intervention schedule
R/engine_solver.R         piecewise integration, simulation, exposure parameters
R/engine_fit.R            parameter estimation
R/engine_import.R         data import, diagram layout, examples
R/engine_batch.R          helpers for analysis scripts (parameter sets, scenarios)
R/mod_model_structure.R   tab 1        R/mod_diagram.R         diagram / editor
R/mod_conditions.R        tab 2        R/mod_element_editor.R  properties panel
R/mod_simulation.R        tab 3        R/utils_ui.R            UI helpers
R/mod_fitting.R           tab 4
www/                      CSS and JavaScript of the editor
inst/examples/            example projects (.rds) and make_examples.R
inst/extdata/             CNBio models (.rds), diclofenac data (CSV)
tools/                    load_engine.R, example_batch_evaporation.R
tests/                    testthat tests, reference outputs of the analysis code, validation_report.R
```

## Licence

Code: GNU General Public License v3.0 or later; data (`inst/extdata/` and the
data embedded in the example projects): CC BY 4.0. See the README and the
LICENSE file at the top of the repository.
