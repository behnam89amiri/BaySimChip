# Data

Processed data used by the scripts in `paper_code/scripts/`. Concentrations are
in µM, times in h, volumes in mL. Entries such as `< 0.05 [LLOQ]` mark values
below the lower limit of quantification of the LC-MS/MS assay.

Licence: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/),
© 2026 Bayer AG. Please cite the paper and the archived version of this
repository (see `CITATION.cff` in the repository root).

## Overview

| Folder / file | Content | Used by |
|---|---|---|
| `models/Kev2DAnalysis.rds` | BaySimChip model of the CN Bio PhysioMimix TL6 gut-liver MPS (gut apical, gut basolateral, liver, circulation) | Figs. 3-5, media-exchange scenario |
| `models/NSBpara2DAnalysis.rds` | the same model with a solid (device material) compartment for non-specific binding | Figs. 4C-7 |
| `fig3/diclofenac_subsets.xlsx` | diclofenac medium concentrations, CN Bio TL6, 13 study arms | Fig. 3 |
| `fig3/literature_comparators.csv` | two literature values of CL_int and P shown in Fig. 3B | Fig. 3 |
| `fig7/Li2007_Papp_CLint.csv` | Caco-2 permeability and hepatocyte CL_int of 29 drugs (literature) | Fig. 7B |
| `si_s5_memantine/` | memantine concentrations, CN Bio TL6 gut-liver MPS | SI Section S5 |
| `empty_chip_nsb/` | diclofenac and FIAU in empty HUMIMIC Chip2 devices (PDMS and non-PDMS) | Fig. 6D, SI Section S6 |
| `si_s7_fiau/` | FIAU concentrations and albumin / ALT, CN Bio PhysioMimix LC12 liver plates | SI Section S7 |

## Details

### `models/`

R objects (lists) read with `load_model()` (`R/engine.R`). Elements: `CompData`
(compartments: name, volume `Vol`, `fu`, clearance `CL`, evaporation rate `Kev`,
dosing and media-exchange settings), `FlowData` (flows `Q` between
compartments), `PermData` (permeable surfaces, `PS`), `MolarMass`, `FinalTime`
and further settings. Compartments: `Comp0` gut apical, `Comp1` gut
basolateral, `Comp2` liver, `Comp3` circulation, `Comp5` solid compartment
(NSB model only). The stored values of CL and PS (0.08 / 0.03 mL/h) are
placeholders; every script sets the values it uses.

### `fig3/diclofenac_subsets.xlsx`

One sheet per study arm. Sheet name `<flow>_<system>_<regimen>_<route>_<dose>`:
flow D = dynamic / S = static; system G-L = gut-liver, G = gut only, L = liver
only; regimen SD = single dose (samples 0.5-72 h), DD = daily dosing (days
7-11, with a daily exchange of one third of the liver-circuit medium); route
and dose of the apical (PO) or liver (IV) dosing. Columns: `Sample_ID` (vial),
`Sample_Name`, `Conc` (µM), `Location` (Apical, Basolateral, Liver),
`Time_h`, `Regimen`, `Flow`, `Config`, `Route`, `Dose_uM`, `key`. In the
single-organ and static arms the sample names say "40 µM DCF" (plate-layout
label); these arms were dosed with 100 µM, as given in `Dose_uM` and in the
sheet names.

### `fig3/literature_comparators.csv`

Columns `SheetName` (row label in Fig. 3B), `Clint` (µL/min/10^6 hepatocytes),
`P` (10^-6 cm/s), `source`.

### `fig7/Li2007_Papp_CLint.csv`

29 drugs from Li et al., J. Biomol. Screen. 2007, 12, 1084-1091: `Compound`,
`Caco2_Papp_nm_per_s`, `CLint_uL_min_per_1e6_cells`, and the same values in
model units, `CL` and `PS` (mL/h; 0.6 x 10^6 hepatocytes per liver compartment,
0.33 cm² gut insert area).

### `si_s5_memantine/`

Files `MNT-<arm>-<location>.xlsx`, arms `PO-1` and `PO-10` (apical dosing,
1 and 10 µM) and `IV-009` (liver dosing, 0.09 µM); location `Ga` = gut apical,
`Lv` = liver. No header; columns: time (h), concentration (µM), compartment
number (1 = gut apical, 3 = liver). One row per sample (replicate chips).

### `empty_chip_nsb/`

`empty_chip_nsb_dcf.xlsx` (diclofenac) and `empty_chip_nsb_fiau.xlsx` (FIAU):
cell-free HUMIMIC Chip2 circuits made of PDMS or of non-PDMS material, two 48-h
dosing phases (`Dose1`, `Dose2`) followed by a 72-h washout with compound-free
medium, nominal concentrations 0 (vehicle), 0.1, 10 and 100 µM. Sheet `README`
explains the sheet names (`<material>_<compound>_<phase>_Conc_<conc>uM`); one
sheet per material, compound, phase and concentration. Columns used: `sample_name`,
`pdms_status`, `compound`, `phase`, `conc_uM` (nominal), `time_h`, `conc`
(measured, µM).

### `si_s7_fiau/`

* `fiau_lc12_lcms_concentrations.csv` - FIAU concentrations (LC-MS/MS) in the
  culture medium of three CN Bio PhysioMimix LC12 liver plates (12 wells each,
  two wells per dose) and in the dosing solutions; one row per sample, in the
  order of the LC-MS report. Columns: `analytical_block` (block of the LC-MS
  report, `Platte1`-`Platte4`; mapped to the plates in `R/fiau_tktd.R`), `vial`,
  `sample_name` (dose, date and study day, time label), `conc_report` (reported
  concentration in µM, or `< 0.005 [LLOQ]`), `conc_high_range` and
  `conc_low_range` (results of the high- and low-range calibration; the
  reported value is one of them).
* `fiau_lc12_albumin_alt.xlsx` - sheet `Albumin + ALT`: one row per well and
  sampling time: date (Excel date serial; 45012 = start of treatment),
  timepoint label (`d0` ... `d10`, `+ 8h` = end of an 8-h pulse), chip
  (= plate), FIAU concentration (or `0.1 % DMSO` vehicle, `empty` = cell-free
  well), ALT (pg/mL), albumin (ng/mL). Sheet `Viability (3D CTG)`: end-point
  viability (relative luminescence; not used by the scripts).
