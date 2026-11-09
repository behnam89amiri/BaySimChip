# BaySimChip

R Shiny app and analysis code for the paper "Quantitative modeling of
platform-specific processes to improve pharmacokinetic studies in
organ-on-a-chip systems" (submitted for publication).

| Folder | Content |
|---|---|
| [`BaySimChipR/`](BaySimChipR/README.md) | The BaySimChip R Shiny app: build a compartmental pharmacokinetic (PK) model of an organ-on-chip (OoC) system, define the experimental conditions (dosing, sampling, media exchange, evaporation, non-specific binding), simulate concentration-time profiles and estimate unknown parameters from measured data. The engine can also be used from R scripts without the app. |
| [`paper_code/`](paper_code/README.md) | R code, processed concentration-time data and results that reproduce the analyses and figures of the paper (Figs. 3-7, SI Figs. S3-S5). |

## Quick start

App (R >= 4.1):

```r
install.packages(c("shiny", "bslib", "DT", "visNetwork", "plotly",
                   "deSolve", "minpack.lm", "jsonlite", "readxl"))
shiny::runApp("BaySimChipR", launch.browser = TRUE)
```

Analyses of the paper (R >= 4.2, from the `paper_code/` folder; the packages are
listed in `paper_code/README.md`):

```
Rscript run_all.R
```

The outputs used for the paper are in `paper_code/output/v1.0.0/`.

## How to cite

Please cite the paper and the archived version of this repository:

* B. Amiri, Ö. Vural, M. Raschke, T. Steger-Hartmann, A. Reichel. Quantitative
  modeling of platform-specific processes to improve pharmacokinetic studies in
  organ-on-a-chip systems. Submitted for publication; DOI to be added.
* BaySimChip v1.0.0, Zenodo. DOI: to be added.

Citation metadata: [`CITATION.cff`](CITATION.cff).

## Licence

Copyright (C) 2026 Bayer AG

The code (the app in `BaySimChipR/` and the scripts in `paper_code/`) is free
software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. It is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License in [`LICENSE`](LICENSE)
for more details.

The data (`paper_code/data/`, `BaySimChipR/inst/extdata/` and the data in the
example projects of the app) are licensed under the Creative Commons
Attribution 4.0 International licence
([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).

### R packages used

The repository contains no third-party code; the packages below are installed
separately by the user (versions used for the paper).

| Package | Version | Licence | Used by |
|---|---|---|---|
| shiny | 1.8.1.1 | GPL-3 | app |
| bslib | 0.7.0 | MIT | app |
| DT | 0.33 | GPL-3 | app |
| visNetwork | 2.1.2 | MIT | app |
| plotly | 4.10.4 | MIT | app |
| htmlwidgets | 1.6.4 | MIT | app |
| jsonlite | 1.8.8 | MIT | app |
| readxl | 1.4.5 | MIT | app, paper code |
| deSolve | 1.40 | GPL (>= 2) | app, paper code |
| minpack.lm | 1.2-4 | GPL-3 | app, paper code |
| ggplot2 | 3.5.2 | MIT | paper code, app batch example |
| patchwork | 1.2.0 | MIT | paper code |
| scales | 1.4.0 | MIT | paper code |
| ragg | 1.3.0 | MIT | paper code |
| dplyr | 1.1.4 | MIT | paper code |
| tidyr | 1.3.1 | MIT | paper code |
| ggrepel | 0.9.5 | GPL-3 | paper code |
| viridisLite | 0.4.2 | MIT | paper code |
| pracma | 2.4.4 | GPL (>= 3) | paper code (one function not used for the figures) |
| testthat | 3.2.1.1 | MIT | app tests |
| callr | 3.7.6 | MIT | app tests |
| curl | 5.2.1 | MIT | app tests |
| httpuv | 1.6.15 | GPL (>= 2) | app tests |

All licences are compatible with GPL-3.0-or-later.

## Authors

Behnam Amiri, Özlem Vural, Marian Raschke, Thomas Steger-Hartmann, Andreas
Reichel (Bayer AG, Berlin, Germany).
