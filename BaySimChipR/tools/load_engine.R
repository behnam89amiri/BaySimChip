# =============================================================================
# load_engine.R - load the BaySimChip engine for use in analysis scripts
#
#   source("path/to/BaySimChipR/tools/load_engine.R")
#
# Sources R/engine_*.R into the global environment (no Shiny needed; needs
# deSolve and minpack.lm; readxl and jsonlite only for data import / JSON).
# Afterwards e.g. bsc_read_project(), bsc_simulate(), bsc_fit(),
# bsc_run_scenarios() and bsc_load_example() are available.
# =============================================================================

local({
  this_dir <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) {
      this_dir <- dirname(normalizePath(of))
      break
    }
  }
  if (is.null(this_dir)) {
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grepl("^--file=", args)])
    this_dir <- if (length(f)) dirname(normalizePath(f[1])) else file.path(getwd(), "tools")
  }
  root <- dirname(this_dir)
  files <- sort(list.files(file.path(root, "R"), pattern = "^engine_.*\\.R$", full.names = TRUE))
  if (!length(files)) stop("BaySimChip engine files not found in ", file.path(root, "R"))
  for (f in files) sys.source(f, envir = globalenv())
  options(baysimchip.app_dir = root)
  invisible(root)
})
