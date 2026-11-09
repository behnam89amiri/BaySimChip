# =============================================================================
# R/setup.R - shared set-up for all paper_code scripts
#
# Every script starts with three lines that find the paper_code root and call
# pc_setup():
#
#   .f    <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
#   .root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
#   source(file.path(.root, "R", "setup.R")); pc_setup(.root)
#
# pc_setup() loads the engine and helper files in R/, and fixes where the
# outputs of this run go:
#   * environment variable PAPER_CODE_RUN_DIR, if set (run_all.R sets it), or
#   * a new folder output/run_YYYY-MM-DD (run_YYYY-MM-DD_2, _3, ... if taken).
# Results are never overwritten: pc_out() refuses to return a path that exists.
# =============================================================================

pc_setup <- function(root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, "run_all.R"))) {
    stop("'", root, "' is not the paper_code root (run_all.R not found). ",
         "Run the scripts from inside paper_code/ or set PAPER_CODE_ROOT.")
  }
  options(paper_code.root = root, stringsAsFactors = FALSE)
  Sys.setenv(PAPER_CODE_ROOT = root)
  for (f in c("engine.R", "scenarios.R", "misspecification.R", "plotting.R")) {
    source(file.path(root, "R", f), local = globalenv())
  }
  invisible(root)
}

pc_root <- function() {
  root <- getOption("paper_code.root")
  if (is.null(root)) stop("Call pc_setup() first.")
  root
}

# Paths inside paper_code/ and paper_code/data/
pc_path <- function(...) file.path(pc_root(), ...)
pc_data <- function(...) {
  p <- file.path(pc_root(), "data", ...)
  if (!file.exists(p)) stop("Input file not found: ", p)
  p
}

# Output folder of the current run (created on first use)
pc_run_dir <- function() {
  rd <- Sys.getenv("PAPER_CODE_RUN_DIR", unset = "")
  if (!nzchar(rd)) {
    base <- file.path(pc_root(), "output", paste0("run_", format(Sys.Date(), "%Y-%m-%d")))
    rd <- base
    k <- 2
    while (dir.exists(rd)) {
      rd <- paste0(base, "_", k)
      k <- k + 1
    }
    Sys.setenv(PAPER_CODE_RUN_DIR = rd)
  }
  if (!dir.exists(rd)) dir.create(rd, recursive = TRUE)
  normalizePath(rd, winslash = "/")
}

# Path of a new output file in the run folder; stops if the file already exists
pc_out <- function(name) {
  p <- file.path(pc_run_dir(), name)
  if (file.exists(p)) {
    stop("Output file already exists (results are never overwritten): ", p,
         "\nUse a new run folder (unset PAPER_CODE_RUN_DIR or point it to an empty folder).")
  }
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  p
}

# Path of an output written earlier in the same run (for scripts that build on
# other scripts, e.g. fig3_plots.R and figures_composite.R)
pc_in_run <- function(name) {
  p <- file.path(pc_run_dir(), name)
  if (!file.exists(p)) {
    stop("Required output of an earlier script not found: ", p,
         "\nRun the producing script first with the same PAPER_CODE_RUN_DIR (see README.md).")
  }
  p
}

# Small logging helper: message with time stamp and script name
pc_log <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
}
