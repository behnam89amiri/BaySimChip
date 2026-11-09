# =============================================================================
# engine_import.R - diagram layout, data files and example projects
#
# (The app reads only the native .rds format and JSON.)
# =============================================================================

#' Default diagram positions for compartments without coordinates
#'
#' Fluid compartments are placed on a circle; solid compartments are placed
#' outside, next to the compartment they exchange with.
bsc_auto_layout <- function(p, radius = 160) {
  cp <- p$compartments
  if (!nrow(cp)) return(p)
  miss <- !is.finite(cp$x) | !is.finite(cp$y)
  if (!any(miss)) return(p)
  fluid <- which(cp$type == "fluid")
  nf <- length(fluid)
  for (k in seq_along(fluid)) {
    i <- fluid[k]
    if (!miss[i]) next
    ang <- base::pi / 2 + 2 * base::pi * (k - 1) / max(nf, 1)
    cp$x[i] <- if (nf == 1) 0 else round(radius * cos(ang))
    cp$y[i] <- if (nf == 1) 0 else round(-radius * sin(ang))
  }
  solid <- which(cp$type == "solid")
  for (k in seq_along(solid)) {
    i <- solid[k]
    if (!miss[i]) next
    nb <- c(p$perms$between2[p$perms$between1 == cp$name[i]],
            p$perms$between1[p$perms$between2 == cp$name[i]])
    j <- match(nb[1], cp$name)
    if (!is.na(j) && is.finite(cp$x[j])) {
      len <- sqrt(cp$x[j]^2 + cp$y[j]^2)
      ux <- if (len > 0) cp$x[j] / len else 1
      uy <- if (len > 0) cp$y[j] / len else 0
      cp$x[i] <- round(cp$x[j] + 130 * ux)
      cp$y[i] <- round(cp$y[j] + 130 * uy)
    } else {
      cp$x[i] <- round(radius * 1.8)
      cp$y[i] <- round(-radius + 90 * k)
    }
  }
  p$compartments <- cp
  p
}

# --- data files -------------------------------------------------------------

#' Sheet names of an Excel file (NULL for other files)
bsc_data_sheets <- function(path, name = path) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx", "xls")) readxl::excel_sheets(path) else NULL
}

#' Read a CSV / TSV / Excel file as a data.frame of character columns
bsc_read_data_file <- function(path, sheet = NULL, name = path) {
  ext <- tolower(tools::file_ext(name))
  df <- switch(ext,
               xlsx = , xls = as.data.frame(readxl::read_excel(path, sheet = sheet %||% 1,
                                                               col_types = "text"),
                                            stringsAsFactors = FALSE),
               csv = utils::read.csv(path, colClasses = "character", check.names = FALSE,
                                     stringsAsFactors = FALSE, strip.white = TRUE),
               tsv = , txt = utils::read.delim(path, colClasses = "character", check.names = FALSE,
                                               stringsAsFactors = FALSE, strip.white = TRUE),
               stop("Unsupported data file type '.", ext, "' (use .csv, .txt, .xlsx or .xls).",
                    call. = FALSE))
  # CSV with semicolons (e.g. German Excel export)
  if (ext == "csv" && ncol(df) == 1 && any(grepl(";", names(df)))) {
    df <- utils::read.csv2(path, colClasses = "character", check.names = FALSE,
                           stringsAsFactors = FALSE, strip.white = TRUE)
  }
  df
}

#' Parse concentrations; "< LLOQ" style entries are below the limit of
#' quantification (BLQ) and are set to 0 (as in the Fig. 3 fits), NA,
#' or half the given limit.
bsc_parse_conc <- function(x, blq = c("zero", "na", "half")) {
  blq <- match.arg(blq)
  x <- trimws(as.character(x))
  is_blq <- grepl("^<", x) | toupper(x) %in% c("BLQ", "BLOQ", "<LLOQ", "LLOQ")
  num <- suppressWarnings(as.numeric(x))
  lim <- suppressWarnings(as.numeric(trimws(sub("^<\\s*=?", "", x))))
  num[is_blq] <- switch(blq, zero = 0, na = NA_real_, half = lim[is_blq] / 2)
  num
}

#' Guess the time, concentration and compartment columns of a data table
bsc_guess_columns <- function(df) {
  nm <- names(df)
  low <- tolower(nm)
  pick <- function(pattern, exclude = character(0)) {
    i <- which(grepl(pattern, low) & !nm %in% exclude)
    if (length(i)) nm[i[1]] else NA_character_
  }
  time <- pick("^time|time_h|^t$|^t \\(|hour")
  conc <- pick("^conc|concentration|^c$|^c \\(|^y$|value", exclude = time)
  comp <- pick("compartment|location|^comp|organ|matrix", exclude = c(time, conc))
  if (is.na(time) && length(nm) >= 1) time <- nm[1]
  if (is.na(conc) && length(nm) >= 2) conc <- setdiff(nm, time)[1]
  list(time = time, conc = conc, compartment = comp)
}

#' Build observed data (compartment, time, conc) from an imported table
#'
#' @param df data.frame read by [bsc_read_data_file()]
#' @param time_col,conc_col column names
#' @param comp_col optional column identifying the compartment; its values
#'   are translated with `comp_map` (named character vector: value in file ->
#'   compartment name; values mapped to "" or NA are dropped).  If `comp_col`
#'   is NULL all rows are assigned to `compartment`.
#' @param blq handling of "< LLOQ" entries, see [bsc_parse_conc()]
bsc_observed_from_table <- function(df, time_col, conc_col, comp_col = NULL, comp_map = NULL,
                                    compartment = NULL, blq = "zero") {
  if (!time_col %in% names(df) || !conc_col %in% names(df)) {
    stop("Time or concentration column not found in the data.", call. = FALSE)
  }
  tt <- suppressWarnings(as.numeric(trimws(df[[time_col]])))
  cc <- bsc_parse_conc(df[[conc_col]], blq = blq)
  if (!is.null(comp_col) && nzchar(comp_col) && comp_col %in% names(df)) {
    raw <- as.character(df[[comp_col]])
    comp <- if (is.null(comp_map)) raw else unname(comp_map[raw])
  } else {
    if (is.null(compartment)) stop("Select the compartment the data belong to.", call. = FALSE)
    comp <- rep(compartment, nrow(df))
  }
  out <- data.frame(compartment = comp, time = tt, conc = cc, stringsAsFactors = FALSE)
  out <- out[!is.na(out$compartment) & nzchar(out$compartment) &
               is.finite(out$time) & is.finite(out$conc), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Read a numeric table (e.g. custom doses or media exchange) from a file;
#' the first `ncol` numeric columns are used.
bsc_read_numeric_table <- function(path, ncol, sheet = NULL, name = path) {
  df <- bsc_read_data_file(path, sheet = sheet, name = name)
  num <- lapply(df, function(x) suppressWarnings(as.numeric(trimws(x))))
  ok <- vapply(num, function(x) any(is.finite(x)), logical(1))
  num <- num[ok]
  if (length(num) < ncol) stop("The file needs at least ", ncol, " numeric columns.", call. = FALSE)
  m <- do.call(cbind, num[seq_len(ncol)])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  unname(m)
}

# --- examples ---------------------------------------------------------------

#' Folder with the example projects
bsc_example_dir <- function() {
  d <- getOption("baysimchip.app_dir", ".")
  file.path(d, "inst", "examples")
}

#' Names of the bundled example projects
bsc_list_examples <- function() {
  f <- list.files(bsc_example_dir(), pattern = "\\.rds$", full.names = FALSE)
  sort(tools::file_path_sans_ext(f))
}

#' Example projects as choices for a select input: label = project name,
#' value = file name (sorted by label)
bsc_example_choices <- function() {
  ex <- bsc_list_examples()
  lab <- vapply(ex, function(e) tryCatch(bsc_load_example(e)$name, error = function(err) ""), "")
  lab <- ifelse(nzchar(lab), lab, ex)
  o <- order(tolower(lab))
  stats::setNames(ex[o], lab[o])
}

#' Load a bundled example project
bsc_load_example <- function(name) {
  path <- file.path(bsc_example_dir(), paste0(name, ".rds"))
  if (!file.exists(path)) stop("Example '", name, "' not found.", call. = FALSE)
  bsc_read_project(path)
}
