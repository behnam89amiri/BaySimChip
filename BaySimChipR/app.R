# =============================================================================
# BaySimChip (R Shiny)
#
# Pharmacokinetic modelling of organ-on-chip (OoC) systems: build the
# compartmental model of a chip, define the experimental conditions (dosing,
# sampling, media exchange, evaporation, non-specific binding), simulate
# concentration-time profiles and estimate unknown parameters from data.
#
# Run from R:       shiny::runApp("path/to/BaySimChipR")
# Engine only:      source the files R/engine_*.R (no Shiny needed)
# =============================================================================

library(shiny)
library(bslib)

# Units such as µM need a UTF-8 locale. R uses one by default (on Windows since
# R 4.2), but an invalid setting (e.g. LC_ALL=C.UTF-8 on Windows) makes R fall
# back to the "C" locale: µ is then garbled and some tables stay empty. In that
# case the character handling is switched to UTF-8 while the app runs.
if (!isTRUE(l10n_info()[["UTF-8"]])) {
  old_ctype <- Sys.getlocale("LC_CTYPE")
  for (loc in c(".UTF-8", "C.UTF-8", "en_US.UTF-8")) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  }
  if (isTRUE(l10n_info()[["UTF-8"]])) {
    onStop(function() Sys.setlocale("LC_CTYPE", old_ctype))
  } else {
    warning("R does not run in a UTF-8 locale: units such as uM may be shown ",
            "wrongly and some tables may stay empty.", call. = FALSE)
  }
}

app_dir <- normalizePath(getwd())   # runApp() sets the working dir to the app folder
options(baysimchip.app_dir = app_dir)

# Shiny sources R/*.R automatically (shiny >= 1.5); fall back to explicit
# sourcing if that is switched off (options(shiny.autoload.r = FALSE)).
if (!exists("bsc_new_project", mode = "function")) {
  for (f in sort(list.files(file.path(app_dir, "R"), pattern = "\\.[Rr]$", full.names = TRUE))) {
    source(f, local = TRUE)
  }
}

bsc_theme <- bs_theme(
  version = 5,
  primary = "#2a78d6",
  secondary = "#52514e",
  "font-size-base" = "0.92rem",
  "card-cap-bg" = "#f6f6f3"
)

ui <- page_navbar(
  title = tags$span(
    tags$strong("BaySimChip"),
    tags$small(class = "text-muted ms-2 d-none d-xl-inline", "PK modelling of organ-on-chip systems")
  ),
  id = "main_tabs",
  theme = bsc_theme,
  window_title = "BaySimChip",
  fillable = FALSE,
  header = tags$head(tags$link(rel = "stylesheet", href = "styles.css"),
                     tags$script(src = "bsc.js")),
  nav_panel("Model structure", icon = icon("diagram-project"), mod_model_structure_ui("structure")),
  nav_panel("Experimental conditions", icon = icon("flask"), mod_conditions_ui("conditions")),
  nav_panel("Simulation", icon = icon("chart-area"), mod_simulation_ui("simulation")),
  nav_panel("Fitting and analysis", icon = icon("chart-line"), mod_fitting_ui("fitting")),
  nav_spacer(),
  nav_item(uiOutput("project_badge"))
)

server <- function(input, output, session) {
  start <- getOption("baysimchip.start_example", "CNBio_gut_liver")
  p0 <- tryCatch(bsc_load_example(start), error = function(e) bsc_new_project())
  state <- bsc_new_state(p0)
  bsc_set_project(state, p0)

  mod_model_structure_server("structure", state)
  mod_conditions_server("conditions", state)
  mod_simulation_server("simulation", state)
  mod_fitting_server("fitting", state)

  output$project_badge <- renderUI({
    p <- state$project
    n_unknown <- sum(bsc_parameter_table(p)$estimate)
    tags$span(class = "bsc-badge",
              icon("cube"), p$name,
              if (n_unknown) tags$span(class = "ms-2 text-warning", sprintf("%d unknown", n_unknown)))
  })
}

# Do not cache the rendered page: rendering also registers the web paths of the
# page's JavaScript and CSS files. The DT and plotly widgets register the jQuery
# of the jquerylib package under the same path as Shiny's jQuery ("jquery-3.6.0")
# but with another file name, so a cached page would ask for a file that is no
# longer served (404; the app then does not react for new visitors).
shinyApp(ui, server)
