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

# Rendering the page (theme and component styles via bslib/sass) takes 1-2 s and
# blocks the R process for every user, so it is rendered once and the HTML is
# reused for later page loads (normal rendering if Shiny's internals change).
ui_cached <- local({
  html <- NULL
  function(req) {
    if (is.null(html)) {
      html <<- tryCatch(shiny:::renderPage(ui, 0, FALSE), error = function(e) NULL)
      if (is.null(html)) return(ui)
    }
    shiny:::httpResponse(200L, content = html)
  }
})

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

shinyApp(ui_cached, server)
