# Smoke test: the app starts in a separate R process and serves its page, and
# every new page load gets all the JavaScript and CSS files it refers to.

test_that("the app starts and serves the four tabs and all page files over HTTP", {
  skip_if_not_installed("callr")
  skip_if_not_installed("curl")
  port <- httpuv::randomPort()
  proc <- callr::r_bg(function(dir, port) {
    app <- shiny::shinyAppDir(dir)
    # Test hook: GET /test-widget-jquery does what the DT and plotly widgets do
    # in a running session: they register the jQuery of the jquerylib package
    # under the web path of Shiny's own jQuery ("jquery-3.6.0"), where Shiny's
    # file name (jquery.min.js) does not exist.
    handler <- app$httpHandler
    app$httpHandler <- function(req) {
      if (identical(req$PATH_INFO, "/test-widget-jquery")) {
        dep <- htmltools::resolveDependencies(list(jquerylib::jquery_core(3)))[[1]]
        shiny::createWebDependency(dep)
        return(list(status = 200L, headers = list("Content-Type" = "text/plain"), body = "ok"))
      }
      handler(req)
    }
    shiny::runApp(app, port = port, launch.browser = FALSE, host = "127.0.0.1")
  }, args = list(dir = bsc_root, port = port), supervise = TRUE)
  on.exit(proc$kill(), add = TRUE)
  url <- sprintf("http://127.0.0.1:%d/", port)
  res <- NULL
  deadline <- Sys.time() + 90
  while (Sys.time() < deadline && proc$is_alive()) {
    res <- tryCatch(curl::curl_fetch_memory(url), error = function(e) NULL)
    if (!is.null(res) && res$status_code == 200) break
    Sys.sleep(0.5)
  }
  if (is.null(res)) fail(paste("App did not start:", paste(proc$read_all_error_lines(), collapse = "\n")))
  expect_equal(res$status_code, 200)
  html <- rawToChar(res$content)
  for (tab in c("Model structure", "Experimental conditions", "Simulation", "Fitting and analysis")) {
    expect_true(grepl(tab, html, fixed = TRUE), info = tab)
  }
  expect_true(grepl("BaySimChip", html, fixed = TRUE))

  # loads the page as a new visitor; returns the scripts and stylesheets it
  # refers to that are not served
  files_not_served <- function() {
    page <- rawToChar(curl::curl_fetch_memory(url)$content)
    refs <- regmatches(page, gregexpr('(src|href)="[^"]+[.](js|css)"', page))[[1]]
    refs <- unique(sub('"$', "", sub('^(src|href)="', "", refs)))
    refs <- refs[!grepl("^[a-z]+:", refs)]
    expect_true(any(grepl("jquery", refs)) && "styles.css" %in% refs)
    status <- vapply(refs, function(r) curl::curl_fetch_memory(paste0(url, r))$status_code, numeric(1))
    refs[status != 200]
  }
  expect_identical(files_not_served(), character(0))
  # Regression test: after a session had shown a table or plot, new visitors
  # got a page whose jQuery file was not found (404), so the app did not react
  # (found before the v1.0.0 release; cause: the rendered page was cached).
  expect_equal(curl::curl_fetch_memory(paste0(url, "test-widget-jquery"))$status_code, 200)
  expect_identical(files_not_served(), character(0))
  expect_true(proc$is_alive())
})
