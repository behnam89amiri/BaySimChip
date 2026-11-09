# Smoke test: the app starts in a separate R process and serves its page.

test_that("the app starts and serves the four tabs over HTTP", {
  skip_if_not_installed("callr")
  skip_if_not_installed("curl")
  port <- httpuv::randomPort()
  proc <- callr::r_bg(function(dir, port) {
    shiny::runApp(dir, port = port, launch.browser = FALSE, host = "127.0.0.1")
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
  # later page loads are served from the cached page: same HTML, fast
  t0 <- Sys.time()
  res2 <- curl::curl_fetch_memory(url)
  expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 0.5)
  expect_identical(rawToChar(res2$content), html)
  # the dependencies referenced by the cached page are still served
  dep <- regmatches(html, regexpr('href="[^"]*bslib-component-css[^"]*\\.css"', html))
  if (length(dep)) {
    dep_res <- curl::curl_fetch_memory(paste0(url, sub('^href="', "", sub('"$', "", dep))))
    expect_equal(dep_res$status_code, 200)
  }
  css <- curl::curl_fetch_memory(paste0(url, "styles.css"))
  expect_equal(css$status_code, 200)
  expect_true(proc$is_alive())
})
