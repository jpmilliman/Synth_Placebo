test_that("West Germany: the documented LTO and placebo results reproduce", {
  skip_if_not_installed("scpi"); skip_on_cran()
  data("scpi_germany", package = "scpi", envir = environment())
  g <- as.data.frame(scpi_germany)[, c("country", "year", "gdp")]
  fit <- dsc_placebo(g, "country", "year", "gdp", "West Germany", c(year = 1991), conf_level = 0.90, diagnostics = FALSE, verbose = FALSE)
  expect_equal(fit$settings$n_donors, 16)
  lto <- dsc_placebo_lto(g, "country", "year", "gdp", "West Germany", c(year = 1991), sensitivity = FALSE, diagnostics = FALSE, warn = FALSE, verbose = FALSE)
  expect_equal(lto$settings$n_matches, 120)
  expect_lt(lto$test$p_lto_naive[1], 0.05)   # 2003 effect: LTO p = 0.017 in the documentation (2 of 120 matches not won)
})
