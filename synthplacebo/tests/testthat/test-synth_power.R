test_that("synth_power() p-values equal dsc_placebo_lto() on data with the effect added", {
  d <- sim_panel(n_units = 8)
  pw <- synth_power(d, "unit", "year", "y", treat_year(), effects = 0.05, effect_scale = "percent", shapes = "constant",
                    alpha = 0.10, verbose = FALSE)
  for (u in c("u2", "u5")) {
    d1 <- d; i <- d1$unit == u & d1$year >= treat_year()[["year"]]; d1$y[i] <- d1$y[i] * 1.05
    fit <- dsc_placebo_lto(d1, "unit", "year", "y", u, treat_year(), sensitivity = FALSE, diagnostics = FALSE, warn = FALSE, verbose = FALSE)
    got <- pw$p_values[pw$p_values$unit == u & pw$p_values$effect == 0.05, ]
    for (tst in c("pointwise", "cumulative", "simultaneous")) {
      row <- match(tst, c("pointwise", "cumulative", "simultaneous"))
      expect_equal(got$p[got$method == "lto" & got$test == tst], fit$test$p_lto_naive[row], tolerance = 1e-10)
      expect_equal(got$p[got$method == "placebo" & got$test == tst], fit$test$p_placebo[row], tolerance = 1e-10)
    }
  }
})

test_that("a placebo-only run gives the same placebo results as the full run", {
  d <- sim_panel(n_units = 8)
  full <- synth_power(d, "unit", "year", "y", treat_year(), effects = c(-0.1, 0.1), shapes = c("constant", "shock"), verbose = FALSE)
  plac <- synth_power(d, "unit", "year", "y", treat_year(), effects = c(-0.1, 0.1), shapes = c("constant", "shock"), methods = "placebo", verbose = FALSE)
  a <- full$p_values[full$p_values$method == "placebo", c("unit", "shape", "effect", "test", "p")]
  b <- plac$p_values[, c("unit", "shape", "effect", "test", "p")]
  m <- merge(a, b, by = c("unit", "shape", "effect", "test"))
  expect_equal(nrow(m), nrow(b)); expect_equal(m$p.x, m$p.y)
})

test_that("synth_power() output: class, false-positive rows, one-sided directions, MDE columns", {
  pw <- synth_power(sim_panel(n_units = 8), "unit", "year", "y", treat_year(), effects = c(-0.2, 0.2), shapes = "constant",
                    alternative = c("two.sided", "one.sided"), verbose = FALSE)
  expect_s3_class(pw, "synth_power")
  z <- pw$power[pw$power$effect == 0, ]
  expect_setequal(unique(z$direction[z$sided == "one.sided"]), c("less", "greater"))
  expect_true(all(c("mde_grid", "mde_interpolated", "max_power_tested", "false_positive") %in% names(pw$mde)))
  expect_true(all(pw$power$power >= 0 & pw$power$power <= 1))
})
