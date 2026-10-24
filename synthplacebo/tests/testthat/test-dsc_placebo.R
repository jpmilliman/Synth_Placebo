test_that("dsc_placebo() returns a classed result with the documented elements", {
  d <- sim_panel()
  fit <- dsc_placebo(d, "unit", "year", "y", "u1", treat_year(), conf_level = 0.80, diagnostics = FALSE, verbose = FALSE)
  expect_s3_class(fit, "dsc_placebo")
  expect_true(all(c("effects", "test", "weights", "placebo", "settings") %in% names(fit)))
  expect_equal(nrow(fit$effects), length(unique(d$year)))
  expect_equal(sum(fit$effects$post), 6)
  expect_true(all(fit$weights >= -1e-8)); expect_equal(sum(fit$weights), 1, tolerance = 1e-6)
  p <- fit$test$p_value[!is.na(fit$test$p_value)]
  expect_true(all(p >= 1 / (fit$settings$n_donors + 1) - 1e-12 & p <= 1 + 1e-12))
})

test_that("a constant effect added after treatment shifts the treated unit's effect exactly", {
  fit0 <- dsc_placebo(sim_panel(), "unit", "year", "y", "u1", treat_year(), conf_level = 0.80, diagnostics = FALSE, verbose = FALSE)
  fit1 <- dsc_placebo(sim_panel(effect = 2), "unit", "year", "y", "u1", treat_year(), conf_level = 0.80, diagnostics = FALSE, verbose = FALSE)
  post <- fit0$effects$post
  expect_equal(fit1$effects$effect[post] - fit0$effects$effect[post], rep(2, sum(post)), tolerance = 1e-5)
  expect_equal(fit1$effects$effect[!post], fit0$effects$effect[!post], tolerance = 1e-6)
})

test_that("too few donors for the level gives infinite intervals, not an error", {
  expect_warning(fit <- dsc_placebo(sim_panel(n_units = 8), "unit", "year", "y", "u1", treat_year(), conf_level = 0.95, diagnostics = FALSE, verbose = FALSE),
                 "intervals are infinite")
  expect_true(all(is.infinite(fit$effects$eff_low_pw[fit$effects$post])))
})

test_that("covariate adjustment and the balance check run (TWFE)", {
  d <- sim_panel(covariate = TRUE)
  fit <- dsc_placebo(d, "unit", "year", "y", "u1", treat_year(), conf_level = 0.80, adjust_covariates = "x", adjust_model = "twfe",
                     diagnostics = FALSE, verbose = FALSE)
  expect_false(is.null(fit$adjustment$gamma))
  expect_equal(fit$balance$covariate, "x")
})

test_that("bad input stops with a message", {
  d <- sim_panel()
  expect_error(dsc_placebo(d, "unit", "year", "y", "nope", treat_year(), diagnostics = FALSE, verbose = FALSE))
  expect_error(dsc_placebo(d[-1, ], "unit", "year", "y", "u1", treat_year(), diagnostics = FALSE, verbose = FALSE))
})
