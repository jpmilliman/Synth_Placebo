test_that("lto_bound() matches Lei and Sudijono's bound at known values", {
  expect_equal(lto_bound(16, 0.05), 1 / 16)
  expect_equal(lto_bound(8, 0.10), 1 / 8)
  expect_true(lto_bound(40, 0.10) >= 0.10)
})

test_that("lto_shift() adds an effect exactly as refitting the shifted data does", {
  d <- sim_panel(n_units = 9)
  Y <- matrix(d$y, ncol = 9); post <- seq_len(nrow(Y)) > 16
  tau <- seq(0.5, 3, length.out = sum(post))
  Y1 <- Y; Y1[post, 1] <- Y1[post, 1] + tau
  c_shift <- synthplacebo:::lto_shift(synthplacebo:::lto_core(Y, 1, post), tau)
  c_refit <- synthplacebo:::lto_core(Y1, 1, post)
  expect_equal(synthplacebo:::lto_pvalues(c_shift, 0), synthplacebo:::lto_pvalues(c_refit, 0), tolerance = 1e-8)
  expect_equal(c_shift$G, c_refit$G, tolerance = 1e-6)
})

test_that("one-sided p-values: mirror symmetry and the placebo count", {
  d <- sim_panel(n_units = 9)
  Y <- matrix(d$y, ncol = 9); post <- seq_len(nrow(Y)) > 16
  a <- synthplacebo:::lto_pvalues_onesided(synthplacebo:::lto_shift(synthplacebo:::lto_core(Y, 2, post), 1), "greater")
  b <- synthplacebo:::lto_pvalues_onesided(synthplacebo:::lto_shift(synthplacebo:::lto_core(-Y, 2, post), -1), "less")
  expect_equal(unlist(a), unlist(b), tolerance = 1e-10)
  core <- synthplacebo:::lto_shift(synthplacebo:::lto_core(Y, 2, post), 1); st <- core$standard
  kd <- st$gaps[post, ] / rep(st$s_donors, each = sum(post)); kt <- st$g / st$s
  manual <- (1 + rowSums(kd >= kt)) / (ncol(kd) + 1)
  expect_equal(synthplacebo:::lto_pvalues_onesided(core, "greater")$p_placebo, manual)
})

test_that("dsc_placebo_lto() returns a classed result with finite intervals for a small pool", {
  fit <- dsc_placebo_lto(sim_panel(n_units = 10), "unit", "year", "y", "u1", treat_year(), conf_level = 0.90,
                         diagnostics = FALSE, warn = FALSE, verbose = FALSE)
  expect_s3_class(fit, "dsc_placebo_lto")
  expect_equal(nrow(fit$test), 4)
  expect_true(all(fit$test$p_lto_naive >= 0 & fit$test$p_lto_naive <= 1))
  expect_true(all(fit$test$gamma_sensitivity >= 1))
  expect_true(all(is.finite(fit$effects$eff_low[fit$effects$post])))
  expect_equal(fit$settings$treated_unit, "u1")
})

test_that("the small-pool warning fires below 10 units with the powered version", {
  expect_warning(dsc_placebo_lto(sim_panel(n_units = 7), "unit", "year", "y", "u1", treat_year(), version = "powered",
                                 sensitivity = FALSE, diagnostics = FALSE, verbose = FALSE), "Small pool|powered version over-rejected|units")
})
