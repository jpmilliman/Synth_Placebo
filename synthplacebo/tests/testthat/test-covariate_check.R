test_that("a covariate that jumps for the treated unit after treatment is flagged; one that does not, is not", {
  d <- sim_panel(n_units = 12, covariate = TRUE)
  d$jump <- d$x; i <- d$unit == "u1" & d$year >= treat_year()[["year"]]; d$jump[i] <- d$jump[i] + 5
  ck <- covariate_trajectory_check(d, "unit", "year", c("x", "jump"), "u1", treat_year())
  expect_s3_class(ck, "covariate_trajectory_check")
  expect_true(ck$summary$flag[ck$summary$covariate == "jump"])
  expect_equal(ck$summary$p_final[ck$summary$covariate == "jump"], 1 / 12)
  expect_true(all(ck$summary$p_final >= ck$summary$smallest_p - 1e-12))
})

test_that("implied_change is -gamma times the covariate's gap; missing covariates are skipped", {
  d <- sim_panel(n_units = 10, covariate = TRUE); d$z <- d$x; d$z[d$unit == "u1"][1] <- NA
  ck <- covariate_trajectory_check(d, "unit", "year", c("x", "z"), "u1", treat_year(), gamma = c(x = 2))
  expect_equal(ck$summary$implied_change_mean, -2 * ck$summary$mean_post_gap)
  expect_length(ck$skipped, 1)
})
