# A small simulated panel: two factors, unit effects, AR(1)-free noise; unit "u1" is the treated unit.
sim_panel <- function(n_units = 10, T0 = 16, T1 = 6, seed = 1, effect = 0, covariate = FALSE) {
  set.seed(seed)
  TT <- T0 + T1
  f <- cbind(cumsum(stats::rnorm(TT, 0.05, 0.3)), sin(seq_len(TT) / 3))
  lam <- matrix(stats::rnorm(n_units * 2), n_units, 2)
  y <- sapply(seq_len(n_units), function(i) 10 + i / 5 + f %*% lam[i, ] + stats::rnorm(TT, 0, 0.15))
  y[(T0 + 1):TT, 1] <- y[(T0 + 1):TT, 1] + effect
  d <- data.frame(unit = rep(paste0("u", seq_len(n_units)), each = TT), year = rep(2000 + seq_len(TT), n_units), y = as.vector(y))
  if (covariate) d$x <- as.vector(sapply(seq_len(n_units), function(i) 5 + 0.3 * f[, 1] * lam[i, 1] + stats::rnorm(TT, 0, 0.1)))
  d
}
treat_year <- function(T0 = 16) c(year = 2000 + T0 + 1)
