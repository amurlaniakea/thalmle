# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Plumbing tests for the options of thal_core() that the other tests do not
# reach. They use a cheap polynomial Q and the oracle g and omega, so they do
# not need hal9001 or haldensify. They check that the code runs and returns
# finite, well-shaped output; they do not validate the estimator.

local_setup_branches <- function(n = 250) {
  set.seed(21)
  d <- simulate_data(n, dgp = 1, treatment = "uniform")
  mkQ <- function(A, W, Y) {
    m <- lm(Y ~ A + I(A^2) + I(A^3) + W)
    function(a, w) as.numeric(predict(m, data.frame(A = a, W = w)))
  }
  list(d = d, mkQ = mkQ, Q = mkQ(d$A, d$W, d$Y),
       g = function(a, w) thalmle:::g0_dgp(a, w, "uniform"),
       om = function(a) thalmle:::omega_oracle(a, "uniform"),
       grid = seq(1, 9, length.out = 7))
}

# The tail of the lambda path (smallest penalties) is ill-conditioned with this
# nearly collinear basis, and glmnet warns that coordinate descent did not
# converge there. Those lambdas are never selected (the size window stops far
# earlier), so the warning is muffled here and nothing else.
quiet_glmnet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("Convergence for", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

check_fit <- function(r, grid) {
  expect_length(r$psi_hat, length(grid))
  expect_true(all(is.finite(r$psi_hat)))
  expect_true(all(is.finite(r$psi_plugin)))
  expect_equal(r$size, length(r$selected))
}

test_that("honest CV, cross-fitting and both plug-in modes run", {
  skip_if_not_installed("glmnet")
  s <- local_setup_branches()
  quiet_glmnet(with(s, {
    check_fit(thal_core(d$A, d$W, d$Y, Q, g, om, grid, fit_Q = mkQ), grid)
    check_fit(thal_core(d$A, d$W, d$Y, Q, g, om, grid, fit_Q = mkQ, crossfit = TRUE), grid)
    check_fit(thal_core(d$A, d$W, d$Y, Q, g, om, grid, fit_Q = mkQ, crossfit = TRUE,
                        cf_plugin = "full"), grid)
  }))
})

test_that("crossfit = TRUE without fit_Q is an error", {
  skip_if_not_installed("glmnet")
  s <- local_setup_branches()
  expect_error(with(s, thal_core(d$A, d$W, d$Y, Q, g, om, grid, crossfit = TRUE)),
               "requires fit_Q")
})

test_that("sigma_tol, k = 0, numeric and obs_max caps, use_cv = FALSE run", {
  skip_if_not_installed("glmnet")
  s <- local_setup_branches()
  quiet_glmnet(with(s, {
    check_fit(thal_core(d$A, d$W, d$Y, Q, g, om, grid, sigma_tol = 1e-8), grid)
    check_fit(thal_core(d$A, d$W, d$Y, Q, g, om, grid, k = 0L), grid)
    check_fit(thal_core(d$A, d$W, d$Y, Q, g, om, grid, w_cap = 3), grid)
    check_fit(thal_core(d$A, d$W, d$Y, Q, g, om, grid, w_cap = "obs_max"), grid)
    r <- thal_core(d$A, d$W, d$Y, Q, g, om, grid, use_cv = FALSE)
    check_fit(r, grid)
    expect_true(is.na(r$cv_size))
  }))
})

test_that("the default cap is sqrt(n) log(n/5) and Inf disables it", {
  skip_if_not_installed("glmnet")
  s <- local_setup_branches()
  # With uniform treatment and oracle nuisances every weight is 1, so any cap
  # >= 1 must give the same curve.
  quiet_glmnet(with(s, {
    set.seed(5); r0 <- thal_core(d$A, d$W, d$Y, Q, g, om, grid)
    set.seed(5); r1 <- thal_core(d$A, d$W, d$Y, Q, g, om, grid, w_cap = Inf)
    expect_equal(r0$psi_hat, r1$psi_hat, tolerance = 1e-8)
  }))
})

test_that("a grid outside range(A) triggers a warning", {
  skip_if_not_installed("glmnet")
  s <- local_setup_branches()
  expect_warning(
    quiet_glmnet(with(s, thal_core(d$A, d$W, d$Y, Q, g, om, c(-1, 5, 11)))),
    "outside range")
})
