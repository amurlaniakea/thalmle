# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Smoke tests for the initial fits. The two HAL wrappers need hal9001 and
# haldensify and are skipped without them (and on CRAN).

test_that("fit_omega_kde returns a positive density on the observed range", {
  set.seed(3)
  A <- runif(200, 0, 10)
  om <- fit_omega_kde(A)
  expect_true(all(om(c(1, 5, 9)) > 0))
  mass <- integrate(om, min(A), max(A))$value
  expect_gt(mass, 0.85)
  expect_lt(mass, 1.05)
})

test_that("fit_Q_hal returns a function of (a, w)", {
  skip_on_cran()
  skip_if_not_installed("hal9001")
  set.seed(1)
  d <- simulate_data(100, dgp = 1, treatment = "uniform")
  Q <- fit_Q_hal(d$A, d$W, d$Y, smoothness_orders = 0L)
  p <- Q(c(2, 5, 8), c(-1, 0, 1))
  expect_type(p, "double")
  expect_length(p, 3)
  expect_true(all(is.finite(p)))
})

test_that("fit_g_haldensify returns a floored density at the training points", {
  skip_on_cran()
  skip_if_not_installed("haldensify")
  set.seed(2)
  d <- simulate_data(200, dgp = 1, treatment = "normal")
  g <- fit_g_haldensify(d$A, d$W)
  gv <- g(d$A, d$W)
  expect_length(gv, 200)
  expect_true(all(gv >= 1e-4))
  expect_error(g(5, 123456), "among the training")
})
