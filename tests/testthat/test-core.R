# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

test_that("true_drc matches Monte Carlo of Q0", {
  set.seed(1)
  w <- rnorm(2e5)
  for (d in 1:3) {
    expect_lt(abs(mean(thalmle:::Q0_dgp(3.7, w, d)) - true_drc(3.7, d)), 0.02)
  }
})

test_that("truncated normal stays in [0, 10]", {
  set.seed(2)
  x <- thalmle:::rtruncnorm01(rnorm(1000, 5, 3), 2, 0, 10)
  expect_true(all(x >= 0 & x <= 10))
})

test_that("pinv_sym is a generalised inverse on rank-deficient input", {
  set.seed(3)
  M <- matrix(rnorm(30), 10, 3) %*% matrix(rnorm(15), 3, 5)
  S <- crossprod(M)
  expect_equal(S %*% pinv_sym(S) %*% S, S, tolerance = 1e-6)
})

test_that("spline_basis has constant first column and right shape", {
  B <- spline_basis(c(1, 2, 3), knots = c(1.5, 2.5), k = 1L)
  expect_equal(dim(B), c(3L, 3L))
  expect_true(all(B[, 1] == 1))
})

test_that("size window selector clamps the CV choice", {
  sizes <- c(0, 2, 5, 10, 15, 22, 28, 35, 50)
  # n = 500, k = 1 -> window [21, 31]
  expect_equal(select_lambda_by_size(sizes, 500, cv_index = 3), 6)
  expect_equal(select_lambda_by_size(sizes, 500, cv_index = 9), 7)
  expect_equal(select_lambda_by_size(sizes, 500, cv_index = 7), 7)
  expect_equal(select_lambda_by_size(sizes, 500), 6)
})

test_that("thal_oracle runs end to end (plumbing only)", {
  skip_if_not_installed("glmnet")
  set.seed(4)
  res <- thal_oracle(n = 300, dgp = 1, treatment = "uniform",
                     q_init = "misspecified")
  expect_length(res$psi_hat, 25)
  expect_true(all(is.finite(res$psi_hat)))
})

test_that("targeting removes most of the bias of a misspecified Q (uniform treatment)", {
  skip_if_not_installed("glmnet")
  set.seed(10)
  r <- thal_oracle(n = 500, dgp = 1, treatment = "uniform", q_init = "misspecified")
  expect_lt(mean(abs(r$psi_hat - r$truth)), 0.5 * mean(abs(r$psi_init - r$truth)))
})

test_that("w_cap = Inf and 'obs_max' are no-ops when every weight equals 1", {
  skip_if_not_installed("glmnet")
  set.seed(10)
  r0 <- thal_oracle(n = 300, dgp = 1, treatment = "uniform", q_init = "misspecified")
  set.seed(10)
  r1 <- thal_oracle(n = 300, dgp = 1, treatment = "uniform", q_init = "misspecified",
                    w_cap = Inf)
  set.seed(10)
  r2 <- thal_oracle(n = 300, dgp = 1, treatment = "uniform", q_init = "misspecified",
                    w_cap = "obs_max")
  expect_equal(r0$psi_hat, r1$psi_hat, tolerance = 1e-8)
  expect_equal(r0$psi_hat, r2$psi_hat, tolerance = 1e-8)
})
