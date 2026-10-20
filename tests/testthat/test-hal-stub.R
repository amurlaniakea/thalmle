test_that("fit_Q_hal returns a working Q function (with hal9001)", {
  A <- runif(30, 0, 10)
  W <- matrix(rnorm(30), 30, 1)
  Y <- rnorm(30)
  Qf <- fit_Q_hal(A, W, Y)
  expect_is(Qf, "function")
  qa <- Qf(5, 1)
  expect_true(is.numeric(qa))
  expect_length(qa, 1)
})

test_that("fit_g_haldensify returns a working g function (haldensify >= 0.2.8)", {
  A <- runif(30, 0, 10)
  W <- matrix(rnorm(30), 30, 1)
  gf <- fit_g_haldensify(A, W)
  expect_is(gf, "function")
  ga <- gf(5, W[1, , drop = FALSE])
  expect_true(is.numeric(ga))
  expect_true(ga > 0)
})

test_that("fit_omega_kde returns a positive density", {
  A <- runif(30, 0, 10)
  omega <- fit_omega_kde(A)
  expect_is(omega, "function")
  expect_true(omega(5) > 0)
})
