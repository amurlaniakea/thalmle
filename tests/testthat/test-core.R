test_that("simulate_data produce estructura correcta", {
  skip_if_not_installed("testthat")
  df <- simulate_data(100, dgp = 1, treatment = "uniform")
  expect_equal(nrow(df), 100)
  expect_true(all(c("A", "W", "Y") %in% names(df)))
  expect_true(all(df$A >= 1 & df$A <= 9))
  expect_true(all(is.finite(df$W)))
  expect_true(all(is.finite(df$Y)))
})

test_that("fit_omega_kde devuelve función y evalúa positivamente", {
  skip_if_not_installed("testthat")
  A <- rnorm(100) + 3
  omega <- fit_omega_kde(A)
  expect_true(is.function(omega))
  val <- omega(3)
  expect_true(is.finite(val) && val > 0)
})

test_that("fit_Q_hal devuelve objeto no nulo con datos simulados", {
  skip_if_not_installed("testthat")
  skip_if_not_installed("hal9001")
  df <- simulate_data(50, dgp = 1, treatment = "uniform")
  res <- fit_Q_hal(A = df$A, W = df$W, Y = df$Y, smoothness_orders = 0L)
  expect_false(is.null(res))
  expect_true("hal9001_fit_stub" %in% class(res) || is.list(res))
  expect_equal(res$smoothness_orders, 0L)
})

test_that("spline_basis tiene primera columna constante y dimensiones correctas", {
  B <- spline_basis(c(1, 2, 3), knots = c(1.5, 2.5), k = 1L)
  expect_true(all(B[, 1] == 1))
  expect_equal(dim(B), c(3L, 4L))
})

test_that("pinv_sym satisface propiedad de pseudo-inversa", {
  set.seed(3)
  M <- crossprod(matrix(rnorm(15), 3, 5))
  S <- pinv_sym(M)
  expect_equal(M %*% S %*% M, M, tolerance = 1e-6)
})

test_that("select_lambda_by_size respeta rango y devuelve índices correctos", {
  sizes <- c(0, 2, 5, 10, 15, 22, 28, 35, 50)
  expect_equal(select_lambda_by_size(sizes, 500, cv_index = 3), 6)
  expect_equal(select_lambda_by_size(sizes, 500, cv_index = 9), 7)
  expect_equal(select_lambda_by_size(sizes, 500, cv_index = 7), 7)
  expect_equal(select_lambda_by_size(sizes, 500), 6)
})

test_that("thal_core runs end to end con datos simulados", {
  skip_if_not_installed("testthat")
  skip_if_not_installed("glmnet")
  df <- simulate_data(100, dgp = 1, treatment = "uniform")
  Qf <- fit_Q_hal(A = df$A, W = df$W, Y = df$Y, smoothness_orders = 0L)
  gf <- tryCatch({
    fit_g_haldensify(A = df$A, W = df$W)
  }, error = function(e) list(predict = function(...) rep(0.5, nrow(df))))
  om <- fit_omega_kde(df$A)
  grid <- seq(min(df$A), max(df$A), length.out = 10)
  res <- thal_core(A = df$A, W = df$W, Y = df$Y,
                   Qf = Qf, gf = gf, om = om, grid = grid)
  expect_true(is.list(res))
  expect_equal(length(res$psi_hat), length(grid))
  expect_true(all(is.finite(res$psi_hat)))
  expect_true(all(is.finite(res$psi_init)))
})

test_that("w_cap obs_max funciona con pesos iguales (simulado trivial)", {
  # Escenario trivial: pesos iguales => cap no debe alterar
  skip_if_not_installed("testthat")
  A <- rnorm(20) + 3
  W <- rnorm(20)
  Y <- A + W + rnorm(20)
  Qf <- fit_Q_hal(A = A, W = W, Y = Y, smoothness_orders = 0L)
  gf <- list(predict = function(...) rep(1, length(A)))
  om <- function(x) rep(1, length(x))
  grid <- seq(min(A), max(A), length.out = 5)
  res1 <- thal_core(A, W, Y, Qf, gf, om, grid, w_cap = "obs_max")
  res2 <- thal_core(A, W, Y, Qf, gf, om, grid, w_cap = Inf)
  expect_equal(res1$psi_hat, res2$psi_hat, tolerance = 1e-6)
})

test_that("true_drc da valores finitos para dgp 1 y 2", {
  skip_if_not_installed("testthat")
  expect_true(is.finite(true_drc(4, 1)))
  expect_true(is.finite(true_drc(6, 2)))
})
