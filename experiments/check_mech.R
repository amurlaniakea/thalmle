# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

pkgload::load_all("thalmle_hermes2", quiet = TRUE)
grid <- seq(1, 9, length.out = 25); tr <- true_drc(grid, 1)
cen <- grid >= 3 & grid <= 7
run <- function(c1, cap, trt = "normal", delta = 0.5, B = 20) {
  set.seed(10)
  res <- replicate(B, {
    d <- simulate_data(500, 1, trt)
    Qf <- function(a, w) Q0_dgp(a, w, 1) + delta * (a - 5)
    r <- suppressWarnings(thal_core(d$A, d$W, d$Y, Qf,
          function(a, w) g0_dgp(a, w, trt), function(a) omega_oracle(a, trt),
          grid, c1 = c1, w_cap = cap))
    c(r$psi_init - tr, r$psi_hat - tr, r$size)
  })
  m <- length(grid)
  E0 <- res[1:m, ]; E1 <- res[(m + 1):(2 * m), ]; sz <- res[2 * m + 1, ]
  f <- function(E) c(sqrt(mean(E^2)), sqrt(mean(E[cen, ]^2)), sqrt(mean(E[!cen, ]^2)))
  a <- f(E1); b <- f(E0)
  cat(sprintf("c1=%-2g cap=%-7s | T-HAL RMSE todo %.2f centro %.2f colas %.2f | inicial todo %.2f centro %.2f colas %.2f | mediana size %g\n",
      c1, if (is.null(cap)) "default" else format(cap), a[1], a[2], a[3], b[1], b[2], b[3], median(sz)))
}
run(6, NULL)
run(0, NULL)
run(6, 30)
run(0, 30)
run(0, Inf)
