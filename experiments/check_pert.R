# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; trt <- args[2]; B <- as.integer(args[3]); delta <- as.numeric(args[4])
pkgload::load_all(pkgdir, quiet = TRUE)
set.seed(10)
grid <- seq(1, 9, length.out = 25)
tr <- true_drc(grid, 1)
res <- replicate(B, {
  d <- simulate_data(500, 1, trt)
  Qf <- function(a, w) Q0_dgp(a, w, 1) + delta * (a - 5)
  r <- suppressWarnings(thal_core(d$A, d$W, d$Y, Qf,
        function(a, w) g0_dgp(a, w, trt), function(a) omega_oracle(a, trt), grid))
  c(r$psi_init - tr, r$psi_plugin - tr, r$psi_hat - tr)
})
m <- nrow(res) / 3
rep1 <- function(tag, x) cat(sprintf("%-24s mean|bias| = %.4f  RMSE = %.4f\n", tag,
  mean(abs(rowMeans(x))), sqrt(mean(x^2))))
cat(sprintf("cfg: trt=%s B=%d delta=%g pkg=%s\n", trt, B, delta, pkgdir))
rep1("initial (sin targeting)", res[1:m, ])
rep1("plug-in dirigido", res[(m + 1):(2 * m), ])
rep1("T-HAL (psi_hat)", res[(2 * m + 1):(3 * m), ])
