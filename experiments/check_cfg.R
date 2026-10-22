# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; dgp <- as.integer(args[2]); trt <- args[3]
qi <- args[4]; B <- as.integer(args[5]); sd0 <- as.integer(args[6])
pkgload::load_all(pkgdir, quiet = TRUE)
set.seed(sd0)
res <- replicate(B, {
  r <- thal_oracle(n = 500, dgp = dgp, treatment = trt, q_init = qi)
  c(r$psi_init - r$truth, r$psi_plugin - r$truth, r$psi_hat - r$truth)
})
m <- nrow(res) / 3
init <- res[1:m, ]; plug <- res[(m + 1):(2 * m), ]; thal <- res[(2 * m + 1):(3 * m), ]
rep1 <- function(tag, x) cat(sprintf("%-24s mean|bias| = %.4f  RMSE = %.4f\n", tag,
  mean(abs(rowMeans(x))), sqrt(mean(x^2))))
cat(sprintf("cfg: dgp=%d trt=%s q_init=%s B=%d seed=%d pkg=%s\n", dgp, trt, qi, B, sd0, pkgdir))
rep1("initial (sin targeting)", init)
rep1("plug-in dirigido", plug)
rep1("T-HAL (psi_hat)", thal)
