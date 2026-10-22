# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

suppressMessages(library(haldensify))
pkgload::load_all("thalmle_hermes2", quiet = TRUE)
lam <- exp(seq(-1, -8, length = 50))
set.seed(1); n <- 300; W <- rnorm(n); A <- runif(n, 0, 10)
hd <- haldensify(A = A, W = W, n_bins = 8, grid_type = "equal_range",
                 lambda_seq = lam, max_degree = 1)
cat("(1) datos de smoke_api.R, en observados A[1:5]:\n")
print(predict(hd, new_A = A[1:5], new_W = W[1:5]))
cat("(1b) mismos datos, en observados A[6:10]:\n")
print(predict(hd, new_A = A[6:10], new_W = W[6:10]))
for (trt in c("uniform", "normal")) {
  set.seed(1); d <- simulate_data(500, 1, trt)
  hd2 <- haldensify(A = d$A, W = d$W, n_bins = 8, grid_type = "equal_range",
                    lambda_seq = lam, max_degree = 1)
  est <- predict(hd2, new_A = d$A, new_W = d$W)
  tru <- g0_dgp(d$A, d$W, trt)
  om <- omega_oracle(d$A, trt)
  cat(sprintf("(2) %-7s en observados: media est=%.3f verdadera=%.3f | min est=%.4f verdadero=%.4f | distintos=%d\n",
      trt, mean(est), mean(tru), min(est), min(tru), length(unique(round(est, 6)))))
  cat(sprintf("    pesos omega/g: verdaderos max=%.1f p99=%.1f | estimados max=%.1f p99=%.1f\n",
      max(om / tru), quantile(om / tru, .99), max(om / est), quantile(om / est, .99)))
}
