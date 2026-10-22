# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

suppressMessages({library(haldensify)})
pkgload::load_all("thalmle_hermes2", quiet = TRUE)
set.seed(1)
ga <- seq(0.05, 9.95, by = 0.1)
for (trt in c("uniform", "normal")) {
  d <- simulate_data(500, 1, trt)
  for (nb in c(8, 20)) {
    hd <- haldensify(A = d$A, W = d$W, n_bins = nb, grid_type = "equal_range",
                     lambda_seq = exp(seq(-1, -8, length = 50)), max_degree = 1)
    for (w in c(-1, 0, 1)) {
      est <- predict(hd, new_A = ga, new_W = rep(w, length(ga)))
      tru <- g0_dgp(ga, rep(w, length(ga)), trt)
      cr <- suppressWarnings(cor(est, tru))
      cat(sprintf("%-7s n_bins=%-2d w=%+d | integral est=%.3f true=%.3f | max est=%.3f true=%.3f | cor=%.3f | valores distintos=%d\n",
          trt, nb, w, sum(est, na.rm = TRUE) * 0.1, sum(tru) * 0.1,
          max(est, na.rm = TRUE), max(tru), cr, length(unique(round(est, 6)))))
    }
  }
}
