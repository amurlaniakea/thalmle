# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

# Quick smoke test with the real nuisance estimators (takes about 3 minutes).
# Requires: install.packages(c("hal9001", "haldensify")) and the installed 'thalmle'.
library(thalmle)
set.seed(1)
d    <- simulate_data(500, dgp = 1, treatment = "normal")
Qf   <- fit_Q_hal(d$A, d$W, d$Y, smoothness_orders = 0L)
gf   <- fit_g_haldensify(d$A, d$W)
om   <- fit_omega_kde(d$A)
grid <- seq(1, 9, length.out = 25)
fit  <- thal_core(d$A, d$W, d$Y, Qf, gf, om, grid,
                  n_proj = length(d$A), w_cap = "obs_max")
out <- data.frame(dose = round(grid, 2), truth = round(true_drc(grid, 1), 3),
                  hal_plugin = round(fit$psi_init, 3), thal = round(fit$psi_hat, 3))
print(out, row.names = FALSE)
cat(sprintf("RMSE  HAL plug-in: %.3f   T-HAL: %.3f   (working-model size: %d)\n",
            sqrt(mean((out$hal_plugin - out$truth)^2)),
            sqrt(mean((out$thal - out$truth)^2)), fit$size))
