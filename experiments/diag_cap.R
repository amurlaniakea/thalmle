# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; trt <- args[2]
reps <- as.integer(strsplit(args[3], ",")[[1]])
pkgload::load_all(pkgdir, quiet = TRUE)
grid <- seq(1, 9, length.out = 25); tr <- true_drc(grid, 1)
rm1 <- function(e) sqrt(mean(e^2))
for (b in reps) {
  set.seed(5000L + b)
  d <- simulate_data(500, 1, trt)
  suppressWarnings(suppressMessages({
    Qf <- fit_Q_hal(d$A, d$W, d$Y)
    gf <- fit_g_haldensify(d$A, d$W, max_degree = 2L)
    om <- fit_omega_kde(d$A)
  }))
  wmax <- max(om(d$A) / gf(d$A, d$W))
  caps <- c(defecto = NA, obs_max = wmax, obs_max_x2 = 2 * wmax, cap10 = 10)
  cat(sprintf("\n===== replica %d | peso maximo observado %.1f =====\n", b, wmax))
  for (nm in names(caps)) {
    set.seed(1000L + b)
    r <- suppressWarnings(suppressMessages(
      thal_core(d$A, d$W, d$Y, Qf, gf, om, grid, n_proj = 500L, c1 = 6,
                w_cap = if (is.na(caps[[nm]])) NULL else caps[[nm]])))
    if (nm == "defecto")
      cat(sprintf("HAL: RMSE %.3f | error a=1 %+.2f, a=9 %+.2f\n", rm1(r$psi_init - tr),
                  (r$psi_init - tr)[1], (r$psi_init - tr)[25]))
    cat(sprintf("%-11s | plug-in dirigido RMSE %.3f | T-HAL RMSE %.3f | error a=1 %+.2f, a=9 %+.2f\n",
                nm, rm1(r$psi_plugin - tr), rm1(r$psi_hat - tr),
                (r$psi_hat - tr)[1], (r$psi_hat - tr)[25]))
  }
  for (a in c(1, 9)) {
    gh <- gf(rep(a, 500), d$W); gt <- g0_dgp(rep(a, 500), d$W, trt)
    cat(sprintf("g en a=%g: mediana(g_est/g_verd) = %.2f | peso max estimado %.1f vs verdadero %.1f\n",
                a, median(gh / gt), max(om(a) / gh), max(omega_oracle(a, trt) / gt)))
  }
}
