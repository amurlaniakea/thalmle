# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; trt <- args[2]; B <- as.integer(args[3])
dgp <- as.integer(args[4]); sb <- as.integer(args[5])
pkgload::load_all(pkgdir, quiet = TRUE)
grid <- seq(1, 9, length.out = 25); tr <- true_drc(grid, dgp)
vars <- list(defecto = list(cap = NULL), obs_max = list(cap = "obs_max"),
             cap10 = list(cap = 10), cap5 = list(cap = 5))
res <- list(); tm <- numeric(0)
for (b in seq_len(B)) {
  t0 <- Sys.time()
  set.seed(sb + b)
  d <- simulate_data(500, dgp, trt)
  r <- tryCatch(suppressWarnings(suppressMessages({
    Qf <- fit_Q_hal(d$A, d$W, d$Y)
    gf <- fit_g_haldensify(d$A, d$W, max_degree = 2L)
    om <- fit_omega_kde(d$A)
    lapply(vars, function(v) {
      set.seed(1000L + b)
      thal_core(d$A, d$W, d$Y, Qf, gf, om, grid, n_proj = 500L, c1 = 6, w_cap = v$cap)
    })
  })), error = function(e) { message("rep ", b, " ERROR: ", conditionMessage(e)); NULL })
  tm <- c(tm, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  message(sprintf("rep %d: %.0f s%s", b, tail(tm, 1), if (is.null(r)) " (FALLO)" else ""))
  if (is.null(r)) next
  res[[length(res) + 1L]] <- list(init = r[[1]]$psi_init - tr,
      thal = lapply(r, function(x) x$psi_hat - tr))
  saveRDS(res, sprintf("capobs2_dgp%d_%s_%d.rds", dgp, trt, sb))
}
stopifnot(length(res) > 0)
rm1 <- function(e) sqrt(mean(e^2))
hal <- sapply(res, function(o) rm1(o$init))
v <- function(k) sapply(res, function(o) rm1(o$thal[[k]]))
pw <- function(f) sqrt(rowMeans(sapply(res, f)^2))
cat(sprintf("cfg: dgp=%d trt=%s n=500 THAL_QSO=%s replicas validas=%d de %d, tiempo medio %.0f s (semillas %d+b)\n",
    dgp, trt, Sys.getenv("THAL_QSO", "(sin definir: orden 1)"), length(res), B, mean(tm), sb))
cat(sprintf("%-10s %9s %16s %9s %9s %6s %6s\n", "vs HAL", "RMSE_pto", "media_dif (EE)",
            "mediana", "reps>0.5", "mejor", "peor"))
cat(sprintf("%-10s %9.3f\n", "HAL", mean(pw(function(o) o$init))))
for (k in names(vars)) {
  dl <- v(k) - hal
  cat(sprintf("%-10s %9.3f %16s %+9.3f %9d %6d %6d\n", k, mean(pw(function(o) o$thal[[k]])),
      sprintf("%+.3f (%.3f)", mean(dl), sd(dl) / sqrt(length(dl))), median(dl),
      sum(dl > 0.5), sum(dl < -1e-12), sum(dl > 1e-12)))
}
cat(sprintf("\n%-10s %16s %9s %6s %6s\n", "vs defecto", "media_dif (EE)", "mediana", "mejor", "peor"))
for (k in names(vars)[-1]) {
  d <- v(k) - v("defecto")
  cat(sprintf("%-10s %16s %+9.3f %6d %6d\n", k,
      sprintf("%+.3f (%.3f)", mean(d), sd(d) / sqrt(length(d))), median(d),
      sum(d < -1e-12), sum(d > 1e-12)))
}
