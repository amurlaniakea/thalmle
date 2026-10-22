# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; trt <- args[2]; B <- as.integer(args[3])
pkgload::load_all(pkgdir, quiet = TRUE)
grid <- seq(1, 9, length.out = 25); tr <- true_drc(grid, 1)
cen <- grid >= 3 & grid <= 7
vars <- list(
  k50_sinSigma  = list(nk = 50L,  tol = NULL),
  k50_tol1e_10  = list(nk = 50L,  tol = 1e-10),
  k50_tol1e_6   = list(nk = 50L,  tol = 1e-6),
  k150_sinSigma = list(nk = 150L, tol = NULL),
  k500_sinSigma = list(nk = 500L, tol = NULL),
  k500_tol1e_8  = list(nk = 500L, tol = 1e-8))
res <- list(); tm <- numeric(0)
for (b in seq_len(B)) {
  t0 <- Sys.time()
  set.seed(5000L + b)
  d <- simulate_data(500, 1, trt)
  r <- tryCatch(suppressWarnings(suppressMessages({
    Qf <- fit_Q_hal(d$A, d$W, d$Y)
    gf <- fit_g_haldensify(d$A, d$W, max_degree = 2L)
    om <- fit_omega_kde(d$A)
    if (b == 1) {
      for (nk in c(50L, 500L)) {
        Phi <- spline_basis(d$A, make_knots(d$A, nk), 1L)
        ev <- eigen(crossprod(Phi) / 500, symmetric = TRUE, only.values = TRUE)$values
        cat(sprintf("diag rep1: %d nudos (N=%d) autovalores de Sigma max %.3g min %.3g | rango con tol 1e-10/1e-8/1e-6/1e-4: %s\n",
            nk, ncol(Phi), max(ev), min(ev),
            paste(sapply(c(1e-10, 1e-8, 1e-6, 1e-4), function(t) sum(ev > t * max(ev))), collapse = "/")))
      }
      wo <- om(d$A) / gf(d$A, d$W)
      cat(sprintf("diag rep1: pesos omega/g en observados: max %.1f p99 %.1f\n", max(wo), quantile(wo, .99)))
      for (a in c(1, 5, 9)) {
        wc <- om(a) / gf(rep(a, 500), d$W)
        cat(sprintf("diag rep1: pesos contrafactuales a=%g: max %.1f p99 %.1f\n", a, max(wc), quantile(wc, .99)))
      }
    }
    lapply(vars, function(v) {
      set.seed(1000L + b)
      thal_core(d$A, d$W, d$Y, Qf, gf, om, grid, n_proj = 500L, c1 = 6,
                n_knots = v$nk, sigma_tol = v$tol)
    })
  })), error = function(e) { message("rep ", b, " ERROR: ", conditionMessage(e)); NULL })
  tm <- c(tm, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  message(sprintf("rep %d: %.0f s%s", b, tail(tm, 1), if (is.null(r)) " (FALLO)" else ""))
  if (is.null(r)) next
  res[[length(res) + 1L]] <- list(init = r[[1]]$psi_init - tr,
      thal = lapply(r, function(x) x$psi_hat - tr),
      size = sapply(r, function(x) x$size), cvsize = sapply(r, function(x) x$cv_size))
  saveRDS(res, sprintf("sigma2_%s.rds", trt))
}
stopifnot(length(res) > 0)
M <- function(f) do.call(cbind, lapply(res, f))
rep1 <- function(tag, E) cat(sprintf("%-26s mean|bias| = %.4f  RMSE = %.4f  (centro %.4f, colas %.4f)\n",
  tag, mean(abs(rowMeans(E))), sqrt(mean(E^2)),
  sqrt(mean(E[cen, , drop = FALSE]^2)), sqrt(mean(E[!cen, , drop = FALSE]^2))))
rm1 <- function(e) sqrt(mean(e^2))
hal <- sapply(res, function(o) rm1(o$init))
cat(sprintf("cfg: trt=%s n=500 n_proj=500 replicas validas=%d de %d, tiempo medio %.0f s\n",
    trt, length(res), B, mean(tm)))
rep1("HAL plug-in (sin targeting)", M(function(o) o$init))
for (k in names(vars)) {
  rep1(paste0("T-HAL ", k), M(function(o) o$thal[[k]]))
  dl <- sapply(res, function(o) rm1(o$thal[[k]])) - hal
  cat(sprintf("   pareado vs HAL: %+.3f (EE %.3f) | mejor %d, empate %d, peor %d de %d | tam. mediana %g | tam. CV sin recortar, mediana %g\n",
      mean(dl), sd(dl) / sqrt(length(dl)), sum(dl < -1e-12), sum(abs(dl) <= 1e-12),
      sum(dl > 1e-12), length(dl), median(sapply(res, function(o) o$size[[k]])),
      median(sapply(res, function(o) o$cvsize[[k]]), na.rm = TRUE)))
}
