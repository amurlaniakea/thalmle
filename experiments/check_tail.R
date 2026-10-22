# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; trt <- args[2]; B <- as.integer(args[3])
pkgload::load_all(pkgdir, quiet = TRUE)
set.seed(10)
grid <- seq(1, 9, length.out = 25); tr <- true_drc(grid, 1)
cen <- grid >= 3 & grid <= 7
vars <- list(
  base   = list(c1 = 6, cap = NULL, oracle = FALSE),
  cap30  = list(c1 = 6, cap = 30,   oracle = FALSE),
  cap10  = list(c1 = 6, cap = 10,   oracle = FALSE),
  c1_0   = list(c1 = 0, cap = NULL, oracle = FALSE),
  g_true = list(c1 = 6, cap = NULL, oracle = TRUE))
res <- list(); tm <- numeric(0)
for (b in seq_len(B)) {
  t0 <- Sys.time()
  d <- simulate_data(500, 1, trt)
  r <- tryCatch(suppressWarnings(suppressMessages({
    Qf <- fit_Q_hal(d$A, d$W, d$Y)
    gf <- fit_g_haldensify(d$A, d$W, max_degree = 2L)
    om <- fit_omega_kde(d$A)
    gt <- function(a, w) g0_dgp(a, w, trt)
    ot <- function(a) omega_oracle(a, trt)
    lapply(vars, function(v) {
      set.seed(1000L + b)
      thal_core(d$A, d$W, d$Y, Qf, if (v$oracle) gt else gf,
                if (v$oracle) ot else om, grid,
                n_proj = 500L, c1 = v$c1, w_cap = v$cap)
    })
  })), error = function(e) { message("rep ", b, " ERROR: ", conditionMessage(e)); NULL })
  tm <- c(tm, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  message(sprintf("rep %d: %.0f s%s", b, tail(tm, 1), if (is.null(r)) " (FALLO)" else ""))
  if (is.null(r)) next
  res[[length(res) + 1L]] <- list(init = r[[1]]$psi_init - tr,
      thal = lapply(r, function(x) x$psi_hat - tr), size = sapply(r, function(x) x$size))
  saveRDS(res, sprintf("tail_%s.rds", trt))
}
stopifnot(length(res) > 0)
M <- function(f) do.call(cbind, lapply(res, f))
rep1 <- function(tag, E) cat(sprintf("%-26s mean|bias| = %.4f  RMSE = %.4f  (centro %.4f, colas %.4f)\n",
  tag, mean(abs(rowMeans(E))), sqrt(mean(E^2)),
  sqrt(mean(E[cen, , drop = FALSE]^2)), sqrt(mean(E[!cen, , drop = FALSE]^2))))
cat(sprintf("cfg: trt=%s n=500 n_proj=500 replicas validas=%d de %d, tiempo medio %.0f s\n",
    trt, length(res), B, mean(tm)))
rep1("HAL plug-in (sin targeting)", M(function(o) o$init))
for (k in seq_along(vars)) rep1(paste0("T-HAL ", names(vars)[k]), M(function(o) o$thal[[k]]))
cat("mediana del tamano:", paste(names(vars), apply(M(function(o) o$size), 1, median), collapse = "  "), "\n")
