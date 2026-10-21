args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; trt <- args[2]; B <- as.integer(args[3])
pkgload::load_all(pkgdir, quiet = TRUE)
set.seed(10)
grid <- seq(1, 9, length.out = 25); tr <- true_drc(grid, 1)
cen <- grid >= 3 & grid <= 7
nps <- c(25L, 200L, 500L)
out <- list(); tm <- numeric(0)
for (b in seq_len(B)) {
  t0 <- Sys.time()
  d <- simulate_data(500, 1, trt)
  r <- tryCatch(suppressWarnings(suppressMessages({
    Qf <- fit_Q_hal(d$A, d$W, d$Y)
    gf <- fit_g_haldensify(d$A, d$W, max_degree = 2L)
    om <- fit_omega_kde(d$A)
    lapply(nps, function(np) {
      set.seed(1000L + b)
      thal_core(d$A, d$W, d$Y, Qf, gf, om, grid, n_proj = np)
    })
  })), error = function(e) { message("rep ", b, " ERROR: ", conditionMessage(e)); NULL })
  tm <- c(tm, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  message(sprintf("rep %d: %.0f s%s", b, tail(tm, 1), if (is.null(r)) " (FALLO)" else ""))
  if (is.null(r)) next
  out[[length(out) + 1L]] <- list(init = r[[1]]$psi_init - tr, plug = r[[1]]$psi_plugin - tr,
        thal = lapply(r, function(x) x$psi_hat - tr), size = sapply(r, function(x) x$size))
}
stopifnot(length(out) > 0)
M <- function(f) do.call(cbind, lapply(out, f))
rep1 <- function(tag, E) cat(sprintf("%-30s mean|bias| = %.4f  RMSE = %.4f  (centro %.4f, colas %.4f)\n",
  tag, mean(abs(rowMeans(E))), sqrt(mean(E^2)),
  sqrt(mean(E[cen, , drop = FALSE]^2)), sqrt(mean(E[!cen, , drop = FALSE]^2))))
cat(sprintf("cfg: trt=%s n=500 replicas validas=%d de %d, tiempo medio %.0f s, pkg=%s\n",
    trt, length(out), B, mean(tm), pkgdir))
rep1("HAL plug-in (sin targeting)", M(function(o) o$init))
rep1("plug-in dirigido", M(function(o) o$plug))
for (k in seq_along(nps)) rep1(sprintf("T-HAL, n_proj=%d", nps[k]), M(function(o) o$thal[[k]]))
cat("mediana del tamano del modelo de trabajo:", median(sapply(out, function(o) o$size[1])), "\n")
