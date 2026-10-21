args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; trt <- args[2]; B <- as.integer(args[3])
dgp   <- if (length(args) >= 4) as.integer(args[4]) else 1L
md    <- if (length(args) >= 5) as.integer(args[5]) else 2L
nproj <- if (length(args) >= 6) as.integer(args[6]) else 25L
pkgload::load_all(pkgdir, quiet = TRUE)
set.seed(10)
grid <- seq(1, 9, length.out = 25); tr <- true_drc(grid, dgp)
cen <- grid >= 3 & grid <= 7
out <- list(); fail <- 0L; tm <- numeric(0)
for (b in seq_len(B)) {
  t0 <- Sys.time()
  d <- simulate_data(500, dgp, trt)
  r <- tryCatch(suppressWarnings(suppressMessages({
    Qf <- fit_Q_hal(d$A, d$W, d$Y)
    gf <- fit_g_haldensify(d$A, d$W, max_degree = md)
    om <- fit_omega_kde(d$A)
    thal_core(d$A, d$W, d$Y, Qf, gf, om, grid, n_proj = nproj)
  })), error = function(e) { message("rep ", b, " ERROR: ", conditionMessage(e)); NULL })
  tm <- c(tm, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  message(sprintf("rep %d: %.0f s%s", b, tail(tm, 1), if (is.null(r)) " (FALLO)" else ""))
  if (is.null(r)) { fail <- fail + 1L; next }
  out[[length(out) + 1L]] <- c(r$psi_init - tr, r$psi_plugin - tr, r$psi_hat - tr)
}
if (!length(out)) stop("todas las replicas fallaron")
res <- do.call(cbind, out); m <- length(grid)
sl <- list(init = 1:m, plug = (m + 1):(2 * m), thal = (2 * m + 1):(3 * m))
rep1 <- function(tag, E) cat(sprintf("%-26s mean|bias| = %.4f  RMSE = %.4f  (centro %.4f, colas %.4f)\n",
  tag, mean(abs(rowMeans(E))), sqrt(mean(E^2)), sqrt(mean(E[cen, ]^2)), sqrt(mean(E[!cen, ]^2))))
cat(sprintf("cfg: dgp=%d trt=%s n=500 max_degree=%d n_proj=%d replicas validas=%d de %d, tiempo medio %.0f s, pkg=%s\n",
    dgp, trt, md, nproj, length(out), B, mean(tm), pkgdir))
rep1("HAL plug-in (sin targeting)", res[sl$init, , drop = FALSE])
rep1("plug-in dirigido", res[sl$plug, , drop = FALSE])
rep1("T-HAL (psi_hat)", res[sl$thal, , drop = FALSE])
