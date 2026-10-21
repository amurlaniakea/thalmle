args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; trt <- args[2]
reps <- as.integer(strsplit(args[3], ",")[[1]])
pkgload::load_all(pkgdir, quiet = TRUE)
grid <- seq(1, 9, length.out = 25); tr <- true_drc(grid, 1)
stored <- readRDS(sprintf("sigma2_%s.rds", trt))
rm1 <- function(e) sqrt(mean(e^2))
for (b in reps) {
  t0 <- Sys.time()
  set.seed(5000L + b)
  d <- simulate_data(500, 1, trt)
  out <- suppressWarnings(suppressMessages({
    Qf <- fit_Q_hal(d$A, d$W, d$Y)
    gf <- fit_g_haldensify(d$A, d$W, max_degree = 2L)
    om <- fit_omega_kde(d$A)
    set.seed(1000L + b)
    list(gf = gf, om = om,
         r = thal_core(d$A, d$W, d$Y, Qf, gf, om, grid, n_proj = 500L, c1 = 6))
  }))
  r <- out$r
  cat(sprintf("\n===== replica %d (%.0f s) =====\n", b,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  cat(sprintf("reproduccion: max |nuevo - guardado| = %.2e\n",
              max(abs((r$psi_hat - tr) - stored[[b]]$thal$k50_sinSigma))))
  cat(sprintf("RMSE: HAL %.3f | plug-in dirigido %.3f | T-HAL %.3f\n",
              rm1(r$psi_init - tr), rm1(r$psi_plugin - tr), rm1(r$psi_hat - tr)))
  cat(sprintf("observaciones con A<2: %d, A>8: %d | rango de A: %.2f a %.2f\n",
              sum(d$A < 2), sum(d$A > 8), min(d$A), max(d$A)))
  wo <- out$om(d$A) / out$gf(d$A, d$W)
  w1 <- out$om(1) / out$gf(rep(1, 500), d$W)
  w9 <- out$om(9) / out$gf(rep(9, 500), d$W)
  cat(sprintf("pesos omega/g: observados max %.1f | contrafactuales a=1 max %.1f | a=9 max %.1f (limite 103)\n",
              max(wo), max(w1), max(w9)))
  cat(sprintf("targeting: tamano %d, tamano CV sin recortar %d, max|eps| = %.3g\n",
              r$size, r$cv_size, max(abs(r$eps))))
  sel <- r$selected; kn <- c(NA, make_knots(d$A, 50L))
  ord <- order(abs(r$eps[sel]), decreasing = TRUE)[seq_len(min(8, length(sel)))]
  print(data.frame(j = sel[ord], nudo = round(kn[sel][ord], 2),
                   eps = signif(r$eps[sel][ord], 3)), row.names = FALSE)
  ix <- c(1, 7, 13, 19, 25)
  E <- rbind(HAL = r$psi_init - tr, plugin_dirigido = r$psi_plugin - tr,
             THAL = r$psi_hat - tr)[, ix]
  colnames(E) <- paste0("a=", round(grid[ix], 1))
  cat("error por dosis:\n"); print(round(E, 2))
}
