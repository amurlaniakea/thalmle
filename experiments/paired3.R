# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

rm1 <- function(e) sqrt(mean(e^2))
for (trt in c("uniform", "normal")) {
  res <- readRDS(sprintf("capobs_%s.rds", trt))
  v <- function(k) sapply(res, function(o) rm1(o$thal[[k]]))
  hal <- sapply(res, function(o) rm1(o$init))
  cat(sprintf("\n== %s (replicas = %d) ==\n", trt, length(res)))
  for (par in list(c("obs_max", "defecto"), c("cap10", "defecto"))) {
    d <- v(par[1]) - v(par[2])
    cat(sprintf("%-8s - %-8s: media %+.3f (EE %.3f) | mediana %+.3f | mejor %d, peor %d, empate %d\n",
        par[1], par[2], mean(d), sd(d) / sqrt(length(d)), median(d),
        sum(d < -1e-12), sum(d > 1e-12), sum(abs(d) <= 1e-12)))
  }
  d <- v("obs_max") - v("defecto")
  o <- order(d)[1:3]; p <- order(d, decreasing = TRUE)[1:3]
  cat("mayores mejoras de obs_max sobre defecto (rep: dif):", paste(sprintf("%d: %+.2f", o, d[o]), collapse = " | "), "\n")
  cat("mayores empeoramientos de obs_max sobre defecto:   ", paste(sprintf("%d: %+.2f", p, d[p]), collapse = " | "), "\n")
  dl <- v("defecto") - hal; q <- order(dl, decreasing = TRUE)[1:3]
  cat("peores replicas de 'defecto' frente al HAL:        ", paste(sprintf("%d: %+.2f", q, dl[q]), collapse = " | "), "\n")
}
