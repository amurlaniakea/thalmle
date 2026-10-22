# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

rm1 <- function(e) sqrt(mean(e^2))
for (trt in c("uniform", "normal")) {
  res <- readRDS(sprintf("sigma2_%s.rds", trt))
  nm <- names(res[[1]]$thal)
  pw <- function(f) sqrt(rowMeans(sapply(res, f)^2))
  hal_rep <- sapply(res, function(o) rm1(o$init))
  cat(sprintf("\n== %s (replicas = %d) ==\n", trt, length(res)))
  cat(sprintf("%-15s %9s %10s %9s %10s %9s\n", "variante", "RMSE_pto", "media_dif", "mediana", "media_10%", "reps>0.5"))
  cat(sprintf("%-15s %9.3f\n", "HAL", mean(pw(function(o) o$init))))
  for (k in nm) {
    dl <- sapply(res, function(o) rm1(o$thal[[k]])) - hal_rep
    cat(sprintf("%-15s %9.3f %+10.3f %+9.3f %+10.3f %9d\n", k,
        mean(pw(function(o) o$thal[[k]])), mean(dl), median(dl), mean(dl, trim = 0.1), sum(dl > 0.5)))
  }
  dl <- sapply(res, function(o) rm1(o$thal[["k50_sinSigma"]])) - hal_rep
  o <- order(dl, decreasing = TRUE)[1:5]
  cat("5 peores replicas de k50_sinSigma (indice: dif):",
      paste(sprintf("%d: %+.2f", o, dl[o]), collapse = " | "), "\n")
}
