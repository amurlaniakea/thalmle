# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

suppressMessages(library(haldensify))
pkgload::load_all("thalmle_real", quiet = TRUE)
tm <- function(label, expr) {
  t0 <- Sys.time(); val <- force(expr)
  cat(sprintf("%-46s %6.1f s\n", label, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(val)
}
pf <- formals(haldensify:::predict.haldensify)
cat("defaults de predict: lambda_select =", deparse(pf$lambda_select),
    " trim =", deparse(pf$trim), "\n")
set.seed(10)
d <- simulate_data(500, 1, "normal")
lam <- exp(seq(-1, -10, length = 100))
fit_hd <- function(nb) tm(sprintf("fit haldensify degree=2 n_bins=%d", nb),
  haldensify(A = d$A, W = d$W, n_bins = nb, grid_type = "equal_range",
             lambda_seq = lam, max_degree = 2))
a <- rep(seq(1, 9, length.out = 25), each = 50)
w <- rep(d$W[1:50], times = 25)          # 1.250 pares
hd10 <- fit_hd(10L)
cat("campos del ajuste:", paste(names(hd10), collapse = ", "), "\n")
p0 <- tm("predict 1.250 pares, por defecto",
         suppressMessages(predict(hd10, new_A = a, new_W = w)))
cat("valores distintos de g (w fijo, 25 dosis):",
    length(unique(round(p0[w == d$W[1]], 8))), "\n")
p1 <- tryCatch(tm("predict 1.250 pares, lambda_select='cv'",
        suppressMessages(predict(hd10, new_A = a, new_W = w, lambda_select = "cv"))),
        error = function(e) { cat("ERROR cv:", conditionMessage(e), "\n"); NULL })
if (!is.null(p1) && length(p1) == length(p0))
  cat("diferencia maxima cv vs defecto:", max(abs(p0 - p1)), "\n")
p2 <- tryCatch(tm("predict 1.250 pares, trim=FALSE",
        suppressMessages(predict(hd10, new_A = a, new_W = w, trim = FALSE))),
        error = function(e) { cat("ERROR trim:", conditionMessage(e), "\n"); NULL })
if (!is.null(p2) && length(p2) == length(p0))
  cat("diferencia maxima trim=FALSE vs defecto:", max(abs(p0 - p2)), "\n")
hd5 <- fit_hd(5L)
tm("predict 1.250 pares, n_bins=5",
   suppressMessages(predict(hd5, new_A = a, new_W = w)))
cat("fin\n")
