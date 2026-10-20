# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' k-th order spline basis in the dose, with a constant column first.
#' k = 0: indicators 1[a >= knot]; k >= 1: (a - knot)_+^k / k!
#' @export
spline_basis <- function(a, knots, k = 1L) {
  a <- as.numeric(a)
  B <- if (k == 0L) {
    outer(a, knots, function(x, s) as.numeric(x >= s))
  } else {
    outer(a, knots, function(x, s) pmax(x - s, 0)^k / factorial(k))
  }
  cbind(1, B)
}

#' Knots from observed doses (quantile grid, capped at n_knots).
#' @export
make_knots <- function(A, n_knots = 50L) {
  unique(quantile(A, probs = seq(0.01, 0.99, length.out = n_knots),
                  names = FALSE))
}

#' Moore-Penrose pseudo-inverse of a symmetric PSD matrix.
#' @export
pinv_sym <- function(S, tol = 1e-10) {
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  keep <- e$values > tol * max(e$values)
  V <- e$vectors[, keep, drop = FALSE]
  V %*% (t(V) / e$values[keep])
}
