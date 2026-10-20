# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Outcome regression: k-th order HAL (hal9001). Returns Q(a, w).
#' @param A Numeric vector of treatment values.
#' @param W Numeric matrix of confounder values.
#' @param Y Numeric vector of outcomes.
#' @param smoothness_orders Integer vector; default from THAL_QSO env.
#' @return A function (a, w) -> Q.
#' @export
fit_Q_hal <- function(A, W, Y, smoothness_orders = as.integer(Sys.getenv("THAL_QSO", "1"))) {
  fit <- hal9001::fit_hal(X = cbind(A = A, W = W), Y = Y, family = "gaussian",
                          smoothness_orders = smoothness_orders)
  function(a, w) as.numeric(predict(fit, new_data = cbind(A = a, W = w)))
}

# Conditional density g(a | w) via haldensify (v0.2.8).
#  * predict() is always called with trim = FALSE: the default trimming floor is
#    1/sqrt(length(new_A)), so results would depend on the size of each call.
#  * haldensify's density is a step function of a (one value per bin). We evaluate
#    it once per bin and per training W_i and answer later queries by table lookup
#    (g is only ever queried at the training W values; any other w is an error).
#  * An explicit, call-size-independent floor (min_density) is applied.
#' Conditional density g(a | w) via haldensify (v0.2.8).
#' @param A Numeric vector of treatment values.
#' @param W Numeric matrix of confounders.
#' @param n_bins Number of bins.
#' @param max_degree Maximum degree.
#' @param min_density Density floor.
#' @return A function (a, w) -> g.
#' @export
fit_g_haldensify <- function(A, W, n_bins = 10L, max_degree = 2L,
                             min_density = 1e-4) {
  hd <- haldensify::haldensify(A = A, W = W, n_bins = n_bins,
          grid_type = "equal_range",
          lambda_seq = exp(seq(-1, -10, length = 100)),
          max_degree = max_degree)
  edges <- c(hd$breaks, hd$range_a[2])          # breaks = lower edges of the bins
  K <- length(edges) - 1L
  mids <- (edges[-1] + edges[-length(edges)]) / 2
  n <- length(W)
  G <- matrix(suppressMessages(as.numeric(predict(hd,
          new_A = rep(mids, each = n), new_W = rep(W, times = K),
          trim = FALSE))), n, K)
  bin_of <- function(a) findInterval(a, edges, rightmost.closed = TRUE,
                                     all.inside = TRUE)
  g <- function(a, w) {
    idx <- match(w, W)
    if (anyNA(idx)) stop("g(a, w): w must be among the training W values")
    pmax(G[cbind(idx, bin_of(a))], min_density)
  }
  if (length(unique(round(g(A, W), 6))) < 5L)
    stop("haldensify returned a degenerate (near-constant) fit")
  g
}

#' Marginal density of A (stabilising weight omega): Sheather-Jones KDE.
#' @param A Numeric vector of treatment values.
#' @return A function a -> omega(a).
#' @export
fit_omega_kde <- function(A) {
  d <- stats::density(A, bw = "SJ", from = min(A), to = max(A), n = 512)
  stats::approxfun(d$x, d$y, rule = 2)
}
