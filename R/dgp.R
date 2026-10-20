# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Data-generating processes of Section 7.1 (arXiv:2607.03824).
# NOTE: the extracted text of the paper is garbled in two places; the
# readings used here are ASSUMPTIONS to check against the PDF:
#   * DGP 2:  g(w) = 4*1[w > -1] - 1[w > 0] + 1[w > 1]
#   * normal treatment: A | W ~ N(5 + 1.5 W, sd = 2), truncated to [0, 10]

g_fun_dgp <- function(w, dgp) {
  switch(as.character(dgp),
    "1" = 5 * (w > 0),
    "2" = 4 * (w > -1) - (w > 0) + (w > 1),
    "3" = 3 * sin(7 * pi * w),
    stop("dgp must be 1, 2 or 3"))
}

Q0_dgp <- function(a, w, dgp) {
  gw <- g_fun_dgp(w, dgp)
  3 * sin(1.5 * a) + 0.5 * a + gw + (a - 5) * gw
}

# E[g(W)], W ~ N(0,1), closed form
Eg_dgp <- function(dgp) {
  switch(as.character(dgp),
    "1" = 5 * 0.5,
    "2" = 4 * pnorm(1) - 0.5 + (1 - pnorm(1)),
    "3" = 0,
    stop("dgp must be 1, 2 or 3"))
}

#' True dose-response curve psi0(a) = E_W[Q0(a, W)]
#' @export
true_drc <- function(a, dgp) {
  3 * sin(1.5 * a) + 0.5 * a + (a - 4) * Eg_dgp(dgp)
}

mu_a <- function(w) 5 + 1.5 * w

rtruncnorm01 <- function(mu, sd, lo, hi) {
  out <- rnorm(length(mu), mu, sd)
  bad <- which(out < lo | out > hi)
  while (length(bad)) {
    out[bad] <- rnorm(length(bad), mu[bad], sd)
    bad <- bad[out[bad] < lo | out[bad] > hi]
  }
  out
}

g0_dgp <- function(a, w, treatment = c("uniform", "normal"),
                   sd = 2, lo = 0, hi = 10) {
  treatment <- match.arg(treatment)
  if (treatment == "uniform") return(ifelse(a >= lo & a <= hi, 1 / (hi - lo), 0))
  mu <- mu_a(w)
  dens <- dnorm(a, mu, sd) / (pnorm(hi, mu, sd) - pnorm(lo, mu, sd))
  ifelse(a >= lo & a <= hi, dens, 0)
}

# Marginal density of A (stabilising weight omega), by quadrature over W
omega_oracle <- function(a, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  if (treatment == "uniform") return(ifelse(a >= 0 & a <= 10, 0.1, 0))
  wq <- qnorm((seq_len(1000) - 0.5) / 1000)
  vapply(a, function(ai) mean(g0_dgp(rep(ai, length(wq)), wq, "normal")),
         numeric(1))
}

#' @export
simulate_data <- function(n, dgp = 1, treatment = c("uniform", "normal")) {
  treatment <- match.arg(treatment)
  W <- rnorm(n)
  A <- if (treatment == "uniform") runif(n, 0, 10) else rtruncnorm01(mu_a(W), 2, 0, 10)
  Y <- Q0_dgp(A, W, dgp) + rnorm(n)
  data.frame(W = W, A = A, Y = Y)
}
