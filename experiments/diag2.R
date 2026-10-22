# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

pkgload::load_all("thalmle", quiet = TRUE)
set.seed(1)
n <- 500
d <- simulate_data(n, 1, "uniform")
A <- d$A; Y <- d$Y
knots <- make_knots(A, 50L)
Phi <- spline_basis(A, knots, 1L)
Q0m <- 3 * sin(1.5 * A) + 0.5 * A
R <- Y - Q0m
S <- crossprod(Phi) / n
ols <- lm.fit(Phi, R)
cat("(a) OLS sin penalizar sobre span(Phi): resid RMSE =", sqrt(mean(ols$residuals^2)), "(esperado ~7.7)\n")
e <- eigen(S, symmetric = TRUE)
Z <- Phi %*% (e$vectors %*% (t(e$vectors) / sqrt(e$values)))
cat("(b) Pn[Z'Z]=I, max desviacion:", max(abs(crossprod(Z)/n - diag(ncol(Z)))), "\n")
show <- function(fit, X, tag) {
  for (i in unique(c(5, 20, 40, length(fit$lambda)))) {
    eps <- as.numeric(fit$beta[, i])
    cat(sprintf("%s idx %3d lambda %.3g df %2d resid RMSE %.3f\n", tag, i,
                fit$lambda[i], fit$df[i], sqrt(mean((R - X %*% eps)^2))))
  }
}
show(glmnet::glmnet(Z, Y, offset = Q0m, intercept = FALSE, standardize = FALSE), Z, "(b) Z ortonormal:")
show(glmnet::glmnet(Phi, Y, offset = Q0m, intercept = FALSE, standardize = FALSE), Phi, "(c) Phi directa:  ")
