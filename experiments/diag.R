pkgload::load_all("thalmle", quiet = TRUE)
set.seed(1)
n <- 500
d <- simulate_data(n, 1, "uniform")
A <- d$A; W <- d$W; Y <- d$Y
knots <- make_knots(A, 50L)
Phi <- spline_basis(A, knots, 1L)
S <- crossprod(Phi) / n
ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
cat("eigenvalues max, min:", max(ev), min(ev), "\n")
cat("eigenvalues kept by pinv_sym:", sum(ev > 1e-10 * max(ev)), "of", length(ev), "\n")
Sp <- pinv_sym(S)
cat("max |Sigma+| =", max(abs(Sp)), "\n")
H <- Phi %*% Sp
cat("H column norms, min and max:", range(sqrt(colSums(H^2))), "\n")
Q0m <- 3 * sin(1.5 * A) + 0.5 * A
R <- Y - Q0m
ols <- lm(R ~ A)
cat("residual on A: intercept, slope =", coef(ols), "(expected about -10 and 2.5)\n")
fit <- glmnet::glmnet(H, Y, offset = Q0m, intercept = FALSE, standardize = FALSE)
cat("path length:", length(fit$lambda), " max df:", max(fit$df), "\n")
idxs <- unique(pmin(c(5, 20, 40, length(fit$lambda)), length(fit$lambda)))
cat("resid RMSE, no targeting:", sqrt(mean(R^2)), "\n")
for (i in idxs) {
  e <- as.numeric(fit$beta[, i])
  cat(sprintf("idx %3d  lambda %.3g  df %3d  resid RMSE %.3f  max|eps| %.3g\n",
              i, fit$lambda[i], fit$df[i], sqrt(mean((R - H %*% e)^2)), max(abs(e))))
}
