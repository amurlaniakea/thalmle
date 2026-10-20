# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' T-HAL-MLE core (Algorithm 1, steps 3-8), continuous outcome, scalar W.
#'
#' The initial fits are passed in as functions so the core can be tested
#' without hal9001/haldensify:
#'   Q_fun(a, w): initial estimate of E[Y | A = a, W = w] (vectorised pairwise)
#'   g_fun(a, w): estimate of the conditional density g(a | w)
#'   omega_fun(a): estimate of the marginal density of A
#'   fit_Q(A, W, Y): optional; refits the outcome regression on a subsample and
#'     returns a function (a, w) -> Q. Needed for honest CV and for cross-fitting.
#'
#' Options:
#'   * fit_Q given, crossfit = FALSE, use_cv = TRUE: honest V-fold CV for the
#'     targeting penalty (Q refitted on each training split; paper Sec. 5.1).
#'   * crossfit = TRUE (requires fit_Q): the offset of the targeting regression is
#'     the OUT-OF-FOLD prediction of Q (K = cf_folds folds), so the residuals Y - Q
#'     are true residuals. cf_plugin = "avg": the plug-in marginal uses the average
#'     of the fold-specific Q fits; "full": it uses the full-sample Q_fun. The
#'     penalty is then chosen with cv.glmnet on those out-of-fold offsets.
#'   * Neither: cv.glmnet reuses the in-sample Q (not honest).
#'
#' Deviations from the paper to review:
#'   * the OLS projection (step 7) uses at most `n_proj` dose points;
#'   * the clever covariate is omega(a)/g(a|W) * phi(a), WITHOUT Sigma^+
#'     (the paper uses Sigma^{-1} phi); the selected working model is phi_j for
#'     nonzero eps_j;
#'   * g and omega are always the ones fitted on the full sample.
#' @param w_cap Cap for counterfactual weights omega/g. NULL uses sqrt(n)*log(n/5). Inf disables.
#' @export
thal_core <- function(A, W, Y, Q_fun, g_fun, omega_fun, grid,
                      k = 1L, n_knots = 50L, c1 = 6, c2 = 9,
                      use_cv = TRUE, n_proj = 200L, nlambda = 100L, w_cap = NULL,
                      fit_Q = NULL, cv_folds = 5L,
                      crossfit = FALSE, cf_folds = 2L, cf_plugin = c("avg", "full"),
                      sigma_tol = NULL) {
  cf_plugin <- match.arg(cf_plugin)
  n <- length(A)
  if (any(grid < min(A, na.rm = TRUE) | grid > max(A, na.rm = TRUE)))
    warning("grid contains values outside range(A); extrapolation will occur without warning.")
  cap <- if (is.null(w_cap)) sqrt(n) * log(n / 5) else if (identical(w_cap, "obs_max")) max(omega_fun(A) / g_fun(A, W)) else as.numeric(w_cap)
  knots <- make_knots(A, n_knots)
  Phi_A <- spline_basis(A, knots, k)
  wA <- pmin(omega_fun(A) / g_fun(A, W), cap)
  Sp <- if (is.null(sigma_tol)) NULL else pinv_sym(crossprod(Phi_A) / n, tol = sigma_tol)
  H <- wA * (if (is.null(Sp)) Phi_A else Phi_A %*% Sp)
  offset <- Q_fun(A, W)
  Q_use <- Q_fun
  if (crossfit) {
    if (is.null(fit_Q)) stop("crossfit = TRUE requires fit_Q")
    fid <- sample(rep_len(seq_len(cf_folds), n))
    Qcf <- lapply(seq_len(cf_folds), function(v) {
      tr <- which(fid != v)
      fit_Q(A[tr], W[tr], Y[tr])
    })
    offset <- numeric(n)
    for (v in seq_len(cf_folds)) {
      va <- which(fid == v)
      offset[va] <- Qcf[[v]](A[va], W[va])
    }
    if (cf_plugin == "avg") {
      Q_use <- function(a, w) Reduce(`+`, lapply(Qcf, function(f) f(a, w))) / cf_folds
    }
  }

  fit <- glmnet::glmnet(H, Y, family = "gaussian", offset = offset,
                        intercept = FALSE, standardize = FALSE,
                        nlambda = nlambda)
  sizes <- fit$df
  cv_index <- NULL
  if (use_cv) {
    if (is.null(fit_Q) || crossfit) {
      cvfit <- glmnet::cv.glmnet(H, Y, family = "gaussian", offset = offset,
                                 intercept = FALSE, standardize = FALSE,
                                 lambda = fit$lambda)
      cv_index <- which.min(cvfit$cvm)
    } else {
      lam <- fit$lambda
      folds <- sample(rep_len(seq_len(cv_folds), n))
      risk <- matrix(NA_real_, cv_folds, length(lam))
      for (v in seq_len(cv_folds)) {
        tr <- which(folds != v); va <- which(folds == v)
        Qv <- fit_Q(A[tr], W[tr], Y[tr])
        fv <- glmnet::glmnet(H[tr, , drop = FALSE], Y[tr], family = "gaussian",
                             offset = Qv(A[tr], W[tr]), intercept = FALSE,
                             standardize = FALSE, lambda = lam)
        pv <- predict(fv, newx = H[va, , drop = FALSE], s = lam,
                      newoffset = Qv(A[va], W[va]))
        risk[v, ] <- colMeans((Y[va] - pv)^2)
      }
      cv_index <- which.min(colMeans(risk))
    }
  }
  cv_size <- if (is.null(cv_index)) NA_integer_ else as.integer(sizes[cv_index])
  idx <- select_lambda_by_size(sizes, n, k, c1, c2, cv_index)
  eps <- as.numeric(fit$beta[, idx])

  # targeted plug-in marginal: mean_i Q*(a, W_i)
  marginal <- function(points) {
    m <- length(points)
    a_rep <- rep(points, each = n)
    w_rep <- rep(W, times = m)
    Qb <- matrix(Q_use(a_rep, w_rep), n, m)
    gb <- matrix(g_fun(a_rep, w_rep), n, m)
    Pm <- spline_basis(points, knots, k)
    if (!is.null(Sp)) Pm <- Pm %*% Sp
    corr <- as.numeric(Pm %*% eps)
    Wmat <- pmin(sweep(1 / gb, 2, omega_fun(points), "*"), cap)
    Qstar <- Qb + sweep(Wmat, 2, corr, "*")
    colMeans(Qstar)
  }

  # baseline: plug-in of the full-sample Q (no targeting, no cross-fitting)
  psi_init <- colMeans(matrix(Q_fun(rep(grid, each = n), rep(W, times = length(grid))),
                              n, length(grid)))
  psi_plugin <- marginal(grid)
  sel <- which(eps != 0)
  if (length(sel) == 0L) {
    warning("targeting step selected no basis functions; returning plug-in")
    return(list(psi_hat = psi_plugin, psi_plugin = psi_plugin, psi_init = psi_init,
                grid = grid, eps = eps, selected = sel, size = 0L,
                cv_size = cv_size, lambda = fit$lambda[idx]))
  }

  # step 7: OLS projection of the plug-in onto the selected basis
  proj_pts <- if (n > n_proj) {
    unique(quantile(A, probs = seq(0, 1, length.out = n_proj), names = FALSE))
  } else sort(A)
  Phi_p <- spline_basis(proj_pts, knots, k)[, sel, drop = FALSE]
  alpha <- pinv_sym(crossprod(Phi_p)) %*% crossprod(Phi_p, marginal(proj_pts))
  psi_hat <- as.numeric(spline_basis(grid, knots, k)[, sel, drop = FALSE] %*% alpha)

  list(psi_hat = psi_hat, psi_plugin = psi_plugin, psi_init = psi_init, grid = grid,
       eps = eps, selected = sel, size = length(sel), cv_size = cv_size,
       lambda = fit$lambda[idx])
}

#' Plumbing check on the paper's simulated DGPs with oracle nuisances.
#' q_init = "misspecified" ignores W in the initial Q, so the plug-in is
#' biased and the targeting step has something to correct.
#' @param w_cap Cap for counterfactual weights omega/g. NULL uses sqrt(n)*log(n/5). Inf disables.
#' @export
thal_oracle <- function(n = 500, dgp = 1, treatment = "uniform",
                        q_init = c("oracle", "misspecified"), k = 1L,
                        grid = seq(1, 9, length.out = 25), ...) {
  q_init <- match.arg(q_init)
  d <- simulate_data(n, dgp, treatment)
  Q_fun <- if (q_init == "oracle") {
    function(a, w) Q0_dgp(a, w, dgp)
  } else {
    function(a, w) 3 * sin(1.5 * a) + 0.5 * a
  }
  res <- thal_core(d$A, d$W, d$Y, Q_fun,
                   g_fun = function(a, w) g0_dgp(a, w, treatment),
                   omega_fun = function(a) omega_oracle(a, treatment),
                   grid = grid, k = k, ...)
  res$truth <- true_drc(grid, dgp)
  res
}
