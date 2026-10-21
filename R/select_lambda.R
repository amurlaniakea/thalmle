# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

#' Rate-optimal window for the working-model size (Section 5.1 reading).
#'
#' `sizes` = number of nonzero coefficients along the lambda path, ordered
#' from the largest penalty to the smallest. Window:
#'   [ceiling(c1 * n^(1/(2k*+1))), floor(c2 * n^(1/(2k*+1)))], k* = k + 1.
#' If `cv_index` is NULL the first index reaching the lower bound is used
#' (the "fixed" selector); otherwise the CV choice is clamped to the window.
#' NOTE: the paper's Algorithm 1 (step 5) has the inequalities reversed
#' with respect to Section 5.1; this follows Section 5.1.
#' @param sizes integer vector: working-model size along the lambda path.
#' @param n sample size.
#' @param k smoothness order of the target function.
#' @param c1,c2 constants of the lower and upper bound of the window.
#' @param cv_index index of the cross-validated lambda, or NULL.
#' @return index into `sizes` of the selected lambda.
#' @export
select_lambda_by_size <- function(sizes, n, k = 1L, c1 = 6, c2 = 9,
                                  cv_index = NULL) {
  kstar <- k + 1L
  rate <- n^(1 / (2 * kstar + 1))
  lo <- ceiling(c1 * rate)
  hi <- floor(c2 * rate)
  first_ge_lo <- function() {
    i <- which(sizes >= lo)
    if (length(i)) min(i) else length(sizes)
  }
  last_le_hi <- function() {
    i <- which(sizes <= hi)
    if (length(i)) max(i) else 1L
  }
  if (is.null(cv_index)) return(first_ge_lo())
  if (sizes[cv_index] < lo) return(first_ge_lo())
  if (sizes[cv_index] > hi) return(last_le_hi())
  cv_index
}
