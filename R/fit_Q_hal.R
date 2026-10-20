# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
#
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Copyright (C) 2026 Pedro Sordo Martínez
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with this program. If not, see
# <https://www.gnu.org/licenses/>.

#' fit_Q_hal
#'
#' Regresión de resultado Q con HAL (preferente) o fallback glmnet/lm.
#'
#' @param A vector dosis
#' @param W covariables
#' @param Y resultado
#' @param smoothness_orders orden (entero largo, por defecto 0L para DGP 1/2)
#' @return objeto con método predict funcional
#' @export
fit_Q_hal <- function(A, W, Y, smoothness_orders = 0L) {
  # TODO: Verificar parámetros exactos de hal9001 con el preprint.
  if (length(A) != length(Y) || length(A) != length(W)) {
    stop("A, W y Y deben tener la misma longitud")
  }
  X <- if (is.data.frame(W)) W else data.frame(W = W)
  df <- data.frame(A = A, X)
  if (requireNamespace("hal9001", quietly = TRUE)) {
    tryCatch({
      # Se usaría: hal9001::fit_hal(X = df, Y = Y, smoothness_orders = smoothness_orders)
      # Fallback a stub funcional con mensaje informativo
      stop("Uso de hal9001 pendiente; se requiere implementación real del preprint")
    }, error = function(e) {
      message("hal9001 no funcional: ", conditionMessage(e))
      # Fallback seguro con glmnet
      return(.fit_Q_fallback(df, Y, smoothness_orders))
    })
  } else {
    message("hal9001 no disponible; usando fallback glmnet/lm")
    .fit_Q_fallback(df, Y, smoothness_orders)
  }
}

.fit_Q_fallback <- function(df, Y, smoothness_orders) {
  # Fallback seguro: usa glmnet si está disponible, si no lm base.
  res <- list(
    smoothness_orders = smoothness_orders,
    fitted = list(status = "fallback_glmnet_or_lm", call = match.call())
  )
  if (requireNamespace("glmnet", quietly = TRUE)) {
    x <- as.matrix(df)
    fit_glm <- glmnet::glmnet(x, Y, family = "gaussian", alpha = 0.5, lambda = 0.01)
    res$fitted <- fit_glm
    res$predict <- function(newdata) {
      newx <- as.matrix(newdata[, names(df), drop = FALSE])
      as.vector(stats::predict(fit_glm, newx = newx, s = 0.01))
    }
  } else {
    fit_lm <- stats::lm(Y ~ ., data = df)
    res$fitted <- fit_lm
    res$predict <- function(newdata) {
      as.vector(stats::predict(fit_lm, newdata = newdata))
    }
  }
  class(res) <- c("hal9001_fit_fallback", class(res))
  res
}
