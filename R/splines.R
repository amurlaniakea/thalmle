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

#' make_knots
#'
#' Genera n puntos de anclaje (knots) basados en cuantiles empíricos.
#'
#' @param x datos observados
#' @param n número de nudos
#' @return vector de nudos
#' @export
make_knots <- function(x, n = 10) {
  ux <- unique(x)
  if (length(ux) < n) n <- length(ux)
  probs <- seq(0, 1, length.out = n + 2)[-c(1, n + 2)]
  stats::quantile(x, probs = probs, names = FALSE)
}

#' spline_basis
#'
#' Base de splines lineales (grado 1) con intercepto constante en primera columna.
#'
#' @param x valores
#' @param knots puntos de anclaje
#' @param k grado (por defecto 1L para lineal)
#' @return matriz base con primera columna constante
#' @export
spline_basis <- function(x, knots, k = 1L) {
  # Base lineal manual: primera columna constante = 1 (intercepto),
  # luego funciones de base por trozos para cada knot.
  B <- cbind(1, x)
  for (kn in knots) {
    B <- cbind(B, pmax(0, x - kn))
  }
  B
}

#' pinv_sym
#'
#' Pseudo-inversa de Moore-Penrose para matriz simétrica mediante SVD.
#'
#' @param M matriz cuadrada
#' @param tol tolerancia para valores singulares (1e-8)
#' @return pseudo-inversa
#' @export
pinv_sym <- function(M, tol = 1e-8) {
  s <- svd(M)
  d <- s$d
  d_inv <- ifelse(d > tol, 1 / d, 0)
  s$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(s$u)
}

#' select_lambda_by_size
#'
#' Selección de lambda por recorte de tamaño del modelo.
#'
#' @param sizes vector de tamaños candidatos
#' @param n tamaño de muestra
#' @param cv_index índice del candidato de CV (opcional)
#' @return índice seleccionado
#' @export
select_lambda_by_size <- function(sizes, n, cv_index = NULL) {
  # TODO: Verificar con la Sección 5.1 del preprint sobre límites exactos.
  L <- ceiling(6 * n^(1/5))
  U <- floor(9 * n^(1/5))
  if (!is.null(cv_index)) {
    val <- sizes[cv_index]
    if (val >= L && val <= U) return(cv_index)
    # Si fuera de rango, devolver el más cercano al rango válido
    if (val < L) {
      idx <- which.min(abs(sizes - L))
    } else {
      idx <- which.min(abs(sizes - U))
    }
    return(idx)
  }
  in_range <- which(sizes >= L & sizes <= U)
  if (length(in_range) > 0) return(in_range[1])
  # Si ninguno en rango, devolver el más cercano al rango
  closest <- which.min(abs(sizes - mean(c(L, U))))
  closest
}
