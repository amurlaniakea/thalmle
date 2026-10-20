# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
#
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Copyright (C) 2026 Pedro Sordo Martinez
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

#' true_drc
#'
#' Funcion verdadera de respuesta a dosis (DGP). Formula placeholder coherente:
#' DGP 1: sigmoide con discontinuidad en a = 5.
#' DGP 2: lineal con pendiente creciente.
#' TODO: Verificar con ecuaciones exactas del preprint arXiv:2607.03824.
#'
#' @param a punto de dosis (escalar)
#' @param d DGP (1 o 2)
#' @return valor escalar
#' @export
true_drc <- function(a, d) {
  # TODO: Verificar con ecuaciones exactas del preprint
  if (d == 1) {
    # Sigmoide centrada con salto en 5
    return(2 * (a - 3) / (1 + exp(-(a - 5))) + 0.5)
  } else if (d == 2) {
    return(0.3 * a + 0.02 * a^2)
  } else {
    return(rep(NA, length(a)))
  }
}

#' simulate_data
#'
#' Datos simulados para DGP 1-3.
#'
#' @param n tamano muestra
#' @param dgp 1, 2 o 3
#' @param treatment "uniform" o "normal"
#' @return data.frame con A, W, Y
#' @export
simulate_data <- function(n = 500, dgp = 1, treatment = "uniform") {
  # TODO: Verificar con ecuaciones exactas del preprint
  W <- rnorm(n, mean = 0, sd = 1)
  if (treatment == "uniform") {
    A <- runif(n, min = 1, max = 9)
  } else {
    A <- pmax(0, pmin(10, rnorm(n, mean = 5, sd = 2)))
  }
  Y <- true_drc(A, dgp) + 0.4 * W + rnorm(n, sd = 0.3)
  data.frame(A = A, W = W, Y = Y)
}
