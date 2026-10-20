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

#' fit_omega_kde
#'
#' Densidad marginal por KDE usando stats::density.
#'
#' @param A vector de dosis
#' @return función (closure) que evalúa la densidad en nuevos puntos
#' @export
fit_omega_kde <- function(A) {
  if (length(A) < 2) stop("A debe tener al menos 2 observaciones")
  dens <- stats::density(A, from = min(A), to = max(A), bw = "nrd0")
  stats::approxfun(dens$x, dens$y, yleft = 0, yright = 0, rule = 2)
}

#' fit_g_haldensify
#'
#' Densidad condicional con haldensify.
#'
#' @param A dosis
#' @param W covariables
#' @return objeto densidad (stub documentado)
#' @export
fit_g_haldensify <- function(A, W) {
  # TODO: Verificar parámetros exactos de haldensify con el preprint.
  # Advertencia (README): predict() debe llamarse con trim = FALSE
  # para evitar floors dependientes del batch.
  if (!requireNamespace("haldensify", quietly = TRUE)) {
    stop("fit_g_haldensify requiere el paquete haldensify (no instalado)")
  }
  stop("Implementación pendiente: fit_g_haldensify — wrapper básico pendiente; usar trim = FALSE en predict() según preprint")
}
