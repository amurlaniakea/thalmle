# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

test_that("DGP 2 is the preprint's g(w) = 4 (1[w>-1] - 1[w>0] + 1[w>1])", {
  w <- c(-2, -1.5, -0.5, 0.5, 1.5, 2)
  expect_equal(thalmle:::g_fun_dgp(w, 2), 4 * ((w > -1) - (w > 0) + (w > 1)))
  expect_equal(thalmle:::Eg_dgp(2), 2)
})

test_that("DGP 1 and 3 are unchanged", {
  w <- c(-1, 0.5, 2)
  expect_equal(thalmle:::g_fun_dgp(w, 1), 5 * (w > 0))
  expect_equal(thalmle:::g_fun_dgp(w, 3), 3 * sin(7 * pi * w))
})
