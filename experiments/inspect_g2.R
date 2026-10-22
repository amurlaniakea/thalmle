# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

suppressMessages(library(haldensify))
pkgload::load_all("thalmle_real", quiet = TRUE)
set.seed(10)
d <- simulate_data(500, 1, "normal"); n <- length(d$A)
hd <- haldensify(A = d$A, W = d$W, n_bins = 10L, grid_type = "equal_range",
                 lambda_seq = exp(seq(-1, -10, length = 100)), max_degree = 2L)
# breaks = bordes inferiores de los 10 bins; el borde superior es range_a[2]
edges <- c(hd$breaks, hd$range_a[2]); K <- length(edges) - 1L
cat("K =", K, " bordes:", round(edges, 3), "\n")
bin_of <- function(a) findInterval(a, edges, rightmost.closed = TRUE, all.inside = TRUE)
mids <- (edges[-1] + edges[-length(edges)]) / 2
sub <- 1:60
G <- matrix(suppressMessages(as.numeric(predict(hd,
        new_A = rep(mids, each = length(sub)), new_W = rep(d$W[sub], times = K)))),
        length(sub), K)
set.seed(2)
ii <- sample(sub, 300, replace = TRUE); aa <- runif(300, 1, 9)
direct <- suppressMessages(as.numeric(predict(hd, new_A = aa, new_W = d$W[ii])))
look <- G[cbind(match(d$W[ii], d$W[sub]), bin_of(aa))]
dd <- abs(direct - look)
cat("(A) 10 bins, bordes corregidos: max dif =", max(dd),
    " pares con dif > 1e-6:", sum(dd > 1e-6), "de 300\n")
bk <- bin_of(aa)
cat("(B) fraccion de desajustes por bin:\n"); print(round(tapply(dd > 1e-6, bk, mean), 2))
rel <- (aa - edges[bk]) / (edges[bk + 1] - edges[bk])
cat("(C) posicion relativa (0 = borde izq, 1 = borde der) de los desajustes:\n")
print(summary(rel[dd > 1e-6]))
a_f <- seq(0.02, 9.93, length.out = 600); w0 <- d$W[1]
g_f <- suppressMessages(as.numeric(predict(hd, new_A = a_f, new_W = rep(w0, 600))))
chg <- which(abs(diff(g_f)) > 1e-10)
cat("(D) saltos empiricos de g (w fijo) en a =", round((a_f[chg] + a_f[chg + 1]) / 2, 3), "\n")
cat("    numero de saltos:", length(chg), "\n")
single <- sapply(1:6, function(k)
  suppressMessages(as.numeric(predict(hd, new_A = aa[k], new_W = d$W[ii[k]]))))
cat("(E) max |vector - par individual| =", max(abs(single - direct[1:6])), "\n")
cat("fin\n")
