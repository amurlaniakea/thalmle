# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

suppressMessages(library(haldensify))
pkgload::load_all("thalmle_real", quiet = TRUE)
set.seed(10)
d <- simulate_data(500, 1, "normal"); n <- length(d$A)
hd <- haldensify(A = d$A, W = d$W, n_bins = 10L, grid_type = "equal_range",
                 lambda_seq = exp(seq(-1, -10, length = 100)), max_degree = 2L)
edges <- c(hd$breaks, hd$range_a[2]); K <- length(edges) - 1L
bin_of <- function(a) findInterval(a, edges, rightmost.closed = TRUE, all.inside = TRUE)
mids <- (edges[-1] + edges[-length(edges)]) / 2
pr <- function(a, w) suppressMessages(as.numeric(predict(hd, new_A = a, new_W = w)))
sub <- 1:60
G <- matrix(pr(rep(mids, each = length(sub)), rep(d$W[sub], times = K)), length(sub), K)
set.seed(2)
ii <- sample(sub, 300, replace = TRUE); aa <- runif(300, 1, 9)
direct <- pr(aa, d$W[ii])
look <- G[cbind(match(d$W[ii], d$W[sub]), bin_of(aa))]
bad <- which(abs(direct - look) > 1e-6)
cat("desajustes tabla vs lote:", length(bad), "de 300\n")
rv <- rev(pr(rev(aa), rev(d$W[ii])))
cat("(1) max |lote - lote con la entrada en orden inverso| =", max(abs(direct - rv)), "\n")
sb <- pr(aa[1:10], d$W[ii[1:10]])
cat("(2) max |10 primeros solos - los mismos dentro del lote| =", max(abs(sb - direct[1:10])), "\n")
jj <- head(bad, 6)
alone <- sapply(jj, function(j) pr(rep(aa[j], 2), rep(d$W[ii[j]], 2))[1])
cat("(3) pares desajustados: valor en el lote, evaluado solo (fila duplicada) y en la tabla\n")
print(round(data.frame(j = jj, bin = bin_of(aa[jj]), a = aa[jj],
                       lote = direct[jj], solo = alone, tabla = look[jj]), 4))
cat("fin\n")
