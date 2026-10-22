# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later

suppressMessages(library(haldensify))
pkgload::load_all("thalmle_real", quiet = TRUE)
set.seed(10)
d <- simulate_data(500, 1, "normal")
n <- length(d$A)
t0 <- Sys.time()
hd <- haldensify(A = d$A, W = d$W, n_bins = 10L, grid_type = "equal_range",
                 lambda_seq = exp(seq(-1, -10, length = 100)), max_degree = 2L)
cat(sprintf("ajuste en %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat("breaks:\n"); print(hd$breaks)
cat("bin_sizes:\n"); print(hd$bin_sizes)
cat("range_a:\n"); print(hd$range_a)
cat("rango observado de A:", range(d$A), "\n")
br <- hd$breaks
if (!is.numeric(br) || length(br) < 3) stop("estructura de breaks inesperada")
bin_of <- function(a) findInterval(a, br, rightmost.closed = TRUE, all.inside = TRUE)
K <- length(br) - 1L
mids <- (br[-1] + br[-length(br)]) / 2
cat("numero de bins K =", K, "\n")
t0 <- Sys.time()
G <- matrix(suppressMessages(as.numeric(predict(hd, new_A = rep(mids, each = n),
                                                new_W = rep(d$W, times = K)))), n, K)
cat(sprintf("tabla %d x %d construida en %.0f s\n", n, K,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
g_lookup <- function(a, w) G[cbind(match(w, d$W), bin_of(a))]
set.seed(2)
ii <- sample(n, 300, replace = TRUE)
aa <- runif(300, 1, 9)
direct <- suppressMessages(as.numeric(predict(hd, new_A = aa, new_W = d$W[ii])))
look <- g_lookup(aa, d$W[ii])
cat("validacion en 300 pares aleatorios (dosis en [1,9], W observados):\n")
cat("  max |directo - tabla| =", max(abs(direct - look)), "\n")
cat("  max directo =", max(direct), "  max tabla =", max(look), "\n")
cat("  pares con diferencia > 1e-6:", sum(abs(direct - look) > 1e-6), "de 300\n")
cat("fin\n")
