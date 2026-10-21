suppressMessages(library(haldensify))
pkgload::load_all("thalmle_real", quiet = TRUE)
set.seed(10)
d <- simulate_data(500, 1, "normal"); n <- length(d$A)
hd <- haldensify(A = d$A, W = d$W, n_bins = 10L, grid_type = "equal_range",
                 lambda_seq = exp(seq(-1, -10, length = 100)), max_degree = 2L)
edges <- c(hd$breaks, hd$range_a[2]); K <- length(edges) - 1L
widths <- diff(edges)
bin_of <- function(a) findInterval(a, edges, rightmost.closed = TRUE, all.inside = TRUE)
mids <- (edges[-1] + edges[-length(edges)]) / 2
pr <- function(a, w) suppressMessages(as.numeric(predict(hd, new_A = a, new_W = w, trim = TRUE)))
sub <- 1:60
G <- matrix(pr(rep(mids, each = length(sub)), rep(d$W[sub], times = K)), length(sub), K)
set.seed(2)
ii <- sample(sub, 300, replace = TRUE); aa <- runif(300, 1, 9)
direct <- pr(aa, d$W[ii])
look <- G[cbind(match(d$W[ii], d$W[sub]), bin_of(aa))]
cat("(1) tabla vs lote (trim = TRUE): max dif =", max(abs(direct - look)),
    " desajustes (>1e-6):", sum(abs(direct - look) > 1e-6), "de 300\n")
sb <- pr(aa[1:10], d$W[ii[1:10]])
cat("(2) 10 primeros solos vs dentro del lote: max dif =", max(abs(sb - direct[1:10])), "\n")
integ <- as.numeric(G %*% widths)
cat("(3) integral de g sobre el soporte, por W (debe ser ~1):\n"); print(summary(integ))
Tm <- sapply(mids, function(m) g0_dgp(rep(m, length(sub)), d$W[sub], "normal"))
cat(sprintf("(4) tabla vs densidad verdadera en los puntos medios: cor = %.3f, RMSE = %.4f\n",
            cor(as.vector(G), as.vector(Tm)), sqrt(mean((G - Tm)^2))))
cat(sprintf("    g minimo estimado = %.5f   g minimo verdadero = %.5f\n", min(G), min(Tm)))
cat("fin\n")
