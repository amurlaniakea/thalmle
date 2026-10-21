pkgload::load_all("thalmle_real", quiet = TRUE)
tm <- function(label, expr) {
  t0 <- Sys.time(); val <- force(expr)
  cat(sprintf("%-42s %6.1f s\n", label, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(val)
}
set.seed(10)
d <- simulate_data(500, 1, "uniform")
Qf <- tm("fit_Q_hal (hal9001, n=500)", fit_Q_hal(d$A, d$W, d$Y))
gf <- NULL
for (md in c(1L, 2L)) {
  g <- tryCatch(tm(sprintf("fit_g_haldensify max_degree=%d", md),
                   fit_g_haldensify(d$A, d$W, max_degree = md)),
                error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(g)) gf <- g
}
grid <- seq(1, 9, length.out = 25)
a1 <- rep(grid, each = 500); w1 <- rep(d$W, times = 25)
tm("Q predict, 12.500 pares (malla)", Qf(a1, w1))
if (!is.null(gf)) tm("g predict, 12.500 pares (malla)", gf(a1, w1))
a2 <- rep(seq(0, 10, length.out = 100), each = 500); w2 <- rep(d$W, times = 100)
tm("Q predict, 50.000 pares", Qf(a2, w2))
if (!is.null(gf)) tm("g predict, 50.000 pares", gf(a2, w2))
cat("fin\n")
