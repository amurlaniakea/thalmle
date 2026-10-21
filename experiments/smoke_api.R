suppressMessages({library(hal9001); library(haldensify)})
cat("hal9001", as.character(packageVersion("hal9001")),
    " haldensify", as.character(packageVersion("haldensify")), "\n")
cat("fit_hal args:\n"); print(names(formals(hal9001::fit_hal)))
cat("predict.hal9001 args:\n"); print(names(formals(hal9001:::predict.hal9001)))
cat("haldensify args:\n"); print(names(formals(haldensify::haldensify)))
cat("predict.haldensify args:\n"); print(names(formals(haldensify:::predict.haldensify)))
set.seed(1)
n <- 300; W <- rnorm(n); A <- runif(n, 0, 10)
Y <- 3 * sin(1.5 * A) + 0.5 * A + rnorm(n)
X <- cbind(A = A, W = W)
t0 <- Sys.time()
fh <- tryCatch(hal9001::fit_hal(X = X, Y = Y, family = "gaussian", smoothness_orders = 1),
               error = function(e) { cat("fit_hal ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(fh)) {
  print(predict(fh, new_data = X[1:5, , drop = FALSE]))
  cat("fit_hal OK en", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
}
t0 <- Sys.time()
hd <- tryCatch(haldensify::haldensify(A = A, W = W, n_bins = 8, grid_type = "equal_range",
                 lambda_seq = exp(seq(-1, -8, length = 50)), max_degree = 1),
               error = function(e) { cat("haldensify ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(hd)) {
  print(predict(hd, new_A = A[1:5], new_W = W[1:5]))
  cat("haldensify OK en", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
}
