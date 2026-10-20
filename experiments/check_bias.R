pkgload::load_all("thalmle", quiet = TRUE)
set.seed(10)
B <- 50
res <- replicate(B, {
  r <- thal_oracle(n = 500, dgp = 1, treatment = "uniform",
                   q_init = "misspecified")
  c(r$psi_plugin - r$truth, r$psi_hat - r$truth)
})
plug <- res[1:25, ]
thal <- res[26:50, ]
cat("plug-in : mean|bias| =", mean(abs(rowMeans(plug))),
    " RMSE =", sqrt(mean(plug^2)), "\n")
cat("T-HAL   : mean|bias| =", mean(abs(rowMeans(thal))),
    " RMSE =", sqrt(mean(thal^2)), "\n")
