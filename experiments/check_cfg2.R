args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; dgp <- as.integer(args[2]); trt <- args[3]
qi <- args[4]; B <- as.integer(args[5]); sd0 <- as.integer(args[6])
cap_arg <- if (length(args) >= 7) args[7] else "default"
pkgload::load_all(pkgdir, quiet = TRUE)
set.seed(sd0)
nwarn <- 0L
res <- withCallingHandlers(replicate(B, {
  r <- if (cap_arg == "default") {
    thal_oracle(n = 500, dgp = dgp, treatment = trt, q_init = qi)
  } else {
    thal_oracle(n = 500, dgp = dgp, treatment = trt, q_init = qi,
                w_cap = as.numeric(cap_arg))
  }
  c(r$psi_init - r$truth, r$psi_plugin - r$truth, r$psi_hat - r$truth)
}), warning = function(w) { nwarn <<- nwarn + 1L; invokeRestart("muffleWarning") })
m <- nrow(res) / 3
init <- res[1:m, ]; plug <- res[(m + 1):(2 * m), ]; thal <- res[(2 * m + 1):(3 * m), ]
rep1 <- function(tag, x) cat(sprintf("%-24s mean|bias| = %.4f  RMSE = %.4f\n", tag,
  mean(abs(rowMeans(x))), sqrt(mean(x^2))))
cat(sprintf("cfg: dgp=%d trt=%s q_init=%s B=%d seed=%d cap=%s pkg=%s\n",
            dgp, trt, qi, B, sd0, cap_arg, pkgdir))
rep1("initial (sin targeting)", init)
rep1("plug-in dirigido", plug)
rep1("T-HAL (psi_hat)", thal)
cat("warnings glmnet/otros:", nwarn, "\n")
