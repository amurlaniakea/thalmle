pkgload::load_all("thalmle_hermes2", quiet = TRUE)
set.seed(10)
B <- 50
grid <- seq(1, 9, length.out = 25)
out <- vector("list", B)
for (b in seq_len(B)) {
  r <- suppressWarnings(thal_oracle(n = 500, dgp = 1, treatment = "normal",
                                    q_init = "misspecified"))
  out[[b]] <- list(err = r$psi_hat - r$truth, err0 = r$psi_init - r$truth,
                   size = r$size, maxeps = max(abs(r$eps)))
}
E  <- sapply(out, function(o) o$err)
E0 <- sapply(out, function(o) o$err0)
cat("RMSE por punto de la malla (T-HAL):\n")
print(round(setNames(sqrt(rowMeans(E^2)), round(grid, 2)), 2))
rr <- sqrt(colMeans(E^2))
cat("RMSE por réplica, cuantiles 0/25/50/75/90/100:\n")
print(round(quantile(rr, c(0, .25, .5, .75, .9, 1)), 2))
cen <- grid >= 3 & grid <= 7
cat(sprintf("T-HAL   global %.2f  centro[3,7] %.2f  colas %.2f\n",
            sqrt(mean(E^2)), sqrt(mean(E[cen, ]^2)), sqrt(mean(E[!cen, ]^2))))
cat(sprintf("inicial global %.2f  centro[3,7] %.2f  colas %.2f\n",
            sqrt(mean(E0^2)), sqrt(mean(E0[cen, ]^2)), sqrt(mean(E0[!cen, ]^2))))
sz <- sapply(out, `[[`, "size"); me <- sapply(out, `[[`, "maxeps")
o <- order(rr, decreasing = TRUE)[1:5]
cat("5 peores réplicas:\n")
print(data.frame(rep = o, rmse = round(rr[o], 1), size = sz[o], maxeps = signif(me[o], 3)))
cat("mediana de size y maxeps en todas:", median(sz), signif(median(me), 3), "\n")
