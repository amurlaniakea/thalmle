pkgload::load_all("thalmle_hermes", quiet = TRUE)
set.seed(10)
n <- 500
d <- simulate_data(n, 1, "normal")
cap <- sqrt(n) * log(n / 5)
w <- omega_oracle(d$A, "normal") / g0_dgp(d$A, d$W, "normal")
cat("cap sqrt(n)*log(n/5) =", cap, "\nobserved omega/g:\n"); print(summary(w))
cat("share of observed weights above cap:", mean(w > cap), "\n")
for (a in c(1, 5, 9)) {
  wa <- omega_oracle(a, "normal") / g0_dgp(rep(a, n), d$W, "normal")
  cat("grid a =", a, " max =", max(wa), " p99 =", quantile(wa, .99),
      " share>cap:", mean(wa > cap), "\n")
}
