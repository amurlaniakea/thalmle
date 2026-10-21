pkgload::load_all("thalmle_hermes2", quiet = TRUE)
for (c1 in c(6, 3, 0)) {
  set.seed(10)
  out <- lapply(1:20, function(b) {
    r <- suppressWarnings(thal_oracle(n = 500, dgp = 1, treatment = "normal",
                                      q_init = "misspecified", c1 = c1))
    list(e = r$psi_hat - r$truth, size = r$size, maxeps = max(abs(r$eps)))
  })
  E <- sapply(out, `[[`, "e")
  cat(sprintf("c1=%g  RMSE=%.2f  mediana size=%g  mediana max|eps|=%.1f\n",
      c1, sqrt(mean(E^2)), median(sapply(out, `[[`, "size")),
      median(sapply(out, `[[`, "maxeps"))))
}
