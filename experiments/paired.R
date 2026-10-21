res <- readRDS("tail_normal.rds")
nm <- names(res[[1]]$thal)
rm <- function(e) sqrt(mean(e^2))
hal <- sapply(res, function(o) rm(o$init))
cat(sprintf("replicas: %d | RMSE medio por replica del HAL plug-in: %.3f\n", length(res), mean(hal)))
for (k in nm) {
  v <- sapply(res, function(o) rm(o$thal[[k]]))
  dl <- v - hal
  cat(sprintf("%-7s RMSE medio %.3f | dif. vs HAL: %+.3f (EE %.3f) | T-HAL mejor en %d de %d replicas\n",
              k, mean(v), mean(dl), sd(dl) / sqrt(length(dl)), sum(dl < 0), length(dl)))
}
