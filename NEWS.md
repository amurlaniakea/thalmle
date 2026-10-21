# thalmle (development)

* DGP 2 now follows the preprint: `g(w) = 4 (1[w > -1] - 1[w > 0] + 1[w > 1])`.
  Earlier versions used `4 1[w > -1] - 1[w > 0] + 1[w > 1]`, so the earlier DGP 2
  figures came from a different DGP. The DGP 2 experiment was re-run on 2026-10-21
  (20 replicates, normal treatment, seeds 7000+b). With the outcome HAL of order 0
  the capped T-HAL variants are slightly better than the HAL plug-in; without `THAL_QSO=0`
  (presumably order 1; the environment was not recorded) they are worse (see the README). The
  earlier files in `results/` carry the suffix `_anterior`, the order-1 run `_qso1`.
* `experiments/check_capobs2.R` prints `THAL_QSO` in its `cfg:` line.
* `thal_core()` warns when `grid` has doses outside `range(A)`.
* Tests for the options of `thal_core()` that were not covered, and smoke tests
  for the wrappers around hal9001 and haldensify (skipped when those packages
  are absent).
* `haldensify (>= 0.2.8)` is now stated in `Suggests`.

# thalmle 0.0.2

* `thal_core()` gains experimental options: `w_cap = "obs_max"` (cap the
  counterfactual weights at the largest weight seen at the observed data; this
  is NOT in the preprint), `sigma_tol` (clever covariate with a truncated
  Moore-Penrose inverse of Sigma), `fit_Q` (honest V-fold CV for the penalty) and
  `crossfit` (out-of-fold offsets for the targeting regression).
* New wrappers `fit_Q_hal()`, `fit_g_haldensify()` and `fit_omega_kde()`.
  `haldensify` predictions are always computed with `trim = FALSE` and looked up
  from a table (the density is a step function of the dose).
* Experiments and logs used to reach the conclusions in the README are in
  `experiments/` and `results/`.

# thalmle 0.0.1

* Skeleton: data-generating processes, spline basis, LASSO targeting step.
