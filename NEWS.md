# thalmle (development)

* DGP 2 now follows the preprint: `g(w) = 4 (1[w > -1] - 1[w > 0] + 1[w > 1])`.
  Earlier versions used `4 1[w > -1] - 1[w > 0] + 1[w > 1]`, so every DGP 2 figure
  in `results/` and in the README was generated with a different DGP and must be
  re-run before it is quoted.
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
