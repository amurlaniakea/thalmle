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
