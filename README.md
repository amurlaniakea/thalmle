# thalmle 🧬

**Experimental, unofficial R implementation of the Targeted HAL-MLE (T-HAL-MLE)
for a continuous-dose dose-response curve.**

> **⚠️ ADVERTENCIA:** Ver `DISCLAIMER.md` antes de usar. Código de investigación; no validado contra el preprint; `w_cap` es heurística no incluida en el artículo.

> **Status: research code.** Point estimates only (no confidence intervals).
> It is *not* validated against the preprint: with uniform treatment it behaves
> as the preprint describes; with normal treatment it does so only in some
> settings and only with a modification that is not in the preprint. Read
> "Evidence" before using it.

It follows the targeting procedure of

> V. Rodriguez, K. Díaz-Ordaz, B. Lehmann, M. J. van der Laan.
> *Targeted Highly Adaptive Lasso Minimum Loss Estimation of Target Functions.*
> arXiv:2607.03824v1 (2026).

The code was written from the preprint. To our knowledge the authors have not
released an implementation (not exhaustively verified). Several details the
preprint does not specify were filled in by us (see "Differences").

## What it does

For observations `(W, A, Y)` (scalar covariate `W`, continuous dose `A`,
continuous outcome `Y`) it estimates `psi(a) = E_W[ E(Y | A = a, W) ]`:

1. outcome regression `Q_n` with a HAL (`hal9001`), `fit_Q_hal()`;
2. conditional density `g_n(a | w)` with `haldensify`, `fit_g_haldensify()`, and
   marginal density `omega` by kernel smoothing, `fit_omega_kde()`;
3. LASSO targeting step (`glmnet`) with `Q_n` as offset and clever covariate
   `omega(A) / g_n(A | W) * phi(A)`, where `phi` is a first-order spline basis;
4. penalty chosen by cross-validation, clipped so that the working model has
   between `ceil(6 n^(1/5))` and `floor(9 n^(1/5))` basis functions;
5. plug-in marginalisation over `W`, then least-squares projection on the
   selected basis functions: `thal_core()`.

## What it does not do

* No confidence intervals, influence-curve variance or bootstrap.
* Scalar `W`, one-dimensional dose, continuous outcome only.
* Only DGP 1 and 2 of the preprint and `n = 500` have been evaluated.

## Quick start

```r
# install.packages(c("glmnet", "hal9001", "haldensify"))
# R CMD INSTALL thalmle_0.0.3.tar.gz
library(thalmle)
set.seed(1)
d    <- simulate_data(500, dgp = 1, treatment = "normal")
Qf   <- fit_Q_hal(d$A, d$W, d$Y, smoothness_orders = 0L)
gf   <- fit_g_haldensify(d$A, d$W)
om   <- fit_omega_kde(d$A)
grid <- seq(1, 9, length.out = 25)
fit  <- thal_core(d$A, d$W, d$Y, Qf, gf, om, grid,
                  n_proj = length(d$A), w_cap = "obs_max")
fit$psi_hat     # T-HAL curve;  fit$psi_init = HAL plug-in
```

A full version is in `experiments/quickstart.R` (about 3 minutes per fit).

## Evidence (read before using)

Setting: DGP of Section 7.1 of the preprint, `n = 500`, outcome HAL of order
zero, `haldensify` with 10 bins and `max_degree = 2`, 50 knots + intercept,
first-order splines for the target function, projection on all observed doses.
Metric: for each replicate, the RMSE over 25 doses in [1, 9]; we report the
**paired difference T-HAL minus HAL plug-in** (negative = T-HAL better) with its
standard error. Samples have 20-30 replicates (the preprint uses 1,000).

> **DGP 2.** Until 2026-10-20 the code used `g(w) = 4 1[w > -1] - 1[w > 0] +
> 1[w > 1]`, which is not the preprint's `g(w) = 4 (1[w > -1] - 1[w > 0] +
> 1[w > 1])`. The DGP 2 rows below come from runs on 2026-10-21 with the
> preprint's DGP 2 (normal treatment, seeds 7000+b, 20 replicates). The one in the
> table used the outcome HAL of order 0 (`THAL_QSO=0`, as in the preprint and in the
> rest of the evidence; commit `d90e02b`, session in
> `results/sessionInfo_dgp2_rerun.txt`). A first run, launched without
> `THAL_QSO=0` (the environment of that shell was not recorded, so it presumably used
> the default, order 1; commit `71cfb72`), is kept as
> `results/capobs2_dgp2_qso1_normal_7000.*` and discussed below. Files with the
> suffix `_anterior` come from the earlier, superseded `g(w)`. DGP 1 is unaffected.

**Uniform treatment, DGP 1, default configuration** (four independent samples of 20):
-0.194 (SE 0.019), -0.162 (0.014), -0.101 (0.023), -0.183 (0.020).
T-HAL is better in 75 of 80 replicates.

**Normal treatment, DGP 1, default configuration** (five independent samples):
+0.201 (0.092), +0.117 (0.111), +0.328 (0.358), +0.027 (0.083), +0.118 (0.083).
T-HAL is better in 59 of 110 replicates: no reliable gain. The mean is driven by
a few replicates with very large errors at the right end of the dose range.

**Normal treatment, `w_cap` variants** (our own modification, see below; samples
with seeds chosen after the variants were devised):

| setting (replicates) | default (cap 103) | `"obs_max"` | cap 10 | cap 5 |
|---|---|---|---|---|
| DGP 1, seeds 7000+b (20) | +0.027 (0.083) | -0.113 (0.039) | -0.054 (0.055) | |
| DGP 1, seeds 9000+b (30) | +0.118 (0.083) | -0.148 (0.040) | -0.104 (0.050) | -0.161 (0.041) |
| DGP 2 (preprint), order 0, seeds 7000+b (20) | -0.044 (0.050) | -0.098 (0.042) | -0.078 (0.044) | -0.096 (0.041) |

Per-dose-point RMSE (the preprint's metric), normal treatment, 30 replicates of
DGP 1: HAL 0.698, default 0.732, `"obs_max"` 0.574. DGP 2 of the preprint (20
replicates, order 0): HAL 0.861, default 0.770, `"obs_max"` 0.735, cap 10 0.748,
cap 5 0.738.

On DGP 2 of the preprint (outcome HAL of order 0) all four T-HAL settings are on
average better than the plain HAL plug-in, but by a small margin. The paired
difference is -0.044 (SE 0.050) for the default cap, which cannot be told apart
from zero, and about two standard errors below zero for `"obs_max"`, cap 10 and
cap 5. T-HAL is better in 13, 12 and 11 of 20 replicates for those three settings
and in 9 for the default. The three capped settings beat the default cap in 15 to
18 of 20 replicates (mean paired differences -0.054, -0.034 and -0.052). The seeds
are the same as in the first DGP 1 sample, so this is not a fresh sample; with 20
replicates it is pilot evidence.

The first DGP 2 run, made without `THAL_QSO=0` and presumably with the outcome HAL of
order 1 (`results/capobs2_dgp2_qso1_*`), went the other way: T-HAL was worse than HAL in 18 to 19 of 20 replicates (paired
differences +0.299, +0.184, +0.229 and +0.179; per-dose RMSE: HAL 0.875, default
1.093, `"obs_max"` 1.002). The setting of the outcome HAL therefore changes the sign of
the comparison, and the conclusions above are specific to order 0.

`w_cap = "obs_max"` caps the counterfactual weights `omega(a) / g_n(a | W_i)` at
the largest weight seen at the observed data. **It is not in the preprint.** It
was devised after looking at failing replicates of DGP 1, then tested on fresh
seeds of DGP 1. On the preprint's DGP 2 (order 0) it beats the default cap in 18 of 20
replicates. A cap of 5 works about as well, so what seems to matter is
not extrapolating the correction to doses that are almost impossible for some
values of `W`; a fixed number such as 5 depends on the scale of the weights and
should not be used as a default.

## Tried, did not help (normal treatment)

* the exact `g` instead of the estimated one (worse);
* honest V-fold cross-validation of the penalty (worse);
* cross-fitting of `Q` with K = 2 and 5 (worse);
* clever covariate with a truncated Moore-Penrose inverse of `Sigma`
  (tolerances 1e-10, 1e-8, 1e-6) and bases with 150 and 500 knots: no gain with
  uniform treatment, no consistent gain with normal treatment.

## Differences from the preprint

| step of Algorithm 1 | this implementation | status |
|---|---|---|
| `Q_n`: HAL, order 0 for DGP 1/2 | `hal9001::fit_hal`, order set with `smoothness_orders` | other HAL settings unknown |
| `g_n`: `haldensify`; `omega`: kernel density | 10 bins, `max_degree = 2`, `trim = FALSE`, table lookup by bin; Sheather-Jones KDE | `haldensify` settings unknown |
| basis `phi_{1,N}`, `N` large | intercept + 50 knots at quantiles of `A` (`N = 51`) | `N` and knots unknown; larger bases were worse |
| `C_g = omega/g * Sigma^{-1} phi` | `omega/g * phi` (no `Sigma`) by default; truncated `Sigma^+` via `sigma_tol` | differs; `Sigma^+` did not help |
| LASSO with `Q_n` as offset | `glmnet`, Gaussian, no intercept, no standardisation | as described |
| penalty: CV clipped to `[6 n^(1/5), 9 n^(1/5)]` | `cv.glmnet` on in-sample `Q_n`, clipped (Section 5.1 version) | Algorithm 1 pseudocode has the inequalities reversed; we follow Section 5.1 |
| weights truncated at `sqrt(n) log(n/5)` | applied to observed and counterfactual weights (about 103 for `n = 500`) | as recommended |
| projection on selected `phi_j` | least squares on all observed doses | as described |
| evaluation | 20-30 replicates, pooled and paired summaries | preprint: 1,000 replicates |

Details the preprint does not specify (HAL settings, basis size and knots, `haldensify`
options, tolerance of the pseudo-inverse) are open questions that remain to be
clarified with the authors.

## Defaults and the settings used in the experiments

| setting | package default | used in the experiments |
|---|---|---|
| `w_cap` (counterfactual weights) | `NULL`: cap `sqrt(n) log(n/5)` as in the preprint | default, `"obs_max"`, 10 or 5 (see Evidence) |
| `n_proj` (doses used by the projection) | 200 | 500, i.e. all observed doses (`n_proj = length(A)`) |
| outcome HAL order (`smoothness_orders`, or `THAL_QSO`) | 1 | 0 |
| `n_knots`, `c1`, `c2`, `k` | 50, 6, 9, 1 | same |
| `fit_g_haldensify()` | 10 bins, `max_degree = 2`, `trim = FALSE` | same |
| penalty | CV clipped to the window | same |

With uniform treatment `w_cap` has no effect (all weights are 1). With normal
treatment the default `w_cap` is the unreliable configuration (see Evidence);
`"obs_max"` is experimental. A projection on only 25 doses was clearly worse
(pooled RMSE 0.760 against 0.380 with all doses, uniform treatment).

## Reproducing the experiments

* R 4.3.3 on Linux/WSL; versions of `glmnet`, `hal9001` and `haldensify` are in
  `results/package_versions.txt`, the session in `results/sessionInfo.txt`.
* Scripts are in `experiments/`, their logs and per-replicate results (`.rds`) in
  `results/`. Most scripts take the package directory as first argument and were
  run with `THAL_QSO=0` (order of the outcome HAL) from the parent directory.
* One replicate takes 3-6 minutes (the `haldensify` fit dominates).
* The logs in `results/` are the record of what was run. Some scripts were run against
  earlier development copies of the package (not included here) and hard-code those
  directory names, so they may need small edits to run against this tree. Early scripts
  also depend on the order in which random numbers are consumed, so rerunning them with a
  different package version can give different samples.

## Known issues / TODO

* `fit_Q_hal()` reads its default smoothness order from the environment variable
  `THAL_QSO` (default 1); it should be an explicit argument. Forgetting `THAL_QSO=0` silently changes the
  experiment (it happened once, see the DGP 2 note under Evidence); `check_capobs2.R`
  now prints it in its `cfg:` line.
* `predict.haldensify()` trims predicted densities at a floor that, according to
  the package NEWS, depends on the number of rows predicted; we saw
  batch-dependent outputs. We always call it with `trim = FALSE`.
* Influence-curve variance and confidence intervals are not implemented.
* The experiments live in the package tree for traceability; they are excluded
  from the built package (`.Rbuildignore`).

## License

AGPL-3.0-or-later. See `LICENSE`.

Copyright (C) 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
