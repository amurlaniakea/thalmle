Subject: Implementation details of the Targeted HAL-MLE (arXiv:2607.03824)

Dear Dr. Rodriguez, Dr. Díaz-Ordaz, Dr. Lehmann and Prof. van der Laan,

I have written an independent R implementation of the T-HAL-MLE for a continuous-dose
dose-response curve, following Algorithm 1 and Section 7, to reproduce the single-
discontinuity DGP (n = 500). With uniform treatment it reproduces the qualitative
ordering of the preprint (T-HAL better than the HAL plug-in in most replicates). With
normal treatment it does not: the targeting step is occasionally unstable at the right
end of the dose range. Before drawing conclusions I would like to be sure I am not
misreading the algorithm. Could you clarify:

1. Outcome HAL: besides the order (0 for DGP 1/2), the values of max_degree,
   num_knots and reduce_basis in hal9001.
2. haldensify: n_bins, grid_type, max_degree, smoothness_orders and the lambda grid; and
   whether predict() was called with its default trimming when g was evaluated at the
   counterfactual pairs (a, W_i) -- the default floor seems to depend on the number of
   rows predicted.
3. The basis phi_{1,N}: its size N, its knots (all observed doses? a grid?) and the
   tolerance of the Moore-Penrose inverse of Sigma.
4. Whether the point estimate uses C_g = (omega/g) Sigma^{-1} phi exactly, and whether
   Q_n is cross-fitted (K = 2) in the point estimate or only in the variance estimator.
5. Whether the weights are truncated at sqrt(n) log(n/5) also at the counterfactual pairs
   in Step 6.
6. In Algorithm 1, Step 5, the inequalities appear reversed with respect to Section 5.1;
   I followed Section 5.1. Could you confirm?
7. Whether code is available.

I am happy to share my code and logs. Thank you for your time.

Best regards,
Pedro Sordo Martínez
[correo]
