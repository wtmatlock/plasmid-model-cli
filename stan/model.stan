data {
  int<lower=1> N;                     // number of tips
  int<lower=1> B;                     // number of blocks
  int<lower=1> K;                     // number of retained eigenvectors (rank)
  array[N,B] int<lower=0,upper=1> Y;  // presence/absence matrix
  matrix[N,K] U;                      // top-K eigenvectors of phylo covariance
  vector[K] lambda;                   // top-K eigenvalues of phylo covariance
}
parameters {
  real mu_bar;                        // overall intercept
  vector<lower=0>[B] sigma_mu;        // scale per-block for intercepts
  cholesky_factor_corr[B] L_block;    // block correlation cholesky (LKJ)
  vector[B] z_mu;                     // standard normals for non-centred mu

  vector<lower=0>[B] sigma_phy;       // phylo sd per block
  vector<lower=0,upper=1>[B] k;       // phylo-structure score per block
  matrix[B,K] z_phy_std;              // standard normals: rows = blocks, cols = eigenvectors
}
transformed parameters {
  // correlated block intercepts (non-centred)
  vector[B] mu = mu_bar + diag_pre_multiply(sigma_mu, L_block) * z_mu;

  // Build the phylo latent S (K x B) such that columns index blocks
  matrix[K,B] S;
  for (k_idx in 1:K) {
    vector[B] z_k = col(z_phy_std, k_idx);     // standard normals across blocks
    vector[B] v_k = L_block * z_k;             // induce block correlation
    for (b in 1:B) {
      real s_bk = sigma_phy[b] * sqrt(k[b] * lambda[k_idx] + (1 - k[b]));
      S[k_idx, b] = s_bk * v_k[b];
    }
  }

  matrix[N,B] phy = U * S;            // phylogenetic contribution at tips

  // Centre phy per block to avoid intercept confounding
  for (b in 1:B) {
    real mn = mean(phy[, b]);
    phy[, b] = phy[, b] - mn;
  }
}
model {
  // Priors
  mu_bar ~ normal(0, 1);
  sigma_mu ~ student_t(3, 0, 0.5);   // half-t
  to_vector(z_mu) ~ normal(0, 1);
  L_block ~ lkj_corr_cholesky(2);    // weakly informative LKJ

  sigma_phy ~ student_t(3, 0, 1);    // half-t
  k ~ beta(1, 1);
  to_vector(z_phy_std) ~ normal(0, 1);

  // Likelihood
  for (b in 1:B)
    Y[, b] ~ bernoulli_logit(mu[b] + phy[, b]);
}
generated quantities {
  matrix[N,B] p;
  for (b in 1:B)
    p[, b] = inv_logit(mu[b] + phy[, b]);

  // full block correlation matrix (from cholesky)
  matrix[B,B] block_corr = multiply_lower_tri_self_transpose(L_block);
}
