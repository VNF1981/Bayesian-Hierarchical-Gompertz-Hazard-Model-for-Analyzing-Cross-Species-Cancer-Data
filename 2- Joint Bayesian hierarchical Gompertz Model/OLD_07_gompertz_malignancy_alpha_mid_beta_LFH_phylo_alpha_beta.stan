functions {
  // Calculate the Gompertz cumulative hazard with a species specific reference age.
  //
  // Hazard:
  // h(t) = alpha_mid * exp(beta * (t - t_ref) + sex_effect)
  //
  // t = animal age in decades
  // t_ref = species specific mid lifespan in decades
  // alpha_mid = species specific malignancy hazard at mid lifespan
  // beta = species specific age related change in malignancy hazard per decade
  real gompertz_cumulative_hazard_midref(
    real t,
    real t_ref,
    real log_alpha_mid,
    real beta,
    real sex_effect
  ) {
    real baseline_scale;
    real cumulative_hazard;

    if (abs(beta) < 1e-6) {
      cumulative_hazard = exp(log_alpha_mid + sex_effect) * t;
    } else {
      baseline_scale = exp(log_alpha_mid + sex_effect - beta * t_ref);
      cumulative_hazard = baseline_scale * expm1(beta * t) / beta;
    }

    return cumulative_hazard;
  }
}

data {
  // Number of individual animals
  int<lower=1> N;

  // Number of species
  int<lower=1> S;

  // Species ID for each animal
  array[N] int<lower=1, upper=S> species_id;

  // Animal age in decades
  vector<lower=0>[N] age_decades;

  // Species specific reference age in decades
  // Here this is mid lifespan:
  // t_ref = 0.5 * max longevity
  vector<lower=0>[S] ref_age_decades;

  // Event indicator
  // 1 = malignancy observed
  // 0 = censored
  array[N] int<lower=0, upper=1> event;

  // Sex indicator
  // 0 = female
  // 1 = male
  vector[N] sex_male;

  // Species level life history predictors
  // These are log transformed and standardized in R.
  vector[S] log_body_mass_z;
  vector[S] log_longevity_z;
  vector[S] log_gestation_z;

  // Lower Cholesky factor of the phylogenetic correlation matrix
  // The R script creates either a Brownian motion matrix or a Pagel lambda matrix.
  matrix[S, S] L_phylo;
}

parameters {
  // Intercept for log malignancy hazard at mid lifespan
  real alpha_intercept;

  // Life history effects on log alpha_mid
  real alpha_body_mass;
  real alpha_longevity;
  real alpha_gestation;

  // Strength of phylogenetic variation in log alpha_mid
  real<lower=1e-6> sigma_phylo_alpha;

  // Strength of non phylogenetic residual variation in log alpha_mid
  real<lower=1e-6> sigma_alpha_resid;

  // Standard normal variables used to create phylogenetic alpha effects
  vector[S] z_phylo_alpha;

  // Standard normal variables used to create non phylogenetic alpha residuals
  vector[S] z_alpha_resid;

  // Intercept for age related malignancy hazard change per decade
  real beta_intercept;

  // Life history effects on beta
  real beta_body_mass;
  real beta_longevity;
  real beta_gestation;

  // Strength of phylogenetic variation in beta
  real<lower=1e-6> sigma_phylo_beta;

  // Strength of non phylogenetic residual variation in beta
  real<lower=1e-6> sigma_beta_resid;

  // Standard normal variables used to create phylogenetic beta effects
  vector[S] z_phylo_beta;

  // Standard normal variables used to create non phylogenetic beta residuals
  vector[S] z_beta_resid;

  // Individual level sex effect on the log hazard scale
  real gamma_sex;
}

transformed parameters {
  // Species specific log alpha at mid lifespan
  vector[S] log_alpha_mid;

  // Species specific beta
  vector[S] beta;

  // Expected log alpha_mid from life history predictors
  vector[S] log_alpha_mid_mean;

  // Expected beta from life history predictors
  vector[S] beta_mean;

  // Phylogenetic and non phylogenetic residuals for alpha_mid
  vector[S] phylo_alpha;
  vector[S] alpha_resid;

  // Phylogenetic and non phylogenetic residuals for beta
  vector[S] phylo_beta;
  vector[S] beta_resid;

  // Create phylogenetically correlated residuals
  phylo_alpha = sigma_phylo_alpha * (L_phylo * z_phylo_alpha);
  phylo_beta = sigma_phylo_beta * (L_phylo * z_phylo_beta);

  // Create independent residuals
  alpha_resid = sigma_alpha_resid * z_alpha_resid;
  beta_resid = sigma_beta_resid * z_beta_resid;

  for (s in 1:S) {
    log_alpha_mid_mean[s] =
      alpha_intercept +
      alpha_body_mass * log_body_mass_z[s] +
      alpha_longevity * log_longevity_z[s] +
      alpha_gestation * log_gestation_z[s];

    beta_mean[s] =
      beta_intercept +
      beta_body_mass * log_body_mass_z[s] +
      beta_longevity * log_longevity_z[s] +
      beta_gestation * log_gestation_z[s];
  }

  // Alpha includes life history effects, a phylogenetic residual, and a non phylogenetic residual.
  log_alpha_mid = log_alpha_mid_mean + phylo_alpha + alpha_resid;

  // Beta also includes life history effects, a phylogenetic residual, and a non phylogenetic residual.
  beta = beta_mean + phylo_beta + beta_resid;
}

model {
  // Priors for life history effects on log alpha_mid
  alpha_intercept ~ normal(-4, 1.5);
  alpha_body_mass ~ normal(0, 1);
  alpha_longevity ~ normal(0, 1);
  alpha_gestation ~ normal(0, 1);

  // Priors for phylogenetic and non phylogenetic alpha variation
  sigma_phylo_alpha ~ exponential(2);
  sigma_alpha_resid ~ exponential(2);

  // Standard normal residuals for alpha
  z_phylo_alpha ~ normal(0, 1);
  z_alpha_resid ~ normal(0, 1);

  // Priors for life history effects on beta
  // beta is per decade and can be positive, near zero, or negative.
  beta_intercept ~ normal(0, 0.5);
  beta_body_mass ~ normal(0, 0.5);
  beta_longevity ~ normal(0, 0.5);
  beta_gestation ~ normal(0, 0.5);

  // Priors for phylogenetic and non phylogenetic beta variation
  // These are more regularizing because beta was harder to estimate.
  sigma_phylo_beta ~ exponential(20);
  sigma_beta_resid ~ exponential(20);

  // Standard normal residuals for beta
  z_phylo_beta ~ normal(0, 1);
  z_beta_resid ~ normal(0, 1);

  // Prior for sex effect
  gamma_sex ~ normal(0, 0.5);

  // Gompertz survival likelihood with right censoring
  for (n in 1:N) {
    int sp;
    real log_hazard_n;
    real cumulative_hazard_n;

    sp = species_id[n];

    // Log hazard at observed age
    log_hazard_n =
      log_alpha_mid[sp] +
      beta[sp] * (age_decades[n] - ref_age_decades[sp]) +
      gamma_sex * sex_male[n];

    // Cumulative hazard from age 0 to observed age
    cumulative_hazard_n = gompertz_cumulative_hazard_midref(
      age_decades[n],
      ref_age_decades[sp],
      log_alpha_mid[sp],
      beta[sp],
      gamma_sex * sex_male[n]
    );

    // Event animals contribute log hazard minus cumulative hazard.
    // Censored animals contribute only minus cumulative hazard.
    if (event[n] == 1) {
      target += log_hazard_n - cumulative_hazard_n;
    } else {
      target += -cumulative_hazard_n;
    }
  }
}

generated quantities {
  // Species specific alpha_mid on the original positive hazard scale
  vector[S] alpha_mid;

  // Proportion of alpha_mid residual variance assigned to phylogeny
  real phylo_fraction_alpha;

  // Proportion of beta residual variance assigned to phylogeny
  real phylo_fraction_beta;

  for (s in 1:S) {
    alpha_mid[s] = exp(log_alpha_mid[s]);
  }

  phylo_fraction_alpha =
    square(sigma_phylo_alpha) /
    (square(sigma_phylo_alpha) + square(sigma_alpha_resid));

  phylo_fraction_beta =
    square(sigma_phylo_beta) /
    (square(sigma_phylo_beta) + square(sigma_beta_resid));
}
