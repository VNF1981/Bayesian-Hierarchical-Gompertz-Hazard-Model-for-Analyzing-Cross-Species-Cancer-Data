###################################################################################################
# 08_Gompertz_Assurance_body_mass_array_worker.R
#
# CLEAN FINAL assurance worker for the Gompertz malignancy model.
#
# Purpose:
#   One SLURM array task = one simulated dataset + one model refit.
#
# Design:
#   5 assumed body-mass effects on log(alpha_mid):
#     0.1, 0.3, 0.5, 0.7, 0.9
#   50 simulations per effect
#   Total tasks = 250
#
# Task mapping:
#   tasks   1-50  -> alpha_body_mass = 0.1
#   tasks  51-100 -> alpha_body_mass = 0.3
#   tasks 101-150 -> alpha_body_mass = 0.5
#   tasks 151-200 -> alpha_body_mass = 0.7
#   tasks 201-250 -> alpha_body_mass = 0.9
#
# This script uses the final model file and the final saved assurance input file.
# It uses exact Stan data names from the final Gompertz model.
###################################################################################################

rm(list = ls(all = TRUE))
ls()

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(readr)
  library(tibble)
})

###################################################################################################
# 0. Fixed project settings
###################################################################################################

project_dir <- "/projects/tollis_lab/TE/mammals/Gompertz_model"
setwd(project_dir)

cmdstanr::set_cmdstan_path("/scratch/vn229/cmdstan-2.38.0")

real_model_dir <- file.path(project_dir, "results_07_full_phylo_lambda")
assurance_out_dir <- file.path(project_dir, "results_08_assurance_body_mass")

dir.create(assurance_out_dir, showWarnings = FALSE, recursive = TRUE)

stan_file <- file.path(
  project_dir,
  "07_gompertz_malignancy_alpha_mid_beta_LFH_phylo_alpha_beta.stan"
)

fit_real_file <- file.path(
  real_model_dir,
  "07_gompertz_malignancy_phylo_alpha_beta_lambda_full_fit.rds"
)

assurance_input_file <- file.path(
  real_model_dir,
  "07_gompertz_malignancy_phylo_alpha_beta_lambda_full_assurance_inputs.rds"
)

if (!file.exists(stan_file)) {
  stop("Stan file not found: ", stan_file)
}

if (!file.exists(fit_real_file)) {
  stop("Full model fit file not found: ", fit_real_file)
}

if (!file.exists(assurance_input_file)) {
  stop("Assurance input file not found: ", assurance_input_file)
}

###################################################################################################
# 1. Assurance settings
###################################################################################################

effect_grid_body_mass <- c(0.1, 0.3, 0.5, 0.7, 0.9)
n_sim_per_effect <- 50
n_total_tasks <- length(effect_grid_body_mass) * n_sim_per_effect

posterior_prob_cutoff <- 0.95

chains <- 2
parallel_chains <- 2
iter_warmup <- 500
iter_sampling <- 500
adapt_delta <- 0.99
max_treedepth <- 15

base_seed <- 20260515

###################################################################################################
# 2. SLURM task mapping
###################################################################################################

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (is.na(task_id)) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) {
    task_id <- as.integer(args[1])
  }
}

if (is.na(task_id)) {
  stop("No SLURM_ARRAY_TASK_ID found and no task id argument was provided.")
}

if (task_id < 1 || task_id > n_total_tasks) {
  stop("Task id must be between 1 and ", n_total_tasks, ". Got: ", task_id)
}

effect_index <- ceiling(task_id / n_sim_per_effect)
sim_i <- task_id - (effect_index - 1) * n_sim_per_effect
assumed_alpha_body_mass <- effect_grid_body_mass[effect_index]

cat("Task id:", task_id, "\n")
cat("Effect index:", effect_index, "\n")
cat("Simulation replicate:", sim_i, "\n")
cat("Assumed alpha_body_mass:", assumed_alpha_body_mass, "\n")

result_file <- file.path(
  assurance_out_dir,
  sprintf(
    "assurance_body_mass_alpha_%0.1f_sim_%03d_task_%03d.csv",
    assumed_alpha_body_mass,
    sim_i,
    task_id
  )
)

if (file.exists(result_file)) {
  cat("Result file already exists. Skipping:", result_file, "\n")
  quit(save = "no", status = 0)
}

###################################################################################################
# 3. Helper functions
###################################################################################################

get_post_median <- function(draws_df, par_name, default_value = NA_real_) {
  if (par_name %in% names(draws_df)) {
    return(stats::median(draws_df[[par_name]], na.rm = TRUE))
  }

  if (is.na(default_value)) {
    stop("Could not find parameter in posterior draws: ", par_name)
  }

  return(default_value)
}

simulate_gompertz_event_times <- function(alpha_mid, beta, t_ref, sex_linear_effect) {
  n <- length(alpha_mid)

  u <- stats::runif(n)
  target_H <- -log(u)

  event_time <- rep(Inf, n)

  beta_small <- abs(beta) < 1e-8

  if (any(beta_small)) {
    rate <- exp(sex_linear_effect[beta_small]) * alpha_mid[beta_small]
    event_time[beta_small] <- target_H[beta_small] / rate
  }

  if (any(!beta_small)) {
    idx <- which(!beta_small)

    b <- beta[idx]
    a <- alpha_mid[idx]
    tref <- t_ref[idx]
    se <- exp(sex_linear_effect[idx])
    Htarget <- target_H[idx]

    rhs <- 1 + Htarget * b / (se * a * exp(-b * tref))

    valid <- rhs > 0 & is.finite(rhs)

    tmp <- rep(Inf, length(idx))
    tmp[valid] <- log(rhs[valid]) / b[valid]

    event_time[idx] <- tmp
  }

  event_time[event_time < 0] <- Inf

  return(event_time)
}

###################################################################################################
# 4. Load final full model and assurance inputs
###################################################################################################

fit_real <- readRDS(fit_real_file)
assurance_inputs <- readRDS(assurance_input_file)

stan_data_real <- assurance_inputs$stan_data

if (is.null(stan_data_real)) {
  stop("Assurance input RDS must contain a list element named stan_data.")
}

###################################################################################################
# 5. Strict checks for final Stan data names
#
# Final model requires exactly:
#   N, S, species_id, age_decades, ref_age_decades, event, sex_male,
#   log_body_mass_z, log_longevity_z, log_gestation_z, L_phylo
###################################################################################################

required_names <- c(
  "N",
  "S",
  "species_id",
  "age_decades",
  "ref_age_decades",
  "event",
  "sex_male",
  "log_body_mass_z",
  "log_longevity_z",
  "log_gestation_z",
  "L_phylo"
)

missing_names <- setdiff(required_names, names(stan_data_real))

if (length(missing_names) > 0) {
  cat("Names in stan_data_real:\n")
  print(names(stan_data_real))
  stop("Missing required final Stan data names: ", paste(missing_names, collapse = ", "))
}

N <- stan_data_real$N
S <- stan_data_real$S

stopifnot(N == length(stan_data_real$species_id))
stopifnot(N == length(stan_data_real$age_decades))
stopifnot(N == length(stan_data_real$event))
stopifnot(N == length(stan_data_real$sex_male))
stopifnot(S == length(stan_data_real$ref_age_decades))
stopifnot(S == length(stan_data_real$log_body_mass_z))
stopifnot(S == length(stan_data_real$log_longevity_z))
stopifnot(S == length(stan_data_real$log_gestation_z))
stopifnot(all(dim(stan_data_real$L_phylo) == c(S, S)))

cat("Loaded assurance Stan data.\n")
cat("Animals:", N, "\n")
cat("Species:", S, "\n")
cat("Observed malignancy events:", sum(stan_data_real$event), "\n")

###################################################################################################
# 6. Extract nuisance parameter medians from real full posterior
###################################################################################################

draws_real <- posterior::as_draws_df(fit_real$draws())

alpha_intercept_real <- get_post_median(draws_real, "alpha_intercept")
alpha_longevity_real <- get_post_median(draws_real, "alpha_longevity", default_value = 0)
alpha_gestation_real <- get_post_median(draws_real, "alpha_gestation", default_value = 0)

sigma_phylo_alpha_real <- get_post_median(draws_real, "sigma_phylo_alpha", default_value = 0)
sigma_alpha_resid_real <- get_post_median(draws_real, "sigma_alpha_resid", default_value = 0)

beta_intercept_real <- get_post_median(draws_real, "beta_intercept")
beta_body_mass_real <- get_post_median(draws_real, "beta_body_mass", default_value = 0)
beta_longevity_real <- get_post_median(draws_real, "beta_longevity", default_value = 0)
beta_gestation_real <- get_post_median(draws_real, "beta_gestation", default_value = 0)

sigma_phylo_beta_real <- get_post_median(draws_real, "sigma_phylo_beta", default_value = 0)
sigma_beta_resid_real <- get_post_median(draws_real, "sigma_beta_resid", default_value = 0)

gamma_sex_real <- get_post_median(draws_real, "gamma_sex", default_value = 0)

###################################################################################################
# 7. Simulate one dataset under the assumed alpha_body_mass effect
###################################################################################################

set.seed(base_seed + task_id + round(assumed_alpha_body_mass * 1000))

species_id <- as.integer(stan_data_real$species_id)

body_mass_species <- as.numeric(stan_data_real$log_body_mass_z)
longevity_species <- as.numeric(stan_data_real$log_longevity_z)
gestation_species <- as.numeric(stan_data_real$log_gestation_z)

z_phylo_alpha <- stats::rnorm(S)
z_alpha_resid <- stats::rnorm(S)
z_phylo_beta <- stats::rnorm(S)
z_beta_resid <- stats::rnorm(S)

phylo_alpha <- sigma_phylo_alpha_real * as.vector(stan_data_real$L_phylo %*% z_phylo_alpha)
alpha_resid <- sigma_alpha_resid_real * z_alpha_resid

phylo_beta <- sigma_phylo_beta_real * as.vector(stan_data_real$L_phylo %*% z_phylo_beta)
beta_resid <- sigma_beta_resid_real * z_beta_resid

log_alpha_mid_species <-
  alpha_intercept_real +
  assumed_alpha_body_mass * body_mass_species +
  alpha_longevity_real * longevity_species +
  alpha_gestation_real * gestation_species +
  phylo_alpha +
  alpha_resid

beta_species <-
  beta_intercept_real +
  beta_body_mass_real * body_mass_species +
  beta_longevity_real * longevity_species +
  beta_gestation_real * gestation_species +
  phylo_beta +
  beta_resid

alpha_mid_ind <- exp(log_alpha_mid_species[species_id])
beta_ind <- beta_species[species_id]
t_ref_ind <- stan_data_real$ref_age_decades[species_id]
sex_linear_effect <- gamma_sex_real * stan_data_real$sex_male

event_time <- simulate_gompertz_event_times(
  alpha_mid = alpha_mid_ind,
  beta = beta_ind,
  t_ref = t_ref_ind,
  sex_linear_effect = sex_linear_effect
)

# Preserve observed follow-up structure.
y_sim <- as.integer(event_time <= stan_data_real$age_decades)

cat("Simulated malignancy events:", sum(y_sim == 1), "\n")
cat("Simulated censored animals:", sum(y_sim == 0), "\n")

stan_data_sim <- stan_data_real
stan_data_sim$event <- y_sim

###################################################################################################
# 8. Refit final model to simulated data
###################################################################################################

gompertz_model <- cmdstanr::cmdstan_model(stan_file)

init_fun <- function() {
  list(
    alpha_intercept = -4,
    alpha_body_mass = 0,
    alpha_longevity = 0,
    alpha_gestation = 0,
    sigma_phylo_alpha = 0.3,
    sigma_alpha_resid = 0.3,
    z_phylo_alpha = stats::rnorm(S, 0, 0.01),
    z_alpha_resid = stats::rnorm(S, 0, 0.01),

    beta_intercept = 0.5,
    beta_body_mass = 0,
    beta_longevity = 0,
    beta_gestation = 0,
    sigma_phylo_beta = 0.02,
    sigma_beta_resid = 0.02,
    z_phylo_beta = stats::rnorm(S, 0, 0.01),
    z_beta_resid = stats::rnorm(S, 0, 0.01),

    gamma_sex = 0
  )
}

fit_sim <- try(
  gompertz_model$sample(
    data = stan_data_sim,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = base_seed + task_id + round(assumed_alpha_body_mass * 1000),
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    init = init_fun,
    refresh = 100,
    output_dir = assurance_out_dir
  ),
  silent = TRUE
)

###################################################################################################
# 9. Save failure row if fit fails
###################################################################################################

if (inherits(fit_sim, "try-error")) {
  result_fail <- tibble(
    task_id = task_id,
    effect_index = effect_index,
    sim_i = sim_i,
    target_parameter = "alpha_body_mass",
    direction = "positive",
    assumed_alpha_body_mass = assumed_alpha_body_mass,
    assumed_HR = exp(assumed_alpha_body_mass),
    n_events_simulated = sum(y_sim == 1),
    n_censored_simulated = sum(y_sim == 0),
    fit_success = FALSE,
    posterior_prob_cutoff = posterior_prob_cutoff,
    posterior_prob_positive = NA_real_,
    median_estimated_alpha_body_mass = NA_real_,
    q5_estimated_alpha_body_mass = NA_real_,
    q95_estimated_alpha_body_mass = NA_real_,
    success_posterior_prob = NA,
    success_q5_positive = NA,
    rhat_alpha_body_mass = NA_real_,
    ess_bulk_alpha_body_mass = NA_real_,
    ess_tail_alpha_body_mass = NA_real_,
    num_divergent = NA_integer_,
    num_max_treedepth = NA_integer_
  )

  readr::write_csv(result_fail, result_file)
  stop("Simulation fit failed. Failure result saved: ", result_file)
}

###################################################################################################
# 10. Extract posterior support for alpha_body_mass
###################################################################################################

draws_sim <- posterior::as_draws_df(fit_sim$draws())

if (!"alpha_body_mass" %in% names(draws_sim)) {
  stop("Could not find alpha_body_mass in posterior draws from simulated fit.")
}

alpha_body_mass_draws <- draws_sim$alpha_body_mass

posterior_prob_positive <- mean(alpha_body_mass_draws > 0, na.rm = TRUE)
median_estimated_alpha_body_mass <- stats::median(alpha_body_mass_draws, na.rm = TRUE)
q5_estimated_alpha_body_mass <- stats::quantile(alpha_body_mass_draws, probs = 0.05, na.rm = TRUE)
q95_estimated_alpha_body_mass <- stats::quantile(alpha_body_mass_draws, probs = 0.95, na.rm = TRUE)

success_posterior_prob <- posterior_prob_positive > posterior_prob_cutoff
success_q5_positive <- q5_estimated_alpha_body_mass > 0

summary_sim <- fit_sim$summary(variables = "alpha_body_mass")
sampler_diag <- fit_sim$sampler_diagnostics()

num_divergent <- sum(sampler_diag[, , "divergent__"])
num_max_treedepth <- sum(sampler_diag[, , "treedepth__"] >= max_treedepth)

result_success <- tibble(
  task_id = task_id,
  effect_index = effect_index,
  sim_i = sim_i,
  target_parameter = "alpha_body_mass",
  direction = "positive",
  assumed_alpha_body_mass = assumed_alpha_body_mass,
  assumed_HR = exp(assumed_alpha_body_mass),
  n_events_simulated = sum(y_sim == 1),
  n_censored_simulated = sum(y_sim == 0),
  fit_success = TRUE,
  posterior_prob_cutoff = posterior_prob_cutoff,
  posterior_prob_positive = posterior_prob_positive,
  median_estimated_alpha_body_mass = median_estimated_alpha_body_mass,
  q5_estimated_alpha_body_mass = as.numeric(q5_estimated_alpha_body_mass),
  q95_estimated_alpha_body_mass = as.numeric(q95_estimated_alpha_body_mass),
  success_posterior_prob = success_posterior_prob,
  success_q5_positive = success_q5_positive,
  rhat_alpha_body_mass = summary_sim$rhat[1],
  ess_bulk_alpha_body_mass = summary_sim$ess_bulk[1],
  ess_tail_alpha_body_mass = summary_sim$ess_tail[1],
  num_divergent = num_divergent,
  num_max_treedepth = num_max_treedepth
)

readr::write_csv(result_success, result_file)

cat("Saved result file:", result_file, "\n")
cat("Success posterior probability:", success_posterior_prob, "\n")
cat("Success q5 positive:", success_q5_positive, "\n")
cat("Posterior probability alpha_body_mass > 0:", posterior_prob_positive, "\n")
cat("Divergences:", num_divergent, "\n")
cat("Max treedepth:", num_max_treedepth, "\n")

###################################################################################################
# End
###################################################################################################
