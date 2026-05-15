######################################################################################################################################
######################## Monsoon version
######################## Bayesian hierarchical phylogenetic Gompertz model for malignancy
######################## Input:  gompertz_model_input.rds and min20Fixed516.nwk
######################## Output: fitted Stan model, posterior summaries, and assurance input RDS
######################## Model: species specific alpha_mid, species specific beta, sex effect, LFH effects, and phylogeny
######################## Phylogeny is added to both alpha_mid and beta
######################## Default phylogenetic mode: Pagel lambda with lambda = 0.46
######################################################################################################################################

rm(list = ls(all = TRUE))
ls()

###################################################################################################
# step 1
# Monsoon project paths
###################################################################################################

project_dir <- "/projects/tollis_lab/TE/mammals/Gompertz_model"
setwd(project_dir)

output_dir <- file.path(project_dir, "results_07_full_phylo_lambda")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

###################################################################################################
# step 2
# Load packages
###################################################################################################

library(cmdstanr)
library(posterior)
library(dplyr)
library(readr)
library(ape)

###################################################################################################
# step 3
# CmdStan path
###################################################################################################

cmdstanr::set_cmdstan_path("/scratch/vn229/cmdstan-2.38.0")

###################################################################################################
# step 4
# Run settings
###################################################################################################

run_mode <- "full"
phylo_mode <- "lambda"
pagel_lambda <- 0.46

input_file <- file.path(project_dir, "gompertz_model_input.rds")
tree_file  <- file.path(project_dir, "min20Fixed516.nwk")
stan_file  <- file.path(project_dir, "07_gompertz_malignancy_alpha_mid_beta_LFH_phylo_alpha_beta.stan")

fit_output_file     <- file.path(output_dir, paste0("07_gompertz_malignancy_phylo_alpha_beta_", phylo_mode, "_", run_mode, "_fit.rds"))
summary_output_file <- file.path(output_dir, paste0("07_gompertz_malignancy_phylo_alpha_beta_", phylo_mode, "_", run_mode, "_summary.csv"))
effects_output_file <- file.path(output_dir, paste0("07_gompertz_malignancy_phylo_alpha_beta_", phylo_mode, "_", run_mode, "_effects.csv"))
species_output_file <- file.path(output_dir, paste0("07_gompertz_malignancy_phylo_alpha_beta_", phylo_mode, "_", run_mode, "_species_summary.csv"))

assurance_input_file <- file.path(
  output_dir,
  "07_gompertz_malignancy_phylo_alpha_beta_lambda_full_assurance_inputs.rds"
)

###################################################################################################
# step 5
# Parallel settings
# This model uses chain level parallelism. Request at least 4 CPUs for 4 parallel chains.
###################################################################################################

chains <- 4
parallel_chains <- 4
iter_warmup <- 1000
iter_sampling <- 1000
adapt_delta <- 0.98
max_treedepth <- 15

options(mc.cores = parallel::detectCores())

###################################################################################################
# step 6
# Read prepared model input
###################################################################################################

gompertz_input <- readRDS(input_file)

individual_data <- gompertz_input$individual_data
species_predictors <- gompertz_input$species_predictors

###################################################################################################
# step 7
# Use malignancy as the event response
# 1 = malignancy observed
# 0 = censored
###################################################################################################

if (!"malignancy" %in% names(individual_data)) {
  cat("Available columns in individual_data:\n")
  print(names(individual_data))
  stop("The column 'malignancy' was not found in individual_data.")
}

individual_data <- individual_data %>%
  mutate(
    event = as.integer(malignancy),
    censored = as.integer(event == 0)
  )

bad_event_values <- setdiff(unique(individual_data$event), c(0, 1))

if (length(bad_event_values) > 0) {
  stop("The malignancy event column must contain only 0 and 1.")
}

cat("Using malignancy as the event response.\n")
cat("Total malignancy events:", sum(individual_data$event), "\n")
cat("Total censored animals:", sum(individual_data$censored), "\n")

###################################################################################################
# step 8
# Clean species names in data and tree style
###################################################################################################

individual_data <- individual_data %>%
  mutate(
    Species = trimws(as.character(Species)),
    Species = gsub(" ", "_", Species)
  )

species_predictors <- species_predictors %>%
  mutate(
    Species = trimws(as.character(Species)),
    Species = gsub(" ", "_", Species)
  )

###################################################################################################
# step 9
# Add species specific reference age
# ref_age_decades = 0.5 * maximum longevity in years / 10
# max_longevity_M is in months, so ref_age_decades = max_longevity_M / 240
###################################################################################################

species_predictors <- species_predictors %>%
  arrange(Species) %>%
  mutate(ref_age_decades = max_longevity_M / 240)

###################################################################################################
# step 10
# Read and clean tree
###################################################################################################

tree <- read.tree(tree_file)
tree$tip.label <- trimws(tree$tip.label)
tree$tip.label <- gsub(" ", "_", tree$tip.label)

###################################################################################################
# step 11
# Match species between data and tree
###################################################################################################

species_in_data <- unique(species_predictors$Species)
species_in_tree <- tree$tip.label
species_keep <- intersect(species_in_data, species_in_tree)

cat("Species in data:", length(species_in_data), "\n")
cat("Species in tree:", length(species_in_tree), "\n")
cat("Overlap species:", length(species_keep), "\n")

if (length(species_keep) < 2) {
  cat("\nExample species in data:\n")
  print(head(sort(species_in_data), 20))
  cat("\nExample species in tree:\n")
  print(head(sort(species_in_tree), 20))
  stop("Fewer than two species overlap between the tree and the dataset after cleaning. Check species names.")
}

species_missing_from_tree <- setdiff(species_in_data, species_in_tree)

if (length(species_missing_from_tree) > 0) {
  write_csv(
    tibble(Species_missing_from_tree = species_missing_from_tree),
    file.path(output_dir, "species_missing_from_tree.csv")
  )
  cat("Species missing from tree saved to species_missing_from_tree.csv\n")
}

###################################################################################################
# step 12
# Prune tree and data to shared species
###################################################################################################

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, species_keep))

species_predictors_phylo <- species_predictors %>%
  filter(Species %in% species_keep)

individual_data_phylo <- individual_data %>%
  filter(Species %in% species_keep)

###################################################################################################
# step 13
# Build phylogenetic correlation matrix
###################################################################################################

A_bm <- ape::vcv(tree_pruned, corr = TRUE)

species_predictors_phylo <- species_predictors_phylo %>%
  slice(match(rownames(A_bm), Species))

if (!all(species_predictors_phylo$Species == rownames(A_bm))) {
  stop("Species order does not match the phylogenetic covariance matrix.")
}

if (phylo_mode == "BM") {
  A_phylo <- A_bm
}

if (phylo_mode == "lambda") {
  A_phylo <- A_bm
  off_diag <- row(A_phylo) != col(A_phylo)
  A_phylo[off_diag] <- pagel_lambda * A_phylo[off_diag]
  diag(A_phylo) <- 1
}

if (!phylo_mode %in% c("BM", "lambda")) {
  stop("phylo_mode must be either 'BM' or 'lambda'.")
}

A_phylo <- A_phylo + diag(1e-6, nrow(A_phylo))

# R's chol() returns upper triangular U, Stan needs lower triangular L.
L_phylo <- t(chol(A_phylo))

###################################################################################################
# step 14
# Create tree aligned species IDs
###################################################################################################

species_lookup_phylo <- species_predictors_phylo %>%
  mutate(species_id_new = row_number()) %>%
  select(Species, species_id_new)

individual_data_phylo <- individual_data_phylo %>%
  select(-species_id) %>%
  left_join(species_lookup_phylo, by = "Species") %>%
  rename(species_id = species_id_new) %>%
  arrange(species_id, individual_id)

species_predictors_phylo <- species_predictors_phylo %>%
  select(-species_id) %>%
  left_join(species_lookup_phylo, by = "Species") %>%
  rename(species_id = species_id_new) %>%
  arrange(species_id)

###################################################################################################
# step 15
# Create Stan data
###################################################################################################

stan_data_use <- list(
  N = nrow(individual_data_phylo),
  S = nrow(species_predictors_phylo),
  species_id = individual_data_phylo$species_id,
  age_decades = individual_data_phylo$age_years / 10,
  ref_age_decades = species_predictors_phylo$ref_age_decades,
  event = individual_data_phylo$event,
  sex_male = individual_data_phylo$sex_male,
  log_body_mass_z = species_predictors_phylo$log_body_mass_z,
  log_longevity_z = species_predictors_phylo$log_longevity_z,
  log_gestation_z = species_predictors_phylo$log_gestation_z,
  L_phylo = L_phylo
)

species_lookup_use <- species_predictors_phylo %>%
  transmute(
    Species = Species,
    species_id = species_id,
    original_species_id = NA_integer_,
    ref_age_decades = ref_age_decades
  )

###################################################################################################
# step 16
# Save assurance input objects
# This file is needed by the Bayesian assurance array workflow.
###################################################################################################

analysis_data <- list(
  individual_data_phylo = individual_data_phylo,
  species_predictors_phylo = species_predictors_phylo,
  species_lookup_use = species_lookup_use,
  phylo_mode = phylo_mode,
  pagel_lambda = pagel_lambda
)

saveRDS(
  list(
    stan_data = stan_data_use,
    A_model = A_phylo,
    A_bm = A_bm,
    analysis_data = analysis_data
  ),
  assurance_input_file
)

cat("Assurance input file saved:", assurance_input_file, "\n")

###################################################################################################
# step 17
# Sanity checks
###################################################################################################

stopifnot(stan_data_use$N == length(stan_data_use$species_id))
stopifnot(stan_data_use$N == length(stan_data_use$age_decades))
stopifnot(stan_data_use$N == length(stan_data_use$event))
stopifnot(stan_data_use$N == length(stan_data_use$sex_male))
stopifnot(stan_data_use$S == length(stan_data_use$ref_age_decades))
stopifnot(stan_data_use$S == length(stan_data_use$log_body_mass_z))
stopifnot(stan_data_use$S == length(stan_data_use$log_longevity_z))
stopifnot(stan_data_use$S == length(stan_data_use$log_gestation_z))
stopifnot(all(dim(stan_data_use$L_phylo) == c(stan_data_use$S, stan_data_use$S)))
stopifnot(min(stan_data_use$species_id) == 1)
stopifnot(max(stan_data_use$species_id) == stan_data_use$S)
stopifnot(!anyNA(stan_data_use$age_decades))
stopifnot(!anyNA(stan_data_use$event))
stopifnot(!anyNA(stan_data_use$sex_male))
stopifnot(!anyNA(stan_data_use$L_phylo))

###################################################################################################
# step 18
# Stable initial values
###################################################################################################

init_fun <- function() {
  list(
    alpha_intercept = -4,
    alpha_body_mass = 0,
    alpha_longevity = 0,
    alpha_gestation = 0,
    sigma_phylo_alpha = 0.3,
    sigma_alpha_resid = 0.3,
    z_phylo_alpha = rnorm(stan_data_use$S, 0, 0.01),
    z_alpha_resid = rnorm(stan_data_use$S, 0, 0.01),
    beta_intercept = 0.5,
    beta_body_mass = 0,
    beta_longevity = 0,
    beta_gestation = 0,
    sigma_phylo_beta = 0.02,
    sigma_beta_resid = 0.02,
    z_phylo_beta = rnorm(stan_data_use$S, 0, 0.01),
    z_beta_resid = rnorm(stan_data_use$S, 0, 0.01),
    gamma_sex = 0
  )
}

###################################################################################################
# step 19
# Print run information
###################################################################################################

cat("Run mode:", run_mode, "\n")
cat("Response:", "malignancy", "\n")
cat("Phylogeny mode:", phylo_mode, "\n")
cat("Pagel lambda:", pagel_lambda, "\n")
cat("Animals:", stan_data_use$N, "\n")
cat("Species:", stan_data_use$S, "\n")
cat("Malignancy events:", sum(stan_data_use$event), "\n")
cat("Censored:", sum(1 - stan_data_use$event), "\n")
cat("Reference age range in decades:", min(stan_data_use$ref_age_decades), "to", max(stan_data_use$ref_age_decades), "\n")
cat("Chains:", chains, "\n")
cat("Parallel chains:", parallel_chains, "\n")
cat("Warmup:", iter_warmup, "\n")
cat("Sampling:", iter_sampling, "\n")
cat("adapt_delta:", adapt_delta, "\n")
cat("max_treedepth:", max_treedepth, "\n")

###################################################################################################
# step 20
# Compile and fit model
###################################################################################################

mod <- cmdstan_model(stan_file)

fit <- mod$sample(
  data = stan_data_use,
  seed = 123,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  init = init_fun,
  refresh = 50,
  output_dir = output_dir
)

###################################################################################################
# step 21
# Save fitted object
###################################################################################################

saveRDS(fit, fit_output_file)

###################################################################################################
# step 22
# Save summaries
###################################################################################################

fit_summary <- fit$summary()
write_csv(fit_summary, summary_output_file)

model_effects <- fit$summary(
  variables = c(
    "alpha_body_mass",
    "alpha_longevity",
    "alpha_gestation",
    "beta_body_mass",
    "beta_longevity",
    "beta_gestation",
    "gamma_sex",
    "sigma_phylo_alpha",
    "sigma_alpha_resid",
    "phylo_fraction_alpha",
    "sigma_phylo_beta",
    "sigma_beta_resid",
    "phylo_fraction_beta"
  )
)

write_csv(model_effects, effects_output_file)

alpha_summary <- fit$summary(variables = "alpha_mid") %>%
  mutate(
    species_id = as.integer(gsub("alpha_mid\\[|\\]", "", variable)),
    parameter = "alpha_mid"
  )

beta_summary <- fit$summary(variables = "beta") %>%
  mutate(
    species_id = as.integer(gsub("beta\\[|\\]", "", variable)),
    parameter = "beta"
  )

species_parameter_summary <- bind_rows(alpha_summary, beta_summary) %>%
  left_join(species_lookup_use, by = "species_id") %>%
  select(
    Species,
    species_id,
    original_species_id,
    ref_age_decades,
    parameter,
    mean,
    median,
    sd,
    q5,
    q95,
    rhat,
    ess_bulk,
    ess_tail
  )

write_csv(species_parameter_summary, species_output_file)

###################################################################################################
# step 23
# Save diagnostics
###################################################################################################

diagnostic_summary <- fit$diagnostic_summary()

# Handle any number of chains robustly.
diagnostic_table <- tibble(
  chain = seq_along(diagnostic_summary$num_divergent),
  num_divergent = diagnostic_summary$num_divergent,
  num_max_treedepth = diagnostic_summary$num_max_treedepth,
  ebfmi = diagnostic_summary$ebfmi
)

write_csv(diagnostic_table, file.path(output_dir, "07_diagnostic_summary.csv"))

sink(file.path(output_dir, "07_cmdstan_diagnose.txt"))
print(fit$cmdstan_diagnose())
sink()

###################################################################################################
# step 24
# Print final summaries
###################################################################################################

cat("Done. Model fitted and outputs saved.\n")
cat("Output directory:", output_dir, "\n")
cat("Fit object:", fit_output_file, "\n")
cat("Full summary:", summary_output_file, "\n")
cat("Effects summary:", effects_output_file, "\n")
cat("Species alpha_mid and beta summary:", species_output_file, "\n")
cat("Assurance input file:", assurance_input_file, "\n\n")

cat("Main parameter summary:\n")
print(
  fit$summary(
    variables = c(
      "alpha_intercept",
      "sigma_phylo_alpha",
      "sigma_alpha_resid",
      "phylo_fraction_alpha",
      "beta_intercept",
      "sigma_phylo_beta",
      "sigma_beta_resid",
      "phylo_fraction_beta",
      "gamma_sex"
    )
  )
)

cat("\nModel effect summary:\n")
print(model_effects)

cat("\nDiagnostic summary:\n")
print(diagnostic_summary)
