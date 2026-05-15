# Gompertz Model for Age-Specific Malignancy Risk

This folder contains a pilot Bayesian hierarchical phylogenetic Gompertz model used to estimate age-specific malignancy risk across mammalian species.

The model uses individual-level simulated records with species identity, age, sex, malignancy status, and species-level life-history predictors.

## Model overview

We modeled malignancy incidence using a Gompertz hazard function. For individual animal `j` from species `i`, the hazard of malignancy at age `t` is:

```text
h_ij(t) = α_mid,i × exp(β_i(t_ij − t_ref,i) + γ sex_ij)
```

where:

i = species
j = individual animal
t_ij = age of animal j in species i, measured in decades
t_ref,i = species-specific reference age, measured in decades
α_mid,i = species-specific malignancy hazard at mid lifespan
β_i = species-specific age-related change in malignancy hazard
γ = sex effect on the log-hazard scale
sex_ij = 1 for male, 0 for female

The cumulative hazard implied by this Gompertz model is used directly in the Stan likelihood, with animals without malignancy treated as right-censored observations.

Age parameterization

Our initial goal was to estimate α_i as the species-specific malignancy hazard at age zero. However, this parameterization was computationally unstable in Stan. Because malignancy is rare at very young ages, estimating the baseline hazard at age zero made α_i and β_i strongly coupled, especially for long-lived species and species with few malignancy events. This produced poor posterior geometry and inefficient sampling.

To improve stability, we reparameterized the model around a biologically meaningful reference age: species-specific mid lifespan.

t_ref,i = 0.5 × maximum longevity_i

Because maximum longevity was recorded in months, the reference age in decades was calculated as:

ref_age_decades = max_longevity_M / 240

The final model therefore estimates:

α_mid,i = malignancy hazard at species-specific mid lifespan

rather than the malignancy hazard at age zero.

This parameterization keeps the Gompertz model intact while improving numerical stability and interpretability. The original age-zero parameterization can be revisited later with additional regularization or alternative model constraints.

Species-level hierarchical model

The final model estimates both species-specific Gompertz parameters jointly.

For malignancy hazard at mid lifespan:

log(α_mid,i) =
α0
+ α_body_mass body_mass_i
+ α_longevity longevity_i
+ α_gestation gestation_i
+ phylogenetic residual_i
+ non-phylogenetic residual_i

For age-related change in malignancy hazard:

β_i =
β0
+ β_body_mass body_mass_i
+ β_longevity longevity_i
+ β_gestation gestation_i
+ phylogenetic residual_i
+ non-phylogenetic residual_i

The species-level predictors are:

adult body mass
maximum longevity
gestation length

Each predictor was log-transformed and standardized before modeling.

Phylogenetic structure

The final model includes phylogenetic covariance in both:

log(α_mid,i)
β_i

The phylogenetic tree used for the mammalian analysis is:

min20Fixed516.nwk

Species names were standardized before matching the data and tree by trimming whitespace and replacing spaces with underscores.

The phylogenetic model used the species present in both the individual-level dataset and the tree:

327 species in the dataset
304 species in the tree
292 overlapping species retained

The default phylogenetic model used Pagel’s λ:

λ = 0.46

A Brownian motion option was kept in the code, but the default analysis used the Pagel’s λ covariance matrix.

Why a joint model?

We used a joint hierarchical model rather than a two-stage model.

In a two-stage approach, species-specific Gompertz parameters would first be estimated separately and then used as outcomes in a second phylogenetic regression. That approach is not ideal here because species differ substantially in sample size and malignancy event count. Some species have few or zero malignancy events, so their species-specific estimates would be highly uncertain.

The joint model estimates all components simultaneously:

individual-level malignancy process
species-specific α_mid
species-specific β
life-history effects
phylogenetic effects
sex effect

This allows uncertainty in species-specific Gompertz parameters to propagate directly into the life-history and phylogenetic effects.

Why a custom Stan model?

A custom Stan model was used because the final model required features that are difficult to combine in standard survival packages:

species-specific α_mid
species-specific β
β allowed to take positive, near-zero, or negative values
individual-level sex effect
right censoring
life-history predictors on both α_mid and β
phylogenetic covariance on both α_mid and β
species-specific reference ages

The Stan model directly implements the Gompertz hazard and cumulative hazard likelihood.

Final full model run

The final full model was run on Monsoon using the phylogeny-matched mammal dataset:

292 species
16,049 individual records
1,032 malignancy events
15,017 censored animals

Main MCMC settings:

chains = 4
warmup = 1000
sampling = 1000
adapt_delta = 0.98
max_treedepth = 15
Pagel λ = 0.46

Model diagnostics from the full run were strong overall:

all 4 chains completed successfully
0 maximum treedepth hits
good E-BFMI in all chains
Rhat values mostly near 1.00
1 divergent transition out of 4000 post-warmup draws

A more conservative rerun with adapt_delta = 0.99 can be used if a fully divergence-free final run is required.

Notes on model development

The final model was developed incrementally, with separate technical checks for data preparation, the Gompertz likelihood, hierarchical species effects, life-history predictors, and phylogenetic covariance. This README describes only the final validated model structure used for analysis.

Main files
07_gompertz_malignancy_alpha_mid_beta_LFH_phylo_alpha_beta.stan

Final Stan model.

07_fit_gompertz_malignancy_alpha_mid_beta_LFH_phylo_alpha_beta_MONSOON.R

Monsoon script for fitting the final model.

min20Fixed516.nwk

Phylogenetic tree used for the mammalian analysis.

gompertz_model_input.rds

Prepared individual-level and species-level input data.

Main outputs

The final model writes outputs to:

results_07_full_phylo_lambda/

Key output files include:

07_gompertz_malignancy_phylo_alpha_beta_lambda_full_fit.rds
07_gompertz_malignancy_phylo_alpha_beta_lambda_full_summary.csv
07_gompertz_malignancy_phylo_alpha_beta_lambda_full_effects.csv
07_gompertz_malignancy_phylo_alpha_beta_lambda_full_species_summary.csv
07_gompertz_malignancy_phylo_alpha_beta_lambda_full_assurance_inputs.rds
07_diagnostic_summary.csv
07_cmdstan_diagnose.txt
species_missing_from_tree.csv
