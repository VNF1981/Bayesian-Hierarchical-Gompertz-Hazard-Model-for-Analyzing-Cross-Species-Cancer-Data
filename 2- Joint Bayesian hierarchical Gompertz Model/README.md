# Gompertz Model for Age-Specific Malignancy Risk

This folder contains a pilot joint Bayesian hierarchical phylogenetic Gompertz model used to estimate age-specific malignancy risk across mammalian species.

The model starts from the prepared input object `gompertz_model_input.rds`, which was generated in Step 1 during the age and sex simulation pipeline. This object contains the individual-level simulated animal data, including age, sex, malignancy status, and species identity, together with species-level life-history predictors. Therefore, users do not need to rerun the age/sex simulation pipeline before fitting the Gompertz model.

## Input data

The main input file is:
```text
gompertz_model_input.rds
```

This file contains:
```text
individual_data
species_predictors
```

The individual-level data include:
```text
species identity
individual animal ID
age at malignancy diagnosis or censoring
sex
malignancy event status
right-censoring status
life-history variables merged to individuals
```

The species-level predictor table includes:
```text
adult body mass
maximum longevity
gestation length
standardized log-transformed life-history predictors
species IDs
```
### *** In the individual-level data, `malignancy = 1` indicates that malignancy was observed at the animal’s recorded age, which is treated as the event time. Animals with `malignancy = 0` had no observed malignancy by their recorded age and are therefore treated as right-censored at that age. ***

## Model overview

We modeled malignancy incidence using a Gompertz hazard function. For individual animal `j` from species `i`, the hazard of malignancy at age `t` is:
```text
h_ij(t) = α_mid,i × exp(β_i(t_ij − t_ref,i) + γ sex_ij)
```

where:

```text
i = species
j = individual animal
t_ij = age of animal j in species i, measured in decades
t_ref,i = species-specific reference age, measured in decades
α_mid,i = species-specific malignancy hazard at mid lifespan
β_i = species-specific age-related change in malignancy hazard
γ = sex effect on the log-hazard scale
sex_ij = 1 for male, 0 for female
```

Animals with malignancy are modeled as events. Animals without malignancy are treated as right-censored observations.

## Age and sex coding

Age was simulated and prepared in the age/sex simulation step before fitting this model. In the Gompertz model, age is converted from years to decades:

```text
age_decades = age_years / 10
```

Sex is encoded as:

```text
sex_male = 0 for female
sex_male = 1 for male
```

This coding is used directly in the hazard model through the term:

```text
γ sex_ij
```

## Age parameterization

The initial modeling goal was to estimate `α_i` as the species-specific malignancy hazard at age zero. However, this parameterization caused computational instability in Stan.

The issue was that malignancy risk is extremely low near age zero, especially for long-lived species. As a result, the baseline hazard at age zero and the Gompertz slope `β_i` became strongly coupled, producing poor posterior geometry and inefficient sampling.

To improve stability, we reparameterized the model around a biologically meaningful reference age: species-specific mid lifespan.

```text
t_ref,i = 0.5 × maximum longevity_i
```

Because maximum longevity was recorded in months, the reference age in decades was calculated as:

```text
ref_age_decades = max_longevity_M / 240
```

The final model therefore estimates:

```text
α_mid,i = malignancy hazard at species-specific mid lifespan
```

rather than malignancy hazard at age zero.

This keeps the Gompertz model structure intact while improving numerical stability and interpretability. The age-zero parameterization can be revisited later with stronger priors or alternative constraints.

## Species-level hierarchical model

The model estimates both species-specific Gompertz parameters jointly.

For malignancy hazard at mid lifespan:

```text
log(α_mid,i) =
α0
+ α_body_mass body_mass_i
+ α_longevity longevity_i
+ α_gestation gestation_i
+ phylogenetic residual_i
+ non-phylogenetic residual_i
```

For age-related change in malignancy hazard:

```text
β_i =
β0
+ β_body_mass body_mass_i
+ β_longevity longevity_i
+ β_gestation gestation_i
+ phylogenetic residual_i
+ non-phylogenetic residual_i
```

The species-level predictors are:

```text
adult body mass
maximum longevity
gestation length
```

Each predictor was log-transformed and standardized before modeling.

## Phylogenetic structure

The final model includes phylogenetic covariance in both:

```text
log(α_mid,i)
β_i
```

The tree file used for the mammalian analysis is:

```text
min20Fixed516.nwk
```

Species names were standardized before matching the data and tree by trimming whitespace and replacing spaces with underscores.

The phylogenetic model used the species present in both the individual-level dataset and the tree:

```text
327 species in the dataset
304 species in the tree
292 overlapping species retained
```

The final model used Pagel’s lambda:

```text
λ = 0.46
```

A Brownian motion option was retained in the code, but the default analysis used the Pagel lambda covariance matrix.

## Why a joint model?

We used a joint hierarchical model rather than a two-stage model.

In a two-stage approach, species-specific Gompertz parameters would first be estimated separately and then used as outcomes in a second phylogenetic regression. That approach is not ideal here because species differ substantially in sample size and malignancy event count. Some species have few or zero malignancy events, so their species-specific estimates would be highly uncertain.

The joint model estimates all components simultaneously:

```text
individual-level malignancy process
species-specific α_mid
species-specific β
life-history effects
phylogenetic effects
sex effect
```

This allows uncertainty in species-specific Gompertz parameters to propagate directly into the life-history and phylogenetic effects.

## Why a custom Stan model?

A custom Stan model was used because the final model required several features that are difficult to combine in standard survival packages:

```text
species-specific α_mid
species-specific β
β allowed to take positive, near-zero, or negative values
individual-level sex effect
right censoring
life-history predictors on both α_mid and β
phylogenetic covariance on both α_mid and β
species-specific reference ages
```

The Stan model directly implements the Gompertz hazard and cumulative hazard likelihood.

## Final full model run

The final full model was run on Monsoon using the phylogeny-matched mammal dataset:

```text
292 species
16,049 individual records
1,032 malignancy events
15,017 censored animals
```

Main MCMC settings:

```text
chains = 4
warmup = 1000
sampling = 1000
adapt_delta = 0.98
max_treedepth = 15
Pagel λ = 0.46
```

Model diagnostics from the full run were strong overall:

```text
all 4 chains completed successfully
0 maximum treedepth hits
good E-BFMI in all chains
Rhat values mostly near 1.00
1 divergent transition out of 4000 post-warmup draws
```

A more conservative rerun with `adapt_delta = 0.99` can be used if a fully divergence-free final run is required.

## Notes on model development

The final model was developed incrementally, with separate technical checks for data preparation, the Gompertz likelihood, hierarchical species effects, life-history predictors, and phylogenetic covariance.

This README describes only the final validated model structure used for analysis.

## Main files

```text
gompertz_model_input.rds
```

Prepared model input containing the individual-level malignancy dataset and species-level life-history predictors.

```text
min20Fixed516.nwk
```

Phylogenetic tree used for the mammalian analysis.

```text
07_gompertz_malignancy_alpha_mid_beta_LFH_phylo_alpha_beta.stan
```

Final Stan model.

```text
07_fit_gompertz_malignancy_alpha_mid_beta_LFH_phylo_alpha_beta_MONSOON_with_assurance_input.R
```

Script used to fit the final model on Monsoon and save the assurance input object.

For a portable GitHub version, users may replace the Monsoon-specific project path with:

```r
project_dir <- getwd()
```

## Main outputs

The final model writes outputs to:

```text
results_07_full_phylo_lambda/
```

Key output files include:

```text
07_gompertz_malignancy_phylo_alpha_beta_lambda_full_fit.rds
07_gompertz_malignancy_phylo_alpha_beta_lambda_full_summary.csv
07_gompertz_malignancy_phylo_alpha_beta_lambda_full_effects.csv
07_gompertz_malignancy_phylo_alpha_beta_lambda_full_species_summary.csv
07_gompertz_malignancy_phylo_alpha_beta_lambda_full_assurance_inputs.rds
07_diagnostic_summary.csv
07_cmdstan_diagnose.txt
species_missing_from_tree.csv
```

## Running the model

Place the following files in the same folder:

```text
gompertz_model_input.rds
min20Fixed516.nwk
07_gompertz_malignancy_alpha_mid_beta_LFH_phylo_alpha_beta.stan
07_fit_gompertz_malignancy_alpha_mid_beta_LFH_phylo_alpha_beta_MONSOON_with_assurance_input.R
```

Then run the fitting script in R:

```r
source("07_fit_gompertz_malignancy_alpha_mid_beta_LFH_phylo_alpha_beta_MONSOON_with_assurance_input.R")
```

For local or non-Monsoon use, update the project path and CmdStan path at the top of the R script.
