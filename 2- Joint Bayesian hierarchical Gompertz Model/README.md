# A Pilot Joint Hierarchical Phylogenetic Gompertz Model for Cross-Species Age-Specific Malignancy Risk

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
standardized log-transformed of these life-history predictors
species IDs
```
## Model overview

We modeled malignancy incidence using a Gompertz hazard function. For individual animal `j` from species `i`, the hazard of malignancy at age `t` is:
```text
h_ij(t) = α_mid,i × exp(β_i(t_ij − t_ref,i) + γ sex_ij)
```

where:

```text
i = species
j = individual animal
t_ij = age of animal j in species i, measured in decades (see below)
t_ref,i = species-specific reference age, measured in decades
α_mid,i = species-specific malignancy hazard at mid lifespan (see explanation below)
β_i = species-specific age-related change in malignancy hazard
γ = sex effect on the log-hazard scale
sex_ij = 1 for male, 0 for female
```

*** Animals with malignancy are modeled as events. Animals without malignancy are treated as right-censored observations. In other words, in the individual-level data, `malignancy = 1` indicates that malignancy was observed at the animal’s recorded age, which is treated as the event time. Animals with `malignancy = 0` had no observed malignancy by their recorded age and are therefore treated as right-censored at that age. ***

*** Age was modeled in decades rather than years or months to improve numerical stability in the Gompertz likelihood. Because the hazard changes exponentially with age, using decades keeps the age scale moderate and makes `β_i` interpretable as the species-specific age-related change in malignancy hazard per decade. ***

If age is measured in months, then β becomes very small because each unit is only one month. If age is measured in years, it is better, but some species still have long lifespans, which can make the exponential term harder for Stan to sample efficiently.

## Age parameterization

The initial modeling goal was to estimate `α_i` as the species-specific malignancy hazard at age zero. However, this parameterization caused computational instability in Stan. The issue was that malignancy risk is extremely low near age zero, especially for long-lived species. As a result, the baseline hazard at age zero and the Gompertz slope `β_i` became strongly coupled, producing poor posterior geometry and inefficient sampling.
To improve stability, I reparameterized the model around a biologically meaningful reference age: species-specific mid lifespan. I will work on this to see how I can include `α_i = malignancy at age zero` for our final model 
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

Species differ substantially in the number of individual animals available for analysis. Some species have many individuals and many malignancy observations, while others have few individuals or few observed malignancy events. The joint hierarchical model addresses this by estimating species-specific `α_mid` and `β` with partial pooling across species. In practice, species with limited data borrow more information from the hierarchical, life-history, and phylogenetic structure of the model, while species with more data have estimates driven more strongly by their own individual-level observations. This helps account for heterogeneous uncertainty due to unequal within-species sample sizes, instead of treating all species-level estimates as equally precise.

The model estimates both species-specific Gompertz parameters jointly. I designed the model so that `α_mid,i` is always positive because it represents a malignancy hazard rate at age zero (here is mid lifespan instead), and hazard rates cannot be negative. In contrast, `β_i` was allowed to take negative, zero, or positive values because it represents how malignancy hazard changes with age. A positive `β_i` indicates increasing hazard with age, a value near zero indicates approximately constant hazard with age, and a negative `β_i` indicates decreasing hazard with age (i.e., for exceptionally cancer resistant species).

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
Ths is the λ reported fro neoplasia and body mass in the big paper (I used a wrong one! It must be for malignancy and mass). The value can be changed in the fitting script by editing:
```text
pagel_lambda <- 0.46
```

A Brownian motion (BM) option was retained in the code. You can change the model to BM by changing this line:
```text
phylo_mode <- "lambda"
to
phylo_mode <- "BM"
```

Briefly, we first construct the Brownian motion phylogenetic variance-covariance matrix from the tree. Pagel’s lambda is then applied by multiplying the off-diagonal elements of this matrix (covariances) by `λ`, while keeping the diagonal elements (variances) equal to 1. This reduces or increases the expected covariance among species according to their shared evolutionary history, while preserving each species’ variance.


## Why a joint model?

I used a joint hierarchical model rather than a two-stage model. In a two-stage approach, species-specific Gompertz parameters would first be estimated separately and then used as outcomes in a second phylogenetic regression. I decided that this approach would not be ideal here because species differ substantially in within-species sample size (i.e., the number of individuals per species) and malignancy event count. Some species have few or zero malignancy events, so their species-specific estimates would be highly uncertain. Thus, it would be difficult to fully propagate uncertainty from the first-stage estimates into the second-stage phylogenetic regression. For this reason, I used a joint hierarchical model that estimates the individual-level Gompertz process, species-level parameters, life-history effects, and phylogenetic effects simultaneously.

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

I used a custom Stan model because the final model structure is not directly supported as a standard off-the-shelf model in common survival modeling packages. The model required a joint Gompertz likelihood with species-specific `α_mid` and `β`, where `β` can take positive, zero, or negative values; individual-level sex effects; right censoring (I need to incorporate left truncation/delayed entries as well); life-history predictors on both Gompertz parameters; and phylogenetic covariance on both `α_mid` and `β`. While some packages can support parts of this structure, implementing the full model in a transparent and flexible way was more straightforward in custom Stan code. 

*** The final model was developed incrementally with AI assistance by adding one component at a time and checking that each version worked correctly before adding the next layer of complexity, including covered data preparation, the Gompertz likelihood, hierarchical species effects, life-history predictors, and phylogenetic covariance. ***

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

## Sex effect in the pilot model

Sex was included as an individual-level covariate in the Gompertz hazard model:
```text
h_ij(t) = α_mid,i × exp(β_i(t_ij − t_ref,i) + γ sex_ij)
```

where:
```text
sex_ij = 0 for females
sex_ij = 1 for males
γ = global sex effect on the log-hazard scale
```

In this parameterization, `γ` is estimated across all species while accounting for age, species identity, life-history predictors, phylogeny, and censoring.

The male-to-female hazard ratio is:
```text
exp(γ)
```

Therefore:
```text
γ > 0  means males have higher malignancy hazard than females
γ < 0  means males have lower malignancy hazard than females
γ ≈ 0  means there is little evidence for an overall sex difference
```

In the current pilot version, sex is modeled as a single global effect. The model does not estimate separate male- and female-specific Gompertz parameters for each species. If the full analysis shows evidence of meaningful sex differences and the data support additional complexity, a future model extension could allow sex effects to vary across species or clades. In the current pilot analysis with simulated individual-level age and sex data, the estimated global sex effect was close to zero, meaning there was no clear evidence for an overall difference in malignancy hazard between males and females in the pilot analysis. However, this result is based on simulated data and should not be interpreted as a biological conclusion about sex-specific cancer risk.


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

adapt_delta is a Stan/HMC tuning parameter that controls how carefully the sampler avoids numerical problems during sampling. Higher adapt_delta values usually mean fewer divergent transitions, more stable sampling, and slower runtime. max_treedepth is another Stan/HMC tuning parameter. It limits how long Stan can explore the posterior during a single MCMC transition.

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

Then run the fitting script in R. 

## *** As a reference, the full model completed in less than one hour on the Monsoon cluster using 1 node, 1 task, 8 CPUs, and 128 GB RAM. ***
