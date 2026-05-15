# A Pilot Joint Bayesian Hierarchical Phylogenetic Gompertz Model for Cross-Species Age-Specific Malignancy Risk

This folder contains a pilot joint Bayesian hierarchical phylogenetic Gompertz model used to estimate age-specific malignancy risk across mammalian species.

The model starts from the prepared input object `gompertz_model_input.rds`, which was generated in Step 1 during the age and sex simulation pipeline. This object contains individual-level simulated animal data, including age, sex, malignancy status, and species identity, together with species-level life-history predictors. Therefore, users do not need to rerun the age/sex simulation pipeline before fitting the Gompertz model.

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
malignancy event indicator
right-censoring indicator
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

In the individual-level data, `malignancy = 1` indicates that malignancy was observed at the animal’s recorded age, which is treated as the event time. Animals with `malignancy = 0` had no observed malignancy by their recorded age and are therefore treated as right-censored at that age.

## Model overview

I modeled malignancy incidence using a Gompertz hazard function. For individual animal `j` from species `i`, the hazard of malignancy at age `t` is:

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

Age was modeled in decades rather than years or months to improve numerical stability in the Gompertz likelihood. Because the hazard changes exponentially with age, using decades keeps the age scale moderate and makes `β_i` interpretable as the species-specific age-related change in malignancy hazard per decade.

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

To improve stability, I reparameterized the model around a biologically meaningful reference age: species-specific mid lifespan.

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

The model estimates both species-specific Gompertz parameters jointly. I constrained `α_mid,i` to be positive because it represents a malignancy hazard rate at species-specific mid lifespan. In contrast, `β_i` was allowed to be negative, zero, or positive because it represents the age-related change in malignancy hazard: positive values indicate increasing hazard with age, values near zero indicate approximately constant hazard, and negative values indicate decreasing hazard with age.

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

This value is the lambda reported for the malignancy and body-mass association in the main comparative analysis. The value can be changed in the fitting script by editing:

```r
pagel_lambda <- 0.46
```

Briefly, I first construct the Brownian motion phylogenetic variance-covariance matrix from the tree. Pagel’s lambda is then applied by multiplying the off-diagonal elements of this matrix by `λ`, while keeping the diagonal elements equal to 1. This changes the expected covariance among species according to their shared evolutionary history, while preserving each species’ variance.

```r
A_phylo <- A_bm
A_phylo[off_diagonal] <- pagel_lambda * A_phylo[off_diagonal]
diag(A_phylo) <- 1
```

A Brownian motion (BM) option was also retained in the code. To use the BM phylogenetic covariance structure instead of the Pagel lambda-transformed covariance matrix, change:

```r
phylo_mode <- "lambda"
```

to:

```r
phylo_mode <- "BM"
```

## Why a joint model?

I used a joint hierarchical model rather than a two-stage model. In a two-stage approach, species-specific Gompertz parameters would first be estimated separately and then used as outcomes in a second phylogenetic regression. I decided that this approach would not be ideal here because species differ substantially in within-species sample size (i.e., the number of individuals per species) and malignancy event count. Some species have few or zero malignancy events, so their species-specific estimates would be highly uncertain. Thus, it would be difficult to fully propagate uncertainty from the first-stage estimates into the second-stage phylogenetic regression. For this reason, I used a joint hierarchical model that estimates the individual-level Gompertz process, species-level parameters, life-history effects, and phylogenetic effects simultaneously.

## Why a custom Stan model?

I used a custom Stan model because the final model structure is not directly supported as a standard off-the-shelf model in common survival modeling packages. The model required a joint Gompertz likelihood with species-specific `α_mid` and `β`, where `β` can take positive, zero, or negative values; individual-level sex effects; right censoring; life-history predictors on both Gompertz parameters; and phylogenetic covariance on both `α_mid` and `β`. While some packages can support parts of this structure, implementing the full model in a transparent and flexible way was more straightforward in custom Stan code. In future extensions, I may also incorporate left truncation/delayed entry if needed.

The Stan model directly implements the Gompertz hazard and cumulative hazard likelihood.

## Model development

The final model was developed incrementally with AI assistance by adding one component at a time and checking that each version worked correctly before adding the next layer of complexity. These technical checks covered data preparation, the Gompertz likelihood, hierarchical species effects, life-history predictors, and phylogenetic covariance.

This README describes only the final validated model structure used for analysis.

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

I used a relatively high `adapt_delta` value to make Stan sample more conservatively and reduce the risk of divergent transitions in this complex hierarchical phylogenetic model.

I also used a relatively high `max_treedepth` value to allow Stan enough trajectory length to explore the complex posterior distribution without frequently hitting the maximum treedepth limit.

Model diagnostics from the full run were strong overall:

```text
all 4 chains completed successfully
0 maximum treedepth hits
good E-BFMI in all chains
Rhat values mostly near 1.00
1 divergent transition out of 4000 post-warmup draws
```

A more conservative rerun with `adapt_delta = 0.99` can be used if a fully divergence-free final run is required.

As a reference, the full model completed in approximately 34 minutes on the Monsoon cluster using 1 node, 1 task, 8 CPUs, and 128 GB RAM. Runtime may vary across systems.

## Body-mass assurance analysis

I also ran a simulation-based assurance analysis to evaluate whether the final joint Gompertz model could detect positive body-mass effects on malignancy hazard at mid lifespan. The target parameter was `alpha_body_mass`, the body-mass effect on `log(α_mid)`.

The analysis tested five assumed positive effect sizes:

```text
alpha_body_mass = 0.1, 0.3, 0.5, 0.7, 0.9
```

These correspond to hazard ratios of approximately:

```text
1.11, 1.35, 1.65, 2.01, 2.46
```

For each assumed effect size, I simulated and refit 50 datasets using the same joint hierarchical phylogenetic Gompertz model structure.

The assurance results were:

| Assumed `alpha_body_mass` | Approx. HR | Number of simulations | Assurance |
|---:|---:|---:|---:|
| 0.1 | 1.11 | 50 | 0.14 |
| 0.3 | 1.35 | 50 | 0.48 |
| 0.5 | 1.65 | 50 | 0.90 |
| 0.7 | 2.01 | 50 | 0.98 |
| 0.9 | 2.46 | 50 | 1.00 |

Here, assurance was defined as the proportion of simulations in which the fitted model recovered evidence for a positive body-mass effect. Specifically, success was evaluated using posterior support for a positive `alpha_body_mass` effect.

The results suggest that the current data and model have low assurance for detecting very small body-mass effects, moderate assurance for an effect of 0.3, and high assurance for effects of 0.5 or larger.

All 250 assurance model fits completed successfully, and no fits hit the maximum treedepth limit. Some divergent transitions occurred across the assurance fits, so these results should be interpreted as a pilot assurance analysis rather than a final power analysis.

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

The assurance analysis writes outputs to:

```text
results_08_assurance_body_mass/
```

Key assurance output files include:

```text
08_gompertz_assurance_body_mass_all_results.csv
08_gompertz_assurance_body_mass_summary.csv
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
