# Gompertz Alpha Mass Assurance Analysis

This folder contains the simulation-based assurance analysis for the final joint hierarchical phylogenetic Gompertz model.

I used the same general assurance logic as in the ZiBBPGLMM analysis: simulate datasets under known effect sizes, refit the full model, and estimate how often the model recovers support for the target effect.

The target parameter was `alpha_body_mass`, which is the body-mass effect on `log(α_mid)`, where `α_mid` is the species-specific malignancy hazard at mid lifespan.

## Simulation design

I tested five assumed positive body-mass effects:

```text
alpha_body_mass = 0.1, 0.3, 0.5, 0.7, 0.9
```

These correspond to approximate hazard ratios of:

```text
HR = 1.11, 1.35, 1.65, 2.01, 2.46
```

For each assumed effect size, I simulated and refit 50 datasets using the final joint hierarchical phylogenetic Gompertz model.

```text
5 effect sizes × 50 simulations = 250 total model refits
```

Each simulated dataset was generated using the fitted full Gompertz model as the baseline, while setting the body-mass effect to the target value.

## Success criterion

A simulation was counted as successful when the refitted model recovered posterior support for a positive body-mass effect:

```text
Pr(alpha_body_mass > 0) > 0.95
```

## Results

| Assumed `alpha_body_mass` | Approx. HR | Simulations | Assurance |
|---:|---:|---:|---:|
| 0.1 | 1.11 | 50 | 0.14 |
| 0.3 | 1.35 | 50 | 0.48 |
| 0.5 | 1.65 | 50 | 0.90 |
| 0.7 | 2.01 | 50 | 0.98 |
| 0.9 | 2.46 | 50 | 1.00 |

The assurance analysis suggests that the current data and model have low ability to reliably detect very small body-mass effects, moderate ability for an effect of `0.3`, and high ability for effects of `0.5` or larger.

All 250 model refits completed successfully. No fits hit the maximum treedepth limit, although some divergent transitions occurred, so these results should be interpreted as a pilot assurance analysis.

## Main files

```text
08_Gompertz_Assurance_body_mass_array_worker_CLEAN_FINAL.R
08_Combine_Gompertz_Assurance_body_mass_array_results_CLEAN_FINAL.R
07_gompertz_malignancy_phylo_alpha_beta_lambda_full_assurance_inputs.rds
```
The last one is the output from the main Gompertz model we performed in the previous step, 2- Joint Bayesian Hierarchical Gompertz Model

The combined outputs are written to:

```text
results_08_assurance_body_mass/
```

The main output files are:

```text
08_gompertz_assurance_body_mass_all_results.csv
08_gompertz_assurance_body_mass_summary.csv
```
