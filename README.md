# Bayesian-Hierarchical-Gompertz-Hazard-Model-for-Analyzing-Cross-Species-Cancer-Data
This repository includes R code for simulating age and sex data for species, fitting a pilot Bayesian hierarchical Gompertz model for malignancy, and developing an assurance analysis framework.

# ACE Individual Level Age and Sex Simulation

This folder contains the first step of the Gompertzian modeling workflow for the ACE comparative oncology dataset.

The goal of this step is to convert the original species level ACE dataset into a simulated individual level dataset. The resulting file includes one row per animal, with simulated age and sex, while preserving the original species level counts for malignancy and necropsy records.

This simulated dataset will later be used to develop and test a Bayesian hierarchical Gompertzian hazard model for malignancy.

## Purpose

The available ACE dataset is organized at the species level. For each species, the dataset includes information such as the number of necropsies, number of benign tumors, number of neoplasia cases, number of malignancy cases, and several life history traits.

However, a Gompertzian hazard model requires individual level information, especially age at diagnosis or age at censoring. Because true age and sex resolved individual records are not currently available, this script simulates age and sex for each animal.

This allows us to build and test the modeling framework before applying it to future data with real age, sex, diagnosis, and denominator information.

## Input file

The script uses the following input file:

```text
Compton_data_plus_simulated_Age_Sex.xlsx
```

This file contains the ACE species level dataset.

The script expects the dataset to include the following types of information:

| Variable type | Example column |
|---|---|
| Species name | `Species` |
| Taxonomic class | `Class` |
| Clade | `Clade` |
| Taxonomic order | `Orders` |
| Family | `Family` |
| Genus | `Genus` |
| Number of necropsies | `Necropsies` |
| Benign tumor count | `Benign` |
| Neoplasia count | `Neoplasia` |
| Malignancy count | `Malignant` |
| Adult body mass | `adult_weight_G` |
| Gestation length | `Gestation_M` |
| Maximum longevity | `max_longevity_M` |

Note that the taxonomic order column in this dataset is named `Orders`.

## Output file

The script creates the following Excel file:

```text
ACE_individual_level_simulated_age_sex.xlsx
```

The output file contains three sheets.

## Output sheets

### 1. `individual_level_simulated`

This is the main output sheet.

It contains one row per simulated animal.

For each species, the number of simulated animals is equal to the number of necropsies in the original ACE dataset.

For example, if a species has 100 necropsies, the output will contain 100 rows for that species.

The malignancy status is also preserved. For example, if a species has 12 malignant cases, then exactly 12 of the 100 simulated animals are assigned malignancy status 1, and the remaining 88 are assigned malignancy status 0.

Main columns include:

| Column | Meaning |
|---|---|
| `Species` | species name |
| `animal_id_within_species` | animal ID within each species |
| `individual_id` | unique simulated individual ID |
| `Class` | taxonomic class |
| `Clade` | clade |
| `Orders` | taxonomic order |
| `Family` | family |
| `Genus` | genus |
| `sex` | simulated sex |
| `age_months` | simulated age in months |
| `age_years` | simulated age in years |
| `malignancy` | simulated individual malignancy status |
| `benign_species_count` | original benign count for the species |
| `neoplasia_species_count` | original neoplasia count for the species |
| `malignant_species_count` | original malignant count for the species |
| `necropsies_species_count` | original necropsy count for the species |
| `adult_weight_G` | observed or filled adult body mass |
| `Gestation_M` | observed or filled gestation length |
| `max_longevity_M` | observed or filled maximum longevity |

### 2. `species_level_imputed`

This sheet contains the species level dataset after cleaning and filling missing life history values.

Missing or invalid values are filled for:

```text
adult_weight_G
Gestation_M
max_longevity_M
```

Invalid values such as 0 or negative values are treated as missing.

This sheet is useful for checking which life history values were observed and which values were filled using taxonomic information.

### 3. `simulation_check`

This sheet checks whether the simulation preserved the original species level information.

It includes:

| Column | Meaning |
|---|---|
| `simulated_animals` | number of simulated individual rows for each species |
| `original_necropsies` | original number of necropsies |
| `simulated_malignant` | number of simulated malignant animals |
| `original_malignant` | original number of malignant cases |
| `necropsy_count_match` | whether simulated animal count matches original necropsy count |
| `malignancy_count_match` | whether simulated malignancy count matches original malignancy count |

Both match columns should be `TRUE` for all species.

## Main assumptions

### One animal per necropsy

The script assumes that the number of individual animals per species is equal to the number of necropsies.

This is a practical starting assumption for developing the model. In the future, this can be replaced with real denominator information from husbandry records or animal years at risk.

### Malignancy count is preserved

For each species, the total number of malignant animals in the simulated individual level dataset is exactly equal to the original species level malignancy count.

This ensures that the simulation does not change the observed malignancy burden in the ACE dataset.

### Sex is simulated using a 50 to 50 probability

Sex is simulated using a simple random draw from a uniform distribution.

The rule is:

```r
sex <- ifelse(runif(n_animals) <= 0.5, "F", "M")
```

This assumes an equal probability of male and female animals.

This is a placeholder assumption and can be replaced later if real sex ratios become available.

### Age is simulated using maximum longevity

Age is simulated using the species maximum longevity value.

Because cancer risk usually increases with age, malignant and non malignant animals are simulated from different age distributions.

Non malignant animals are simulated with a broader age distribution:

```r
rbeta(n_animals, shape1 = 2.0, shape2 = 2.5)
```

Malignant animals are simulated with an older shifted age distribution:

```r
rbeta(n_animals, shape1 = 4.0, shape2 = 1.8)
```

The simulated age fraction is then multiplied by the filled species maximum longevity:

```r
age_months = age_fraction * max_longevity_M
```

This means that malignant animals tend to be older on average than non malignant animals, which is biologically consistent with the general expectation that malignancy risk increases with age.

## Filling missing life history values

Some species have missing or invalid life history values.

The script treats values of 0 or negative values as missing for the life history traits.

Missing values are filled using a taxonomic hierarchy from the most specific level to the broadest level.

The hierarchy is:

```text
Genus
Family
Orders
Clade
Class
Global median
```

For example, if a species is missing maximum longevity, the script first tries to use the median longevity of species in the same genus. If that is not available, it tries family, then order, then clade, then class, and finally the global median.

This approach allows the script to keep all species in the dataset while using the closest available taxonomic information to fill missing values.

## Reproducibility

The script uses a fixed random seed:

```r
set.seed(123)
```

This means that the same input file should produce the same simulated ages and sexes each time the script is run.

## Required R packages

The script uses the following R packages:

```r
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
```

If needed, they can be installed using:

```r
install.packages(c("readxl", "dplyr", "tidyr", "writexl"))
```

## How to run the script

Place the R script and the input Excel file in the same folder.

Then run:

```r
source("simulate_ACE_age_sex.R")
```

or run from the terminal using:

```bash
Rscript simulate_ACE_age_sex.R
```

After the script finishes, it will create:

```text
ACE_individual_level_simulated_age_sex.xlsx
```

## Current stage of the workflow

This script only performs the data simulation step.

It does not fit the Gompertzian model yet.

The next step will be to use the simulated individual level dataset to fit a Bayesian hierarchical Gompertzian hazard model for malignancy.

## Planned next step

The next script will use the individual level dataset generated here to model malignancy risk as a function of age, sex, species level life history traits, and species level random effects.

The initial model will focus on malignancy only.

Later, the same framework can be extended to:

```text
Benign neoplasia
Overall neoplasia
Cancer mortality
```

## Notes

This simulation is intended as a prototype for model development and assurance analysis.

The simulated age and sex values should not be interpreted as real observations.

The purpose is to build a working statistical framework that can later be applied to real age and sex resolved data when those data become available.
