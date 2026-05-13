######################################################################################################################################
######################## This code is designed for simulating age and sex for animals in our ACE dataset  
######################## The output is an Excel sheet including benign, neoplasia, malignancy, three life history
######################## traits, plus age and sex for various number of animals in each species, 
######################## for 327 speices in ACE dataset
######################################################################################################################################
rm(list = ls(all = TRUE))
ls()

# step 1
# Load packages
# Install these first if needed:
# install.packages(c("readxl", "dplyr", "tidyr", "writexl"))
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

# step 2
# Set input and output file names
input_file  <- "Compton_data.xlsx"
output_file <- "ACE_individual_level_simulated_age_sex.xlsx"

# step 3
# Read the ACE species level dataset
ace_species <- read_excel(input_file)

# step 4
# Column names expected by this script
# If your file uses slightly different names, change only this section
species_col    <- "Species"
class_col      <- "Class"
clade_col      <- "Clade"
order_col      <- "Orders"
family_col     <- "Family"
genus_col      <- "Genus"
necropsy_col   <- "Necropsies"
benign_col     <- "Benign"
neoplasia_col  <- "Neoplasia"
malignant_col  <- "Malignant"
body_mass_col  <- "adult_weight_G"
gestation_col  <- "Gestation_M"
longevity_col  <- "max_longevity_M"

# step 5
# Check that all required columns exist
required_cols <- c(
  species_col,
  class_col,
  clade_col,
  order_col,
  family_col,
  genus_col,
  necropsy_col,
  benign_col,
  neoplasia_col,
  malignant_col,
  body_mass_col,
  gestation_col,
  longevity_col
)

missing_cols <- setdiff(required_cols, names(ace_species))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "These required columns are missing from the input file: ",
      paste(missing_cols, collapse = ", "),
      "\nPlease edit the column name section at the top of the script."
    )
  )
}

# step 6
# Set random seed so the simulation can be repeated
set.seed(123)

# step 7
# Clean species level data
# Values of -1 and 0 are treated as missing for life history traits
ace_clean <- ace_species %>%
  mutate(
    across(
      all_of(c(necropsy_col, benign_col, neoplasia_col, malignant_col,
               body_mass_col, gestation_col, longevity_col)),
      as.numeric
    ),
    across(
      all_of(c(body_mass_col, gestation_col, longevity_col)),
      ~ ifelse(.x <= 0, NA_real_, .x)
    ),
    across(
      all_of(c(necropsy_col, benign_col, neoplasia_col, malignant_col)),
      ~ ifelse(is.na(.x), 0, .x)
    ),
    !!necropsy_col := round(.data[[necropsy_col]]),
    !!benign_col := round(.data[[benign_col]]),
    !!neoplasia_col := round(.data[[neoplasia_col]]),
    !!malignant_col := round(.data[[malignant_col]])
  ) %>%
  filter(.data[[necropsy_col]] > 0)

# step 8
# Helper function for filling missing numeric values by taxonomic hierarchy
fill_by_taxonomy <- function(data, value_col, hierarchy_cols) {
  filled_col <- paste0(value_col, "_filled")
  source_col <- paste0(value_col, "_fill_source")
  out <- data
  out[[filled_col]] <- out[[value_col]]
  out[[source_col]] <- ifelse(is.na(out[[value_col]]), NA_character_, "observed")
  
  for (tax_col in hierarchy_cols) {
    medians <- out %>%
      group_by(.data[[tax_col]]) %>%
      summarise(group_median = median(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
      mutate(group_median = ifelse(is.infinite(group_median), NA_real_, group_median))
    
    out <- out %>%
      left_join(medians, by = setNames(tax_col, tax_col)) %>%
      mutate(
        !!filled_col := ifelse(is.na(.data[[filled_col]]) & !is.na(group_median), group_median, .data[[filled_col]]),
        !!source_col := ifelse(is.na(.data[[source_col]]) & !is.na(group_median), tax_col, .data[[source_col]])
      ) %>%
      select(-group_median)
  }
  
  global_median <- median(out[[value_col]], na.rm = TRUE)
  
  if (is.na(global_median) || is.infinite(global_median)) {
    stop(paste0("No valid values found for ", value_col, ". Cannot fill missing values."))
  }
  
  out <- out %>%
    mutate(
      !!filled_col := ifelse(is.na(.data[[filled_col]]), global_median, .data[[filled_col]]),
      !!source_col := ifelse(is.na(.data[[source_col]]), "global_median", .data[[source_col]])
    )
  
  return(out)
}

# step 9
# Fill missing longevity using the most specific available taxonomic levels first
# The order here means we first try genus, then family, then order, then clade, then class, then global median
taxonomy_hierarchy <- c(genus_col, family_col, order_col, clade_col, class_col)

ace_filled <- ace_clean %>%
  fill_by_taxonomy(longevity_col, taxonomy_hierarchy) %>%
  fill_by_taxonomy(body_mass_col, taxonomy_hierarchy) %>%
  fill_by_taxonomy(gestation_col, taxonomy_hierarchy)

# step 10
# Make sure malignant case counts do not exceed necropsy counts
# This protects the simulation from impossible records
ace_filled <- ace_filled %>%
  mutate(
    !!malignant_col := pmin(.data[[malignant_col]], .data[[necropsy_col]]),
    !!benign_col := pmin(.data[[benign_col]], .data[[necropsy_col]]),
    !!neoplasia_col := pmin(.data[[neoplasia_col]], .data[[necropsy_col]])
  )

# step 11
# Function to simulate individual animals for one species
simulate_species_animals <- function(row_data) {
  n_animals <- as.integer(row_data[[necropsy_col]])
  n_malignant <- as.integer(row_data[[malignant_col]])
  
  malignancy_status <- c(rep(1, n_malignant), rep(0, n_animals - n_malignant))
  malignancy_status <- sample(malignancy_status, size = n_animals, replace = FALSE)
  
  sex <- ifelse(runif(n_animals) <= 0.5, "F", "M")
  
  max_life <- as.numeric(row_data[[paste0(longevity_col, "_filled")]])
  
  age_fraction <- ifelse(
    malignancy_status == 1,
    rbeta(n_animals, shape1 = 4.0, shape2 = 1.8),
    rbeta(n_animals, shape1 = 2.0, shape2 = 2.5)
  )
  
  age_months <- age_fraction * max_life
  age_years <- age_months / 12
  
  tibble(
    Species = row_data[[species_col]],
    animal_id_within_species = seq_len(n_animals),
    individual_id = paste0(row_data[[species_col]], "_", seq_len(n_animals)),
    Class = row_data[[class_col]],
    Clade = row_data[[clade_col]],
    Orders = row_data[[order_col]],
    Family = row_data[[family_col]],
    Genus = row_data[[genus_col]],
    sex = sex,
    age_months = age_months,
    age_years = age_years,
    malignancy = malignancy_status,
    benign_species_count = row_data[[benign_col]],
    neoplasia_species_count = row_data[[neoplasia_col]],
    malignant_species_count = row_data[[malignant_col]],
    necropsies_species_count = row_data[[necropsy_col]],
    adult_weight_G = row_data[[paste0(body_mass_col, "_filled")]],
    Gestation_M = row_data[[paste0(gestation_col, "_filled")]],
    max_longevity_M = row_data[[paste0(longevity_col, "_filled")]],
    adult_weight_G_fill_source = row_data[[paste0(body_mass_col, "_fill_source")]],
    Gestation_M_fill_source = row_data[[paste0(gestation_col, "_fill_source")]],
    max_longevity_M_fill_source = row_data[[paste0(longevity_col, "_fill_source")]]
  )
}

# step 12
# Create individual level dataset
individual_data <- bind_rows(
  lapply(seq_len(nrow(ace_filled)), function(i) {
    simulate_species_animals(ace_filled[i, ])
  })
)

# step 13
# Create a species level summary to check that the simulation preserved the malignancy counts
simulation_check <- individual_data %>%
  group_by(Species) %>%
  summarise(
    simulated_animals = n(),
    simulated_malignant = sum(malignancy),
    simulated_female = sum(sex == "F"),
    simulated_male = sum(sex == "M"),
    mean_age_months = mean(age_months),
    median_age_months = median(age_months),
    min_age_months = min(age_months),
    max_age_months = max(age_months),
    .groups = "drop"
  ) %>%
  left_join(
    ace_filled %>%
      select(
        Species = all_of(species_col),
        original_necropsies = all_of(necropsy_col),
        original_malignant = all_of(malignant_col),
        original_benign = all_of(benign_col),
        original_neoplasia = all_of(neoplasia_col),
        max_longevity_M_filled = all_of(paste0(longevity_col, "_filled")),
        max_longevity_M_fill_source = all_of(paste0(longevity_col, "_fill_source"))
      ),
    by = "Species"
  ) %>%
  mutate(
    necropsy_count_match = simulated_animals == original_necropsies,
    malignancy_count_match = simulated_malignant == original_malignant
  )

# step 14
# Create an imputed species level table for record keeping
species_level_imputed <- ace_filled

# step 15
# Write output Excel file
# Only the individual level simulated sheet is saved because this is the sheet needed for the Gompertzian model
write_xlsx(
  list(
    individual_level_simulated = individual_data
  ),
  output_file
)

# step 16
# Print basic checks
cat("Done. Output file created:", output_file, "\n")
cat("Number of species:", n_distinct(individual_data$Species), "\n")
cat("Number of simulated animals:", nrow(individual_data), "\n")
cat("Number of malignant animals:", sum(individual_data$malignancy), "\n")
cat("All necropsy counts preserved:", all(simulation_check$necropsy_count_match), "\n")
cat("All malignancy counts preserved:", all(simulation_check$malignancy_count_match), "\n")
