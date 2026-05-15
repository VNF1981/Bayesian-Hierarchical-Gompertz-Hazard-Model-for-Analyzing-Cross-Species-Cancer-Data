###################################################################################################
# 08_Combine_Gompertz_Assurance_body_mass_array_results.R
#
# CLEAN FINAL combine script for the Gompertz body-mass assurance run.
#
# This script matches:
#   08_Gompertz_Assurance_body_mass_array_worker.R
#
# Expected input directory:
#   results_08_assurance_body_mass
#
# Expected input files:
#   assurance_body_mass_alpha_0.1_sim_001_task_001.csv
#
# Outputs:
#   08_gompertz_assurance_body_mass_all_results.csv
#   08_gompertz_assurance_body_mass_summary.csv
###################################################################################################

rm(list = ls(all = TRUE))
ls()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

project_dir <- "/projects/tollis_lab/TE/mammals/Gompertz_model"
setwd(project_dir)

assurance_out_dir <- file.path(project_dir, "results_08_assurance_body_mass")

if (!dir.exists(assurance_out_dir)) {
  stop("Could not find assurance output directory: ", assurance_out_dir)
}

result_files <- list.files(
  assurance_out_dir,
  pattern = "^assurance_body_mass_alpha_[0-9]+\\.[0-9]_sim_[0-9]+_task_[0-9]+\\.csv$",
  full.names = TRUE
)

cat("Assurance output directory:", assurance_out_dir, "\n")
cat("Found result files:", length(result_files), "\n")

if (length(result_files) == 0) {
  cat("First 30 files in directory:\n")
  print(head(list.files(assurance_out_dir), 30))
  stop("No matching assurance CSV result files found.")
}

all_results <- result_files %>%
  lapply(readr::read_csv, show_col_types = FALSE) %>%
  dplyr::bind_rows() %>%
  arrange(assumed_alpha_body_mass, sim_i, task_id)

cat("\nCounts by assumed_alpha_body_mass:\n")
print(table(all_results$assumed_alpha_body_mass))

cat("\nCounts by effect_index:\n")
print(table(all_results$effect_index))

cat("\nTask ID range:\n")
print(range(all_results$task_id, na.rm = TRUE))

cat("\nTotal result rows:", nrow(all_results), "\n")

summary_results <- all_results %>%
  group_by(target_parameter, direction, assumed_alpha_body_mass, assumed_HR) %>%
  summarise(
    n_sim_requested_or_found = n(),
    n_fit_success = sum(fit_success == TRUE, na.rm = TRUE),

    n_valid_posterior_prob = sum(!is.na(success_posterior_prob)),
    n_success_posterior_prob = sum(success_posterior_prob == TRUE, na.rm = TRUE),
    assurance_posterior_prob = mean(success_posterior_prob, na.rm = TRUE),

    n_valid_q5 = sum(!is.na(success_q5_positive)),
    n_success_q5 = sum(success_q5_positive == TRUE, na.rm = TRUE),
    assurance_q5_positive = mean(success_q5_positive, na.rm = TRUE),

    median_estimated_alpha_body_mass = median(median_estimated_alpha_body_mass, na.rm = TRUE),
    median_q5_estimated_alpha_body_mass = median(q5_estimated_alpha_body_mass, na.rm = TRUE),
    median_q95_estimated_alpha_body_mass = median(q95_estimated_alpha_body_mass, na.rm = TRUE),

    median_posterior_prob_positive = median(posterior_prob_positive, na.rm = TRUE),

    mean_events_simulated = mean(n_events_simulated, na.rm = TRUE),
    min_events_simulated = min(n_events_simulated, na.rm = TRUE),
    max_events_simulated = max(n_events_simulated, na.rm = TRUE),

    total_divergent = sum(num_divergent, na.rm = TRUE),
    total_max_treedepth = sum(num_max_treedepth, na.rm = TRUE),

    median_rhat_alpha_body_mass = median(rhat_alpha_body_mass, na.rm = TRUE),
    max_rhat_alpha_body_mass = max(rhat_alpha_body_mass, na.rm = TRUE),
    median_ess_bulk_alpha_body_mass = median(ess_bulk_alpha_body_mass, na.rm = TRUE),
    min_ess_bulk_alpha_body_mass = min(ess_bulk_alpha_body_mass, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(assumed_alpha_body_mass)

raw_file <- file.path(
  assurance_out_dir,
  "08_gompertz_assurance_body_mass_all_results.csv"
)

summary_file <- file.path(
  assurance_out_dir,
  "08_gompertz_assurance_body_mass_summary.csv"
)

readr::write_csv(all_results, raw_file)
readr::write_csv(summary_results, summary_file)

cat("\nAssurance summary:\n")
print(summary_results)

cat("\nSaved raw combined results to:\n")
cat(raw_file, "\n")

cat("\nSaved summary results to:\n")
cat(summary_file, "\n")
