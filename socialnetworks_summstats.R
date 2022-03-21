########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Load data, Summary Statistics & Database check
########################### Written by: Nicolas Fajardo
########################### Last modified: 24/06/2020

## Check existence of all necessary .csv files
csv_filenames <- c("network_data.csv", "covariates_data.csv", "perceptions_idalter.csv", "coincidences_nominations.csv")
if (!all(file.exists(csv_filenames))) {
  message("Some (all) .csv files are missing") # Any of the csv_filenames is missing
  if (!all(file.exists(csv_filenames) %>% head(2))) { # If first two are missing
    message("Network and covariates data was not found (Are you in the right directory?)")
  }
  if (!all(file.exists(csv_filenames) %>% tail(2))) { # If last two are missing
    if (file.exists(csv_filenames)[[3]]) { # If only perceptions_idalter is not missing, calculate coincidences
      calculate_coincidences <- TRUE
      source("socialnetworks_perceptions.R")
    } else if (file.exists(csv_filenames)[[4]]) { # If only coincidences_nominations is not missing, calculate perceptions
      calculate_perceptions <- TRUE
      source("socialnetworks_perceptions.R")
    } else { # If both perceptions_idalter and coincidences_nominations are missing, calculate them
      source("socialnetworks_perceptions.R")
    }
  }
}

## Load network data
if (class(data)[[1]] == "function") {
  data <- read_csv("network_data.csv", # Read .csv on directory
    col_types = as.list(c(rep("f", 4), rep("i", 42), rep("f", 4), rep("i", 60)))
  ) %>% # Set column types
    rename_all(~ str_replace_all(.x, "lider", "leader")) %>% # Replace some typos
    rename_all(~ str_replace_all(.x, "skill", "academic"))
  loaded_data <- TRUE # Loaded data signal
}

## Load perceptions data
if (!exists("perceptions_data")) {
  perceptions_data <- read_csv("perceptions_idalter.csv", # Read .csv on directory
    col_types = as.list(c(rep("f", 4), rep("n", 5)))
  ) %>% # Set column types
    rename_if(is.numeric, ~ paste0("n_", .)) # Change nomination variables names (each starts with "n_")
  loaded_perceptions_data <- TRUE # Loaded data signal
}

## Load covariates data
if (!exists("covariates_data")) {
  covariates_data <- read_csv("covariates_data.csv", # Read .csv on directory
    col_types = as.list(c(rep("f", 15), rep("n", 7)))
  ) %>% # Set column types
    left_join(perceptions_data %>% # Sum perceptions by idalter
      pivot_wider(
        id_cols = c(coar, cohort, idalter), # Reshape tibble to match dimensions of covariates data
        names_from = characteristic_t,
        values_from = c(n_friendly, n_leader, n_popular, n_shy, n_academic)
      ) %>%
      rename(id = idalter), by = c("coar", "cohort", "id")) # Join number of perceptions of each variable as covariates
  loaded_covariates_data <- TRUE # Loaded data signal
}

## Load coincidences data
if (!exists("coincidences_data")) {
  coincidences_data <- read_csv("coincidences_nominations.csv", # Read .csv on directory
    col_types = as.list(c(rep("f", 4), rep("n", 10)))
  ) # Set column types
  loaded_coincidences_data <- TRUE # Loaded data signal
}

## Summary
if (display_summary == TRUE) {
  skim_without_charts(data) %>% print() # Summary of all variables
  cat(paste0("\n", "Repeated observations:"))
  pmap_dfr(expand_grid(variable = c("idego", "idalter")),
    check_repeated,
    data = data,
    group = c("coar", "cohort")
  ) %>%
    print() # Check repeated values in idalter and idego grouping by "coar" and "cohort"
  skim_without_charts(covariates_data) %>% print() # Summary of all variables
  skim_without_charts(perceptions_data) %>% print() # Summary of all variables
  skim_without_charts(coincidences_data) %>% print() # Summary of all variables
}
