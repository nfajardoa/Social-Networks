########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Join Community Membership with Perception Data
########################### Written by: Nicolas Fajardo
########################### Last modified: 23/05/2020
########################### Comments:
###########################

## Check if data is loaded, if not, load it.
load_datasrc("socialnetworks_summstats.R")
load_lastdataset("weighted_distance", module = "wdistance")
load_lastdataset("s_communities_membership", module = "communities")

## Create file
if (dir.exists("regressions/")) {
  setwd("regressions/")
} else {
  dir.create("regressions/")
  setwd("regressions/")
}

## Set parameters
keep_covariates <- c("male") # Which covariates will be regressed (All regressions)

## Generate variables for regressions
merged_data <- data %>%
  left_join(covariates_data %>%
    dplyr::select(coar, cohort, id, !!!syms(keep_covariates)) %>%
    rename(idego = id), by = c("coar", "cohort", "idego")) %>%
  left_join(covariates_data %>%
    dplyr::select(coar, cohort, id, !!!syms(keep_covariates)) %>%
    rename(idalter = id), by = c("coar", "cohort", "idalter"), suffix = c(".idego", ".idalter")) %>%
  mutate(
    idego = as.character(idego),
    idalter = as.character(idalter),
    s_sex = if_else(male.idego == male.idalter, 1, 0),
    s_sex = factor(s_sex, levels = 0:1, labels = c("diff", "same")),
    i_sex = case_when(
      male.idego == 0 & male.idalter == 0 ~ 1,
      male.idego == 1 & male.idalter == 1 ~ 2,
      male.idego == 0 & male.idalter == 1 ~ 0,
      male.idego == 1 & male.idalter == 0 ~ 0
    ),
    i_sex = factor(i_sex, levels = 0:2, labels = c("other", "both girls", "both boys")),
    network = paste(coar, cohort, sep = "_"),
    network_id = as.integer(factor(network, levels = unique(network)))
  ) %>%
  left_join(coincidences_data, by = c("coar", "cohort", "idego", "idalter")) %>%
  mutate_at(vars(matches("coincidences")), replace_na, 0) %>%
  filter(idego != idalter) # Delete same idego and same idalter observations

if (save_csv == TRUE) {
  merged_data %>%
    mutate_all(as.character) %>%
    mutate_all(replace_na, replace = ".") %>%
    write_csv(path = paste0("regressions_data_", Sys.Date(), ".csv"))
}

if (use_db == TRUE) { ## Load Database
  database <- src_postgres(dbname = "socialnetworks", user = "test", password = "12345", host = "localhost", port = "5432")
  database_connection <- dbConnect(RPostgres::Postgres(), dbname = "socialnetworks", host = "localhost", port = 5432, user = "test", password = "12345")
  weighted_distance %<>% mutate(data = map(name, tbl, src = database)) ## Renew connections
}

## Extract Variable names
wdistance_names <- names(data) %>%
  str_subset("contiguo") %>%
  paste0("w", .) %>%
  paste(collapse = " + ")
distance_names <- names(data) %>%
  str_subset("contiguo") %>%
  paste(collapse = " + ")

#### Network ~ (Weighted) Distance Models
## Clustering: network_id
## Calculated factors: s_sex or i_sex
## Absorbed factors (fixed effects): idego + idalter
# Comments: By design, the element weighted_distance only contains the
# identification variables [idvars] (coar, cohort, idego, idalter), and the w- and distance
# variables, in order to add other variables to the regressions is necessary to create
# a dataframe which will be joined before the regression is executed (must have the idvars).
# The aforementioned join, preserves only the observations in the added data, while joining 
# both data frames columns (see right_join, where add_data is the right tibble).

add_data <- merged_data %>%
  dplyr::select(coar, cohort, idego, idalter, i_sex)

# List of formulas of desired models (see felm)
weighted_distance_formulas <- list(
  paste0(distance_names, "| idego + idalter + factor(i_sex) | 0 | network_id"),
  paste0(wdistance_names, "| idego + idalter + factor(i_sex) | 0 | network_id"),
  paste0(wdistance_names, " + ", distance_names, "| idego + idalter + factor(i_sex) | 0 | network_id")
)

# Run regressions
weighted_distance_models <-
  regress_network(
    weighted_distance,
    weighted_distance_formulas,
    response_col = "name",
    nested_data = "data",
    add_data = add_data,
    baseline_data = merged_data,
    control_baseline = TRUE,
    cmethod = "reghdfe"
  )

# Reshape and print LaTeX tables of estimated models
weighted_distance_models %<>%
  unnest(models) %>%
  group_by(variable, selector) %>%
  group_split() %>%
  walk(sink_latex_output, filename = "weighted_distance_regressions.txt") %>%
  reduce(bind_rows)

# Save models object on file
save(weighted_distance_models, file = paste0("weighted_distance_models_", Sys.Date(), ".RData"))

if (save_memory == TRUE) {
  rm(weighted_distance_models)
}

#### Coincidences ~ Distance Models
# Clustering: network_id
# Calculated factors: none
# Absorbed factors (fixed effects): idego + idalter + factor(i_sex)

coincidences_distance_formula <- list(
  paste0(distance_names, "| idego + idalter + factor(i_sex) | 0 | network_id"),
  "ingroup | idego + idalter  + factor(i_sex) | 0 | network_id"
)

# Run regressions
coincidences_distance_models <-
  regress_network(
    expand_grid(characteristic = characteristics, time = characteristics_t) %>%
      mutate(variable = paste("coincidences", characteristic, time, sep = "_")),
    coincidences_distance_formula,
    response = "variable",
    use_data = merged_data,
    control_baseline = TRUE,
    cmethod = "reghdfe"
  )

# Save models object on file
save(coincidences_distance_models, file = paste0("coincidences_distance_models_", Sys.Date(), ".RData"))

# Print models
coincidences_distance_models %>%
  group_by(characteristic) %>%
  group_split() %>%
  walk(sink_latex_output, filename = "coincidences_distance_regressions.txt")

if (save_memory == TRUE) {
  rm(coincidences_distance_models)
}

#### Nominations ~ Distance Models
# Clustering: network_id
# Calculated factors: none
# Absorbed factors (fixed effects): idego + idalter + factor(i_sex)

nominations_distance_formula <- list(
  paste0(distance_names, " | idego + idalter + factor(i_sex) | 0 | network_id"),
  "ingroup | idego + idalter + factor(i_sex) | 0 | network_id"
)

# Run regressions
nominations_distance_models <-
  regress_network(
    expand_grid(characteristic = characteristics, time = characteristics_t) %>%
      mutate(variable = paste(characteristic, time, sep = "_")),
    nominations_distance_formula,
    response = "variable",
    use_data = merged_data,
    control_baseline = TRUE,
    cmethod = "reghdfe"
  )

# Save models object on file
save(nominations_distance_models, file = paste0("nominations_distance_models_", Sys.Date(), ".RData"))

# Print models
nominations_distance_models %>%
  group_by(characteristic) %>%
  group_split() %>%
  walk(sink_latex_output, filename = "nominations_distance_regressions.txt")

if (save_memory == TRUE) {
  rm(nominations_distance_models)
}

## Exit
setwd("..")
