########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Estimate Heterogeneous Effects with ML Methods
########################### (Chernozhukov, Demirer, Duflo, Fernandez-Val)
########################### Written by: Nicolas Fajardo
########################### Last modified: 04/07/2020
########################### Comments:
###########################

## Check if data is loaded, if not, load it.
load_datasrc("socialnetworks_summstats.R")

## Prepare data
initial_data <- data %>%
  left_join(covariates_data %>% rename(idego = id), by = c("coar", "cohort", "idego")) %>%
  left_join(covariates_data %>% rename(idalter = id), by = c("coar", "cohort", "idalter"), suffix = c(".idego", ".idalter")) %>%
  mutate(
    idego = as.character(idego),
    idalter = as.character(idalter),
    soft_skls = if_else((friendly_end | popular_end | leader_end) == 1, 1, 0),
    soft_base_skls = if_else((friendly_base | popular_base | leader_base) == 1, 1, 0),
    hard_skls = academic_end,
    hard_base_skls = academic_base,
    ingroup_5 = if_else((contiguo1 | contiguo2 | contiguo3 | contiguo4 | contiguo5) == 1, 1, 0),
    ingroup_3 = if_else((contiguo1 | contiguo2 | contiguo3) == 1, 1, 0),
    any_base_coll_ingroup = any_base_coll * ingroup,
    any_base_coll_ingroup_5 = any_base_coll * ingroup_5,
    any_base_coll_ingroup_3 = any_base_coll * ingroup_3,
    any_base_mut_ingroup = any_base_mut * ingroup,
    any_base_mut_ingroup_5 = any_base_mut * ingroup_5,
    any_base_mut_ingroup_3 = any_base_mut * ingroup_3,
    any_base_dir_ingroup = any_base_dir * ingroup,
    any_base_dir_ingroup_5 = any_base_dir * ingroup_5,
    any_base_dir_ingroup_3 = any_base_dir * ingroup_3,
    soft_base_skls_ingroup = soft_base_skls * ingroup,
    soft_base_skls_ingroup_5 = soft_base_skls * ingroup_5,
    soft_base_skls_ingroup_3 = soft_base_skls * ingroup_3,
    hard_base_skls_ingroup = hard_base_skls * ingroup,
    hard_base_skls_ingroup_5 = hard_base_skls * ingroup_5,
    hard_base_skls_ingroup_3 = hard_base_skls * ingroup_3,
    network = paste(coar, cohort, sep = "_"),
    network_id = as.integer(factor(network, levels = unique(network))),
  )

h_controls <- c(
  # "male", Removed due to collinearity
  "highc",
  "highs",
  # "Hsocial_exp",
  # "Hcognitive_exp",
  "rural",
  "vraem",
  "poor_score",
  "educ_padre",
  "educ_madre"
)

same_vars <- initial_data %>%
  dplyr::select(coar, cohort, idego, idalter, matches(h_controls)) %>%
  pivot_longer(
    cols = matches(h_controls),
    names_to = c("variable", "id_type"),
    names_pattern = "(.*)\\.(.*)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = id_type, values_from = value, names_prefix = "value_") %>%
  mutate(across(matches("^s_"), as.character)) %>%
  mutate(s = if_else(value_idego == value_idalter, 1, 0)) %>%
  pivot_wider(
    id_cols = c(coar, cohort, idego, idalter),
    names_from = variable,
    names_prefix = "s_",
    values_from = s
  ) %>%
  mutate(across(matches("^s_"), ~ factor(.x, levels = 0:1, labels = c("diff", "same"))))

initial_data %<>% left_join(same_vars, by = c("coar", "cohort", "idego", "idalter"))

nominations <- names(initial_data) %>%
  str_subset("^n_") %>%
  str_subset("end", negate = TRUE)
homophily <- names(initial_data) %>% str_subset("^s_")

### PARAMETERS SET-UP
# vector of outcome variables
response_vars <- c(
  "any_coll",
  "any_mut",
  "any_dir",
  "soft_skls",
  "hard_skls"
)

# vector of treatment variables
treatment_vars <- c("ingroup", "ingroup_3", "ingroup_5", "x_ingroup", "x_ingroup_3", "x_ingroup_5")
# Clustered errors for regressions (if no cluster use cluster <- "0")
clusters <- "network_id"
# Fixed effects for regressions (if no fixed_effect use fixed_effect <- "0")
fixed_effects <- "network_id"

# create a vector of control variables
controls <- c(
  "male",
  "highc",
  "highs",
  "Hsocial_exp",
  "Hcognitive_exp",
  "rural",
  "vraem",
  "poor_score",
  "educ_padre",
  "educ_madre",
  "admission_test",
  "mate_baseline",
  "liter_baseline",
  "centrality_base"
)

controls <- c(paste0(controls, ".idego"), paste0(controls, ".idalter"), nominations, homophily)

# characteristics for most/least affected analysis
affected <- c(
  "male",
  "highc",
  "highs",
  "rural",
  "vraem",
  "poor_score"
)

affected <- c(paste0(affected, ".idego"), paste0(affected, ".idalter"), nominations, homophily)

# Join parameters
analysis_parameters <-
  expand_grid(
    response = response_vars,
    controls = controls %>% list(),
    treatment = treatment_vars,
    cluster = clusters,
    fixed_effect = fixed_effects,
    affected = affected %>% list()
  ) %>%
  mutate(
    treatment = map2_chr(response, treatment, add_interaction2treatment),
    controls = map2(response, controls, add_baseline2parameters),
    affected = map2(response, affected, add_baseline2parameters)
  )

### METHODS SET-UP
## Elastic Net

elastic_net <- tibble(
  method_name = "elastic_net",
  recipe = quote(
    step_meanimpute(all_numeric()) %>%
      step_modeimpute(all_nominal()) %>%
      step_dummy(all_predictors(), -all_numeric()) %>%
      step_center(all_predictors()) %>%
      step_scale(all_predictors()) # Steps to modify data
  ) %>% list(),
  model = quote(
    parsnip::linear_reg(penalty = tune::tune(), mixture = tune::tune()) %>%
      parsnip::set_engine("glmnet")
  ) %>% list(),
  tuner = TRUE
  #   quote(
  #   process %>% # Default argument (Do not change)
  #     tune::tune_grid(
  #       resamples = data_cv, #data_cv is the resampling object by default, don't change it
  #       grid = 100,
  #       metrics = yardstick::metric_set(yardstick::rmse)
  #     )
  # ) %>% list()
)

## Random Forest

random_forest <- tibble(
  method_name = "random_forest",
  recipe = quote(
    step_meanimpute(all_numeric()) %>%
      step_modeimpute(all_nominal()) %>%
      step_dummy(all_predictors(), -all_numeric()) %>%
      step_center(all_predictors()) %>%
      step_scale(all_predictors()) # Steps to modify data
  ) %>% list(),
  model = quote(
    parsnip::rand_forest(mtry = 5, mode = "regression") %>%
      parsnip::set_engine("ranger")
  ) %>% list(),
  tuner = NA
)

## Join especified methods
methods <- bind_rows(random_forest, elastic_net)

### Join set-ups
setup <- expand_grid(
  parameters = analysis_parameters %>% list(),
  method = methods %>% list()
) %>%
  unnest(cols = parameters) %>%
  unnest(cols = method)

filtered_data <- initial_data %>%
  filter(cohort == "2015" | cohort == "2016") %>%
  select_data4anaylsis(analysis_parameters)

# Remove unused objects
rm(
  data, initial_data, coincidences_data, covariates_data, perceptions_data,
  elastic_net, random_forest, methods, same_vars, analysis_parameters,
  affected, clusters, controls, fixed_effects, h_controls, homophily,
  treatment_vars, nominations, response_vars, csv_filenames
)

# Remove problematic model
setup %<>% filter(response != "hard_skls" | treatment != "hard_base_skls_ingroup_5" | method_name != "elastic_net") %>%
  filter(response != "hard_skls" | treatment != "hard_base_skls_ingroup_3" | method_name != "elastic_net") %>%
  filter(response != "hard_skls" | treatment != "hard_base_skls_ingroup" | method_name != "elastic_net")

save.image("Setup.RData")

tic("Random Forest and Tuned Elastic Net: 100 Simulations (30 models)")
simulations <-
  simulate_heterogeneity_analysis(
    move_seed = 64,
    simulations = 1,
    data = filtered_data,
    stratum = "network_id",
    control = setup,
    folds = 2,
    repetitions = 2,
    groups = 5,
    affected_quantile = 0.25,
    significance = 0.05,
    save_data = FALSE
  )
closeAllConnections()

results <- read_n_merge_simulations() %>%
  correct_simulation_tests() %>%
  filter(simulation >= 100)
check_results <- read_n_merge_simulations(check_seeds = TRUE) %>%
  correct_simulation_tests() %>%
  filter(simulation >= 100)

## Print tables
best <- print_results(results, "best")
blp <- print_results(results, "blp")
gate <- print_results(results, "gate")
clan <- print_results(results, "clan") %>% arrange(`CLAN Variable`)

## Print images
ggsave("HET_groups.png", plot_results(results), width = 15, height = 60, limitsize = FALSE)
ggsave("Nominations.png", plot_clan(results, filter = "^n_", modify_vars = TRUE, use_grouping = TRUE), width = 15, height = 20, limitsize = FALSE)
ggsave("Homophily.png", plot_clan(results, filter = "(highc)|(highs)|(poor_score)|(rural)|(vraem)", modify_vars = TRUE, use_grouping = TRUE), width = 15, height = 25, limitsize = FALSE)
ggsave("Other.png", plot_clan(results, filter = "(_base_)|(educ_padre)|(educ_madre)", modify_vars = TRUE, use_grouping = FALSE), width = 15, height = 7.5, limitsize = FALSE)

## Save tables
sink("Results_2015-2016_Best.txt")
gt(best) %>%
  fmt_markdown(columns = everything()) %>%
  as_latex() %>%
  cat()

sink("Results_2015-2016_BLP.txt")
gt(blp) %>%
  fmt_markdown(columns = everything()) %>%
  as_latex() %>%
  cat()

sink("Results_2015-2016_GATE.txt")
gt(gate) %>%
  fmt_markdown(columns = everything()) %>%
  as_latex() %>%
  cat()

sink("Results_2015-2016_CLAN.txt")
gt(clan) %>%
  fmt_markdown(columns = everything()) %>%
  as_latex() %>%
  cat()
closeAllConnections()

## Save results
save.image(file = "Results_2015-2016.RData", version = 3)
