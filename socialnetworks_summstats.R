########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Load data, Summary Statistics & Database check
########################### Written by: Nicolas Fajardo
########################### Last modified: 23/05/2020

## Load data
data <- read_csv("network_data.csv",
                 col_types = as.list(c(rep("f", 4), rep("i", 46-5+1), rep("f", 4), rep("i", 60)))) %>%
  rename_all(~str_replace_all(.x, "lider", "leader")) # Replace some typos
loaded_data <- TRUE

## ID Summary
#summary_id <- pmap_dfr(expand_grid(variable = c("idego", "idalter")),
#                       check_repeated, data = data, group = c("coar", "cohort"))