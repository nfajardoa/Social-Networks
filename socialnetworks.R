########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Written by: Nicolas Fajardo
########################### Last modified: 23/05/2020

## Set Up
set.seed(5596)
setwd("F:/Social Networks") # Set working directory
source("socialnetworks_packages.R") # Load packages
source("socialnetworks_functions.R") # Load functions
plan(multiprocess,
  workers = 57 # Enable multiprocessing (parallelization)
)              # Be aware of the resources of your computer while setting the workers up,
               # each of them can consume up to 3GB of RAM.

display_summary <- FALSE # Display summary statistics for the data (TRUE/FALSE)

## Load data, Summary Statistics & Database check
# inputs: "network_data.csv" "covariates_data.csv" "perceptions_idalter.csv" "coincidences_nominations.csv"
# output objects: data, covariates_data, coincidences_data
# output values: loaded_data, loaded_covariates_data, loaded_perceptions_data
# comments: covariates_data contains also sum of perceptions
#           if the .csv files are not found, module:perceptions is automatically loaded
source("socialnetworks_summstats.R")

## Set parameters (Submodule WDistance)
save_csv <- FALSE # if TRUE, saves cvs files of finished computations.
save_progress <- TRUE # if TRUE, saves the environments after one computation is finished
save_memory <- FALSE # if TRUE, deletes objects after their use
use_db <- FALSE # Use a database environment (TRUR/FALSE). Using databases greatly 
                # reduces the RAM usage, however, they have to be set up independently and
                # no guide is provided. 

## Analysis variables (Submodules Communities and WDistance)
# Networks
selectors <- c("coll", "mut") # Which networks to use, collapsed or mutual
variables <- c("friend", "study", "social", "any") # Which variables to analyze
times <- c("none", "base", "mid", "end") # Which cross-sections to include 
coars <- distinct(data, coar) %>% pull() # Which schools
cohorts <- distinct(data, cohort) %>% pull() # Which cohorts

# Community identification (Submodule communities)
community_algorithms <- c(
  "edge_betweenness",
  "fast_greedy",
  "label_prop",
  "infomap",
  "leading_eigen",
  "louvain",
  "walktrap"
) # Optimal and spinglass were deleted

# Perceptions
characteristics <- c("friendly", "leader", "popular", "shy", "academic")
characteristics_t <- c("base", "end")

#### MAIN ANALYSIS ####

## Identify Communities, join Community Membership with Perception Data and summarize [NOT SUPPORTED]
# inputs: data, loaded_data
# outputs: communities, communities_membership, joint_communities_membership, membership_perceptions
#          membership_perceptions_summary, n_communities_summary, s_communities_membership
# file outputs: "communities_(date).RData" "communities_membership_(date).csv"
#               "membership_perceptions_summary_(date).csv" "n_communities_summary_(date).csv"
#               "adjacency_matrices_summary" "communities_summary"
source("socialnetworks_communities.R")

## Calculate Weighted Distance (network's adjacency matrix x distance matrix)
# inputs: "data, loaded_data
# outputs: distance_variables, adjacency_matrices, distance_matrices, distance_summary, id_cases, weighted_distance
# file outputs: "adjacency_matrices_(date).RData" "weighted_distance_(date).RData"
source("socialnetworks_wdistance.R")

## Regressions
# inputs: data, covariates_data, coincidences_data,
#         "weighted_distance_(date).RData" "s_communities_membership_(date).csv"
# outputs:
source("socialnetworks_regressions.R")

## Heterogeneous Treatment Effects
# inputs:
# outputs:
source("socialnetworks_heteffects.R")

#### SIDE CALCULATIONS ####

## Generate perceptions or coincidences data from the original network file
## and store each on separate .csv files
# inputs: "network_data.csv"
# outputs: perception_conditionals, perceptions_data, perceptions_school, coincidences
# file outputs: "perceptions_idalter.csv", "perceptions_school.csv", "coincidences_nominations.csv"
source("socialnetworks_perceptions.R")

## Generate plots and histrograms from communities data (Only one coar at a time) [NOT REFACTORED]
# inputs: communities_membership,
# outputs:
source("socialnetworks_plots.R")
