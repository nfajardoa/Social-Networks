########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Written by: Nicolas Fajardo
########################### Last modified: 23/05/2020

## Set Up
set.seed(5596)
setwd("Documents/BRIQ/Roman/Data Analysis/") # Set working directory
source("socialnetworks_packages.R") # Load packages
source("socialnetworks_functions.R") # Load functions
#plan(multiprocess) # Enable multiprocessing 
    
## Start

## Load data, Summary Statistics & Database check
# inputs: "network_data.csv"
# outputs: data, loaded_data, summary_id
source("socialnetworks_summstats.R")

## Set parameters
selectors <- c("coll", "mut")
variables <- c("friend", "study" ,"social", "any")
times <- c("none", "base", "mid", "end")
coar <- distinct(data, coar) %>% pull()
cohort <- distinct(data, cohort) %>% pull()
community_algorithms <- c("edge_betweenness", "fast_greedy", "label_prop", "infomap", "leading_eigen", "louvain", "walktrap") # Optimal and spinglass were deleted

## Identifiying Communities
# inputs: data, loaded_data
# outputs: selectors, variables, times, coar, cohort, community_algorithms, communities, communities_membership
# database output: "communities_date.RData" "communities_membership_date.RData"
source("socialnetworks_communities.R")

## Summarize Perception Data
# inputs: data, loaded_data
# outputs: perceptions_idalter, perceptions_school
source("socialnetworks_perceptions.R")

## Join Community Membership with Perception Data
# inputs: communities_membership, perceptions_idalter, perceptions_school
# outputs: 
source("socialnetworks_jointdatabase.R")

## Plots and Histograms
# inputs: communities_membership, 
# outputs:
source("socialnetworks_plots.R")

## Calculate mean, community mean, within and between variance of perceptions x algorithms
# inputs:
# outputs:
source("socialnetworks_communitysumm.R")

## Covariates
# inputs:
# outputs:
source("socialnetworks_covariates.R")

## Count Coincidences
# inputs:
# outputs:
source("socialnetworks_coincidences.R")

## Calculate Weighted Distance
# inputs:
# outputs:
source("socialnetworks_wdistance.R")

## Regressions
# inputs:
# outputs:
source("socialnetworks_regressions.R")

## Heterogeneous Treatment Effects
# inputs:
# outputs:
source("socialnetworks_hetteffects.R")
