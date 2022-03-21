# Social-Networks
Implementation of “Generic Machine Learning Inference on Heterogeneous Treatment Effects in Randomized Experiments” from Chernozhukov et al. (2018)
Written by: Nicolas Fajardo
Last modified: 23/05/2020

## Set Up

See socialnetworks.R, socialnetworks_functions.R contains all self-made functions.

## Load data, Summary Statistics & Database check 
*source("socialnetworks_summstats.R")*
inputs: "network_data.csv" "covariates_data.csv" "perceptions_idalter.csv" "coincidences_nominations.csv"
output objects: data, covariates_data, coincidences_data
output values: loaded_data, loaded_covariates_data, loaded_perceptions_data
comments: covariates_data contains also sum of perceptions if the .csv files are not found, module:perceptions is automatically loaded

## Identify Communities, join Community Membership with Perception Data and summarize [NOT SUPPORTED]
*source("socialnetworks_communities.R")*
inputs: data, loaded_data
outputs: communities, communities_membership, joint_communities_membership, membership_perceptions membership_perceptions_summary, n_communities_summary, s_communities_membership
file outputs: "communities_(date).RData" "communities_membership_(date).csv" "membership_perceptions_summary_(date).csv" "n_communities_summary_(date).csv" "adjacency_matrices_summary" "communities_summary"

## Calculate Weighted Distance (network's adjacency matrix x distance matrix)
*source("socialnetworks_wdistance.R")*
inputs: "data, loaded_data
outputs: distance_variables, adjacency_matrices, distance_matrices, distance_summary, id_cases, weighted_distance
file outputs: "adjacency_matrices_(date).RData" "weighted_distance_(date).RData"

## Regressions
*source("socialnetworks_regressions.R")*
inputs: data, covariates_data, coincidences_data, "weighted_distance_(date).RData" "s_communities_membership_(date).csv"
outputs: LATEX tables

## Heterogeneous Treatment Effects
Implement “Generic Machine Learning Inference on Heterogeneous Treatment Effects in Randomized Experiments” from Chernozhukov et al. (2018)
*source("socialnetworks_heteffects.R")*

## Sida Calculations

*source("socialnetworks_perceptions.R")*
Generate perceptions or coincidences data from the original network file and store each on separate .csv files
inputs: "network_data.csv"
outputs: perception_conditionals, perceptions_data, perceptions_school, coincidences
file outputs: "perceptions_idalter.csv", "perceptions_school.csv", "coincidences_nominations.csv"

*source("socialnetworks_plots.R")*
Generate plots and histrograms from communities data (Only one coar at a time) [NOT REFACTORED]
inputs: communities_membership,
outputs: .pdf images
