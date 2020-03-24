########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Written by: Nicolas Fajardo
########################### Last modified: 21/03/2020

##### LOADING ####

# install.packages("tidygraph")
# install.packages("ggraph")
# install.packages("tictoc")

library(tictoc)
library(tidyverse)
library(magrittr)
library(tidygraph)
library(ggraph)
library(igraph)
library(reshape2)

##### START #####

## Setting working directory
setwd("Documents/BRIQ/")
## Sourcing functions
source("network_functions.R")
## Load data
data <- read_csv("network_data.csv")
data$idalter %<>% as.character
names(data) %<>% str_replace_all("lider", "leader") # Replace some typos

## Summary Statistics & Database check

# ID's consistency

arguments <- expand_grid(variable = c("idego", "idalter"))
summary_id <- pmap_df(arguments, check_repeated, data = data, group = c("coar", "cohort"))

##### IDENTIFYING COMMUNITIES #####

## UNDIRECTED GRAPHS
# MUTUAL RELATIONS ONLY

## Set parameters

# Test arguments
variables <- c("any")
times <- c("base", "end") # mid was deleted
coar <- "lima"
cohort <- c(2015, 2016)
community_algorithms <- c("edge_betweenness", "fast_greedy", "label_prop", "infomap", "leading_eigen", "louvain", "walktrap") # Optimal and spinglass were deleted

# variables <- c("friend", "study" ,"social", "any")
# times <- c("none", "base", "end") #mid was deleted
# coar <- distinct(data, coar) %>% pull()
# coar <- coar[!coar %in% c("cajamarca","limaprovincias","tumbes")] #Deleting no running regions
# cohort <- distinct(data, cohort) %>% pull()

communities <- pmap_df(
  expand_grid(school = coar, cohort = cohort, variable = variables, time = times),
  wrapped_undirected_adjacency_matrix, data
) %>%
  expand_grid(algorithm = community_algorithms) %>%
  mutate(
    communities_data = map2(algorithm, network_data, wrapped_communities),
    communities_n = map_dbl(communities_data, wrapped_communities_n)
  ) %>%
  mutate(membership = map2(communities_data, algorithm, enframe_membership))

## DIRECTED (not yet implemented)

##### CHARACTERISTICS #####

## Set parameters

characteristics <- c("friendly", "leader", "popular", "shy")
characteristics_t <- c("base", "end")

## Summarize information of perceptions at the individual level

arguments <- expand_grid(characteristic = characteristics, time = characteristics_t) %>% mutate(filtering = paste(characteristic, time, sep = "_"))
characteristics_idalter <- pmap_df(arguments, function_one, data = data, groups = c("coar", "cohort", "idalter")) %>%
  group_by(coar, cohort, time, idalter) %>%
  summarise_all(sum, na.rm = TRUE)
names(characteristics_idalter)[names(characteristics_idalter) == "coar"] <- "school" # Change column name

## Summarize information of identified individuals according to the threshold (friendly, shy, popular, leader)
## School x Cohort

## WARNING: Threshold is set at 1 by default

arguments <- expand_grid(characteristic = characteristics)
characteristics_school <- pmap_df(arguments, function_two, data = characteristics_idalter, groups = c("school", "cohort", "time"), threshold = 1) %>%
  group_by(school, cohort, time) %>%
  summarise_all(sum, na.rm = TRUE)

## Summarize of Variable x Cohort x School x Time identifying membership by algorithm per idalter

communities_membership <- communities %>%
  group_by(variable, school, cohort, time, .drop = TRUE) %>%
  summarize(groups = merge_n_wrap(membership)) %>%
  unnest(cols = c(groups))

## JOIN MEMBERSHIP WITH PERCEPTION DATA

network_data_merged <- left_join(communities_membership, characteristics_idalter) %>% replace(is.na(.), 0)

# Calculate mean, community_mean, within and between variance of all characteristics and algorithms

arguments <- expand_grid(characteristic = characteristics, algorithm = community_algorithms)
network_summary <- pmap_df(arguments, summarize_by_algorithm, data = network_data_merged, threshold = 1, group = c("variable", "school", "cohort", "time"))

## STORING DATA

write_csv(summary_id %>% select_if(negate(is.list)), path = "summary_id.csv")
write_csv(network_summary %>% select_if(negate(is.list)), path = "network_summary.csv")
write_csv(characteristics_school, path = "characteristics_school.csv")
write_csv(characteristics_idalter, path = "characteristics_idalter.csv")
write_csv(network_data_merged, path = "network_data_merged_lima_any.csv")

save.image(file='social_networks.RData')

  ## PLOTS
# Only as an example
# WARNING: Implementation of automatized process is pending

test <- undirected_adjacency_matrix(data, school = "amazonas", cohort = 2015, variable = "friend", time = "end") %>%
  graph_from_adjacency_matrix() %>%
  as.undirected() %>%
  as_tbl_graph()

test_communities <- cluster_edge_betweenness(test)

test %<>% activate(nodes) %>%
  mutate(community = as.integer(membership(test_communities))) %>%
  mutate(community = as.factor(community))

# Plot with igraph
plot(test, vertex.label = NA, vertex.palette = diverging_pal, vertex.size = 6, vertex.frame.color = NA) # Network
plot(test, mark.groups = test_communities, vertex.label = NA, vertex.palette = diverging_pal, vertex.size = 6, vertex.frame.color = NA) # Communities

# Plot with ggraph
ggraph(test) + geom_edge_link(color = "gray") + geom_node_point(color = "orange", size = 3) + theme_graph() # Network
ggraph(test) + geom_edge_link(color = "gray") + geom_node_point(aes(color = community), size = 3) + theme_graph() # Communities

#################
