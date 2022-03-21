########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Identifiying Communities
########################### Written by: Nicolas Fajardo
########################### Last modified: 26/06/2020
########################### Comments:
########################### UNDIRECTED GRAPHS ONLY

## Check if data is loaded, if not, load it.
load_datasrc("socialnetworks_summstats.R")

## Create file
if(dir.exists("communities/")) {
  setwd("communities/")
}else {
  dir.create("communities/")
  setwd("communities/")}

## Generate Graphs
communities_data <- pmap_dfr(
  expand_grid(coar = coars, selector = selectors, cohort = cohorts, variable = variables, time = times), #Generate all possible combinations
  generate_undirected_adjacency_matrix, #Generate the undirected adjacency matrix (Custom function)
  quietly = TRUE, #All details are printed on adjacency_matrices_summary.txt
  data = data #Use data
  ) %>% #Nest process summary data
  nest(adjacency_matrix_summary = c(data_percentage,           #Percentage of original data remaining by idego and idalter
                                    drop_singles,              #Should singles be dropped? (TRUE by default)
                                    dropped_values,            #Character vector of dropped values
                                    matrix_symmetry,           #Is the matrix symmetrical?
                                    forced_diagonal,           #Is the matrix diagonal converted to 0's? (TRUE by default)
                                    corrected_diagonal_values) #Which diagonal values were corrected (NULL if no value was forced)
  ) %>%  
  mutate(variable = as.factor(variable), #Transform to factor
         selector = as.factor(selector), #Transform to factor
         time = as.factor(time),         #Transform to factor
         nodes_n = as.integer(nodes_n)   #Transform to integer
  ) %>% 
  filter(nodes_n != 0L) #Drop empty adjacency matrices

closeAllConnections() #Stop writting on process_summary.txt

# Save generated graphs
if (save_progress == TRUE) {
save(communities_data, file = paste0("networks_graphs_", Sys.Date()))
}

## Identify Communities
communities_data %<>%
  expand_grid(algorithm = community_algorithms %>% as.factor()) %>% #Add all the specified algorithms
  mutate(name = paste(variable, selector, coar, cohort, time, sep = " "), #Generate a name variable
         community = pmap(list(algorithm = algorithm, network_data = network_data, name = name), #Generate all possible combinations
                          wrap_communities #Calculate communities (communities-object from igraph)
                          ) 
         ) %>%
  dplyr::select(-name) %>% #Drop name variable
  mutate(
  communities_n = map_int(community, wrap_communities_n), #Calculate number of communities for each network
  membership = map2(community, algorithm, enframe_membership, name = "idego"), #Generate a nested tibble of membership data
  network_data = map2(network_data, community, pipe_addcomm2graph) #Add the communities information to the generated graphs
  ) %>% drop_na(communities_n) #Drop empty networks

## Summarize number of communities by network and algorithm
n_communities_summary <- communities_data %>% 
  dplyr::select(variable, time, selector, algorithm, coar, cohort, communities_n) %>%
  group_by(variable, time, selector, algorithm) %>%
  dplyr::summarize(
    mean = mean(communities_n),
    min = min(communities_n),
    q25 = quantile(communities_n, probs = 0.25),
    q75 = quantile(communities_n, probs = 0.75),
    max = max(communities_n)
  ) %>% arrange(variable, selector, time)

## Summarize Variable x Cohort x School x Time identifying membership by algorithm per idalter
communities_membership <- communities_data %>%
  group_by(variable, selector, coar, cohort, time, .drop = TRUE) %>%
  summarize(groups = merge_n_wrap(membership, variable = "idego")) %>%
  unnest(cols = c(groups)) %>% 
  ungroup() %>% 
  mutate(idego = idego %>% as.factor())

## Summarize number of characteristics per community
membership_perceptions <- left_join(communities_membership, 
                                    perceptions_data %>% 
                                    rename(idego = idalter), # idalter was renamed as idego (for merging purposes)
                                    by = c("coar", "cohort", "idego"))

membership_perceptions_summary <- pmap_dfr(expand_grid(characteristic = paste("n", characteristics, sep = "_"), algorithm = community_algorithms), 
                                                  summarize_by_algorithm, 
                                                  data = membership_perceptions, 
                                                  threshold = 1, 
                                                  variable = "idego", 
                                                  group = c("variable", "selector", "coar", "cohort", "time", "characteristic_t"))

## Save Data
if(save_csv == TRUE) {
  write_csv(communities_membership, path = paste0("communities_membership_", Sys.Date(), ".csv"))
  write_csv(membership_perceptions_summary %>% dplyr::select(-community_mean), path = paste0("membership_perceptions_summary_", Sys.Date(), ".csv"))
  write_csv(n_communities_summary, path = paste0("n_communities_summary_", Sys.Date(), ".csv"))
}
if (save_progress == TRUE) save.image(file = paste0("communities_", Sys.Date(), ".RData"))
if (save_memory == TRUE) rm(communities, n_communities_summary, membership_perceptions_summary)

## Create same community variable by community algorithm (s_*_network)
joint_communities_membership <- communities_membership %>%
  pivot_wider(id_cols = c(coar, cohort, idego),
              names_from = c(variable, time, selector), 
              values_from = contains(community_algorithms)
  ) %>% 
  rename_at(vars(matches("none")), str_remove_all, "none_") %>%
  full_join(., rename(., idalter = idego), 
            by = c("coar", "cohort"), 
            suffix = c(".idego", ".idalter"))

community_variables <- names(joint_communities_membership) %>% 
  str_subset(c(".idego", ".idalter")) %>% 
  str_remove_all("(.idego)|(.idalter)") %>% 
  unique()

s_communities_membership <- map(community_variables, mutate_by_algorithm, data = joint_communities_membership) %>%
  reduce(inner_join, by = c("coar", "cohort", "idalter", "idego"))

## Save Data
if (save_memory == TRUE) rm(communities_membership, joint_communities_membership)
if (save_csv == TRUE) write_csv(s_communities_membership, path = paste0("s_communities_membership_", Sys.Date(), ".csv"))

## Exit
setwd("..")
