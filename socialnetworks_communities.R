########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Identifiying Communities
########################### Written by: Nicolas Fajardo
########################### Last modified: 23/05/2020
########################### Comments:
########################### UNDIRECTED GRAPHS ONLY

## Check if data is loaded, if not, load it.
load_datasrc("socialnetworks_summstats.R", main = TRUE)

## Generate Graphs
tic("Graph generation")
communities <- future_pmap_dfr(
  expand_grid(school = coar, selector = selectors, cohort = cohort, variable = variables, time = times) %>% head(2),
  generate_undirected_adjacency_matrix,
  quietly = TRUE,
  data = data
) %>% nest(adjacency_matrix_summary = c(data_percentage, drop_singles, dropped_values, matrix_symmetry, forced_diagonal, corrected_diagonal_values)) 
closeAllConnections()
communities %<>% mutate(variable = as.factor(variable),
                        selector = as.factor(selector),
                        time = as.factor(time),
                        nodes_n = as.integer(nodes_n)) %>% 
  filter(nodes_n != 0L)
toc()

## Identify Communities
tic("Community identification")
communities %<>%
  expand_grid(algorithm = community_algorithms %>% as.factor()) %>%
  mutate(name = paste(variable, selector, school, cohort, time, sep = " "),
         communities_data = pmap(list(algorithm = algorithm, network_data = network_data, name = name), wrap_communities)) %>%
  dplyr::select(-name)
toc()

communities %<>% mutate(
  communities_n = map_int(communities_data, wrap_communities_n),
  membership = map2(communities_data, algorithm, enframe_membership, name = "idego"),
  network_data = map2(network_data, communities_data, pipe_addcomm2graph)
) %>% drop_na(communities_n)

## Save Data
save.image(file = paste0("communities_", Sys.Date(), ".RData"))

## Summarize of Variable x Cohort x School x Time identifying membership by algorithm per idalter
communities_membership <- communities %>%
  group_by(variable, selector, school, cohort, time, .drop = TRUE) %>%
  summarize(groups = merge_n_wrap(membership, variable = "idego")) %>%
  unnest(cols = c(groups)) %>% ungroup() %>% mutate(idego = idego %>% as.factor())

## Save Data
save(communities_membership, file = paste0("communities_membership_", Sys.Date(), ".RData"))
