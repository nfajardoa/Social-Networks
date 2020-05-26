########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Calculate Weighted Distance
########################### Written by: Nicolas Fajardo
########################### Last modified: 23/05/2020
########################### Comments:
########################### 

## Retrieve variable names
distance_variables <- names(data) %>% str_subset("contiguo")

## Generate Adjacency Matrices
tic("Adjacency Matrices generation")
adjacency_matrices <- pmap_dfr(
  expand_grid(school = coar, selector = selectors, cohort = cohort, variable = variables, time = times),
  generate_undirected_adjacency_matrix,
  quietly = TRUE,
  data = data,
  matrix2graph = FALSE,
  drop_singles = FALSE,
  force_diagonal = FALSE
) %>%
  nest(adjacency_matrix_summary = c(data_percentage, drop_singles, dropped_values, matrix_symmetry, forced_diagonal, corrected_diagonal_values)) %>%
  rename(adjacency_matrix = network_data)
closeAllConnections()
adjacency_matrices %<>% mutate(variable = as.factor(variable),
                        selector = as.factor(selector),
                        time = as.factor(time),
                        nodes_n = as.integer(nodes_n)) %>% 
  filter(nodes_n != 0L)
toc()

## Calculate Weighted Distance Matrices
tic("Weighted Distance Matrices Generation")
distance_matrices <- pmap_dfr(
  expand_grid(school = coar, selector = "none", cohort = cohort, variable = distance_variables, time = "none"),
  generate_undirected_adjacency_matrix,
  quietly = TRUE,
  data = data,
  matrix2graph = FALSE,
  drop_singles = FALSE,
  force_diagonal = FALSE
) %>%
  nest(adjacency_matrix_summary = c(data_percentage, drop_singles, dropped_values, matrix_symmetry, forced_diagonal, corrected_diagonal_values)) %>%
  dplyr::select(-c(selector, time)) %>%
  rename(distance_variable = variable, distance_adjacency_matrix = network_data, distance_nodes_n = nodes_n, distance_adjacency_matrix_summary = adjacency_matrix_summary)
closeAllConnections()
distance_matrices %<>% mutate(distance_nodes_n = as.integer(distance_nodes_n)) %>% filter(distance_nodes_n != 0L)
toc()

distance_summary <- data %>%
  group_by(coar, cohort) %>%
  dplyr::select(coar, cohort, idego, idalter, starts_with("contiguo")) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE) %>%
  pivot_longer(cols = starts_with("contiguo"), names_to = "distance_variable", values_to = "distance_edge_sum") %>%
  rename(school = coar) %>%
  right_join(distance_matrices) %>%
  drop_na(distance_edge_sum) %>%
  nest(distance_data = c(distance_variable, distance_adjacency_matrix, distance_nodes_n, distance_adjacency_matrix_summary, distance_edge_sum))

adjacency_matrices %<>% filter(nodes_n != 0) %>%
  left_join(distance_summary, by = c("school", "cohort")) %>%
  mutate(distance_data = map2(adjacency_matrix, distance_data, ~ add_column(.y, distance_weighted_network = calculate_distance(.x, .y)))) %>%
  mutate(distance = map(distance_data, ~ .x$distance_weighted_network %>%
                          map2(.x$distance_variable, pipe_matrix2tibble) %>%
                          reduce(full_join, by = c("idego", "idalter"))))

rm(list = c("distance_matrices", "distance_summary"))
save(adjacency_matrices, file = paste0("weighted_distance_", Sys.Date(), ".RData"))