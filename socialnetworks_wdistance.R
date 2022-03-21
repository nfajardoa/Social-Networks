########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Calculate Weighted Distance
########################### Written by: Nicolas Fajardo
########################### Last modified: 26/06/2020
########################### Comments:
###########################

## Check if data is loaded, if not, load it.
load_datasrc("socialnetworks_summstats.R")

## Create file
if(dir.exists("wdistance/")) {
  setwd("wdistance/")
}else {
  dir.create("wdistance/")
  setwd("wdistance/")}

## Retrieve variable names
distance_variables <- names(data) %>% str_subset("contiguo")

## Generate Adjacency Matrices
adjacency_matrices <- pmap_dfr(
  expand_grid(coar = coars, selector = selectors, cohort = cohorts, variable = variables, time = times), #Generate all possible combinations
  generate_undirected_adjacency_matrix, #Generate the undirected adjacency matrix (Custom function)
  quietly = TRUE, #All details are printed on adjacency_matrices_summary.txt
  data = data, #Use data
  matrix2graph = FALSE, #Do not turn the adjacency matrix to igraph object
  drop_singles = FALSE, #Do not drop single nodes
  force_diagonal = FALSE #Do not force diagonal values if different from 0
  ) %>%
  nest(adjacency_matrix_summary = c(data_percentage,           #Percentage of original data remaining by idego and idalter
                                    drop_singles,              #Should singles be dropped? (TRUE by default)
                                    dropped_values,            #Character vector of dropped values
                                    matrix_symmetry,           #Is the matrix symmetrical?
                                    forced_diagonal,           #Is the matrix diagonal converted to 0's? (TRUE by default)
                                    corrected_diagonal_values) #Which diagonal values were corrected (NULL if no value was forced)
  ) %>%
  rename(adjacency_matrix = network_data) %>% #Rename variable
  filter(nodes_n != 0L) # Filter empty values

closeAllConnections()

## Calculate Weighted Distance Matrices
distance_matrices <- pmap_dfr(
  expand_grid(coar = coars, selector = "none", cohort = cohorts, variable = distance_variables, time = "none"), #Generate all possible combinations
  generate_undirected_adjacency_matrix, #Generate the undirected adjacency matrix (Custom function)
  quietly = TRUE, #All details are printed on adjacency_matrices_summary.txt
  data = data, #Use data
  matrix2graph = FALSE, #Do not turn the adjacency matrix to igraph object
  drop_singles = FALSE, #Do not drop single nodes
  force_diagonal = FALSE #Do not force diagonal values if different from 0
  ) %>%
  nest(adjacency_matrix_summary = c(data_percentage,           #Percentage of original data remaining by idego and idalter
                                    drop_singles,              #Should singles be dropped? (TRUE by default)
                                    dropped_values,            #Character vector of dropped values
                                    matrix_symmetry,           #Is the matrix symmetrical?
                                    forced_diagonal,           #Is the matrix diagonal converted to 0's? (TRUE by default)
                                    corrected_diagonal_values) #Which diagonal values were corrected (NULL if no value was forced)
  ) %>%
  dplyr::select(-c(selector, time)) %>% #Drop selector and time variables, since the distance matrices do not depend on them
  rename(distance_variable = variable, #Rename all columns (in order to merge data)
         distance_adjacency_matrix = network_data, 
         distance_nodes_n = nodes_n, 
         distance_adjacency_matrix_summary = adjacency_matrix_summary
  ) %>% 
  mutate(distance_nodes_n = as.integer(distance_nodes_n)) %>% # Change variable type
  filter(distance_nodes_n != 0L) # Filter empty values

closeAllConnections()

# Calculate the sum of edges by distance variable (Summary variable)
distance_summary <- data %>%
  group_by(coar, cohort) %>% #Group data
  dplyr::select(coar, cohort, idego, idalter, starts_with("contiguo")) %>% #Select only necessary columns
  summarise_if(is.numeric, sum, na.rm = TRUE) %>% #Sum all numeric variables
  pivot_longer(cols = starts_with("contiguo"), names_to = "distance_variable", values_to = "distance_edge_sum") %>% #Reshape dataframe
  right_join(distance_matrices) %>% #Join summary with distance matrices
  drop_na(distance_edge_sum) %>% #Drop NA values
  nest(distance_data = c(distance_variable,                 #Name of the variable
                         distance_adjacency_matrix,         #Matrix object
                         distance_nodes_n,                  #NUmber of nodes
                         distance_adjacency_matrix_summary, #Summary of the generating process
                         distance_edge_sum))                #Sum of edges (easy check for sparse matrices)

# Calculate the matrix multiplication (Weighted Distance)
adjacency_matrices %<>% left_join(distance_summary, by = c("coar", "cohort")) %>% # Add distance matrices data (with summary)
  mutate(distance_data = map2(adjacency_matrix, distance_data, #First 2 arguments to the function
                              function(adjacency_matrix, distance_data) {
                                add_column(distance_data, #Add column to distance_data
                                           distance_weighted_network = calculate_distance(adjacency_matrix, distance_data) #Calculate new column (custom function)
                                           )
                            })
  ) %>%
  mutate(distance = map(distance_data, #Only argument to the function
                        function(distance_data) {
                          distance_data$distance_weighted_network %>% #Transform all the calculated matrices into tibble form
                            map2(distance_data$distance_variable, pipe_matrix2tibble) %>%
                            reduce(full_join, by = c("idego", "idalter")) #Merge all tibbles (Bind database again)
  }))

if (save_memory == TRUE) {
  rm(list = c("distance_matrices", "distance_summary"))
}
if (save_progress == TRUE) {
  save(adjacency_matrices, file = paste0("adjacency_matrices_", Sys.Date(), ".RData"))
}

### WARNING, use_db option was created in order to free low RAM devices (less than 40GB) 
### in order to make all necessary calculations, if this argument is set to true, the regressions
### module also uses it, and expects the database object to exist. In general, before being able to
### use this option is mandatory to set up a Postgres SQL database. Options can be directly changed
### in this code, but no database creation step is provided. I recommend using DBeaver community in 
### order to easily perform such a step. 

if (use_db == TRUE) { ## Load Database
  database <- src_postgres(dbname = "socialnetworks", user = "test", password = "12345", host = "localhost", port = "5432")
  database_connection <- dbConnect(RPostgres::Postgres(), dbname = "socialnetworks", host = "localhost", port = 5432, user = "test", password = "12345")
}

## Prepare Data of Weighted distance
id_cases <- quo(data %>% #Quoted expression of id_cases (This is made to reduce the generated data to the actual observed values, since the adjacency matrix generation step attributes 0's to missing values automatically)
                  dplyr::select(coar, cohort, idego, idalter, matches("contiguo")) %>% 
                  mutate_if(is.factor, as.character)) 

weighted_distance <- adjacency_matrices %>%
  dplyr::select(variable, selector, coar, cohort, time, adjacency_matrix, distance) %>% #Select columns
  rename(wdistance = distance) %>% #Rename variable
  group_by(variable, selector, time) %>% #Group data
  group_nest() %>% #Nest data by groups (see: tidyr)
  mutate(
    name = paste(variable, time, selector, sep = "_") %>% str_remove_all("_none"), #Generate a name variable
    data = map2(data, name, #First two arguments of the function
                merge_groups, #Merge tibbles by groups (Custom function)
                id_cases = id_cases, #Filter observed cases
                groups = c("coar", "cohort", "idego", "idalter"), #Group data by these variables
                .use_database = use_db) #Use database option
  )

if (save_memory == TRUE) {
  rm(list = c("adjacency_matrices"))
}
if (save_progress == TRUE) {
  save(weighted_distance, file = paste0("weighted_distance_", Sys.Date(), ".RData"))
}
if (use_db == TRUE) {
  dbDisconnect(database_connection)
} # Disconnect Database

## Exit
setwd("..")
