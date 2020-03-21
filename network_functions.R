########################### NETWORKS FUNCTIONS ###########################
########################### Written by: Nicolas Fajardo
########################### Last modified: 21/03/2020

undirected_graph_generator <- function() {
  cat("\n Variable:", i, " School:", j, " Cohort:", k, " Time:", t, "\n") # Information about process

  network <- undirected_adjacency_matrix(data, school = j, cohort = k, variable = i, time = t) %>%
    graph_from_adjacency_matrix() %>%
    as.undirected() %>%
    as_tbl_graph() # Using undirected_adjacency_matrix(.) to create adjacency matrix
  # Also graph is read as undirected and converted to a tidygraph object
  nodes_n <- network %>%
    unlist(recursive = FALSE) %>%
    `[[`(1)

  summary_networks %<>% add_row(
    variable = i, school = j, cohort = k, time = t, nodes_n = nodes_n,
    network_data = list(network)
  ) # Storing information on data frame for replication purposes
}

#### GENERATE UNDIRECTED ADJACENCY MATRIX FOR THE DATASET
## Arguments:
# data := Dataset
# school := School to filter
# cohort := Cohort to filter
# variable := Variable to filter ("friend", "study" ,"social", "any")
# time := Time to filter ("none", "base", "end", "mid")
# delete_missing := Delete missing observations
# selector := Which columns must be used ("mut", "any")
# replace_na := Replace missing values with 0s
# drop_singles := Drop unconnected nodes
# force_diagonal := Delete loops inside data [Setting all diagonal elements to 0]

undirected_adjacency_matrix <- function(data, school, cohort, variable, time = NULL, delete_missing = TRUE,
                                        selector = "mut", replace_na = TRUE, drop_singles = TRUE, force_diagonal = TRUE) {

  # Create variable to filter using the function arguments
  if (time == "none") {
    variable <- c(variable, selector) %>%
      paste(collapse = "_") %>%
      paste(collapse = "") # If none, only the variable name
    # and the selector is concatenated
  } else {
    variable <- c(variable, time, selector) %>%
      paste(collapse = "_") %>%
      paste(collapse = "")
  }

  # Filter the desired data
  data_group <- filter(data, coar == school, cohort == cohort) %>%
    select(-c(coar, cohort)) %>%
    select(idego, idalter, !!variable) %>%
    arrange(idego, idalter)

  # If no filtered data exists, exit the function
  if (is_empty(data_group)) {
    return(NULL)
  }

  # Perform missing deletion
  if (delete_missing == TRUE) {
    data_group %<>% drop_na()
  }

  # Reshaping data frame to wide form (almost as a matrix)
  data_group %<>% spread(key = idalter, value = !!variable)

  # Perform NA replacement
  if (replace_na == TRUE) {
    data_group %<>% replace(is.na(.), 0)
  }

  # Renaming rows with each node id and converting the data into a matrix object
  data_group %<>% column_to_rownames("idego") %>% as.matrix()

  # Identify which ids are not both in rows and columns (Conformability)
  difference <- union(setdiff(rownames(data_group), colnames(data_group)), setdiff(colnames(data_group), rownames(data_group)))
  intersection <- intersect(rownames(data_group), colnames(data_group))
  diffrows <- setdiff(rownames(data_group), intersection)
  diffcols <- setdiff(colnames(data_group), intersection)

  # Delete the outlier ids
  if (!is_empty(difference)) {
    cat("Eliminating id(s)", difference, "for ensuring complete pairwise relations\n")
    cat("Row id(s):", diffrows, " Column id(s):", diffcols, "\n")
    data_group <- data_group[rownames(data_group) %in% intersection, colnames(data_group) %in% intersection]
  }

  # Perform singles dropping
  if (drop_singles == TRUE) {
    empty_rows <- rownames(data_group[rowSums(data_group) == 0, ])
    empty_cols <- colnames(data_group[, colSums(data_group) == 0])
    singles <- intersect(empty_rows, empty_cols)

    # If the all network is empty, exit the function
    if (all(!(rownames(data_group) %in% singles)) & all(!(colnames(data_group) %in% singles))) {
      return(NULL)
    }

    # Identify which ids are not both in rows and columns (Conformability)
    data_group <- data_group[!(rownames(data_group) %in% singles), !(colnames(data_group) %in% singles)]
    difference <- union(setdiff(rownames(data_group), colnames(data_group)), setdiff(colnames(data_group), rownames(data_group)))
    intersection <- intersect(rownames(data_group), colnames(data_group))
    diffrows <- setdiff(rownames(data_group), intersection)
    diffcols <- setdiff(colnames(data_group), intersection)

    # Delete the outlier ids
    if (!is_empty(difference)) {
      cat("Eliminating id(s)", difference, "for ensuring complete pairwise relations\n")
      cat("Row id(s):", diffrows, " Column id(s):", diffcols, "\n")
      data_group <- data_group[rownames(data_group) %in% intersection, colnames(data_group) %in% intersection]
    }
  }

  # CHECKS

  check_symmetry <- all(data_group == t(data_group))
  cat("Is adjacency Matrix symmetric? ", check_symmetry, "\n")
  check_diagonal <- all(diag(data_group) == 0)
  cat("Are all the elements of the diagonal of the adjacency matrix equal to 0? ", check_diagonal, "\n")

  # Force all diagonal elements to be 0
  if (force_diagonal == TRUE) {
    diag(data_group) <- 0
  }

  return(data_group)
}

#### MISCELLANEOUS FUNCTIONS FOR SUMMARIZING

### Calculate membership and store it in a proper tibble
## Arguments:
# x := Dataset
# y := Value of variable

enframe_membership <- function(x, y) {
  membership(x) %>% enframe(name = "idalter", value = y)
}

### Merge every membership?s tibble
## Arguments:
# x := Tibble list

merge_n_wrap <- function(x) {
  reduce(x, full_join, by = "idalter") %>% list()
}

### Calculate overall and community mean, within and between variance for one algorithm x one characteristic
## Arguments:
# x := Dataset
# algorithm := Algorithm to filter ("edge_betweenness", "fast_greedy", "label_prop", "infomap", "leading_eigen", "louvain", "walktrap")
# characteristic := Characteristic to filter ("friendly", "leader", "popular", "shy")
# threshold := Minimum counting threshold
# group := Grouping variables ("variable", "school", "cohort", "time") [Variable, School, Cohort, Time]

summarize_by_algorithm <- function(data, algorithm, characteristic, threshold, group) {
  group_comm <- append(group, algorithm) # Group by pre-defined grouping variable plus algorithm (Operate by community)
  group_result <- c(group, "overall_mean", "within_var", "between_var") # Group by everything (To return a clean outcome)

  data %<>% select(group, idalter, !!algorithm, !!characteristic) %>% # Select appropiate columns
    mutate(!!(characteristic) := ifelse((!!rlang::sym(characteristic) >= threshold), 1, 0)) %>% # Turn the characteristic variable to 1, if greater or equal to threshold, 0, otherwise
    ungroup() %>%
    group_by_at(group) %>%
    mutate("overall_mean" := mean(!!rlang::sym(characteristic))) %>% # Remove and add appropriate grouping and calculate overall mean
    ungroup() %>%
    group_by_at(group_comm) %>%
    mutate("community_mean" := mean(!!rlang::sym(characteristic))) %>% # Remove and add appropriate grouping and calculate community mean
    mutate(
      "within_var" := (!!rlang::sym(characteristic) - community_mean)^2,
      "between_var" := (community_mean - overall_mean)^2
    ) %>% # Using individual data, and the overall and community mean, calculate squared within and between variation
    ungroup() %>%
    group_by_at(group) %>%
    mutate(
      within_var := 1 / (length(!!rlang::sym(characteristic)) - 1) * sum(within_var),
      between_var := 1 / (length(!!rlang::sym(characteristic)) - 1) * sum(between_var)
    ) %>% # Remove and add appropriate grouping and
    # calculate variance (1/n * sum(corresponding variations))
    select(-idalter, -!!characteristic) %>%
    unique() %>%
    ungroup() %>%
    group_by_at(group_result) %>% # Delete idalter and characteristic columns, summarize by grouping and group again to
    # deliver proper result tibble
    nest(community_mean := c(!!rlang::sym(algorithm), community_mean)) %>% # Nest the community data to facilitate visualization
    mutate("characteristic" := characteristic, "algorithm" := algorithm) # Add identification variables

  return(data)
}

#### FUNCTION WRAPPERS FOR USE INSIDE RECURSION

### Wrap and storage function of undirected_adjacency_matrix()
## Arguments:
# ... := Everything passed to undirected_adjacency_matrix()

wrapped_undirected_adjacency_matrix <- function(...) {
  args <- list(...)[c("variable", "school", "cohort", "time")] %>% as_tibble() #select proper arguments to pass to undirected_adjacency_matrix()
  print(args %>% as.character()) # Print information about process

  network <- undirected_adjacency_matrix(...) %>% #Set network
    graph_from_adjacency_matrix() %>%             #Create igraph object from adjacency matrix
    as.undirected() %>%                           #Specify undirectedness
    as_tbl_graph()                                #Transform it to tbl_graph object
  nodes_n <- network %>%                          #Calculate number of nodes for summarizing
    unlist(recursive = FALSE) %>%
    `[[`(1)
  summary_networks <- args %>% bind_cols(tibble(nodes_n = nodes_n, network_data = list(network))) #Storing data frame
  return(summary_networks)
}

### Use community detection algorithms and store output in a proper tibble
## Arguments:
# algorithm := Community detection Algorithm to be used
# network_data := tbl_graph object storing the network data

wrapped_communities <- function(algorithm, network_data) {
  tic(algorithm) # Time starts
  print(algorithm)
  data <- tryCatch(get(paste0("cluster_", algorithm))(network_data), error = list("error"))
  toc() # Time ends
  return(data)
}

### Extract the number of communities identified by algorithm
## Arguments:
# communities_data := communities object storing the communities data

wrapped_communities_n <- function(communities_data) {
  tryCatch(communities_data %>% pluck(membership) %>% max(), error = "error") # Extracting number of communities identified
}

### Count number of nominations for each id by groups
## Arguments:
# data := raw database
# characteristic := chracteristic to count (friendly, shy, leader, popular)
# time := survey data to use (suffixes as base or end)
# filtering := filtering variable at database (column name)
# groups := grouping options for counting

function_one <- function(data, characteristic, time, filtering, groups) {
  grouped_data <- data %>%
    group_by_at(groups) %>%
    count(wt = get(filtering), name = characteristic) %>%
    mutate(time = time)
}

### Count the number of diferent nominated ids above certain threshold by groups and characteristics
## Arguments:
# data := Data of the number of nominations for each id by groups and characteristics
# characteristic := chracteristic to count (friendly, shy, leader, popular)
# groups := grouping options for counting
# threshold := Minimum Threshold to count the data

function_two <- function(data, characteristic, groups, threshold) {
  filtered_data <- data %>%
    filter(get(characteristic) >= threshold) %>%
    group_by_at(groups) %>%
    count(name = paste0(characteristic, "_n"))
}

### Check repeated observations for a variable for some grouping
## Arguments:
# data := Raw database
# variable := Variable to check
# groups := grouping options for checking

check_repeated <- function(data, variable, group) {
  variable_n <- data %>%
    distinct(!!rlang::sym(variable)) %>%
    nrow()
  variable_grouped <- data %>% distinct(!!!rlang::syms(c(group, variable)))
  variable_grouped_n <- variable_grouped %>% nrow()

  repeated_values <- variable_grouped %>%
    group_by_at(variable) %>%
    filter(n() > 1)

  arguments <- tibble(group = group)
  repeated_values_grouped <- pmap(arguments, group_n_count, data = repeated_values, variable = variable)
  names(repeated_values_grouped) <- paste0("by_", group)
  repeated_values_grouped_n <- repeated_values_grouped %>% map_df(nrow)
  repeated_values_grouped %<>% reduce(full_join, by = variable) %>% list()
  
  result <- tibble(
    variable = variable,
    total_obs = variable_n,
    grouped_obs = variable_grouped_n,
    repeated_values = repeated_values
  ) %>%
    nest(repeated_values = starts_with("repeated_values")) %>%
    mutate(repeated_obs = map_dbl(repeated_values, nrow)) %>%
    bind_cols(repeated_values_grouped_n) %>%
    add_column(repeated_values_grouped)

  return(result)
}

### Group a database and count number of repetitions by a variable (group)
## Arguments:
# data := Raw database
# variable := Variable to count
# groups := grouping options for counting

group_n_count <- function(data, variable, group) {
  grouped_data <- data %>%
    distinct(!!!rlang::syms(c(group, variable))) %>%
    #group_by_at(c(group, variable)) %>%
    count() %>%
    filter(n > 1) %>%
    #group_by_at(variable) %>%
    #ungroup() %>%
    #select(-!!sym(group)) %>%
    rename(!!sym(paste0("by_", group)) := n) #%>%
    #summarise(!!sym(paste0("mean_by_", group)) := sum(n))
}
