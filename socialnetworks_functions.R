########################### NETWORKS FUNCTIONS ###########################
########################### Written by: Nicolas Fajardo
########################### Last modified: 25/03/2020

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

generate_undirected_adjacency_matrix_ <- function(data, school, cohort, variable, time, selector,
                                                  delete_missing = TRUE, replace_na = TRUE,
                                                  drop_singles = TRUE, force_diagonal = TRUE) {
  args <- c(selector, variable, school, cohort, time)
  cat("\n------", args %>% as.character() %>% str_to_upper(), "------ \n") # Print information about process
  
  # Create variable to filter using the function arguments
  if (time == "none" & selector == "none") {
    variable2filter <- c(variable) # If none, only the variable name
    # and the selector is concatenated
  } else if (time == "none") {
    variable2filter <- c(variable, selector) %>%
      paste(collapse = "_") %>%
      paste(collapse = "") # If none, only the variable name
    # and the selector is concatenated
  } else {
    variable2filter <- c(variable, time, selector) %>%
      paste(collapse = "_") %>%
      paste(collapse = "")
  }

  # Filter the desired data
  data_group <- dplyr::filter(data, coar == !!school, cohort == !!cohort) %>%
    dplyr::select(-c(coar, cohort)) %>%
    dplyr::select(idego, idalter, !!variable2filter) %>%
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
  data_group %<>% spread(key = idalter, value = !!variable2filter)

  # Perform NA replacement
  if (replace_na == TRUE) {
    data_group %<>% replace(is.na(.), 0)
  }

  # Renaming rows with each node id and converting the data into a matrix object
  data_group %<>% column_to_rownames("idego") %>% as.matrix()

  cat("\n --- Conform square matrix ---\n")

  # Identify which ids are not both in rows and columns (Conformability)
  intersection <- intersect(rownames(data_group), colnames(data_group))
  diffrows <- setdiff(rownames(data_group), intersection)
  diffcols <- setdiff(colnames(data_group), intersection)
  difference <- union(diffrows, diffcols)
  original_dim <- data_group %>% dim()

  # Delete the outlier ids
  if (!is_empty(intersection)) {
    cat(" Eliminate id(s) for ensuring complete pairwise relations\n")
    cat(" Number of row id(s) dropped [idego]:", diffrows %>% length(), "\n")
    cat(" Number of column id(s) dropped [idalter]:", diffcols %>% length(), "\n")
    cat(" Total number of dropped id(s):", difference %>% length(), "\n")
    data_group <- data_group[rownames(data_group) %in% intersection, colnames(data_group) %in% intersection]
    info_percentage <- original_dim %>%
      `/`(data_group %>% dim()) %>%
      `^`(-1) %>%
      `*`(100)
    names(info_percentage) <- c("idego", "idalter")
    cat(" Remaining percentage of data [idego] [idalter]: ", info_percentage, "\n --- \n\n")
  }

  # Perform singles dropping
  if (drop_singles == TRUE) {
    empty_rows <- rownames(data_group[rowSums(data_group) == 0, ])
    empty_cols <- colnames(data_group[, colSums(data_group) == 0])
    singles <- intersect(empty_rows, empty_cols)

    cat(" --- Drop Singles ---\n")
    cat(" Eliminate id(s) having no links\n")

    # If the all network is empty, exit the function
    # if (all((rownames(data_group) %in% singles)) & all((colnames(data_group) %in% singles))) {
    #  return(NULL)
    # }

    data_group <- data_group[!rownames(data_group) %in% singles, !colnames(data_group) %in% singles]
    cat(" Number of single edges dropped: ", singles %>% length(), "\n")
    info_percentage <- original_dim %>%
      `/`(data_group %>% dim()) %>%
      `^`(-1) %>%
      `*`(100)
    names(info_percentage) <- c("idego", "idalter")
    cat(" Remaining percentage of data [idego] [idalter]: ", info_percentage, "\n --- \n\n")
  }

  # CHECKS

  cat(" --- Adjacency Matrix Checks --- \n")
  check_symmetry <- all(data_group == t(data_group))
  cat(" Is adjacency Matrix symmetric? ", check_symmetry, "\n")
  check_diagonal <- all(diag(data_group) == 0)
  cat(" Are all the elements of the diagonal of the adjacency matrix equal to 0? ", check_diagonal, "\n")

  # Force all diagonal elements to be 0
  if (force_diagonal == TRUE) {
    corrected_diagonal <- diag(data_group)[diag(data_group) != 0]
    cat(" Number of corrected diagonal values: ", corrected_diagonal %>% length(), "\n --- \n\n")
    diag(data_group) <- 0
  }

  if (!exists("info_percentage")) {
    info_percentage <- 100
  }

  result <- tibble(
    variable = variable,
    selector = selector,
    school = school,
    cohort = cohort,
    time = time,
    network_data = data_group %>% list(),
    data_percentage = info_percentage %>% enframe(name = "variable", "percentage") %>% list(),
    drop_singles = drop_singles,
    dropped_values = if (exists("singles")) {
      union(difference, singles) %>% list()
    } else {
      difference %>% list()
    },
    matrix_symmetry = check_symmetry,
    forced_diagonal = force_diagonal,
    corrected_diagonal_values = if (exists("corrected_diagonal")) {
      corrected_diagonal %>%
        names() %>%
        list()
    } else {
      NULL %>% list()
    }
  )

  return(result)
}

#### GENERATE GRAPH FROM MATRIX
## Arguments:
# x := Matrix
# ... := Arguments to graph_from_adjacency_matrix(), as.undirected(), as_tbl_graph() [igraph]

pipe_matrix2graph <- function(x, ...) {
  if (is_null(x)) {
    x %<>% make_empty_graph(n = 0, directed = FALSE) %>%
      as_tbl_graph(...)
  } else {
    x %<>% graph_from_adjacency_matrix(...) %>% # Create igraph object from adjacency matrix
      as.undirected(...) %>% # Specify undirectedness
      as_tbl_graph(...) # Transform it to tbl_graph object
  }
}

#### ADD COMUNITIES TO GRAPH NODES
## Arguments:
# x := Graph object
# y := Community object

pipe_addcomm2graph <- function(x, y) {
  if (is_null(y) | is(y, "warning")) {
    return(x)
  }
  else {
    x %<>% activate(nodes) %>%
      mutate(community = as.integer(membership(y))) %>%
      mutate(community = as.factor(community))
    return(x)
  }
}

#### GENERATE TIBBLE FROM MATRIX
## Arguments:
# matrix := Matrix
# name := Name of values column
# row_id_name := Name of rownames column (idego)
# col_id_name := Name of colnames column (idalter)

pipe_matrix2tibble <- function(matrix, name, row_id_name = "idego", col_id_name = "idalter") {
  matrix %>% as_tibble(rownames = row_id_name) %>% pivot_longer(cols = !matches("idego"), names_to = col_id_name, values_to = name)
}

#### MISCELLANEOUS FUNCTIONS FOR SUMMARIZING

### Calculate membership and store it in a proper tibble
## Arguments:
# x := Dataset
# y := Value of variable

enframe_membership <- function(x, y, name) {
  if (is(x, "warning")) {
    NA
  }
  else if (is_null(x)) {
    NULL
  }
  else {
    membership(x) %>% enframe(name = name, value = as.character(y)) %>% mutate_if(is.numeric, as.integer)
  }
}

### Merge every membership's tibble
## Arguments:
# x := Tibble list

merge_n_wrap <- function(x, variable) {
  reduce(x, full_join, by = variable) %>% list()
}

### Calculate overall and community mean, within and between variance for one algorithm x one characteristic
## Arguments:
# x := Dataset
# algorithm := Algorithm to filter ("edge_betweenness", "fast_greedy", "label_prop", "infomap", "leading_eigen", "louvain", "walktrap")
# characteristic := Characteristic to filter ("friendly", "leader", "popular", "shy")
# threshold := Minimum counting threshold
# group := Grouping variables ("variable", "school", "cohort", "time") [Variable, School, Cohort, Time]

summarize_by_algorithm <- function(data, algorithm, variable, characteristic, threshold, group) {
  group_comm <- append(group, algorithm) # Group by pre-defined grouping variable plus algorithm (Operate by community)
  group_result <- c(group, "overall_mean", "within_var", "between_var") # Group by everything (To return a clean outcome)

  data %<>% dplyr::select(group, !!variable, !!algorithm, !!characteristic) %>% # Select appropiate columns
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
    dplyr::select(-!!variable, -!!characteristic) %>%
    unique() %>%
    ungroup() %>%
    group_by_at(group_result) %>% # Delete variable and characteristic columns, summarize by grouping and group again to
    # deliver proper result tibble
    nest(community_mean := c(!!rlang::sym(algorithm), community_mean)) %>% # Nest the community data to facilitate visualization
    mutate("characteristic" := characteristic, "algorithm" := algorithm) # Add identification variables

  return(data)
}

cross_ifself <- function(x) {
  x %<>% pull(idego)
  tidyr::crossing(student_a = x, student_b = x) %>% filter(student_a != student_b)
}

count_coincidences <- function(variable, data) {
  data %>%
    dplyr::select(coar, cohort, idego, idalter, !!sym(variable)) %>%
    filter(!!sym(variable) > 0) %>%
    dplyr::select(-!!sym(variable)) %>%
    nest(data = idego) %>%
    mutate(data = map(data, cross_ifself)) %>%
    unnest(cols = data) %>%
    group_by(coar, cohort, student_a, student_b) %>%
    count(name = paste0("coincidences_", variable))
}

mutate_by_algorithm <- function(name, data) {
  idego_sel <- paste(name, "idego", sep = ".")
  idalter_sel <- paste(name, "idalter", sep = ".")
  data %<>% mutate(!!paste0("same_", name) := ifelse(!!sym(idego_sel) == !!sym(idalter_sel), 1, 0)) %>%
    dplyr::select(coar, cohort, idego, idalter, contains("same"))
}

conform_matrices_distance2adj <- function(distance_matrices, adjacency_matrix) {
  distance_matrices %<>% map(~ .[rownames(.) %in% rownames(adjacency_matrix), colnames(.) %in% colnames(adjacency_matrix)])
  return(distance_matrices)
}

conform_matrices_adj2distance <- function(adjacency_matrix, distance_matrices) {
  columns <- distance_matrices %>%
    map(colnames) %>%
    reduce(intersect)
  rows <- distance_matrices %>%
    map(rownames) %>%
    reduce(intersect)
  adjacency_matrix <- adjacency_matrix[rownames(adjacency_matrix) %in% rows, colnames(adjacency_matrix) %in% columns]
  return(adjacency_matrix)
}

calculate_distance <- function(adjacency_matrix, distance) {
  distance_matrices <- distance$distance_adjacency_matrix
  check_distance_conformability <- distance_matrices %>%
    map(dim) %>%
    map(~ . >= dim(adjacency_matrix)) %>%
    unlist() %>%
    all()
  
  if (check_distance_conformability == TRUE) {
    distance_matrices %<>% conform_matrices_distance2adj(adjacency_matrix)
    adjacency_matrix %<>% conform_matrices_adj2distance(distance_matrices)
    distance_matrices %<>% map(~ . %*% adjacency_matrix)
    return(distance_matrices)
  } else if (check_distance_conformability == FALSE) {
    adjacency_matrix %<>% conform_matrices_adj2distance(distance_matrices)
    distance_matrices %<>% conform_matrices_distance2adj(adjacency_matrix)
    distance_matrices %<>% map(~ . %*% adjacency_matrix)
    return(distance_matrices)
  } else {
    return(list())
  }
}

#### FUNCTION WRAPPERS FOR USE INSIDE RECURSION

### Wrap and storage function of generate_undirected_adjacency_matrix()
## Arguments:
# ... := Everything passed to generate_undirected_adjacency_matrix()

generate_undirected_adjacency_matrix <- function(quietly = TRUE, matrix2graph = TRUE, ...) {
  if (quietly == TRUE) {
    sink("process_summary.txt", append = TRUE)
    if (matrix2graph == TRUE) {
    output <- generate_undirected_adjacency_matrix_(...) %>% # Generate adjacency matrix from data
      mutate(network_data = map(network_data, pipe_matrix2graph)) %>% # Set network
      mutate(nodes_n = map_dbl(network_data, ~ unlist(.x, recursive = FALSE) %>% `[[`(1))) 
    } else if (matrix2graph == FALSE) {
    output <- generate_undirected_adjacency_matrix_(...) %>%
      mutate(nodes_n = map_dbl(network_data, ~ unlist(.x, recursive = FALSE) %>% dim() %>% `[[`(1)))
    }
    sink()
  } else {
    if (matrix2graph == TRUE) {
      output <- generate_undirected_adjacency_matrix_(...) %>% # Generate adjacency matrix from data
        mutate(network_data = map(network_data, pipe_matrix2graph)) %>% # Set network
        mutate(nodes_n = map_dbl(network_data, ~ unlist(.x, recursive = FALSE) %>% `[[`(1))) 
    } else if (matrix2graph == FALSE) {
      output <- generate_undirected_adjacency_matrix_(...) %>%
        mutate(nodes_n = map_dbl(network_data, ~ unlist(.x, recursive = FALSE) %>% dim() %>% `[[`(1))) 
    }
  }
  return(output)
}

### Use community detection algorithms and store output in a proper tibble
## Arguments:
# algorithm := Community detection Algorithm to be used
# network_data := tbl_graph object storing the network data

wrap_communities <- function(algorithm, network_data, name) {
  sink("communities_summary", append = TRUE)
  cat(paste0("\n", name, "\n"))
  if ((network_data %>% gorder() == 0) & (network_data %>% gsize() == 0)) {
    data <- NULL
  }
  else {
    tic(algorithm) # Time starts
    data <- tryCatch(get(paste0("cluster_", algorithm))(network_data), error = function(e) {
      print(e)
      return(NULL)
    }, warning = function(w) {
      print(w)
    })
    toc() # Time ends
  }
  sink()
  return(data)
}

### Extract the number of communities identified by algorithm
## Arguments:
# communities_data := communities object storing the communities data

wrap_communities_n <- function(communities_data) {
  if (is(communities_data, "warning")) {
    NA
  }
  else if (is_null(communities_data)) {
    0L
  }
  else {
    tryCatch(communities_data %>% pluck(membership) %>% max() %>% as.integer, error = "error")
  } # Extracting number of communities identified
}

### Count number of nominations for each id by groups
## Arguments:
# data := raw database
# characteristic := chracteristic to count (friendly, shy, leader, popular)
# time := survey data to use (suffixes as base or end)
# filtering := filtering variable at database (column name)
# groups := grouping options for counting

count_nominations <- function(data, characteristic, time, groups) {
  name = paste(characteristic, time, sep = "_")
  grouped_data <- data %>%
    group_by_at(groups) %>%
    dplyr::count(wt = get(name), name = characteristic) %>%
    mutate(charc_time = time)
}

### Count the number of diferent nominated ids above certain threshold by groups and characteristics
## Arguments:
# data := Data of the number of nominations for each id by groups and characteristics
# characteristic := chracteristic to count (friendly, shy, leader, popular)
# groups := grouping options for counting
# threshold := Minimum Threshold to count the data

count_nominations_by_thld <- function(data, characteristic, groups, threshold) {
  filtered_data <- data %>%
    filter(get(characteristic) >= threshold) %>%
    group_by_at(groups) %>%
    dplyr::count(name = paste0(characteristic, "_n"))
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
    # group_by_at(c(group, variable)) %>%
    count() %>%
    filter(n > 1) %>%
    # group_by_at(variable) %>%
    # ungroup() %>%
    # select(-!!sym(group)) %>%
    rename(!!sym(paste0("by_", group)) := n) # %>%
  # summarise(!!sym(paste0("mean_by_", group)) := sum(n))
}

#### PLOTTING FUNCTIONS

### Plot Graph highlighting communities
## Arguments:
# name := Plot's name
# network_data := Graph object [igraph]
# algorithm := Name of algorithm used for community identification
# communities_n := Number of identified communities

plot_communities <- function(name, network_data, algorithm, communities_n) {
  if ((network_data %>% gorder() == 0) & (network_data %>% gsize() == 0)) {
    return(ggplot())
  }

  if (is.na(communities_n)) {
    return(ggplot())
  }
  else if (communities_n <= 0) {
    return(ggplot())
  }
  else if (communities_n > 1) {
    colors <- qualpal(communities_n)$hex
  }
  else {
    colors <- "4e79a7"
  }

  ggraph(network_data, layout = "stress") + ggtitle(
    label = name,
    subtitle = algorithm
  ) + labs(color = "Community") +
    geom_edge_link(color = "gray", edge_width = 0.01, edge_alpha = 0.1) +
    geom_node_point(aes(color = community), size = 1) + theme_graph() +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(override.aes = list(size = 1), nrow = 4)) +
    theme(
      text = element_text(size = 5, family = "DejaVu Sans Light"),
      legend.key.size = unit(0.25, "cm"),
      legend.position = "bottom",
      legend.justification = "left",
      legend.direction = "horizontal",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}

### Save multiple table_grob objects to multiple .pdf pages
## Arguments:
# name := .pdf file name
# table_grob := table_grob object [grobs]
# path := Name of the path to store the individual .pdf files
# ... := Arguments passed to ggsave() [ggplot]

ggsave_marrangeGrob <- function(name, table_grob, path = "Community_Plots", ...) {
  main_directory <- getwd()
  dir.create(file.path(main_directory, path))
  setwd(file.path(main_directory, path))
  filename <- paste0(name, ".pdf")
  ggsave(filename = filename, plot = table_grob, ...)
  setwd(main_directory)
}

### Plot centrality histograms of perceptions by communities
## Arguments:
# name := Plot's name
# network_data := Graph object [igraph]
# algorithm := Name of algorithm used for community identification
# communities_n := Number of identified communities

plot_perceptions_by_communities <- function(name, data) {
  ggplot(data, aes(deviation, color = community_algorithm, fill = community_algorithm)) +
    geom_density(alpha = 0.1) +
    facet_grid(rows = vars(characteristic), cols = vars(charc_time), scales = "free") +
    scale_color_manual(values = qualpal(7)$hex) +
    scale_fill_manual(values = qualpal(7)$hex) +
    labs(x = "Community", fill = "Community Algorithm", color = "Community Algorithm", title = name) +
    scale_y_continuous(breaks = extended_range_breaks(), labels = scales::number_format(accuracy = 0.01), expand = c(0, 0)) +
    scale_x_continuous(breaks = extended_range_breaks(), labels = scales::number_format(accuracy = 0.01), expand = c(0, 0)) +
    theme_light() +
    theme(
      panel.spacing.x = unit(1.5, "lines"),
      panel.spacing.y = unit(1.5, "lines"),
      axis.title = element_blank(),
      text = element_text(size = 8, family = "DejaVu Sans Light"),
      legend.key.size = unit(0.25, "cm"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.position = "bottom",
      legend.justification = "left",
      legend.direction = "horizontal"
    )
}

load_lastdataset <- function(name) {
  if (!exists(name)) {
    position <- file.info(list.files() %>% str_subset(name)) %>% `[[` ("mtime") %>% which.max()
    full_name <- list.files() %>% str_subset(name) %>% `[[` (position)
    load(full_name, envir = globalenv())
  }
}

load_datasrc <- function(src_name, main = FALSE) {
  if (!exists("loaded_data") | main == FALSE) {
    source(src_name) } else if (loaded_data != TRUE) {
      source(src_name) }
}
