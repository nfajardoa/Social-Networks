########################### NETWORKS FUNCTIONS ###########################
########################### Written by: Nicolas Fajardo
########################### Last modified: 24/06/2020
########################### Comments: This document groups the functions based on the module
###########################           in which they are used. If two or more modules use the
###########################           function, it will appear on the first module which uses it
###########################           according to the order displayed in "socialnetworks.R"

############## MAIN MODULES

###### MODULE: Load data, Summary Statistics & Database check (summstats)

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
  data %>%
    distinct(!!!rlang::syms(c(group, variable))) %>%
    dplyr::count() %>%
    filter(n > 1) %>%
    rename(!!sym(paste0("by_", group)) := n)
}

###### MODULE: Identify Communities (communities)

### Wrap and storage function of generate_undirected_adjacency_matrix()
## Arguments:
# ... := Everything passed to generate_undirected_adjacency_matrix()

generate_undirected_adjacency_matrix <- function(quietly = TRUE, matrix2graph = TRUE, ...) {
  if (quietly == TRUE) {
    sink("adjacency_matrices_summary.txt", append = TRUE)
    if (matrix2graph == TRUE) {
      output <- generate_undirected_adjacency_matrix_(...) %>% # Generate adjacency matrix from data
        mutate(network_data = map(network_data, pipe_matrix2graph)) %>% # Set network
        mutate(nodes_n = map_dbl(network_data, ~ unlist(.x, recursive = FALSE) %>% `[[`(1)))
    } else if (matrix2graph == FALSE) {
      output <- generate_undirected_adjacency_matrix_(...) %>%
        mutate(nodes_n = map_dbl(network_data, ~ unlist(.x, recursive = FALSE) %>%
          dim() %>%
          `[[`(1)))
    }
    sink()
  } else {
    if (matrix2graph == TRUE) {
      output <- generate_undirected_adjacency_matrix_(...) %>% # Generate adjacency matrix from data
        mutate(network_data = map(network_data, pipe_matrix2graph)) %>% # Set network
        mutate(nodes_n = map_dbl(network_data, ~ unlist(.x, recursive = FALSE) %>% `[[`(1)))
    } else if (matrix2graph == FALSE) {
      output <- generate_undirected_adjacency_matrix_(...) %>%
        mutate(nodes_n = map_dbl(network_data, ~ unlist(.x, recursive = FALSE) %>%
          dim() %>%
          `[[`(1)))
    }
  }
  return(output)
}

### Generate undirected adjacency matrix for the Dataset
## Arguments:
# data := Dataset
# coar := School to filter
# cohort := Cohort to filter
# variable := Variable to filter ("friend", "study" ,"social", "any")
# time := Time to filter ("none", "base", "end", "mid")
# delete_missing := Delete missing observations
# selector := Which columns must be used ("mut", "any")
# replace_na := Replace missing values with 0s
# drop_singles := Drop unconnected nodes
# force_diagonal := Delete loops inside data [Setting all diagonal elements to 0]

generate_undirected_adjacency_matrix_ <- function(data, coar, cohort, variable, time, selector,
                                                  delete_missing = TRUE, replace_na = TRUE,
                                                  drop_singles = TRUE, force_diagonal = TRUE) {
  args <- c(selector, variable, coar %>% as.character(), cohort %>% as.character(), time)
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
  data_group <- dplyr::filter(data, coar == !!coar, cohort == !!cohort) %>%
    dplyr::select(-c(coar, cohort)) %>%
    dplyr::select(idego, idalter, !!variable2filter) %>%
    dplyr::arrange(idego, idalter)

  # If no filtered data exists, exit the function
  if (is_empty(data_group)) {
    return(NULL)
  }

  # Perform missing deletion
  if (delete_missing == TRUE) {
    data_group %<>% drop_na()
  }

  # Reshaping data frame to wide form (almost as a matrix)
  data_group %<>% tidyr::pivot_wider(id_cols = idego, names_from = idalter, values_from = !!sym(variable2filter))

  # Perform NA replacement
  if (replace_na == TRUE) {
    data_group %<>% replace(is.na(.), 0)
  }

  # Renaming rows with each node id and converting the data into a matrix object
  data_group %<>% tibble::column_to_rownames("idego") %>% as.matrix()

  cat("\n --- Conform square matrix ---\n")

  # Identify which ids are not both in rows and columns (Conformability)
  intersection <- dplyr::intersect(rownames(data_group), colnames(data_group))
  diffrows <- dplyr::setdiff(rownames(data_group), intersection)
  diffcols <- dplyr::setdiff(colnames(data_group), intersection)
  difference <- dplyr::union(diffrows, diffcols)
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
    coar = coar,
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

#### Generate graph from matrix
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

### Use community detection algorithms and store output in a proper tibble
## Arguments:
# algorithm := Community detection Algorithm to be used
# network_data := tbl_graph object storing the network data

wrap_communities <- function(algorithm, network_data, name) {
  sink("communities_summary", append = TRUE)
  cat(paste0("\n", name, "\n"))
  if ((network_data %>% igraph::gorder() == 0) & (network_data %>% igraph::gsize() == 0)) {
    data <- NULL
  }
  else {
    tic(algorithm) # Time starts
    data <- tryCatch(get(paste0("cluster_", algorithm), envir = asNamespace("igraph"))(network_data), error = function(e) {
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
    tryCatch(communities_data %>% pluck(membership) %>% max() %>% as.integer(), error = "error")
  } # Extracting number of communities identified
}

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
    membership(x) %>%
      enframe(name = name, value = as.character(y)) %>%
      mutate_if(is.numeric, as.integer)
  }
}

#### Add communities information to graph nodes
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
# group := Grouping variables ("variable", "coar", "cohort", "time") [Variable, coar, Cohort, Time]

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

### Calculate overall and community mean, within and between variance for one algorithm x one characteristic
## Arguments:
# x := Dataset
# algorithm := Algorithm to filter ("edge_betweenness", "fast_greedy", "label_prop", "infomap", "leading_eigen", "louvain", "walktrap")
# characteristic := Characteristic to filter ("friendly", "leader", "popular", "shy")
# threshold := Minimum counting threshold
# group := Grouping variables ("variable", "coar", "cohort", "time") [Variable, coar, Cohort, Time]

mutate_by_algorithm <- function(name, data) {
  idego_sel <- paste(name, "idego", sep = ".")
  idalter_sel <- paste(name, "idalter", sep = ".")
  data %<>% mutate(!!paste0("s_", name) := ifelse(!!sym(idego_sel) == !!sym(idalter_sel), 1, 0)) %>%
    dplyr::select(coar, cohort, idego, idalter, matches("^s_"))
}

###### MODULE: Calculate Weighted Distance (wdistance)

#### Generate Tibble from Matrix
## Arguments:
# matrix := Matrix
# name := Name of values column
# row_id_name := Name of rownames column (idego)
# col_id_name := Name of colnames column (idalter)

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

#### Generate Tibble from Matrix
## Arguments:
# matrix := Matrix
# name := Name of values column
# row_id_name := Name of rownames column (idego)
# col_id_name := Name of colnames column (idalter)

pipe_matrix2tibble <- function(matrix, name, row_id_name = "idego", col_id_name = "idalter") {
  matrix %>%
    as_tibble(rownames = row_id_name) %>%
    pivot_longer(cols = !matches("idego"), names_to = col_id_name, values_to = name)
}

#### WARNING: UGLY FUNCTION ####
#### Generate Tibble from Matrix
## Arguments:
# matrix := Matrix
# name := Name of values column
# row_id_name := Name of rownames column (idego)
# col_id_name := Name of colnames column (idalter)

merge_groups <- function(data, name, id_cases = tibble(), groups = NULL, .use_database) {
  data %<>% mutate(wdistance = map2(data$wdistance, data$adjacency_matrix, pm2t_join, name = name)) %>%
    dplyr::select(-adjacency_matrix) %>%
    unnest(cols = wdistance) %>%
    rename_at(vars(matches("contiguo")), ~ paste0("w", .x)) %>%
    mutate_if(is.factor, as.character) %>%
    right_join(rlang::eval_tidy(id_cases), by = groups) %>%
    mutate(
      network = paste(coar, cohort, sep = "_"),
      network_id = as.integer(factor(network, levels = unique(network)))
    )
  if (.use_database == TRUE) {
    dbWriteTable(conn = database_connection, name = name, value = data, overwrite = TRUE, temporary = FALSE)
    link <- tbl(database, name)
    return(link)
  } else if (.use_database == FALSE) {
    return(data)
  } else {
    (stop("No condition was set"))
  }
}

###### MODULE: Regressions (regressions)

multiple_reformulate <- function(formulas, response) {
  formulas %<>% as_tibble_col(column_name = "model_formula") %>%
    mutate(
      model_formula = map(model_formula, ~ reformulate(.x, response = response)),
      model_type = 1:n()
    )
}

sink_latex_output <- function(object, column = "model", filename, append = TRUE, ...) {
  if (missing(filename)) {
    filename <- as.character(object)
  }
  sink(filename, append = append, ...)
  object %>%
    `[[`(column) %>%
    stargazer()
  closeAllConnections()
}

add_baseline <- function(data, formula, id_cols, baseline_data) {
  response <- terms(formula) %>%
    as.character() %>%
    `[[`(2)
  if (str_detect(response, "_base")) {
    return(data)
  } else {
    response <- if_else(response %>% str_detect("(_mid)|(_end)"),
      response %>% str_replace("(_mid)|(_end)", "_base"),
      response %>% str_replace("_", "_base_")
    )
    baseline_data %<>% dplyr::select(all_of(id_cols), !!sym(response))
    return(collect(data) %>% right_join(baseline_data, by = id_cols))
  }
}

update_formula <- function(formula, baseline = TRUE) {
  response <- terms(formula) %>%
    as.character() %>%
    `[[`(2)
  if (baseline == TRUE) {
    if (str_detect(response, "_base")) {
      return(formula)
    } else {
      terms <- terms(formula) %>%
        as.character() %>%
        `[[`(3) %>%
        paste0(
          if_else(response %>% str_detect("(_mid)|(_end)"),
            response %>% str_replace("(_mid)|(_end)", "_base"),
            response %>% str_replace("_", "_base_")
          ),
          " + ",
          .
        )
      return(reformulate(terms, response = response))
    }
  } else {
    warning("No formula was updated (No baseline control)")
    return(formula)
  }
}

regress_network <- function(regressions_object, formulas,
                            response_col = "name", engine = "felm",
                            nested_data = NULL,
                            use_data, add_data,
                            baseline_data,
                            control_baseline = FALSE,
                            id_cols = c("coar", "cohort", "idego", "idalter"),
                            ...) {
  regressions_object %<>%
    mutate(models = map(!!sym(response_col), multiple_reformulate, formulas = formulas))
  if (!is.null(nested_data)) {
    regressions_object %>%
      mutate(models = map2(models, !!sym(nested_data), function(models, data) {
        mutate(models,
          model_formula = map(model_formula, update_formula, baseline = control_baseline),
          model = map(model_formula, function(model_formula) {
            felm(
              formula = model_formula,
              data = join_base_n_data(data,
                add_data,
                id_cols,
                formula = model_formula,
                baseline_data = baseline_data
              ),
              ...
            )
          })
        )
      })) %>%
      dplyr::select(-!!sym(nested_data))
  } else {
    if (missing(use_data)) {
      stop("No data provided for regressions")
    } else {
      unnest(regressions_object, models) %>%
        mutate(
          model_formula = map(model_formula, update_formula, baseline = control_baseline),
          model = map(model_formula, felm, data = use_data, ...)
        )
    }
  }
}

join_base_n_data <- function(data, add_data, id_cols, ...) {
  data <- collect(data) %>%
    add_baseline(id_cols = id_cols, ...)

  if (missing(add_data)) {
    return(data)
  } else if (!missing(id_cols)) {
    return(data %>% right_join(add_data, by = id_cols))
  } else {
    stop("No id columns for adding data were provided")
  }
}

###### MODULE: Heterogeneous Effects (heteffects)

# This function unifies splitting, estimation and inference frameworks using all specified methods
# seed := seed to fix the split
# data := data frame containing relevant variables
# stratum := character vector of variable to perform the statified split
# folds := number of folds if a tuning operation is performed (otherwise is ignored)
# repetitions := number of repetitions if a tuning operation is performed (otherwise is ignored)
# groups := number of desired groups (could be automatically lowered depending on data density)
# affected_quantile := lowest and highest quantile to perform CLAN analysis
# significance := significance level to calculate the confidence intervals
# save_splits := return data object with the initial split? (FALSE)
# save_model := return regression object with the result? (FALSE)
# save_data := return data object with the result? (FALSE)
# save_control := return control with the result? (TRUE)
# save_parameters := return parameters with the result? (TRUE)
# ... := options to stratified split

analise_heterogeneity <- function(seed, data, stratum, control, folds,
                                  repetitions = 2, groups, affected_quantile, significance,
                                  save_splits = FALSE, save_model = FALSE,
                                  save_data = FALSE, save_parameters = TRUE, ...) {

  # Generate Split
  split <- stratified_split(data, stratum, seed, ...)
  # Add Controls (Formula and keep only used variables)
  control %<>% generate_controls.ml()

  # Define Globals
  global_objects <- create_named_list(
    estimate_treatment_effect, select_split, run_workflow, cut_quantiles, expand_quantiles,
    check_n_jitter, create_named_list, regress_blp, regress_gate, regress_clan, estimate_heteffects,
    calculate_accuracy, correct_pvalue, calculate_correct_pvalue, generate_coefficients, generate_hettest,
    analysis_clan
  )
  cluster_packages <- c(
    "tidyverse", "magrittr", "broom", "lfe", "multcomp", "modelr", "tidymodels", "ranger", "glmnet",
    "dplyr", "parsnip", "recipes", "tune", "yardstick", "fastDummies"
  )

  # Directory to save
  dir <- file.path(getwd(), "heteffects_simulations")
  dir.create(dir, showWarnings = FALSE)

  # Options' Hash
  H <- new.env()
  H[["stratum"]] <- stratum
  H[["control"]] <- control
  H[["folds"]] <- folds
  H[["repetitions"]] <- repetitions
  H[["groups"]] <- groups
  H[["affected_quantile"]] <- affected_quantile
  H[["significance"]] <- significance

  # Save Objects
  hash <- .robustDigest(H)
  rm(data, H)

  # Run Analysis
  analysis_result <- future_pmap_dfr(
    .l = control,
    .f = analise_heterogeneity_,
    seed = seed,
    hash = hash,
    split = split,
    folds = folds,
    repetitions = repetitions,
    groups = groups,
    affected_quantile = affected_quantile,
    significance = significance,
    save_model = save_model,
    save_data = save_data,
    save_parameters = save_parameters,
    .options = future_options(
      packages = cluster_packages,
      scheduling = Inf,
      globals = global_objects,
      seed = TRUE
    )
  )

  # # Wrap-up results
  # final_result <- tibble(
  #   seed = seed,
  #   heterogeneity_analysis = analysis_result %>% list(),
  # )
  #
  # if (save_splits == TRUE) {
  #   final_result %<>% add_column(.split = split %>% list())
  # }

  rm(list = ls() %>% str_subset("^analysis_result$", negate = TRUE), envir = environment())
  gc(reset = TRUE, verbose = FALSE)
  return(analysis_result)
}

# This function makes the initial split by strata and saves the result as a split object (see rsample).
# data := data frame to be splitted
# stratum := variable containing the strata
# seed := seed to fix the random selection of observations
# ... := more options to createDataPartition()

stratified_split <- function(data, stratum, seed, ...) {
  if (!is_tibble(data) & !is.matrix(data)) {
    stop("`data` must be a tibble.", call. = FALSE)
  }

  set.seed(seed)

  in_id <- data %>%
    dplyr::select(all_of(stratum)) %>%
    pull() %>%
    createDataPartition(list = FALSE, ...) %>%
    as.integer()

  out_id <- setdiff(1:nrow(data), in_id)

  if (!is.integer(in_id) | any(in_id < 1)) {
    stop("`in_id` must be a positive integer vector.", call. = FALSE)
  }

  if (!all(is.na(out_id))) {
    if (!is.integer(out_id) | any(out_id < 1)) {
      stop("`out_id` must be a positive integer vector.", call. = FALSE)
    }
  }

  if (length(in_id) == 0) {
    stop("At least one row should be selected for the analysis set.",
      call. = FALSE
    )
  }

  structure(
    list(
      data = data,
      in_id = in_id,
      out_id = out_id
    ),
    class = "rsplit"
  )
}

# This function generates a formula to calculate the ML model, given a set of controls
# controls := analysis_parameters data frame

generate_controls.ml <- function(controls) {
  controls %<>%
    mutate(
      use_vars = pmap(list(response, controls, treatment, cluster, fixed_effect, affected), c),
      use_vars = map(use_vars, str_subset, "^0$", negate = TRUE),
      ml_formula = map2(controls, response, ~ reformulate(paste(.x, collapse = " + "), .y))
    )
}

# This function unifies both estimation and inference frameworks given a method
# method_name := character vector of a ML method name
# ... := options to estimate_treatment_effect() and estimate_heteffects(...)

analise_heterogeneity_ <- function(..., hash, seed, response, treatment, method_name) {

  # Execute Framework
  result <- estimate_treatment_effect(treatment = treatment, ...) %>%
    check_n_jitter() %>%
    estimate_heteffects(response = response, ...) %>%
    add_column(
      simulation = seed,
      method_name = all_of(method_name),
      hash %>% as_tibble_row(),
      .before = 1
    )

  save(result,
    file = paste0(
      "heteffects_simulations/",
      response, "-",
      treatment, "-",
      method_name, "-",
      hash$control, "-",
      as.numeric(Sys.time()) %>% `*`(1000) %>% round(0), ".RData"
    )
  )

  rm(list = ls() %>% str_subset("^seed$", negate = TRUE), envir = environment()) # Save memory
  gc(reset = TRUE, verbose = FALSE) # Collect garbage
  tibble(seed = seed, result = TRUE, time = Sys.time())
}

# This function runs the estimation framework to estimate the sample treatment effect (by simulation)
# split := split object containing the data (training/test) or (analysis/assessment)
# treatment := character of treatment variable name
# use_vars := variables used in the estimation (to reduce memory usage)
# groups := number of groups to estimate heterogeneity (quantiles)
# affected_quantile := lowest and highest quantile to perform CLAN analysis
# ... := options to run_workflow()

estimate_treatment_effect <- function(split, treatment, use_vars, groups = 5, affected_quantile = 0.2, ...) {
  if (!missing(use_vars)) split %<>% select_split(use_vars) # Keep only necessary columns from split
  assessment_data <- assessment(split) # Save assessment data
  propensity_score <- (split %>% analysis() %>% filter(!!sym(treatment) == 1) %>% nrow() +
    split %>% assessment() %>% filter(!!sym(treatment) == 1) %>% nrow()) /
    nrow(split$data) # Calculation of propensity score from sample

  yz0.fit <- run_workflow(
    data_analysis = rsample::analysis(split) %>% filter(!!sym(treatment) == 0),
    ...
  ) # Estimate B(Z)

  yz1.fit <- run_workflow(
    data_analysis = rsample::analysis(split) %>% filter(!!sym(treatment) == 1),
    ...
  ) # Estimate S(Z) + B(Z)

  # Predict proxies
  B <- predict(yz0.fit, new_data = assessment_data) %>% pull()
  S <- predict(yz1.fit, new_data = assessment_data) %>%
    pull() %>%
    `-`(B)

  if ((assessment_data %>% nrow()) > B %>% length()) {
    assessment_data %>% drop_na()
  }

  assessment_data %<>%
    add_column(
      B = all_of(B), # Save B(Z) by observation
      S = all_of(S) # Save S(Z) by observation
    ) %>%
    mutate(
      D_ort := as.numeric(!!sym(treatment)) - all_of(propensity_score), # Calculate the Orthogonalized D values
      S_ort := D_ort * (S - mean(S)), # Calculate the Orthogonalized S values
      S_group := cut_quantiles(S + runif(1, 0, 0.00001), probs = seq(0, 1, length.out = groups + 1)), # Calculate groups of S(Z)
      model_matrix(., ~ -1 + S_group), # Transform the groups into dummies
      most_affected := if_else(S > quantile(S, 1 - affected_quantile), 1, 0), # Mark the sample least affected units
      least_affected := if_else(S < quantile(S, affected_quantile), 1, 0), # Mark the sample least affected units
      weight = all_of(as.numeric((propensity_score * (1 - propensity_score))**-1)) # Calculate weights to regress
    ) %>%
    mutate(across(matches("S_groupQ"), ~ .x * D_ort)) # Orthogonalizes group dummies relative to all other regressors that are functions of Z

  attr(assessment_data, "treatment_var") <- treatment
  rm(list = ls() %>% str_subset("^assessment_data$", negate = TRUE), envir = environment())
  gc(reset = TRUE, verbose = FALSE)

  return(assessment_data)
}

# This function select some variables (vars) in a split object (see rsample).
# split := split object
# vars := character vector of variable names

select_split <- function(split, vars) {
  vars %<>% str_subset("^0$", negate = TRUE)
  split$data %<>%
    dplyr::select(!!!syms(vars)) %>%
    drop_na(!!!syms(vars))
  split
}

# This function runs the workflow of the Machine Learning environment (tidymodels)
# data_analysis := data to estimate the model
# ml_formula := formula specifying the predictors and the target variable
# folds := number of folds if a tuning operation is performed (otherwise is ignored)
# repetitions := number of repetitions if a tuning operation is performed (otherwise is ignored)
# recipe := recipe CALL to pre-process the data (see recipes)
# model := model CALL to lay the model out (see parsnip)
# tuner := tuner CALL to tune the model (see parsnip)

run_workflow <- function(data_analysis, ml_formula, folds, repetitions,
                         recipe, model, tuner, ...) {

  # Transform data into data frame
  data_analysis %<>% as.data.frame()

  ## Evaluate the quoted arguments
  recipe <- recipes::recipe(ml_formula, data = data_analysis) %>% # Pre-processing Recipe
    step_meanimpute(all_numeric()) %>% # Add processing steps
    step_modeimpute(all_nominal()) %>%
    step_zv(all_predictors()) %>%
    step_dummy(all_predictors(), -all_numeric()) %>%
    step_normalize(all_predictors())

  model <- rlang::eval_tidy(model) # Model

  ## Generate Workflow (Pre-processing, Model)
  process <-
    workflows::workflow() %>%
    workflows::add_recipe(recipe) %>%
    workflows::add_model(model)

  ## Generate Cross-Validation Samples
  # if (!missing(tuner)) {
  #  data_cv <- rsample::vfold_cv(data_analysis, repeats = repetitions, v = folds)
  #  tuner <- rlang::eval_tidy(tuner) # Tuner
  # }

  gc(reset = TRUE, verbose = FALSE)

  if (missing(tuner)) {
    return(fit(process, data = data_analysis)) # Return the fit object
  } else if (all(is.na(tuner))) {
    return(fit(process, data = data_analysis)) # Return the fit object
  } else {
    data_cv <- rsample::vfold_cv(data_analysis, repeats = repetitions, v = folds)
    tuner <- process %>%
      tune::tune_grid(
        resamples = data_cv, # data_cv is the resampling object by default, don't change it
        grid = 50,
        metrics = yardstick::metric_set(yardstick::rmse)
      )
    best_model <- tuner %>% select_best()
    return(finalize_workflow(process, best_model) %>% fit(data = data_analysis)) # Return the fit object with the best model by tuning
  }
}

# This function expands the calculated quantiles of a numeric vector
# ... := options to quantile()
# amount := amount to expand the generated breaks

expand_quantiles <- function(..., amount = 0.001) {
  breaks <- quantile(..., include.lowest = TRUE, na.rm = TRUE)
  breaks[1] <- breaks[1] - amount
  breaks[length(breaks)] <- breaks[length(breaks)] + amount
  breaks
}

# This function cut the quantiles from a variable and automatically
# lowers the group numbers if two quantiles have the same density
# variable := variable to calculate and cut the quantiles
# ... := options to expand_quantiles()

cut_quantiles <- function(variable, ...) {
  breaks <- expand_quantiles(x = variable, ...)
  labels <- 1:(length(breaks) - 1) %>% paste0("Q", .)
  tryCatch(
    {
      cut(
        x = variable,
        breaks = breaks, labels = labels
      )
    },
    error = function(e) {
      breaks <- expand_quantiles(x = variable, ...) %>% unique()
      labels <- 1:(length(breaks) - 1) %>% paste0("Q", .)
      warning(paste("Decreasing number of groups to", length(breaks) - 1))
      cut(
        x = variable,
        breaks = breaks, labels = labels
      )
    }
  )
}

# This function adds the jittering as prefered by Chernozhukov et al. (2018) to solve weak identification
# data := data frame from estimate_treatment_effect()

check_n_jitter <- function(data) {
  if (var(data$B) == 0) { # Variance of Baseline
    data %<>% mutate(B = B + rnorm(1, mean = 0, sd = 0.1))
  }
  if (var(data$S) == 0) { # Variance of Treatment Effect
    data %<>% mutate(S = S + rnorm(1, mean = 0, sd = 0.1))
  }
  if (var(data$most_affected) == 0) { # Variance of Most Affected Observations
    data %<>% mutate(most_affected = n() %>% runif() %>% `<`(0.1) %>% as.numeric())
  }
  if (var(data$least_affected) == 0) { # Variance of Least Affected Observations
    data %<>% mutate(least_affected = n() %>% runif() %>% `<`(0.1) %>% as.numeric())
  }

  rm(list = ls() %>% str_subset("^data$", negate = TRUE), envir = environment())
  gc(reset = TRUE, verbose = FALSE)

  return(data)
}

# This function unifies the Chernozhukov et al. (2018) framework, calculating BLP, GATES and CLAN on
# a assessment data frame from estimate_treatment_effect()
# assessment_data := assessment data frame from estimate_treatment_effect()
# response := character vector of response variable (Y)
# fixed_effect := character of absorbed fixed effects, separated by "+", as used for a felm() estimation (see lfe)
# cluster := character of cluster variables, separated by "+", as used for a felm() estimation (see lfe)
# affected := length of quantile to calculate most and least affected units
# significance := significance level to calculate the confidence intervals
# save_data := return data object with the result? (FALSE)
# save_parameters := return parameters with the result? (TRUE)
# ... := non-supported

estimate_heteffects <- function(assessment_data, response, fixed_effect,
                                cluster, affected, significance, save_model,
                                save_data = FALSE, save_parameters = TRUE, ...) {
  result <- bind_cols(
    regress_blp(
      response = response,
      fixed_effect = fixed_effect,
      cluster = cluster,
      significance = significance,
      data = assessment_data,
      save_model = save_model
    ),
    regress_gate(
      response = response,
      fixed_effect = fixed_effect,
      cluster = cluster,
      significance = significance,
      data = assessment_data,
      save_model = save_model
    )
  ) %>% add_column(
    clan = regress_clan(
      variables = affected,
      data = assessment_data,
      significance = significance,
      save_model = save_model
    ) %>% list()
  )

  treatment <- attr(assessment_data, "treatment_var")
  observations <- nrow(assessment_data)

  if (save_data == TRUE) {
    result %<>% mutate(assessment_data = assessment_data %>% list(), .before = 1)
  }
  if (save_parameters == TRUE) {
    parameters <- create_named_list(response, treatment, fixed_effect, cluster, observations) %>%
      as_tibble_row()
    result %<>% add_column(parameters, .before = 1)
  }

  rm(list = ls() %>% str_subset("^result$", negate = TRUE), envir = environment())
  gc(reset = TRUE, verbose = FALSE)

  return(result)
}

# Helper function to store data
create_named_list <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}

# This function performs the BLP regression specified by Chernozhukov et al. (2018) as the first approach to
# identify the coefficients from a weighted linear projection
# response := character of response variable (Y)
# fixed_effect := character of absorbed fixed effects, separated by "+", as used for a felm() estimation (see lfe)
# cluster := character of cluster variables, separated by "+", as used for a felm() estimation (see lfe)
# significance := significance level to calculate the confidence intervals
# data := data frame from estimate_treatment_effect()
# save_model := return regression object with the result? (FALSE)

regress_blp <- function(response, fixed_effect, cluster, significance, data, save_model = FALSE) {

  # Y = alpha_0 + alpha_B*B + beta1*(D − p(Z)) + beta2*(D − p(Z))(S − ES) + epsilon
  # S was added to the code as in https://github.com/demirermert/MLInference/blob/master/Heterogeneity/EL1.R
  blp_formula <- as.formula(paste(response, "~", "B + S + D_ort + S_ort | ", fixed_effect, "| 0 |", cluster, sep = ""))
  blp_regression <- tryCatch(
    {
      felm(blp_formula, data = data, weights = data$weight)
    },
    error = function(e) {
      cat(" Error in Best Linear Predictor Regression \n Regressing only orthogonalized values \n")
      blp_formula <- as.formula(paste(response, "~", "D_ort + S_ort| ", fixed_effect, "| 0 |", cluster, sep = ""))
      blp_regression <- felm(blp_formula, data = data, weights = data$weight)
      return(blp_regression)
    },
    warning = function(w) {
      cat(" Warning in Best Linear Predictor Regression \n Regressing only orthogonalized values \n")
      blp_formula <- as.formula(paste(response, "~", "D_ort + S_ort| ", fixed_effect, "| 0 |", cluster, sep = ""))
      blp_regression <- felm(blp_formula, data = data, weights = data$weight)
      return(blp_regression)
    }
  )

  # Save the coefficients in a data frame (a la broom)
  coefficients <- generate_coefficients(blp_regression, significance)
  # Save the BLP accuracy measure
  accuracy <- calculate_accuracy("S_ort", blp_regression, data = data)

  result <- tibble(
    blp_coefficients = coefficients %>% list(),
    blp_accuracy = accuracy
  )

  if (save_model == TRUE) {
    result %<>% add_column(blp_regression = blp_regression %>% list(), .before = 1)
  }

  rm(list = ls() %>% str_subset("^result$", negate = TRUE), envir = environment()) # Save memory
  gc(reset = TRUE, verbose = FALSE) # Collect garbage

  return(result)
}

# This function calculates the accuracy of the Machine Learning Methods (Best BLP)
# as proposed by Chernozhukov et al. (2018), maximizing the correlation between the ML proxy predictor S(Z)
# and the true score s0(Z), or equivalent to maximizing the R^2 in the regression of s0 (Z) on S(Z)
# name := name of the variable to calculate (S_ort)
# regression := regression object
# data := data frame to calculate the sample variance of S

calculate_accuracy <- function(name, regression, data) {
  estimate <- regression %>%
    tidy() %>%
    filter(term == name) %>%
    pull(estimate)
  term <- str_remove(name, "_ort")
  variance <- data %>%
    pull(!!sym(term)) %>%
    var()
  return(abs(estimate) * sqrt(variance))
}

# This function performs the GATES regression specified by Chernozhukov et al. (2018) as the first approach to
# identify the coefficients from a weighted linear projection
# response := character of response variable (Y)
# fixed_effect := character of absorbed fixed effects, separated by "+", as used for a felm() estimation (see lfe)
# cluster := character of cluster variables, separated by "+", as used for a felm() estimation (see lfe)
# significance := significance level to calculate the confidence intervals
# data := data frame from estimate_treatment_effect()
# save_model := return regression object with the result? (FALSE)

regress_gate <- function(response, fixed_effect, cluster, significance, data, save_model = FALSE) {
  # Save group dummies name
  group_variables <- names(data) %>% str_subset("S_groupQ")
  # Find the number of groups
  groups <- group_variables %>%
    str_remove_all("S_groupQ") %>%
    as.numeric() %>%
    max()

  # Y = alpha_0 + alpha_b*B + alpha_s*S + sum(k=1, K)gamma_k*(D − p(Z))*1(G_k) + ν, k is the number of groups and 1(G_k) is the group dummy
  gate_formula <- as.formula(paste(response, "~", "B + S + ", paste(group_variables, collapse = " + "), " | ", fixed_effect, "| 0 |", cluster, sep = ""))

  gate_regression <- tryCatch(
    {
      felm(gate_formula, data = data, weights = data$weight)
    },
    error = function(e) {
      cat(" Error in GATE Regression \n Regressing only group values \n")
      gate_formula <- as.formula(paste(response, "~", paste(group_variables, collapse = " + "), " | ", fixed_effect, "| 0 |", cluster, sep = ""))
      gate_regression <- felm(gate_formula, data = data, weights = data$weight)
      return(gate_regression)
    },
    warning = function(w) {
      cat(" Warning in GATE Regression \n Regressing only group values \n")
      gate_formula <- as.formula(paste(response, "~", paste(group_variables, collapse = " + "), " | ", fixed_effect, "| 0 |", cluster, sep = ""))
      gate_regression <- felm(gate_formula, data = data, weights = data$weight)
      return(gate_regression)
    }
  )

  # Save the coefficients in a data frame (a la broom)
  coefficients <- generate_coefficients(gate_regression, significance)
  # Sort the groups coefficients given the estimate
  sorted_groups <- coefficients %>%
    filter(str_detect(term, "S_group")) %>%
    arrange(estimate) %>%
    mutate(term = paste0("S_groupQ", row_number()))
  # Save final sorted coefficients and save corresponding confidence intervals (corrected by group estimation)
  coefficients %<>% filter(str_detect(term, "S_group", negate = TRUE)) %>%
    bind_rows(sorted_groups) %>%
    mutate(
      mod.conf.low = ifelse(str_detect(term, "^S_groupQ"), estimate - qnorm(1 - significance / (groups)) * std.error, NA),
      mod.conf.high = ifelse(str_detect(term, "^S_groupQ"), estimate + qnorm(1 - significance / (groups)) * std.error, NA)
    )

  # Perform the heterogeneity test as depicted by Chernozhukov et al. (2018)
  # higher quantile - lower quantile == 0, with a significance level of "significance"
  heterogeneity_test <- generate_hettest(gate_regression, c(paste0(group_variables %>% tail(1), "-S_groupQ1==0")), significance)

  # Calculate accuracy of GATES by the part of variation of s0(Z) explained by the grouped S(Z),
  # maximizing the R^2 in the regression (without a constant), assuming equal groups
  accuracy <- coefficients %>%
    filter(str_detect(term, "^S_groupQ")) %>%
    dplyr::select(estimate) %>%
    pull() %>%
    `^`(2) %>%
    sum() %>%
    `/`(groups)

  result <- tibble(
    gate_coefficients = coefficients %>% list(),
    gate_hettest = heterogeneity_test %>% list(),
    gate_accuracy = accuracy
  )

  if (save_model == TRUE) {
    result %<>% add_column(gate_regression = gate_regression %>% list(), .before = 1)
  }

  rm(list = ls() %>% str_subset("^result$", negate = TRUE), envir = environment()) # Save memory
  gc(reset = TRUE, verbose = FALSE) # Collect garbage

  return(result)
}

# This function perform the CLAN analysis for a group of variables
# variables := character vector of variables to analyse
# ... := options to analysis_clan()

regress_clan <- function(variables, ...) {
  variables %>%
    as_tibble_col(column_name = "clan_variable") %>%
    mutate(clan = map(clan_variable, analysis_clan, ...)) %>%
    unnest(clan)
}

# This function performs the CLAN analysis specified by Chernozhukov et al. (2018), simply calculating the means and its difference
# of the most and least affected units
# clan_variable := character of variable g(Y, Z)
# significance := significance level to calculate the confidence intervals
# data := data frame from estimate_treatment_effect()
# save_model := return regression object with the result? (FALSE)

analysis_clan <- function(clan_variable, data, significance, save_model = FALSE) {
  # Keep only most and least affected units
  data %<>% filter((most_affected == 1) | (least_affected == 1)) %>%
    dplyr::select(sym(!!clan_variable), least_affected, most_affected)

  # Transform the clan_variable in dummies (if factor)
  if (data %>% pull(!!sym(clan_variable)) %>% vctrs::vec_ptype_abbr() %>% `==`("fct")) {
    tryCatch(
      {
        data %<>% dummy_cols(
          remove_first_dummy = TRUE,
          remove_selected_columns = TRUE
        )
      },
      error = function(e) {
        warning(paste(clan_variable, " in CLAN analysis is not a codifiable factor, mean is not intrepretable \n"))
        result <- tibble(
          clan_coefficients = NULL %>% list(),
          clan_hettest = NULL %>% list(),
        )
        if (save_model == TRUE) {
          result %<>% add_column(clan_regression = NULL %>% list(), .before = 1)
        }

        rm(list = ls() %>% str_subset("^result$", negate = TRUE), envir = environment())
        gc(reset = TRUE, verbose = FALSE)

        return(result)
      }
    )
  }

  var_names <- names(data) %>% str_subset("(^most_affected$)|(^least_affected$)", negate = TRUE)

  tryCatch(
    {
      clan_formula <- reformulate("most_affected + least_affected", response = paste0("cbind(", paste(var_names, collapse = ", "), ")"), intercept = FALSE)
      clan_regression <- lm(clan_formula, data = data)
      coefficients <- generate_coefficients(clan_regression, significance)
      hettest_specification <- rep("most_affected - least_affected == 0", length(var_names)) %>% c()
      names(hettest_specification) <- var_names
      heterogeneity_test <- generate_hettest(clan_regression, hettest_specification, significance)
      result <- tibble(
        clan_coefficients = coefficients %>% list(),
        clan_hettest = heterogeneity_test %>% list(),
      )

      if (save_model == TRUE) {
        result %<>% add_column(clan_regression = clan_regression %>% list(), .before = 1)
      }

      rm(list = ls() %>% str_subset("^result$", negate = TRUE), envir = environment())
      gc(reset = TRUE, verbose = FALSE)

      return(result)
    },
    error = function(e) {
      cat(e %>% as.character())
      cat(" Error in CLAN analysis: Using means without inference instead \n", paste0("(", clan_variable, ")"), "\n")
      clan_regression <- "ERROR" %>% list()

      coefficients <- data %>%
        group_by(most_affected, least_affected) %>%
        summarise(across(everything(), mean), .groups = "drop") %>%
        pivot_longer(cols = all_of(var_names), names_to = "response", values_to = "estimate") %>%
        pivot_longer(cols = c(most_affected, least_affected), names_to = "term") %>%
        filter(value != 0) %>%
        dplyr::select(response, term, estimate, -value) %>%
        mutate(
          std.error = NA,
          statistic = NA,
          p.value = NA,
          conf.low = estimate,
          conf.high = estimate,
          modconf.low = 0.5,
          modconf.high = 0.5,
        )

      # Check if there is a significant difference in means
      heterogeneity_test <-
        tibble(
          contrast = "most_affected - least_affected",
          null.value = 0,
          estimate = coefficients$estimate[1] - coefficients$estimate[2],
          std.error = c(NA),
          statistic = c(NA),
          p.value = c(NA),
          conf.low = estimate,
          conf.high = estimate,
          modconf.low = c(0.5),
          modconf.high = c(0.5)
        )

      result <- tibble(
        clan_coefficients = coefficients %>% list(),
        clan_hettest = heterogeneity_test %>% list(),
      )

      if (save_model == TRUE) {
        result %<>% add_column(clan_regression = NULL %>% list(), .before = 1)
      }

      rm(list = ls() %>% str_subset("^result$", negate = TRUE), envir = environment())
      gc(reset = TRUE, verbose = FALSE)

      return(result)
    }
  )
}

# This function tidies the regression object coefficients (see broom) and calculates
# the corresponding the standard confidence interval
# regression := regression object
# significance := significance level to calculate the confidence intervals

generate_coefficients <- function(regression, significance) {
  tidy(regression, conf.int = TRUE, conf.level = 1 - significance) %>%
    correct_pvalue()
}

# This function calculates a test and delivers a tidied data frame (see broom), and calculates
# the corresponding the standard confidence intervals
# regression := regression object
# hypothesis := character of the linear hypothesis to test (see multcomp, glht())
# significance := significance level to calculate the confidence intervals

generate_hettest <- function(regression, hypothesis, significance) {
  if (length(hypothesis) > 1) {
    stop(" Multilevel factors comparison is not yet supported \n")
  } else {
    test <- glht(regression, linfct = hypothesis)
    test %<>%
      summary() %>%
      tidy() %>%
      full_join(confint(test, level = 1 - significance) %>% tidy(), by = c("contrast", "estimate")) %>%
      correct_pvalue()
    return(test)
  }
}

# This function corrects the p values given the Chernozhukov et al. (2018) inference
# framework, using calculate_correct_pvalue()
# tidied_regression := tidied regression object

correct_pvalue <- function(tidied_regression) {
  if ("p.value" %in% names(tidied_regression)) {
    tidied_regression %<>% mutate(
      mod.p.value.low = map2_dbl(estimate, p.value, calculate_correct_pvalue),
      mod.p.value.high = map2_dbl(estimate, p.value, calculate_correct_pvalue, low = FALSE)
    )
  } else if ("adj.p.value" %in% names(tidied_regression)) {
    tidied_regression %<>% mutate(
      mod.p.value.low = map2_dbl(estimate, adj.p.value, calculate_correct_pvalue),
      mod.p.value.high = map2_dbl(estimate, adj.p.value, calculate_correct_pvalue, low = FALSE)
    )
  }
  return(tidied_regression)
}

# This function calculates the corrected p values given the Chernozhukov et al. (2018) inference
# framework, depending on the estimate sign
# estimate := numeric value
# p.value := standard calculated p.value
# low := should the lower ci be estimated? (TRUE)

calculate_correct_pvalue <- function(estimate, p.value, low = TRUE) {
  if (is.na(estimate)) {
    NA
  }
  else if (low == TRUE) {
    if (estimate < 0) {
      (1 - p.value / 2)
    } else {
      (p.value / 2)
    }
  } else {
    if (estimate < 0) {
      (p.value / 2)
    } else {
      (1 - p.value / 2)
    }
  }
}

# This function perform multiple simulations of analise_heterogeneity() with a parallelized map (see furrr and future)
# WARNING: A packages vector is needed to load every worker
# simulations := number of desired simulations
# ... := options to analise_heterogeneity()
# move_seed = add a number to every element of the sequence of seeds (1:simulations by default)

simulate_heterogeneity_analysis <- function(simulations, ..., move_seed = 0) {
  seq(move_seed, simulations + move_seed - 1) %>%
    map_dfr(analise_heterogeneity, ...)
}

# This function calculates the (mean, length, median) of a simulation object
# from simulate_heterogeneity_analysis()
# simulations := simulation object
# check := calculate length (FALSE)
# cov_prob := calculate mean (FALSE)
# WARNING: if both parameters are set to TRUE, the mean is calculated

collect_simulations <- function(simulations, check = FALSE, cov_prob = FALSE, check_seeds = FALSE, .function) {
  if (missing(.function)) {
    if (check == TRUE) {
      .function <- "length"
    } else {
      .function <- "median"
    }
    if (cov_prob == TRUE) {
      .function <- "length"
    } else {
      .function <- "median"
    }
  }
  seeds <- simulations %>%
    pull(simulation) %>%
    sort()
  independent_seeds <- unique(seeds) %>% length()
  diff_seed <- abs(diff(seeds))
  non_sequential_seeds <- (abs(diff(seeds)) > 1)
  repeating_seeds <- (abs(diff(seeds)) < 1)
  simulations %<>%
    dplyr::select(-simulation) %>%
    group_by(method_name, response, treatment, fixed_effect, cluster) %>%
    summarize(across(matches("(observations)|(blp)|(gate)|(clan)"), ~ summarize_simulations(.x, cur_column(), .function = .function)), .groups = "keep") %>%
    arrange(response, treatment, fixed_effect, cluster)
  if (check_seeds == TRUE) {
    simulations %<>%
      mutate(
        n_seeds = all_of(independent_seeds),
        rep_seeds = all_of(seeds[repeating_seeds]) %>% list(),
        non_seq_seeds = all_of(seeds[non_sequential_seeds]) %>% `+`(all_of(diff_seed[non_sequential_seeds]) / 2) %>% list()
      )
  }
  return(simulations)
}

# Helper function to calculate the .function within a list of nested data frames or values
# data := a list of data frames or values from analise_heterogeneity()
# name := name of the analysis (blp, gates, clan)

summarize_simulations <- function(data, name, .function = "median") {
  if (str_detect(name, "clan")) {
    data %>%
      reduce(bind_rows) %>%
      group_by(clan_variable) %>%
      summarise(across(everything(), ~ calculate_summary(.x, cur_column(), .function = .function)), .groups = "keep") %>%
      list()
  } else {
    calculate_summary(data, name, .function = .function)
  }
}

# Calculate the .function of a given variable
# data := data frame resulting from an analysis function or statistical test (BLP, GATES, CLAN)
# name := name of the analysis (blp, gates, clan)
# .function := function to use
# .clean_output = eliminate variables not needed for the main analysis? (TRUE)

calculate_summary <- function(data, name, .function = "median", .clean_output = TRUE) {
  selection <- list("term", c("contrast", "null.value"))[[1 + str_detect(name, "_hettest")]] %>% unlist()
  group_selection <- list(c("conf.low", "conf.high"), c("mod.conf.low", "mod.conf.high"))[[1 + str_detect(name, "gate_coefficients")]] %>% unlist()
  if (is_list(data)) {
    result <- data %>%
      reduce(bind_rows) %>%
      group_by(!!!syms(selection)) %>%
      summarise(across(everything(), !!sym(.function)), .groups = "rowwise") %>%
      mutate(mod.p.value = min(1, 4 * min(mod.p.value.low, mod.p.value.high))) %>% # Uniform Validity of Variational P-Value
      ungroup() %>%
      dplyr::select(-c(mod.p.value.low, mod.p.value.high))
    if (.clean_output == TRUE) {
      if (str_detect(name, "(blp_coefficients)|(gate_coefficients)")) {
        tryCatch(
          {
            result %<>% filter(str_detect(term, "(^B$)|(^S$)", negate = TRUE))
          },
          error = function(e) {
            warning("B or S variables were not found")
            result
          }
        )
      }
      result %<>% dplyr::select(!!!syms(selection), estimate, !!!syms(group_selection), mod.p.value)
    }
    result %<>% list()
  } else {
    rlang::exec(.function, data)
  }
}

# Print a latex table of the results from a collect_simulation() object
# results := a collect_simulation() object
# table := name of the analysis (blp, gates, clan)
# .filter := delete other groups information

print_results <- function(results, table, .filter = TRUE) {
  if (table == "best") {
    results %<>%
      dplyr::select(`method_name`, `response`, `treatment`, `fixed_effect`, `cluster`, blp_accuracy, gate_accuracy) %>%
      mutate(
        across(where(is.character), ~ map(.x, change_name)),
        across(where(is.numeric), round, digits = 4)
      ) %>%
      rename_with(~ .x %>% map_chr(change_name), dplyr::everything()) %>%
      group_by(response, treatment, fixed_effect, cluster)
    return(results)
  }

  tbl_selector <- list("blp_coefficients", "gate_coefficients", "clan")[[1 + str_detect(table, "gate") + 2 * str_detect(table, "clan")]] %>% unlist()
  group_selector <- list(c("response", "treatment", "fixed_effect", "cluster"), c("response", "treatment", "fixed_effect", "cluster", "clan_variable"))[[1 + str_detect(table, "clan")]] %>% unlist()
  ci <- c("conf.low", "conf.high")

  results_table <- results %>%
    dplyr::select(`method_name`, `response`, `treatment`, `fixed_effect`, `cluster`, !!sym(tbl_selector)) %>%
    unnest(!!sym(tbl_selector))

  if (table == "clan") {
    results %<>% unnest(clan)
    results_table %<>%
      dplyr::select(method_name, !!!syms(group_selector), clan_coefficients) %>%
      unnest(clan_coefficients)
  }

  results_table %<>%
    mutate(term = map_chr(term, change_name)) %>%
    group_by(!!!syms(group_selector)) %>%
    rename_with(~ .x %>% str_remove_all("mod."), matches("mod.conf"))

  if (table == "gate") {
    groups <- results_table$term %>%
      str_remove_all("Group ") %>%
      as.numeric() %>%
      max(na.rm = TRUE)
    if (.filter == TRUE) {
      results_table %<>% filter(term == "Group 1" | term == paste("Group", groups))
    }
  }

  if (str_detect(table, "(gate)|(clan)")) results_table %<>% bind_rows(print_test(results, table))

  results_table %<>%
    mutate(
      estimate = pmap_chr(
        list(estimate, !!!syms(ci), mod.p.value),
        function(estimate, conf.low, conf.high, mod.p.value) {
          paste(round(estimate, digits = 4),
            paste0("(", round(conf.low, digits = 4), ", ", round(conf.high, digits = 4), ")"),
            paste0("[", print_pvalue(mod.p.value), "]"),
            sep = "<br>"
          )
        }
      ),
      method_name = map_chr(method_name, change_name)
    ) %>%
    dplyr::select(method_name, !!!syms(group_selector), term, estimate, -c(!!!syms(ci), mod.p.value)) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    group_by(response, treatment, fixed_effect, cluster)
  names(results_table) %<>% map_chr(change_name)
  return(results_table)
}

# Helper function to print the results of a statistical test from a collect_simulation() object (BLP, CLAN)
# results := a collect_simulation() object
# table := name of the analysis (blp, gates)

print_test <- function(results, table) {
  tbl_selector <- list("gate_hettest", c("clan_hettest", "clan_variable"))[[1 + str_detect(table, "clan")]] %>% unlist()
  var_selector <- list("estimate", c("clan_variable", "estimate"))[[1 + str_detect(table, "clan")]] %>% unlist()
  results_table <- results %>%
    dplyr::select(`method_name`, `response`, `treatment`, `fixed_effect`, `cluster`, !!!syms(tbl_selector)) %>%
    unnest(!!sym(tbl_selector[1])) %>%
    dplyr::select(`method_name`, `response`, `treatment`, `fixed_effect`, `cluster`, contrast, !!!syms(var_selector), conf.low, conf.high, mod.p.value) %>%
    rename(term = contrast) %>%
    mutate(term = "Difference") %>%
    group_by(response, treatment, fixed_effect, cluster, !!sym(var_selector[1]))
  return(results_table)
}

# Helper function to print a nice p.value if it is too small
# p.value := a p.value (numeric)

print_pvalue <- function(p.value) {
  if (p.value < 0.001) {
    p.value <- "< 0.001"
  } else {
    p.value %<>% round(digits = 3)
  }
  return(p.value)
}

# Dictionary to print translated results for presentation
# name := a name

change_name <- function(name) {
  if (name == "D_ort") name %<>% str_replace("D_ort", "ATE")
  if (name == "S_ort") name %<>% str_replace("S_ort", "HET")
  if (name == "method_name") name %<>% str_replace("method_name", "ML Method")
  if (name == "estimate") name %<>% str_replace("estimate", "Coefficient")
  if (name == "conf.low") name %<>% str_replace("conf.low", "Lower Band (CI)")
  if (name == "conf.high") name %<>% str_replace("conf.high", "Upper Band (CI)")
  if (name == "mod.conf.low") name %<>% str_replace("mod.conf.low", "Lower Band (CI)")
  if (name == "mod.conf.high") name %<>% str_replace("mod.conf.high", "Upper Band (CI)")
  if (name == "mod.p.value") name %<>% str_replace("mod.p.value", "P value")
  if (str_detect(name, "S_groupQ.")) name %<>% str_replace("S_groupQ", "Group ")
  if (str_detect(name, "random_forest")) name %<>% str_replace("random_forest", "Random Forest ") %>% str_remove("_")
  if (str_detect(name, "elastic_net")) name %<>% str_replace("elastic_net", "Elastic Net ") %>% str_remove("_")
  if (name == "blp_accuracy") name %<>% str_replace("blp_accuracy", "Best BLP")
  if (name == "gate_accuracy") name %<>% str_replace("gate_accuracy", "Best GATES")
  if (name == "most_affected") name %<>% str_replace("most_affected", "Most Affected")
  if (name == "least_affected") name %<>% str_replace("least_affected", "Least Affected")
  if (name == "clan_variable") name %<>% str_replace("clan_variable", "CLAN Variable")
  if (name == "n_leader_base") name %<>% str_replace("n_leader_base", "Nominations Leader")
  if (name == "n_popular_base") name %<>% str_replace("n_popular_base", "Nominations Popular")
  if (name == "n_friendly_base") name %<>% str_replace("n_friendly_base", "Nominations Friendly")
  if (name == "n_shy_base") name %<>% str_replace("n_shy_base", "Nominations Shy")
  if (name == "n_academic_base") name %<>% str_replace("n_academic_base", "Nominations Academic")
  return(name)
}

# Plot the results of the heteorgeneity analysis a la Chernozhukov et al. (2018) [Figures 4-6]
# results := a collect_simulations() object, using median as a summarizing function
# WARNING: No checks yet available, do not use it with a collect_simulations() object
# with other summarizing functions different than "median"

plot_results <- function(results) {
  data_treatment <- results %>%
    dplyr::select(`method_name`, `treatment`, `response`, `fixed_effect`, `cluster`, blp_coefficients) %>%
    unnest(blp_coefficients) %>%
    filter(term == "D_ort") %>%
    mutate(
      method_name = map_chr(method_name, change_name),
      treatment = str_replace_all(treatment, "ingroup$", "Distance (10)"),
      treatment = str_replace_all(treatment, "ingroup_3$", "Distance (3)"),
      treatment = str_replace_all(treatment, "ingroup_5$", "Distance (5)"),
      treatment = str_replace_all(treatment, "(.*)_base_(.*)_", "Interaction with "),
      treatment = str_to_sentence(treatment)
    )
  data_groups <- results %>%
    dplyr::select(`method_name`, `treatment`, `response`, `fixed_effect`, `cluster`, gate_coefficients) %>%
    unnest(gate_coefficients) %>%
    filter(str_detect(term, "S_group")) %>%
    mutate(
      term = map_chr(term, str_replace, "S_groupQ", "Group "),
      group = map_dbl(term, ~ str_remove_all(.x, "Group ") %>% as.numeric()),
      method_name = map_chr(method_name, change_name),
      treatment = str_replace_all(treatment, "ingroup$", "Distance (10)"),
      treatment = str_replace_all(treatment, "ingroup_3$", "Distance (3)"),
      treatment = str_replace_all(treatment, "ingroup_5$", "Distance (5)"),
      treatment = str_replace_all(treatment, "(.*)_base_(.*)_", "Interaction with "),
      treatment = str_to_sentence(treatment)
    )

  ggplot(data_treatment) +
    theme_gray(base_size = 14) +
    facet_grid(cols = vars(response), rows = vars(treatment, method_name), scales = "free_y") +
    geom_point(data = data_groups, aes(y = estimate, x = group), size = 3) +
    geom_errorbar(data = data_groups, aes(ymax = mod.conf.high, ymin = mod.conf.low, x = group, y = estimate, width = 0.7), show.legend = TRUE) +
    geom_hline(data = data_treatment, aes(yintercept = estimate), linetype = 2) +
    geom_hline(data = data_treatment, aes(yintercept = conf.low), linetype = 2, color = "blue") +
    geom_hline(data = data_treatment, aes(yintercept = conf.high), linetype = 2, color = "blue") +
    labs(title = "GATES", y = "Treatment Effect", x = "Group by Het Score")
}

plot_clan <- function(results, filter = ".", modify_vars = FALSE, use_grouping = TRUE) {
  data_clan <- results %>%
    dplyr::select(`method_name`, `treatment`, `response`, `fixed_effect`, `cluster`, clan) %>%
    unnest(clan) %>%
    filter(str_detect(clan_variable, filter)) %>%
    mutate(
      clan_hettest = map(clan_hettest, ~ dplyr::select(.x, -null.value) %>%
        mutate(contrast = "difference") %>%
        rename(term = contrast)),
      clan_data = map2(clan_coefficients, clan_hettest, bind_rows)
    ) %>%
    dplyr::select(-c(clan_coefficients, clan_hettest)) %>%
    unnest(clan_data)

  if (modify_vars == TRUE) {
    data_clan %<>%
      mutate(
        clan_variable = if_else(str_detect(clan_variable, "_base_"), "Baseline", clan_variable),
        treatment = str_replace_all(treatment, "ingroup$", "Distance (10)"),
        treatment = str_replace_all(treatment, "ingroup_3$", "Distance (3)"),
        treatment = str_replace_all(treatment, "ingroup_5$", "Distance (5)"),
        treatment = str_replace_all(treatment, "(.*)_base_(.*)_", "Interaction with "),
        treatment = str_to_sentence(treatment),
        grouping = case_when(
          str_detect(clan_variable, "idego") ~ "Ego",
          str_detect(clan_variable, "idalter") ~ "Alter",
          TRUE ~ "Homophily"
        ),
        clan_variable = str_remove_all(clan_variable, "(s_)|(.idego)|(.idalter)"),
        across(where(is.character), ~ .x %>% map_chr(change_name)),
        mod.p.value = if_else(mod.p.value < 0.05, TRUE, FALSE)
      ) %>%
      pivot_wider(values_from = c(estimate, conf.low, conf.high, mod.p.value), names_from = term)
    names(data_clan) %<>% str_to_lower() %>% str_replace_all(" ", "_")
  }

  if (use_grouping == TRUE) {
    plot <- ggplot(data_clan) +
      geom_linerange(aes(
        y = grouping,
        xmin = estimate_least_affected,
        xmax = estimate_most_affected,
        group = treatment,
        color = treatment,
        linetype = mod.p.value_difference
      ),
      position = position_dodgev(0.75), size = 1
      ) +
      geom_point(aes(
        y = grouping,
        x = estimate_least_affected,
        group = treatment,
        color = treatment,
        shape = mod.p.value_least_affected
      ),
      position = position_dodgev(0.75), size = 2
      ) +
      geom_point(aes(
        y = grouping,
        x = estimate_most_affected,
        group = treatment,
        color = treatment,
        shape = mod.p.value_most_affected
      ),
      position = position_dodgev(0.75), size = 2
      ) +
      geom_errorbarh(aes(
        y = grouping,
        xmin = conf.low_least_affected,
        xmax = conf.high_least_affected,
        group = treatment,
        color = treatment
      ),
      height = 0.2,
      position = position_dodgev(0.75)
      ) +
      geom_errorbarh(aes(
        y = grouping,
        xmin = conf.low_most_affected,
        xmax = conf.high_most_affected,
        group = treatment,
        color = treatment
      ),
      height = 0.2,
      position = position_dodgev(0.75)
      ) +
      facet_grid(clan_variable + method_name ~ response, scales = "free_x") +
      scale_linetype_manual(values = c("dotted", "solid")) +
      scale_shape_manual(values = c(17, 19)) +
      labs(
        shape = "Is mean significant at 95%?",
        linetype = "Is difference significant at 95%?",
        color = "Treatment",
        x = "CLAN Means",
        y = ""
      )
  } else {
    plot <- ggplot(data_clan) +
      geom_linerange(aes(
        y = clan_variable,
        xmin = estimate_least_affected,
        xmax = estimate_most_affected,
        group = treatment,
        color = treatment,
        linetype = mod.p.value_difference
      ),
      position = position_dodgev(0.75), size = 1
      ) +
      geom_point(aes(
        y = clan_variable,
        x = estimate_least_affected,
        group = treatment,
        color = treatment,
        shape = mod.p.value_least_affected
      ),
      position = position_dodgev(0.75), size = 2
      ) +
      geom_point(aes(
        y = clan_variable,
        x = estimate_most_affected,
        group = treatment,
        color = treatment,
        shape = mod.p.value_most_affected
      ),
      position = position_dodgev(0.75), size = 2
      ) +
      geom_errorbarh(aes(
        y = clan_variable,
        xmin = conf.low_least_affected,
        xmax = conf.high_least_affected,
        group = treatment,
        color = treatment
      ),
      height = 0.2,
      position = position_dodgev(0.75)
      ) +
      geom_errorbarh(aes(
        y = clan_variable,
        xmin = conf.low_most_affected,
        xmax = conf.high_most_affected,
        group = treatment,
        color = treatment
      ),
      height = 0.2,
      position = position_dodgev(0.75)
      ) +
      facet_grid(method_name ~ response, scales = "free_x") +
      scale_linetype_manual(values = c("dotted", "solid")) +
      scale_shape_manual(values = c(17, 19)) +
      labs(
        shape = "Is mean significant at 95%?",
        linetype = "Is difference significant at 95%?",
        color = "Treatment",
        x = "CLAN Means",
        y = ""
      )
  }
  return(plot)
}

# Helper function to add baseline variables depending on the selected response variable and the controls vector
add_baseline2parameters <- function(response, vector) {
  baseline <- str_replace(response, "_", "_base_")
  vector %<>% c(baseline)
}

# Helper function to add baseline variables depending on the selected response variable and the controls vector
add_interaction2treatment <- function(response, treatment) {
  baseline <- str_replace(response, "_", "_base_")
  treatment %<>% str_replace_all("^x", baseline)
  return(treatment)
}

# Helper function to reduce memory usage
select_data4anaylsis <- function(data, parameters) {
  keep_variables <- parameters %>%
    mutate(use_vars = pmap(.l = list(response, controls, treatment, cluster, fixed_effect, affected), c)) %>%
    pull(use_vars) %>%
    reduce(c) %>%
    unique()
  data %<>% dplyr::select(!!!syms(keep_variables))
  return(data)
}

read_n_merge_simulations <- function(dir = "heteffects_simulations", ...) {
  files <- list.files(dir) %>%
    str_remove_all(".RData") %>%
    str_split_fixed("-", 5)
  colnames(files) <- c("response", "treatment", "method_name", "hash", "time")
  files %<>% as_tibble() %>%
    group_by(response, treatment, method_name) %>%
    group_nest(.key = "simulations") %>%
    mutate(
      simulations = map_int(simulations, nrow),
      key = paste(response, treatment, method_name, sep = "-"),
      key = map(key, ~ list.files(path = dir, pattern = .x, full.names = TRUE)),
      results = map(key, read_by_name, ...),
      check = map_int(results, nrow) == 1
    )
  result <- files %>% filter(check == TRUE)
  dropped <- files %>% filter(check == FALSE)
  if (nrow(dropped) != 0) {
    warning(paste("Inconsistent simulations were found at",
                  dropped$key,
                  "Those results are automatically dropped.",
                  "Please check the files or re-run the simulations!",
                  sep = "\n"
    ))
  }
  result %<>% dplyr::select(simulations, results) %>%
    unnest(results)
}

read_by_name <- function(names, ...) {
  names %>%
    map(function(file) {
      load(file)
      return(result)
    }) %>%
    reduce(bind_rows) %>%
    group_by(affected_quantile, folds, groups, repetitions, significance, stratum) %>%
    group_split() %>%
    map(collect_simulations, ...) %>%
    reduce(bind_rows)
}

correct_simulation_tests <- function(results) {
  results %>%
    mutate(gate_hettest = map(gate_hettest, merge_hettests))
}

merge_hettests <- function(gate_hettest) {
  gate_hettest %>%
    summarise(
      across(where(is.character), tail, n = 1),
      across(where(is.numeric), median)
    )
}

############## SIDE MODULES

###### MODULE: Generate perceptions or coincidences data from the original network file (perceptions)

### Count number of nominations for each id by groups
## Arguments:
# data := raw database
# characteristic := chracteristic to count (friendly, shy, leader, popular)
# time := survey data to use (suffixes as base or end)
# filtering := filtering variable at database (column name)
# groups := grouping options for counting

count_nominations <- function(data, characteristic, time, groups) {
  name <- paste(characteristic, time, sep = "_")
  grouped_data <- data %>%
    group_by_at(groups) %>%
    dplyr::count(wt = get(name), name = characteristic) %>%
    mutate(characteristic_t = time)
}

### Count the number of diferent nominated ids above certain threshold by groups and characteristics
## Arguments:
# data := Data of the number of nominations for each id by groups and characteristics
# characteristic := chracteristic to count (friendly, shy, leader, popular)
# groups := grouping options for counting
# threshold := Minimum Threshold to count the data

count_nominations_by_thld <- function(data, characteristic, groups, threshold) {
  filtered_data <- data %>%
    filter(!!sym(characteristic) >= threshold) %>%
    group_by_at(groups) %>%
    dplyr::count(name = paste0(characteristic, "_n"))
}

### Count number of nominations for each id by groups
## Arguments:
# data := raw database
# characteristic := chracteristic to count (friendly, shy, leader, popular)
# time := survey data to use (suffixes as base or end)
# filtering := filtering variable at database (column name)
# groups := grouping options for counting

count_coincidences <- function(variable, data) {
  data %>%
    dplyr::select(coar, cohort, idego, idalter, !!sym(variable)) %>%
    filter(!!sym(variable) > 0) %>%
    dplyr::select(-!!sym(variable)) %>%
    nest(data = idego) %>%
    mutate(data = map(data, cross_itself)) %>%
    unnest(cols = data) %>%
    group_by(coar, cohort, student_a, student_b) %>%
    dplyr::count(name = paste0("coincidences_", variable))
}

###### MODULE: Generate plots and histrograms from communities data (plots)

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
    facet_grid(rows = vars(characteristic), cols = vars(characteristic_t), scales = "free") +
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

############## BASE FUNCTIONS (They are used inside other functions)

cross_itself <- function(x) {
  x %<>% pull(idego)
  if (is.factor(x)) {
    x %<>% as.character()
  }
  tidyr::crossing(student_a = x, student_b = x) %>% filter(student_a != student_b)
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

pm2t_join <- function(tibble, matrix, name, groups = c("idego", "idalter")) {
  pipe_matrix2tibble(matrix, name = name) %>% right_join(tibble, by = groups)
}

############## LOADING FUNCTIONS (They are used in all modules)

load_lastdataset <- function(name, module) {
  setwd(paste0(module, "/"))
  position <- file.info(list.files() %>% str_subset(name)) %>%
    `[[`("mtime") %>%
    which.max()
  if (length(position) != 0) {
    full_name <- list.files() %>%
      str_subset(name) %>%
      `[[`(position)
    extension <- full_name %>% str_extract("\\..*")
    if (extension == ".csv") {
      assign(name, read_csv(full_name), envir = globalenv())
    } else if (extension == ".RData") {
      load(full_name, envir = globalenv())
    } else {
      setwd("..")
      stop(paste("Non suported extension:", full_name))
    }
    setwd("..")
  } else {
    setwd("..")
    stop(paste("File does not exist:", name, "(Run appropiate module)"))
  }
}

load_datasrc <- function(src_name) {
  load_parameters <- c(
    "loaded_data",
    "loaded_perceptions_data",
    "loaded_covariates_data",
    "loaded_coincidences_data"
  )
  if (!exists(load_parameters, envir = globalenv())) {
    source(src_name)
  } else if (any(map_lgl(load_parameters, ~ get(.x, envir = globalenv()) == FALSE))) {
    source(src_name)
  }
}
