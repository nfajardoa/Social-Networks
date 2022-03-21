########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Summarize Perception Data
########################### Written by: Nicolas Fajardo
########################### Last modified: 24/06/2020
########################### Comments:
########################### 

## Load network data
data <- read_csv("network_data.csv", #Read .csv on directory
                 col_types = as.list(c(rep("f", 4), rep("i", 42), rep("f", 4), rep("i", 60)))) %>% #Set column types
  rename_all(~str_replace_all(.x, "lider", "leader")) %>% # Replace some typos
  rename_all(~str_replace_all(.x, "skill", "academic"))
loaded_data <- TRUE

## Arguments
characteristics <- c("friendly", "leader", "popular", "shy", "academic")
characteristics_t <- c("base", "end")

## Conditionals for running
perception_conditionals <- c("calculate_perceptions", "calculate_perceptions_by_school", "calculate_coincidences") #Conditionals to calculate
if(all(!map_lgl(perception_conditionals, exists)) == TRUE) {map_lgl(perception_conditionals, function(cond) assign(cond, TRUE, envir = globalenv())) #If no conditional exists, all objects are calculated
}else {map_lgl(perception_conditionals[!map_lgl(perception_conditionals, exists)], function(cond) assign(cond, FALSE, envir = globalenv()))} #If at least one conditional exists, all the non-existent are NOT calculated

## Summarize information of perceptions at the individual level
if(calculate_perceptions == TRUE) {
  
  perceptions_data <- pmap_dfr(expand_grid(characteristic = characteristics, time = characteristics_t), #Generate all combinations of characteristics and characteristics_t
                                      count_nominations, data = data, groups = c("coar", "cohort", "idalter")) %>% #Count nominations by groups (Custom function)
    mutate(characteristic_t = characteristic_t %>% as.factor()) %>% #Convert "characteristic_t" to factor
    group_by(coar, cohort, characteristic_t, idalter) %>% #Group tibble
    summarise_all(sum, na.rm = TRUE) 
  
  # Write .csv file
  write_csv(perceptions_data, path = "perceptions_idalter.csv")
  
  ## Summarize information of identified individuals according to the threshold (friendly, shy, popular, leader)
  ## School x Cohort | WARNING: Threshold is set at 1 by default
  if(calculate_perceptions_by_school == TRUE) {
    perceptions_school <- pmap_dfr(tibble(characteristic = characteristics), #Create tibble with "characteristics" name
                                          count_nominations_by_thld, #Count nominations by threshold (Custom function)
                                          data = perceptions_data, #Data from perceptions_data
                                          groups = c("coar", "cohort", "characteristic_t"), 
                                          threshold = 1) %>%
      group_by(coar, cohort, characteristic_t) %>% #Group tibble
      summarise_all(sum, na.rm = TRUE) #Collapse values
    
    # Write .csv file
    write_csv(perceptions_school, path = "perceptions_school.csv")
  }
}

## Calculate coincidences on perpeceptions
if(calculate_coincidences == TRUE) {
  coincidences_data <- expand_grid(characteristics, characteristics_t) %>% #Generate all combinations of characteristics and characteristics_t
    mutate(variable = paste(characteristics, characteristics_t, sep = "_")) %>% #Generate a variable pasting both columns
    dplyr::select(variable) %>% #Select only "variable" column
    mutate(coincidences_by_var = map(variable, count_coincidences, data = data)) %>% #Calculate coincidences (Custom function)
    summarize(data = reduce(coincidences_by_var, full_join, by = c("coar", "cohort", "student_a", "student_b")) %>% list()) %>% #Merge observations by idego and idalter 
    unnest(cols = data) %>% #Unnest calculations
    replace(is.na(.), 0) %>% #Replace missing values
    rename(idego = student_a, idalter = student_b) #Rename variables
  
  # Write .csv file
  write_csv(coincidences_data, path = "coincidences_nominations.csv")
}