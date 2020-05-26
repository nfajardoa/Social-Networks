########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Join Community Membership with Perception Data
########################### Written by: Nicolas Fajardo
########################### Last modified: 23/05/2020
########################### Comments:
########################### 

## Check if data is loaded, if not, load it.
load_datasrc("socialnetworks_perceptions.R")
load_lastdataset("communities_membership")

## 
jointdatabase <- left_join(communities_membership, perceptions_idalter %>% rename(idego = idalter)) %>%
  mutate(idego = idego %>% as_factor()) %>%
  complete(charc_time, nesting(variable, selector, school, cohort, time, idego, !!!rlang::syms(community_algorithms)), fill = setNames(rep(times = length(characteristics), x = 0L), characteristics) %>% as.list()) %>%
  filter(!is.na(charc_time)) %>%
  dplyr::select(variable, selector, school, cohort, time, idego, charc_time, !!!rlang::syms(community_algorithms), !!!rlang::syms(characteristics)) # idalter was renamed as idego (for merging purposes)

communities_membership_wide <- full_join(communities_membership, communities_membership %>% rename(idalter = idego), by = c("variable", "selector", "school", "cohort", "time"), suffix = c(".idego", ".idalter"))
communities_membership_wide <- communities_membership %>% 
  pivot_wider(id_cols = c(school, cohort, idego), names_from = c(variable, time, selector), values_from = contains(community_algorithms)) %>%
  rename(coar = school)
communities_membership_wide %<>% full_join(communities_membership_wide %>% rename(idalter = idego), by = c("coar", "cohort"), suffix = c(".idego", ".idalter"))
arguments <- expand_grid(algorithm = community_algorithms, variable = variables, time = times, selector = selectors) %>% 
  mutate(name = paste(algorithm, variable, time, selector, sep="_")) %>% dplyr::select(name)
communities_membership_same <- pmap(arguments %>% as.list, mutate_by_algorithm, data = communities_membership_wide) %>%
  reduce(inner_join, by = c("coar", "cohort", "idalter", "idego"))
communities_wide <- inner_join(communities_membership_wide, communities_membership_same, by = c("coar", "cohort", "idego", "idalter")) %>%
  right_join(data, by = c("coar", "cohort", "idego", "idalter"))

write_csv(communities_wide, path = paste0("network_data_regressions_", Sys.Date(), ".csv"))
#rm(data, arguments, communities, communities_membership, communities_membership_same, communities_membership_wide, summary_id)
save.image(file = paste0("network_data_regressions_", Sys.Date(), ".RData"))