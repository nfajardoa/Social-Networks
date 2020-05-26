########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Summarize Perception Data
########################### Written by: Nicolas Fajardo
########################### Last modified: 23/05/2020
########################### Comments:
########################### 

## Check if data is loaded, if not, load it.
load_datasrc("socialnetworks_summstats.R", main = TRUE)

## Set parameters
characteristics <- c("friendly", "leader", "popular", "shy")
characteristics_t <- c("base", "end")

## Summarize information of perceptions at the individual level
perceptions_idalter <- future_pmap_dfr(expand_grid(characteristic = characteristics, time = characteristics_t), 
                                    count_nominations, data = data, groups = c("coar", "cohort", "idalter")) %>%
  mutate(charc_time = charc_time %>% as.factor()) %>% group_by(coar, cohort, charc_time, idalter) %>%
  summarise_all(sum, na.rm = TRUE) %>% rename(school = coar)

## Summarize information of identified individuals according to the threshold (friendly, shy, popular, leader)
## School x Cohort | WARNING: Threshold is set at 1 by default
perceptions_school <- future_pmap_dfr(expand_grid(characteristic = characteristics), 
                                      count_nominations_by_thld, data = perceptions_idalter, groups = c("school", "cohort", "charc_time"), threshold = 1) %>%
  group_by(school, cohort, charc_time) %>%
  summarise_all(sum, na.rm = TRUE)