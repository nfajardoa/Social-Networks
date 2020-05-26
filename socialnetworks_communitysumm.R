# arguments <- expand_grid(characteristic = characteristics, algorithm = community_algorithms)
# network_summary <- future_pmap_dfr(arguments, summarize_by_algorithm, data = network_data_merged, threshold = 1, variable = "idego", group = c("variable", "selector", "school", "cohort", "time"))
# network_summary2 <- future_pmap_dfr(arguments, summarize_by_algorithm, data = network_data_merged, threshold = 1, variable = "idego", group = c("variable", "selector", "school", "cohort", "time", "charc_time"))

## STORING DATA

closeAllConnections()
# write_csv(summary_id %>% select_if(negate(is.list)), path = "summary_id.csv")
# write_csv(network_summary %>% select_if(negate(is.list)), path = "network_summary.csv")
# write_csv(network_summary2 %>% select_if(negate(is.list)), path = "network_summary2.csv")
# write_csv(characteristics_school, path = "characteristics_school.csv")
# write_csv(characteristics_idalter, path = "characteristics_idalter.csv")
# write_csv(network_data_merged, path = "network_data_merged.csv")
write_csv(communities_membership, path = paste0("network_data_by_community_", Sys.Date(), ".csv"))

load("communities_membership_2020-04-02.RData")
source("network_functions.R")

communities %<>% select(algorithm, variable, time, selector, school, cohort, communities_n) %>%
  group_by(algorithm, variable, time, selector) %>% summarize(mean = mean(communities_n),
                                                              min = min(communities_n),
                                                              q1 = quantile(communities_n, probs = 0.25),
                                                              q3 = quantile(communities_n, probs = 0.75),
                                                              max = max(communities_n)) %>%
  mutate(name = paste(algorithm, variable, time, selector, sep="_")) %>% ungroup() %>%
  select(name, mean, min, q1, q3, max) %>% mutate(name = str_replace_all(name, "none_", ""))

write_csv(communities, path = paste0("communities_summary_", Sys.Date(), ".csv"))

communities %<>% filter_at(vars(starts_with("name")), any_vars(!(str_detect(.,"mid") | str_detect(.,"end")))) 

write_csv(communities, path = paste0("communities_short_summary_", Sys.Date(), ".csv"))