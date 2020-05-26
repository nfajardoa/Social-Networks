coincidences <- expand_grid(characteristics, characteristics_t) %>% 
  mutate(variable = paste(characteristics, characteristics_t, sep = "_")) %>% 
  select(variable) %>% 
  mutate(coincidences_by_var = map(variable, count_coincidences, data = data)) %>% 
  summarize(data = reduce(coincidences_by_var, full_join, by = c("coar", "cohort", "student_a", "student_b")) %>% list()) %>% 
  unnest(cols = data) %>%
  replace(is.na(.), 0)

write_csv(coincidences, path = paste0("coincidences_", Sys.Date(), ".csv"))
