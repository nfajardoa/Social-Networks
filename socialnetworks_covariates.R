covariates_data <- read_csv("covariates_data.csv") %>% mutate(network = paste(coar, cohort, sep = "_")) %>%
  select(-c(coar, cohort))
regressions_data <- read_csv("regressions_data_2020-04-18.csv") 
coincidences_data <- read_csv("coincidences_2020-04-06.csv") %>% rename(idego = student_a, idalter = student_b) %>%
  mutate(network = paste(coar, cohort, sep = "_")) %>% select(-c(coar, cohort)) %>%
  rename_all(~ str_replace_all(., pattern = "friendly", "frdly")) %>%
  rename_all(~ str_replace_all(., pattern = "leader", "ldr")) %>%
  rename_all(~ str_replace_all(., pattern = "popular", "pop")) %>%
  rename_all(~ str_replace_all(., pattern = "coincidences", "coincs"))

merged_data <- regressions_data %>% 
  left_join(covariates_data %>% rename(idego = id), by = c("idego", "network")) %>%
  left_join(covariates_data %>% rename(idalter = id), by = c("idalter", "network"), suffix = c("_idego", "_idalter")) %>%
  left_join(coincidences_data, by = c("idalter", "idego", "network"))

merged_data$idalter %<>% as.character
merged_data$idego %<>% as.character

data <- read_csv("network_data.csv")
data$idalter %<>% as.character
data$idego %<>% as.character
names(data) %<>% str_replace_all("lider", "leader")
data %<>% select(idalter, idego, ingroup)

merged_data %<>% left_join(data, by = c("idego", "idalter")) %>%
  select_if(!(stringr::str_detect(names(.), "end") & stringr::str_detect(names(.), "coincs", negate = TRUE)) | stringr::str_detect(names(.), "mid")) %>%
  mutate(s_sex = if_else(male_idego == male_idalter, "1", "0", ".")) %>%
  mutate_all(replace_na, replace = ".")

write_csv(merged_data, path = paste0("regressions_data_", Sys.Date(), ".csv"))
