########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Join Community Membership with Perception Data
########################### Written by: Nicolas Fajardo
########################### Last modified: 23/05/2020
########################### Comments:
########################### 

## Load Database
database <- src_postgres(dbname = "socialnetworks", user = "test", password = "perras", host = "localhost", port = "5432")
database_connection <- dbConnect(RPostgres::Postgres(), dbname = 'socialnetworks', host = 'localhost', port = 5432, user = 'test', password = 'perras')

## Check if data is loaded, if not, load it.
load_datasrc("socialnetworks_summstats.R", main = TRUE)
#load_lastdataset("weighted_distance")
load("weighted_distance_2020-05-27.RData")
wdistance_names <- names(data) %>% str_subset("contiguo") %>% paste0("w", .) %>% paste(collapse = " + ")
distance_names <- names(data) %>% str_subset("contiguo") %>% paste(collapse = " + ")
regression_formulas <- list(paste0(wdistance_names, " | idego + idalter +  | 0 | network_id"),
                 paste0(wdistance_names, " + ", distance_names, " | idego + idalter +  | 0 | network_id"))

## Renew Connections & Run Models
weighted_distance_models <- weighted_distance %>% mutate(data = map(name, tbl, src = database)) %>% #Renew Connections
  mutate(models = map(name, ~ multiple_reformulate(regression_formulas, response = .x))) %>%
  unnest(models) %>% mutate(model = map2(model_formula, data, ~ felm(formula = .x, data = collect(.y), cmethod = "reghdfe")))
save(weighted_distance_models, file = paste0("weighted_distance_models_", Sys.Date(), ".RData"))

## Print models
weighted_distance_models %>% group_by(variable, selector) %>% 
  group_split() %>% walk(sink_latex_output, filename = "weighted_distance_regressions_output.txt")

##### Normal distance

test_data <- read_csv("regressions_data_2020-05-04.csv", na = ".") 
test_data$idalter %<>% as.character
test_data$idego %<>% as.character

est <- felm(coincs_frdly_base ~ 1 + ingroup + factor(s_sex) | idego + idalter | 0 | network_id, data=test_data, cmethod ='reghdfe')
est_2 <- felm(coincs_frdly_end ~ 1 + ingroup + factor(s_sex) | idego + idalter | 0 | network_id, data=data, cmethod ='reghdfe')

est_a <- felm(coincs_frdly_end ~ ingroup + factor(s_sex) + s_sex*male_idego | idego + idalter | 0 | network_id, data=data, cmethod ='reghdfe')

coincidences_names <- names(data) %>% str_subset("coincs")
contiguo_names <- names(data) %>% str_subset("contiguo") %>% paste(collapse = " + ")
variables <- c(contiguo_names, "ingroup")

models <- expand_grid(response = coincidences_names, variables) %>% 
  mutate(variables = map(variables, paste0, values = " + factor(s_sex) | idego + idalter | 0 | network_id"),
         formula = map2(variables, response, reformulate),
         model = map(formula, felm, data = data, cmethod ='reghdfe'),
         coefficients = map(model, tidy),
         summary = map(model, glance)) 

models %>% unnest(coefficients) %>% dplyr::select(response, formula, term, estimate, std.error, statistic, p.value) %>%
  pivot_longer(cols = c(estimate, std.error, statistic, p.value), names_to = "type") 

%>% 
  pivot_wider(id_cols = c(formula, term, type), values_from = value)

pivot_wider(id_cols = c("term"))

communities_wide %<>% rename_all(~ str_replace_all(., pattern = "none_", "")) %>% 
  rename_all(~ str_replace_all(., pattern = "betweenness", "btwn")) %>%
  rename_all(~ str_replace_all(., pattern = "greedy", "grdy")) %>%
  rename_all(~ str_replace_all(., pattern = "leading", "lding")) %>%
  rename_all(~ str_replace_all(., pattern = "eigen", "eign")) %>%
  rename_all(~ str_replace_all(., pattern = "label", "lbl")) %>%
  rename_all(~ str_replace_all(., pattern = "infomap", "infmp")) %>%
  rename_all(~ str_replace_all(., pattern = "louvain", "lvain")) %>%
  rename_all(~ str_replace_all(., pattern = "walktrap", "wktrp")) %>%
  rename_all(~ str_replace_all(., pattern = "same_", "s_")) %>%
  rename_all(~ str_replace_all(., pattern = "friend", "frd")) %>%
  rename_all(~ str_replace_all(., pattern = "study", "stdy")) %>%
  rename_all(~ str_replace_all(., pattern = "social", "soc")) %>%
  dplyr::select(coar, cohort, idego, idalter, matches("^s_"), contains("contiguo")) %>%
  mutate(network = paste(coar, cohort, sep="_")) %>% 
  ungroup() %>%
  dplyr::select(-c(coar, cohort)) %>%
  mutate_all(replace_na, replace = ".") %>%
  mutate(idego = factor(idego, levels = unique(idego)),
         idalter = factor(idalter, levels = unique(idalter)), 
         network_id = as.numeric(factor(network, levels = unique(network))))

test <- communities_wide %>% dplyr::select(coar, cohort, idego, idalter, matches("(any_mut$)|(any_coll$)"), contains("contiguo")) %>%
  filter(idego != idalter) %>%
  mutate(network = paste(coar, cohort, sep="_")) %>% dplyr::select(-c(coar, cohort)) %>%
  drop_na() %>% mutate(idego = factor(idego, levels = unique(idego)),
                       idalter = factor(idalter, levels = unique(idalter)), 
                       network = as.numeric(factor(network, levels = unique(network))))

test %>% mutate(contiguo = reduce(dplyr::select(., starts_with("contiguo")), `+`)) %>% dplyr::select(contiguo) %>% View

lfe.index <- felm(same_edge_betweenness_any_coll ~ contiguo1 + contiguo2 + contiguo3 + contiguo4 + contiguo5 + contiguo6 + contiguo7 + contiguo8 + contiguo9 | index | (0) | network,
                  data=test) 

est <- felm(same_edge_betweenness_any_coll ~ contiguo1 + contiguo2 | index | 0 | network, data = test)

plm.index <- plm(formula= same_edge_betweenness_any_coll ~ factor(contiguo1) + factor(contiguo2) + factor(contiguo3) + factor(contiguo4) + factor(contiguo5) + factor(contiguo6) + factor(contiguo7) + factor(contiguo8) + factor(contiguo9),
                 data=test, model="within", index=c("index"), effect = "individual") 

lm.index <- lm(formula= same_edge_betweenness_any_mut ~ factor(contiguo1) + factor(contiguo2) + factor(contiguo3) + factor(contiguo4) + factor(contiguo5) + factor(contiguo6) + factor(contiguo7) + factor(contiguo8) + factor(contiguo9), data=test)

write_csv(test, path = paste0("test_regressions_", Sys.Date(), ".csv"))
write_csv(communities_wide, path = paste0("regressions_data_", Sys.Date(), ".csv"))
