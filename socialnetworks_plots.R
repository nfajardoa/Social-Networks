########################### NETWORKS AND COMMUNITY IDENTIFICATION ###########################
########################### Plots and histrograms from communities data
########################### Written by: Nicolas Fajardo
########################### Last modified: 26/06/2020

characteristics_summary <- characteristics_idalter %>% pivot_longer(cols = contains(characteristics), names_to = "characteristic", values_to = "nominations") %>%
  group_by(characteristic, add = TRUE) %>% summarize(q3 = quantile(nominations, probs = 0.75))

load_lastdataset()

load("communities/communities_2020-06-26.RData")


network_data_summary_by_community <- communities_membership %>% 
  pivot_longer(cols = all_of(community_algorithms), names_to = "community_algorithm", values_to = "community") %>%
  group_by(time, charc_time, community_algorithm, community, add = TRUE) %>%
  pivot_longer(cols = contains(characteristics), names_to = "characteristic", values_to = "nominations") %>%
  left_join(characteristics_summary, by = c("school", "cohort", "charc_time", "characteristic"))

network_data_summary_by_community_1 <- network_data_summary_by_community %>%
  mutate(perception = if_else(nominations >= 1, 1, 0)) %>% group_by(characteristic, add = TRUE) %>%
  summarize(count = sum(perception)) %>% group_by(variable, selector, time, charc_time, community_algorithm, characteristic) %>%
  mutate(centrality = mean(count), deviation = count - centrality, type = "1")

network_data_summary_by_community_75 <- network_data_summary_by_community %>%
  mutate(perception = if_else(nominations >= q3, 1, 0)) %>% group_by(characteristic, add = TRUE) %>%
  summarize(count = sum(perception)) %>% group_by(variable, selector, time, charc_time, community_algorithm, characteristic) %>%
  mutate(centrality = mean(count), deviation = count - centrality, type = "q3")

network_data_summary_by_community_plot <- bind_rows(network_data_summary_by_community_1, network_data_summary_by_community_75) %>%
  group_by(type, variable, selector, time) %>% 
  mutate(name = paste(variable, selector, time, type, sep = "_")) %>%
  select(name, charc_time, characteristic, community_algorithm, deviation) %>% 
  nest(data = c(charc_time, characteristic, community_algorithm, deviation)) %>%
  ungroup() %>% select(name, data) %>% mutate(plot = map2(name, data, plot_characteristics_by_communities))

cairo_pdf("density.pdf",onefile=T)
network_data_summary_by_community_plot %>% pull(plot) %>% print()
dev.off()

## PLOTS      

communities_lima <- communities %>%
  filter(school == "lima") %>%
  mutate(
    name = paste(variable, selector, school, cohort, time, sep = " "),
    plot = pmap(list(name, network_data, algorithm, communities_n), plot_communities)
  )
null_lima <- expand_grid(distinct(communities_lima, variable, selector, school, cohort, time), algorithm = c("none", "leading_eigen")) %>% 
  anti_join(communities_lima %>% select(variable, selector, school, cohort, time, algorithm) %>% filter(algorithm == "leading_eigen")) %>% 
  add_column(plot = nullGrob() %>% list())

communities_lima %<>% group_by(variable, selector, school, cohort, time) %>% bind_rows(null_lima) %>% arrange(variable, selector, school, cohort, time)
communities_lima_plots <- summarize(communities_lima, table_grob = marrangeGrob(grob = plot, nrow = 2, ncol = 4, top = NULL)) %>% 
  mutate(name = paste(variable, selector, school, cohort, time, sep = "_")) %>% ungroup() %>%
  select(name, table_grob)

# null <- expand_grid(distinct(communities, variable, selector, school, cohort, time), algorithm = c("none", "leading_eigen")) %>% 
#   anti_join(communities %>% select(variable, selector, school, cohort, time, algorithm) %>% filter(algorithm == "leading_eigen")) %>% 
#   add_column(plot = nullGrob() %>% list())
# 
# communities %<>% group_by(variable, selector, school, cohort, time) %>% bind_rows(null) %>% arrange(variable, selector, school, cohort, time)
# communities_plots <- summarize(communities_lima, plot = marrangeGrob(grob = plot, nrow = 2, ncol = 4, top = NULL))

pwalk(communities_lima_plots, ggsave_marrangeGrob, width = 32.385, height = 21.59, device = cairo_pdf, units = "cm")
save.image(file = paste0("social_networks_", Sys.Date(), ".RData"))