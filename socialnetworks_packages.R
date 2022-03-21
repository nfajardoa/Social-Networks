### PACKAGES LOADING AND INSTALLATION

# Package names

packages <- c("furrr", "tictoc", "tidyverse", "magrittr", "tidygraph", "ggraph", "igraph", "reproducible",
              "reshape2", "Rglpk", "grid", "gridExtra", "ggplot2", "extrafont", "qualpalr", "ggthemes",
              "broom", "estimatr", "lfe", "stargazer", "foreign", "quantreg", "dbplyr", "RPostgreSQL", "future.apply",
              "RPostgres", "pryr", "skimr", "multcomp", "caret", "modelr", "tidymodels", "ranger", "glmnet", "ggstance",
              "dplyr", "parsnip", "recipes", "tune", "yardstick", "fastDummies", "gt", "future", "furrr", "scales", "ggalt")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Options
options(future.globals.maxSize = 1.5*1024^3, future.debug = FALSE)
