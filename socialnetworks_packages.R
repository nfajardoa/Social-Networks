### PACKAGES LOADING AND INSTALLATION

# Package names

packages <- c("furrr", "tictoc", "tidyverse", "magrittr", "tidygraph", "ggraph", "igraph", 
              "reshape2", "Rglpk", "grid", "gridExtra", "ggplot2", "extrafont", "qualpalr", "ggthemes",
              "broom", "estimatr", "lfe", "stargazer", "foreign", "quantreg", "gbm", "glmnet",
              "MASS", "rpart", "doParallel", "sandwich", "randomForest",
              "nnet", "matrixStats", "xtable", "readstata13", "car", "lfe", "doParallel",
              "caret", "foreach", "multcomp","cowplot", "neuralnet", "dbplyr", "RPostgreSQL", "RPostgres", "pryr") 

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Options
options(future.globals.maxSize = 891289600)

