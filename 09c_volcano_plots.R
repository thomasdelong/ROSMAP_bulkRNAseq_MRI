library(yaml)
library(tidyverse)
library(gprofiler2)
library(ggrepel)
library(eulerr)
library(ComplexUpset)
env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

example_for_names <- read.csv(paste0(scratch_path, "/processed_data/toptables/dkt/dorsolateralprefrontalcortex/caudalanteriorcingulate.CurvInd.csv"))
genenames <- gconvert(example_for_names$X, filter_na = FALSE) %>% dplyr::select(target, name)

