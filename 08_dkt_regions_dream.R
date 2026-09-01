library(limma)
library(tidyverse)
library(variancePartition)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
rna_tissue <- args[1]
region_feature <- args[2]

env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

combined_data <- read.csv(paste0(scratch_path, "/processed_data/07_combined_wholebrain_metadata.csv"), check.names = FALSE)
rnaseq_subset_filtered <- read.csv(paste0(scratch_path, "/processed_data/07_rnaseq_subset.csv"), row.names = 1, check.names = FALSE)
dkt <- read.csv(paste0(scratch_path, "/processed_data/DKT_stats_combat_averaged.csv")) %>%
    rename(projid = Subject) %>%
    group_by(projid) %>%
    filter(Session == max(Session)) %>%
    ungroup()


data <- combined_data[, -1] %>% # the [,-1] is just indices, it won't filter with that column since its unnamed
  filter(gsub(" ", "", tissue) == rna_tissue) %>%
  left_join(dkt[, c('projid', region_feature)], by = 'projid') %>%
  column_to_rownames('specimenID')

rna_ordered <- rnaseq_subset_filtered[, rownames(data)] 

formula <- as.formula(paste0("~ `", region_feature, "` + age_diff + age_death + master_brain_side + msex"))
fit <- dream(rna_ordered, formula, data)
fit <- eBayes(fit)

tt <- topTable(fit, coef = paste0(region_feature), number = Inf)
out_dir <- paste0(scratch_path, "/processed_data/toptables/dkt/", rna_tissue, "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(tt, paste0(out_dir, region_feature, ".csv"))