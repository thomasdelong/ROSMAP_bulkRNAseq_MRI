library(tidyverse)
library(variancePartition)
library(limma)
library(here)
library(yaml)
library(gprofiler2)


env <- read_yaml(here::here("01_environment.yml"))
data_path <- env$paths$data
scratch_path <- env$paths$scratch

combined_data <- read.csv(paste0(scratch_path, '/processed_data/07_combined_wholebrain_metadata.csv'), check.names = FALSE)
rnaseq_subset_filtered <- read.csv(paste0(scratch_path, '/processed_data/07_rnaseq_subset.csv'), check.names = FALSE, row.names = 1)

# data <- combined_data %>%
#   dplyr::select(MeanThickness, tissue, projid, age_diff, age_death, master_brain_side, msex, specimenID) %>%
#   column_to_rownames('specimenID')
# data$tissue <- factor(gsub(" ", "", data$tissue))
# rnaseq_all_tissue <- rnaseq_subset_filtered[, rownames(data)]
# formula_mixed <- as.formula("~ MeanThickness + age_diff + age_death + master_brain_side + msex + tissue + (1|projid)")
# fit <- dream(rnaseq_all_tissue, formula_mixed, data)
# fit <- eBayes(fit)

# MeanThickness_alltissue <- topTable(fit, coef = 'MeanThickness', number = Inf)
# write.csv(MeanThickness_alltissue, paste0(scratch_path, "/processed_data/toptables/MeanThickness_alltissue.csv"))

# data <- combined_data[,c('cogdx', 'tissue', 'projid', 'age_diff', 'age_death', 'master_brain_side', 'msex', 'specimenID')] %>%
#   column_to_rownames('specimenID') 
# data$tissue <- gsub(" ", "", data$tissue)
# data$cogdx <- factor(data$cogdx, levels = c(1, 2, 4))
# rnaseq_cogdx <- rnaseq_subset_filtered[, rownames(data)]
# formula <- as.formula(paste0("~ 0 + tissue + cogdx + (1 | projid) + age_diff + age_death + master_brain_side + msex"))
# print("fitting cogdx")
# fit <- dream(rnaseq_cogdx, formula, data)
# fit <- eBayes(fit)

# toptable_AD <- topTable(fit, coef = 'cogdx4', number = Inf)
# write.csv(toptable_AD, paste0(scratch_path, "/processed_data/toptables/AD_all.csv"))

data <- combined_data %>%
  dplyr::select(MeanThickness, tissue, projid, age_diff, age_death, master_brain_side, msex, specimenID) %>%
  column_to_rownames('specimenID') 
data <- data %>%
  filter(tissue != 'Head of caudate nucleus')
data$tissue <- factor(gsub(" ", "", data$tissue))
rnaseq_all_tissue <- rnaseq_subset_filtered[, rownames(data)]
formula_mixed <- as.formula("~ MeanThickness + age_diff + age_death + master_brain_side + msex + tissue + (1|projid)")
fit <- dream(rnaseq_all_tissue, formula_mixed, data)
fit <- eBayes(fit)

MeanThickness_alltissue <- topTable(fit, coef = 'MeanThickness', number = Inf)
write.csv(MeanThickness_alltissue, paste0(scratch_path, "/processed_data/toptables/MeanThickness_cortex.csv"))

data <- combined_data[,c('cogdx', 'tissue', 'projid', 'age_diff', 'age_death', 'master_brain_side', 'msex', 'specimenID')] %>%
  column_to_rownames('specimenID')
data <- data %>%
  filter(tissue != 'Head of caudate nucleus')
data$tissue <- factor(gsub(" ", "", data$tissue))
data$cogdx <- factor(data$cogdx, levels = c(1, 2, 4))
rnaseq_cogdx <- rnaseq_subset_filtered[, rownames(data)]
formula <- as.formula(paste0("~ 0 + tissue + cogdx + (1 | projid) + age_diff + age_death + master_brain_side + msex"))
print("fitting cogdx")
fit <- dream(rnaseq_cogdx, formula, data)
fit <- eBayes(fit)

toptable_AD <- topTable(fit, coef = 'cogdx4', number = Inf)
write.csv(toptable_AD, paste0(scratch_path, "/processed_data/toptables/AD_cortex.csv"))