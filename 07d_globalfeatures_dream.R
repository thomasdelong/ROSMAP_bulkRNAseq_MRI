library(tidyverse)
library(variancePartition)
library(limma)
library(edgeR)
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

counts_subset <- read.csv(paste0(scratch_path, '/processed_data/07_rnaseq_counts_subset.csv'), check.names = FALSE, row.names = 1)

# Cortex-wide fits with donor as a random effect. Sequencing batch goes in the model as a fixed effect and voom weights every
# observation (voomWithDreamWeights, from the design WITHOUT the term being tested), instead of removing the batch effect from
# a VST matrix first. Library prep is left out because every sequencing batch used a single prep (checked in fit_cortex), so
# batch already carries it.
#
# Batch and tissue still overlap: batches 15 and 16 hold only frontal and temporal cortex samples, and those two tissues were
# sequenced only in batches 15 and 16. Knowing a sample is in batch 15 or 16 therefore already says it is frontal or temporal
# cortex, so one column of the fixed-effect matrix repeats information in the others. dream stops on such a column.
# drop_redundant() removes columns that are exact combinations of earlier ones (the tested term is placed first, so it is
# never the one removed) and prints what it removed.
drop_redundant <- function(X, keep) {
  X <- X[, c("(Intercept)", keep, setdiff(colnames(X), c("(Intercept)", keep))), drop = FALSE]
  q <- qr(X)
  redundant <- colnames(X)[q$pivot[-seq_len(q$rank)]]
  if (length(redundant) > 0) message("redundant fixed-effect columns dropped: ", paste(redundant, collapse = ", "))
  X[, setdiff(colnames(X), c("(Intercept)", redundant)), drop = FALSE]
}
fit_cortex <- function(term) {
  data <- combined_data[, names(combined_data) != ""] %>% filter(tissue != 'Head of caudate nucleus') %>% column_to_rownames('specimenID')   # drop 07's unnamed index column
  data$tissue <- factor(gsub(" ", "", data$tissue)); data$sequencingBatch <- factor(data$sequencingBatch)
  stopifnot("a sequencing batch used more than one library prep" = all(count(distinct(data, sequencingBatch, libraryPrep), sequencingBatch)$n == 1))
  if (term == "cogdx") data$cogdx <- factor(data$cogdx, levels = c(1, 2, 4))
  data <- data[!is.na(data[[term]]), ]
  Xf <- model.matrix(as.formula(paste("~", term, "+ age_diff + age_death + master_brain_side + msex + tissue + sequencingBatch")), data)
  keep <- grep(paste0("^", term), colnames(Xf), value = TRUE); Xf <- drop_redundant(Xf, keep)
  colnames(Xf) <- make.names(colnames(Xf)); fixed <- cbind(data[, "projid", drop = FALSE], Xf)
  tested <- make.names(keep); null_terms <- setdiff(colnames(Xf), tested)
  dge <- calcNormFactors(DGEList(as.matrix(counts_subset[, rownames(fixed)])))
  vobj <- voomWithDreamWeights(dge, as.formula(paste("~", paste(c(null_terms, "(1|projid)"), collapse = " + "))), fixed)
  fit <- eBayes(dream(vobj, as.formula(paste("~", paste(c(colnames(Xf), "(1|projid)"), collapse = " + "))), fixed))
  list(fit = fit, coef = tested)
}

f <- fit_cortex("MeanThickness")
write.csv(topTable(f$fit, coef = "MeanThickness", number = Inf), paste0(scratch_path, "/processed_data/toptables/MeanThickness_cortex.csv"))
f <- fit_cortex("CortexVol_norm")
write.csv(topTable(f$fit, coef = "CortexVol_norm", number = Inf), paste0(scratch_path, "/processed_data/toptables/CortexVol_norm_cortex.csv"))
print("fitting cogdx")
f <- fit_cortex("cogdx")
write.csv(topTable(f$fit, coef = "cogdx4", number = Inf), paste0(scratch_path, "/processed_data/toptables/AD_cortex.csv"))
