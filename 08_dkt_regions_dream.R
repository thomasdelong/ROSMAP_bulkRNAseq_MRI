library(limma)
library(tidyverse)
library(variancePartition)
library(edgeR)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
rna_tissue <- args[1]
region_feature <- args[2]
n_perm <- if (length(args) >= 3) as.integer(args[3]) else 0

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

counts <- read.csv(paste0(scratch_path, "/processed_data/07_rnaseq_counts_subset.csv"), row.names = 1, check.names = FALSE)
data <- data[!is.na(data[[region_feature]]), ]
data$x <- data[[region_feature]]
data$sequencingBatch <- factor(data$sequencingBatch); data$libraryPrep <- factor(data$libraryPrep)

# Sequencing batch goes in the model when this tissue has more than one. Library prep does not: every sequencing batch
# used a single library prep, so once batch is in the model, prep adds nothing (its effect is already part of the batch
# effect). The check below stops the script if that ever stops being true.
batches_with_two_preps <- data %>% distinct(sequencingBatch, libraryPrep) %>% count(sequencingBatch) %>% filter(n > 1)
stopifnot("a sequencing batch used more than one library prep; add libraryPrep to covs" = nrow(batches_with_two_preps) == 0)
covs <- c("age_diff", "age_death", "master_brain_side", "msex")
if (n_distinct(data$sequencingBatch) > 1) covs <- c(covs, "sequencingBatch")

# voom precision weights with one quality weight per specimen, estimated from the design WITHOUT the MRI term. The weights
# therefore cannot adapt to the MRI feature, they are the same for every feature of this tissue, and the same weights
# can be reused for permutations of the MRI table. RIN is not in the model: it tracks the neuronal content of the sample
# (the signal itself), so adjusting for it removes part of what the MRI association measures.
null_design <- model.matrix(reformulate(covs), data)          # covariates only
full_design <- model.matrix(reformulate(c("x", covs)), data)  # covariates plus the MRI feature
dge <- calcNormFactors(DGEList(as.matrix(counts[, rownames(data)])))
v <- voomWithQualityWeights(dge, null_design)
fit <- eBayes(lmFit(v, full_design))
tt <- topTable(fit, coef = "x", number = Inf)
out_dir <- paste0(scratch_path, "/processed_data/toptables/dkt/", rna_tissue, "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(tt, paste0(out_dir, region_feature, ".csv"))
if (n_perm > 0) {
  # Shuffle MRI rows across donors. Built from the same sorted donor list and seed in every job,
  # so permutation b is identical for every tissue and every MRI feature (shared nulls).
  donors <- sort(unique(combined_data$projid))
  set.seed(6543)
  perms <- replicate(n_perm, sample(donors))   # perms[i, b]: donor whose MRI row donor i receives
  rows <- match(data$projid, donors)
  null_pi1 <- vapply(seq_len(n_perm), function(b) {
    xp <- dkt[[region_feature]][match(perms[rows, b], dkt$projid)]
    ok <- !is.na(xp)
    d <- full_design[ok, , drop = FALSE]; d[, "x"] <- xp[ok]
    1 - propTrueNull(eBayes(lmFit(v[, ok], d))$p.value[, "x"]) 
  }, numeric(1))
  dir.create(paste0(out_dir, "perm"), showWarnings = FALSE)
  write.csv(data.frame(perm = seq_len(n_perm), pi1 = null_pi1),
            paste0(out_dir, "perm/", region_feature, ".csv"), row.names = FALSE)
}