library(tidyverse)
library(variancePartition)
library(limma)
library(edgeR)
library(here)
library(yaml)

env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

combined_data <- read.csv(paste0(scratch_path, "/processed_data/07_combined_wholebrain_metadata.csv"), check.names = FALSE)
counts_subset <- read.csv(paste0(scratch_path, "/processed_data/07_rnaseq_counts_subset.csv"), check.names = FALSE, row.names = 1)
dkt <- read.csv(paste0(scratch_path, "/processed_data/DKT_stats_combat_averaged.csv")) %>%
  rename(projid = Subject) %>%
  group_by(projid) %>%
  filter(Session == max(Session)) %>%
  ungroup()

data <- combined_data[, names(combined_data) != ""] %>%
  filter(tissue != "Head of caudate nucleus") %>%
  left_join(dkt[, c("projid", "middletemporal.GrayVol_norm")], by = "projid") %>%
  filter(!is.na(middletemporal.GrayVol_norm), !is.na(cogn_globaln_lv)) %>%
  column_to_rownames("specimenID")
data$tissue <- factor(gsub(" ", "", data$tissue))
data$sequencingBatch <- factor(data$sequencingBatch)

# batches 15/16 overlap with tissue
X <- model.matrix(~ age_diff + age_death + master_brain_side + msex + tissue + sequencingBatch, data)
q <- qr(X)
X <- X[, sort(q$pivot[seq_len(q$rank)])]
X <- X[, colnames(X) != "(Intercept)"]
fixed <- data.frame(projid = data$projid, vol = data$middletemporal.GrayVol_norm, cog = data$cogn_globaln_lv, X, row.names = rownames(data))
covs <- setdiff(names(fixed), c("projid", "vol", "cog"))

dge <- calcNormFactors(DGEList(as.matrix(counts_subset[, rownames(data)])))
vobj <- voomWithDreamWeights(dge, reformulate(c(covs, "(1|projid)")), fixed)

vol <- topTable(eBayes(dream(vobj, reformulate(c("vol", covs, "(1|projid)")), fixed)), coef = "vol", number = Inf, sort.by = "none")
cog <- topTable(eBayes(dream(vobj, reformulate(c("cog", covs, "(1|projid)")), fixed)), coef = "cog", number = Inf, sort.by = "none")

genes <- tibble(t_vol = vol$t, t_cog = cog$t, sig = vol$adj.P.Val < 0.05)
rho <- cor(genes$t_vol, genes$t_cog, method = "spearman")

p <- ggplot(genes, aes(t_vol, t_cog)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point_rast(aes(color = "#33538d"), size = 0.4, alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "#33538d", linewidth = 0.5, linetype = 'dashed') +
  #guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, label = paste0("Spearman's rho = ", rho)) +
  labs(x = "t, middle temporal volume", y = "t, global cognition (last visit)") +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(file.path(scratch_path, "figures", "main", "07f_tc_volume_vs_cog_genes.pdf"), p, width = 5, height = 4)