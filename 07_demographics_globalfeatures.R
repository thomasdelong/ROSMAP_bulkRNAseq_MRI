library(tidyverse)
library(variancePartition)
library(limma)
library(here)
library(yaml)
library(gprofiler2)
library(ggrepel)
library(eulerr)
#library(ComplexUpset)
#library(ComplexHeatmap)

env <- read_yaml(here::here("01_environment.yml"))
data_path <- env$paths$data
scratch_path <- env$paths$scratch

# load in rna data and merge 
dlpfc <- read.csv(paste0(scratch_path, "/processed_data/dlpfc_RNAseq_counts_filtered_corrected.csv"), 
                        row.names = 1, check.names = FALSE)
tc <- read.csv(paste0(scratch_path, "/processed_data/tc_RNAseq_counts_filtered_corrected.csv"), 
                        row.names = 1, check.names = FALSE)
fc <- read.csv(paste0(scratch_path, "/processed_data/fc_RNAseq_counts_filtered_corrected.csv"), 
                        row.names = 1, check.names = FALSE)
pcc <- read.csv(paste0(scratch_path, "/processed_data/pcc_RNAseq_counts_filtered_corrected.csv"), 
                        row.names = 1, check.names = FALSE)
hcn <- read.csv(paste0(scratch_path, "/processed_data/hcn_RNAseq_counts_filtered_corrected.csv"), 
                        row.names = 1, check.names = FALSE)

rnaseq_df_unfiltered <- merge(dlpfc, tc, by=0) %>% 
    column_to_rownames('Row.names') %>%
    merge(fc, by=0) %>% 
    column_to_rownames('Row.names')%>%
    merge(pcc, by=0) %>% 
    column_to_rownames('Row.names')%>%
    merge(hcn, by=0) %>%
    column_to_rownames('Row.names')

print('merged?')
ncol(rnaseq_df_unfiltered) == ncol(dlpfc)+ncol(tc)+ncol(fc)+ncol(pcc)+ncol(hcn)
# loading just the mri data with rna, then subsetting rna to people with mris
projid_mapping <- read.csv(paste0(data_path, "/ROSMAP_RNA/Metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv"))
projid_mapping <- subset(projid_mapping, assay == 'rnaSeq' & specimenID %in% colnames(rnaseq_df_unfiltered))[c('specimenID','projid', 'tissue')] 
rosmaster_df <- read.csv(paste0(data_path, "/ROSMAP_RNA/Metadata/rosmaster.csv"))
side <- read.csv(paste0(data_path, "/ROSMAP_RNA/Metadata/rosmapAD_data_phenotypes_masterBrainSide.csv"))
confounds_df <- read.csv(paste0(data_path, "/ROSMAP_MRI/rosmapAD_data_imaging_confounds_ageScan.csv"))
wholebrain_all <- read.csv(paste0(scratch_path, "/processed_data/DKT_stats_combat_global.csv")) %>%
    rename(projid = Subject)
projid_mapping_subset <- subset(projid_mapping, projid %in% wholebrain_all$projid)

rnaseq_df <- rnaseq_df_unfiltered[, colnames(rnaseq_df_unfiltered) %in% projid_mapping_subset$specimenID]

print('numbers correct?')
paste0('sample ids: ', length(unique(projid_mapping_subset$specimenID)))
paste0('sample ids in RNA data: ', length(unique(colnames(rnaseq_df))))

# Calculate the age difference
merged_age <- merge(rosmaster_df[, c("projid", "age_death", "msex", "cogdx", "gpath", "cogn_global_random_slope", "educ")],
                    confounds_df[, c("projid", "age_scan")], by = "projid") %>%
              merge(side[,c('projid', 'master_brain_side')], by = 'projid')
merged_age <- merged_age %>%
  group_by(projid) %>% 
  filter(age_scan == max(age_scan)) %>%
  ungroup() %>%
  filter(projid %in% wholebrain_all$projid)
merged_age$age_diff <- merged_age$age_death - merged_age$age_scan
wholebrain <- subset(wholebrain_all, projid %in% projid_mapping_subset$projid) %>%
  group_by(projid) %>%
  filter(Session == max(Session)) %>%
  ungroup()
combined_data <- merge(merged_age, wholebrain, by = "projid") %>%
  right_join(projid_mapping_subset, by = 'projid') %>%
  drop_na(msex) %>%
  subset(cogdx == 1 | cogdx == 2 | cogdx == 4)

# Reorder the RNAseq columns to match the order in combined_data for dream
rnaseq_subset_filtered <- rnaseq_df[, combined_data$specimenID]

# Get all features from the radiomics data
features <- colnames(wholebrain)

print('total RNA samples:')
length(projid_mapping$specimenID)
print('total individuals (RNA batch correction):')
length(unique(projid_mapping$projid))
print('total individuals:')
print(length(unique(combined_data$projid)))
print('samples per tissue:')
table(combined_data$tissue)
tissues <- c('dorsolateral prefrontal cortex', 'frontal cortex', 'temporal cortex', 'posterior cingulate cortex', 'Head of caudate nucleus')
for (i in tissues){
    print(paste(i, 'hemisphere proportions'))
    print(prop.table(table(subset(combined_data, tissue == paste0(i))$master_brain_side)))
    print("______________________________________________________________________________________________________________________")
}
projid_unique <- combined_data %>% group_by(projid) %>% slice_max(tissue)
print('mean age difference:')
mean(projid_unique$age_diff)
print('sd age difference:')
sd(projid_unique$age_diff)
print('mean age death:')
mean(projid_unique$age_death)
print('sd age death:')
sd(projid_unique$age_death)
print('sex:')
table(projid_unique$msex)
print('cogdx:')
table(projid_unique$cogdx)

p <- ggplot(projid_unique, aes(x = age_death, fill = factor(msex), group = msex)) +
    geom_histogram(binwidth = 2, position = 'stack')+
    theme_minimal()+
    labs(x = "Age at Death", y = paste0("Count (Total = ",length(projid_unique$projid)," Individuals)"), fill = "Sex")+
    scale_fill_manual(values = c('#f23b26', '#3339e8')) 

ggsave(file=paste0(scratch_path,"/figures/main/07_agedeath_sex.svg"), plot=p, width=10, height=4, device = svg)

#https://www.radc.rush.edu/docs/var/detail.htm?category=Clinical+Diagnosis&subcategory=Final+consensus+diagnosis&variable=cogdx
cogdx <- as.data.frame(table(projid_unique$cogdx)) %>%
  rename(group_code = 1, count = 2) %>%
  mutate(group = case_when(group_code == 1 ~ "Control", group_code == 2 ~ "MCI", group_code == 4 ~ "AD"))

p <- ggplot(cogdx, aes(x = '', y = count, fill = group)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar(theta = "y", start=0) +
  theme_minimal()+
  scale_fill_manual(values = c('#D93911','#011638','#FFBD23')) #ad, control, mci

ggsave(file=paste0(scratch_path,"/figures/main/07_diagnosis.svg"), plot=p, width=4, height=4, device = svg)

# # RUN DREAM FOR VOLUME, THICKNESS, AD
# rownames(combined_data) <- NULL
# # Volume
# data <- combined_data[,c('CortexVol_norm', 'tissue', 'projid', 'age_diff', 'age_death', 'master_brain_side', 'msex', 'specimenID')] %>%
#       column_to_rownames('specimenID') %>%
#       filter(tissue == 'dorsolateral prefrontal cortex')
# data$tissue <- gsub(" ", "", data$tissue)
# rnaseq_dlpfc <- rnaseq_subset_filtered[, rownames(data)]
# #define seperately, model.matrix has issues with the random effect
# formula <- as.formula(paste0("~ CortexVol_norm + age_diff + age_death + master_brain_side + msex"))
# # Fit the linear model using limma
# fit <- dream(rnaseq_dlpfc, formula, data)
# fit <- eBayes(fit)
# print('CortexVol_norm term, no interaction')
# tt_volume <- topTable(fit, coef = 'CortexVol_norm', number = Inf)
# write.csv(tt_volume, paste0(scratch_path, "/processed_data/toptables/CortexVol_norm.csv"))

# # Thickness 
# data <- combined_data[,c('MeanThickness', 'tissue', 'projid', 'age_diff', 'age_death', 'master_brain_side', 'msex', 'specimenID')] %>%
#   column_to_rownames('specimenID') %>%
#   filter(tissue == 'dorsolateral prefrontal cortex')
# data$tissue <- gsub(" ", "", data$tissue)

# formula <- as.formula("~ MeanThickness + age_diff + age_death + master_brain_side + msex")
# fit <- dream(rnaseq_dlpfc, formula, data)
# fit <- eBayes(fit)

# tt_thickness <- topTable(fit, coef = 'MeanThickness', number = Inf)
# write.csv(tt_thickness, paste0(scratch_path, "/processed_data/toptables/MeanThickness.csv"))

# #AD
# data <- combined_data[,c('cogdx', 'tissue', 'projid', 'age_diff', 'age_death', 'master_brain_side', 'msex', 'specimenID')] %>%
#   column_to_rownames('specimenID') %>%
#   filter(tissue == 'dorsolateral prefrontal cortex')
# data$tissue <- gsub(" ", "", data$tissue)
# data$cogdx <- factor(data$cogdx, levels = c(1, 2, 4))

# formula_cogdx <- as.formula("~ cogdx + age_diff + age_death + master_brain_side + msex")
# fit_cogdx <- dream(rnaseq_dlpfc, formula_cogdx, data)
# fit_cogdx <- eBayes(fit_cogdx)

# toptable_AD <- topTable(fit_cogdx, coef = 'cogdx4', number = Inf)
# write.csv(toptable_AD, paste0(scratch_path, "/processed_data/toptables/AD.csv"))


# volume_gene <- "FUZ" 
# gene_lookup <- gconvert(volume_gene, organism = "hsapiens", target = "ENSG")
# gene_id <- gene_lookup$target[1]
# toptable_volume_labeled <- toptable_volume %>%
#   rownames_to_column('input') %>%
#   mutate(label = ifelse(input == gene_id, volume_gene, NA))

# p <- ggplot(toptable_volume_labeled, aes(x = logFC, y = -log10(P.Value), color = BH_1)) +
#   geom_point() +
#   ggtitle("") +
#   labs(color = 'pFDR < 0.1') +
#   scale_color_manual(values = c('gray', '#1683A6')) +
#   ylim(0, 8) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         text = element_text(family = "roboto", size = 12)) +
#   geom_text_repel(aes(label = label), show.legend = FALSE, box.padding = 0.75,
#                    max.overlaps = Inf, na.rm = TRUE) +
#   ggtitle('Volume')
# p

# ggsave(file = paste0(scratch_path, '/figures/main/07_volume_volcano.svg'), plot = p, width = 4, height = 4, device = svg)

# thickness_gene <- "HES5" 
# gene_lookup <- gconvert(thickness_gene, organism = "hsapiens", target = "ENSG")
# gene_id <- gene_lookup$target[1]
# toptable_thickness_labeled <- toptable_thickness %>%
#   rownames_to_column('input') %>%
#   mutate(label = ifelse(input == gene_id, thickness_gene, NA))

# p <- ggplot(toptable_thickness_labeled, aes(x = logFC, y = -log10(P.Value), color = BH_1)) +
#   geom_point() +
#   ggtitle("") +
#   labs(color = 'pFDR < 0.1') +
#   scale_color_manual(values = c('gray', '#2D6E2E')) +
#   ylim(0, 8) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         text = element_text(family = "roboto", size = 12)) +
#   geom_text_repel(aes(label = label), show.legend = FALSE, box.padding = 0.75,
#                    max.overlaps = Inf, na.rm = TRUE) +
#   ggtitle('Global Cortical Thickness')
# p

# ggsave(file = paste0(scratch_path, '/figures/main/07_thickness_volcano.svg'), plot = p, width = 4, height = 4, device = svg)
data_all_tissue <- combined_data %>%
  dplyr::select(CortexVol_norm, tissue, projid, age_diff, age_death, master_brain_side, msex, specimenID) %>%
  column_to_rownames('specimenID')
data_all_tissue$tissue <- factor(gsub(" ", "", data_all_tissue$tissue))

# rnaseq_subset_filtered was already built from all tissues earlier in the script -
# just reuse it directly instead of re-subsetting to DLPFC
rnaseq_all_tissue <- rnaseq_subset_filtered[, rownames(data_all_tissue)]
formula_mixed <- as.formula("~ CortexVol_norm + age_diff + age_death + master_brain_side + msex + tissue + (1|projid)")

fit <- dream(rnaseq_all_tissue, formula_mixed, data_all_tissue)
fit <- eBayes(fit)

CortexVol_norm_alltissue <- topTable(fit, coef = 'CortexVol_norm', number = Inf)
write.csv(CortexVol_norm_alltissue, paste0(scratch_path, "/processed_data/toptables/CortexVol_norm_alltissue.csv"))