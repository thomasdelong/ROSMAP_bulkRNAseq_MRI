library(tidyverse)
library(variancePartition)
library(limma)
library(here)
library(yaml)
library(gprofiler2)

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

# save data to run dream with separate file (job submission since its hungry for memory)
rownames(combined_data) <- NULL
write.csv(combined_data, paste0(scratch_path, '/processed_data/07_combined_wholebrain_metadata.csv'))
write.csv(rnaseq_subset_filtered, paste0(scratch_path, '/processed_data/07_rnaseq_subset.csv'))

