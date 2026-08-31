# Takes the raw count matrices and metadata from ROSMAP and applies normalization, technical covariate/batch correction
library(tidyverse)
library(DESeq2)
library(edgeR)
library(yaml)

env <- read_yaml(here::here("MRI_RNA_Code/01_environment.yml"))
data_path <- env$paths$data
scratch_path <- env$paths$scratch
pdf(file = paste0(scratch_path, '/figures/supplemental/RNASeq_prepocessing.pdf'))

# metadata loading
synapsetospecimen <- read_csv(paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch1_Stranded/ROSMAP_batch1_provenance.csv"))
batch1 <- read.delim(paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch1_Stranded/ROSMAP_batch1_gene_all_counts_matrix.txt")) %>%
  dplyr::slice(-(1:4)) %>%
  dplyr::select(-(48)) %>% #duplicate ID in synapsetospecimen
  rename_with(~ifelse(is.na(match(., synapsetospecimen$id)), ., synapsetospecimen$specimenID[match(., synapsetospecimen$id)]), .cols = everything())
synapsetospecimen <- read_csv(paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch2_Stranded/ROSMAP_batch2_provenance.csv"))
batch2 <- read.delim(paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch2_Stranded/ROSMAP_batch2_gene_all_counts_matrix.txt")) %>%
  dplyr::slice(-(1:4)) %>%
  rename_with(~ifelse(is.na(match(., synapsetospecimen$id)), ., synapsetospecimen$specimenID[match(., synapsetospecimen$id)]), .cols = everything())
rna_combined <- left_join(batch1, batch2, by="feature")
synapsetospecimen <- read_csv(paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch3_Stranded/ROSMAP_batch3_provenance.csv"))
batch3 <- read.delim(paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch3_Stranded/ROSMAP_batch3_gene_all_counts_matrix.txt")) %>%
  dplyr::slice(-(1:4)) %>%
  rename_with(~ifelse(is.na(match(., synapsetospecimen$id)), ., synapsetospecimen$specimenID[match(., synapsetospecimen$id)]), .cols = everything())
rna_combined <- left_join(rna_combined, batch3, by="feature")
synapsetospecimen <- read_csv(paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch4_Stranded/ROSMAP_batch4_provenance.csv"))
batch4 <- read.delim(paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch4_Stranded/ROSMAP_batch4_gene_all_counts_matrix.txt")) %>%
  dplyr::slice(-(1:4)) %>%
  rename_with(~ifelse(is.na(match(., synapsetospecimen$id)), ., synapsetospecimen$specimenID[match(., synapsetospecimen$id)]), .cols = everything())


rna_combined <- left_join(rna_combined, batch4, by="feature")

metadata <- read_csv(paste0(data_path, "/ROSMAP_RNA/Metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv"))
metadata_rnaSeq <- metadata %>%
  filter(assay == "rnaSeq")

# defining looping variables
tissue_list <- c(
  'dlpfc' = 'dorsolateral prefrontal cortex',
  'tc' = 'temporal cortex',
  'fc' = 'frontal cortex',
  'pcc' = 'posterior cingulate cortex',
  'hcn' = 'Head of caudate nucleus'
)

# this loop makes the initial raw count matrices
for (tissue in names(tissue_list)){
  tissue_to_use <- tissue_list[tissue]
  tissue_abbreviation <- tissue
  
  # Process Batch 1
  star_log_file1 <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch1_Stranded/ROSMAP_batch1_Star_Log_Merged.txt")
  provenance_file1 <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch1_Stranded/ROSMAP_batch1_provenance.csv")
  allmetrics1file <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch1_Stranded/ROSMAP_batch1_Study_all_metrics_matrix_clean.txt")
  star_log_df1 <- read.table(star_log_file1, sep = "\t", header = TRUE)
  provenance_df1 <- read.csv(provenance_file1)
  allmetrics1 <- read.delim(allmetrics1file)
  merged_df1 <- merge(star_log_df1, provenance_df1, by.x = "Sample", by.y = "id")%>%
  dplyr::slice(-47) #same duplicate from earlier
  merged_df1$Sample <- NULL

  # Process Batch 2
  star_log_file2 <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch2_Stranded/ROSMAP_batch2_Star_Log_Merged.txt")
  provenance_file2 <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch2_Stranded/ROSMAP_batch2_provenance.csv")
  allmetrics2file <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch2_Stranded/ROSMAP_batch2_Study_all_metrics_matrix_clean.txt")
  star_log_df2 <- read.table(star_log_file2, sep = "\t", header = TRUE)
  provenance_df2 <- read.csv(provenance_file2)
  allmetrics2 <- read.delim(allmetrics2file)
  merged_df2 <- merge(star_log_df2, provenance_df2, by.x = "Sample", by.y = "id")
  all1_2 <- rbind(allmetrics1, allmetrics2)
  merged_df2$Sample <- NULL

  # Process Batch 3
  star_log_file3 <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch3_Stranded/ROSMAP_batch3_Star_Log_Merged.txt")
  provenance_file3 <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch3_Stranded/ROSMAP_batch3_provenance.csv")
  allmetrics3file <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch3_Stranded/ROSMAP_batch3_Study_all_metrics_matrix_clean.txt")
  star_log_df3 <- read.table(star_log_file3, sep = "\t", header = TRUE)
  provenance_df3 <- read.csv(provenance_file3)
  allmetrics3 <- read.delim(allmetrics3file)
  merged_df3 <- merge(star_log_df3, provenance_df3, by.x = "Sample", by.y = "id")
  all1_2_3 <- rbind(all1_2, allmetrics3)
  merged_df3$Sample <- NULL

  # Process Batch 4
  star_log_file4 <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch4_Stranded/ROSMAP_batch4_Star_Log_Merged.txt")
  provenance_file4 <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch4_Stranded/ROSMAP_batch4_provenance.csv")
  allmetrics4file <- paste0(data_path, "/ROSMAP_RNA/Rosmap_Gene_Quantification/Rosmap_Batch4_Stranded/ROSMAP_batch4_Study_all_metrics_matrix_clean.txt")
  star_log_df4 <- read.table(star_log_file4, sep = "\t", header = TRUE)
  provenance_df4 <- read.csv(provenance_file4)
  allmetrics4 <- read.delim(allmetrics4file)
  merged_df4 <- merge(star_log_df4, provenance_df4, by.x = "Sample", by.y = "id")
  all1_2_3_4 <- rbind(all1_2_3, allmetrics4)
  merged_df4$Sample <- NULL

  # Combine all batches into a single data frame
  star_combined <- rbind(merged_df1, merged_df2, merged_df3, merged_df4)
  star_combined <- merge(star_combined, all1_2_3_4, by.x = "specimenID", by.y = "sample")
  metadata_rnaSeq_full <- left_join(metadata_rnaSeq, star_combined, by = 'specimenID')



  metadata_filtered <- as.data.frame(metadata_rnaSeq_full) %>%
  filter(tissue == tissue_to_use) %>% 
  column_to_rownames(var = "specimenID")
  rna_filtered <- rna_combined[, colnames(rna_combined) %in% rownames(metadata_filtered) | colnames(rna_combined) == "feature"] %>% 
  column_to_rownames(var = "feature")
  # Remove version number from ensembl gene ids
  rownames(rna_filtered) = rownames(rna_filtered) %>% substr(., 1, 15) %>% make.names(., unique = T)

  message(
  paste(
    'processed', tissue_to_use, ':\n',
    'Metadata dimensions before: ', paste(dim(metadata_rnaSeq_full), collapse = " x "), '\n',
    'RNASeq dimensions before: ', paste(dim(rna_combined), collapse = " x "), '\n',
    'Metadata dimensions after: ', paste(dim(metadata_filtered), collapse = " x "), '\n',
    'RNASeq dimensions after: ', paste(dim(rna_filtered), collapse = " x ")
  )
  )
  head(metadata_filtered[1:5])
  head(rna_filtered[1:5])
  length(unique(rownames(metadata_filtered)))
  length(rownames(metadata_filtered))
  length(unique(colnames(rna_filtered)))
  length(colnames(rna_filtered))

  write.csv(rna_filtered, paste0(scratch_path, "/processed_data/", tissue_abbreviation, "_RNAseq_counts.csv"))
  write.csv(metadata_filtered, paste0(scratch_path, "/processed_data/", tissue_abbreviation, "_RNAseq_metadata.csv")) 

}

# now that all the data has been combined and formatted, filter out genes that don't meet our cpm threshold in at least three datasets
threshold_mean <- -1 

files <- c(
  pcc   = paste0(scratch_path, "/processed_data/pcc_RNAseq_counts.csv"),
  hcn   = paste0(scratch_path, "/processed_data/hcn_RNAseq_counts.csv"),
  fc    = paste0(scratch_path, "/processed_data/fc_RNAseq_counts.csv"),
  dlpfc = paste0(scratch_path, "/processed_data/dlpfc_RNAseq_counts.csv"),
  tc    = paste0(scratch_path, "/processed_data/tc_RNAseq_counts.csv")
)

compute_gene_stats <- function(path) {
  rnaseq_data <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  y <- DGEList(counts = rnaseq_data)
  logCPM <- cpm(y, log = TRUE)
  tibble(
    gene = rownames(logCPM),
    Mean = rowMeans(logCPM),
    #Variance = apply(logCPM, 1, var),
    #GenePresence = rowSums(rnaseq_data > 0) / ncol(rnaseq_data)
  )
}

stats_list <- imap(files, ~{
  df <- compute_gene_stats(.x)
  df$dataset <- .y
  df
})

all_genes <- stats_list %>% map(~.$gene) %>% purrr::reduce(union)

pass_mat <- sapply(stats_list, function(df) {
  v <- setNames(df$Mean > threshold_mean, df$gene)
  v[all_genes] %>% replace_na(FALSE)
})
rownames(pass_mat) <- all_genes

genes_pass_all <- all_genes[rowSums(pass_mat) == ncol(pass_mat)]

exclusive_by_dataset <- lapply(colnames(pass_mat), function(d) {
  all_genes[pass_mat[, d] & rowSums(pass_mat) == 1]
}) %>% set_names(colnames(pass_mat))

n_all <- length(genes_pass_all)
n_exclusive <- sapply(exclusive_by_dataset, length)
message("Pass in all datasets: ", n_all, '\n Exclusive:')
print(n_exclusive)

pass_count <- rowSums(pass_mat)

# genes_fail_any <- all_genes[rowSums(pass_mat) < ncol(pass_mat)]
# genes_fail_some <- all_genes[rowSums(pass_mat) > 0 & rowSums(pass_mat) < ncol(pass_mat)]
# genes_pass_one <- all_genes[pass_count == 1]
# genes_pass_three <- all_genes[pass_count >= 3]

write.csv(data.frame(gene = genes_pass_all), paste0(scratch_path, "/processed_data/genes_pass_in_all_tissues.csv"), row.names = FALSE)


for (tissue in names(tissue_list)){
  tissue_to_use <- tissue_list[tissue]
  tissue_abbreviation <- tissue
  rna_filtered <- read.csv(paste0(scratch_path, "/processed_data/", tissue_abbreviation, "_RNAseq_counts.csv"), check.names = FALSE, row.names = 1)
  metadata_filtered <- read.csv(paste0(scratch_path, "/processed_data/", tissue_abbreviation, "_RNAseq_metadata.csv"), check.names = FALSE, row.names = 1)

  tech_covar_list <- colnames(metadata_filtered)
  tech_covar_list <- tech_covar_list[!tech_covar_list %in% c("ethnicity", "specimenID", "id","individualID","tissue","organ", "sex", "race", "ageDeath", "diagnosis", "apoeGenotype", "synapseID", "Age_norm",
  "Started.job.on","Started.mapping.on","Finished.on", "thal", "Braak", "CREAD", "CDR", "plaqueMean", "cogdx", "dcfdx_lv", "ceradsc", "braaksc", "msex", "educ", "spanish","apoe_genotype","age_at_visit_max",
  "age_first_ad_dx","age_death","cts_mmse30_first_ad_dx","cts_mmse30_lv", "exclude","excludeReason","isPostMortem","X","BrodmannArea","libraryPreparationMethod","readStrandOrigin", "platform")] 
  print(tech_covar_list) 

  convert_percent_to_numeric <- function(x) {
    as.numeric(gsub("%", "", x)) / 100
  }

  tech_covariates <- data.frame(metadata_filtered)[,tech_covar_list]

  # Define which columns need conversions
  percentage_columns <- c("Uniquely.mapped.reads..", "Mismatch.rate.per.base...", 
                          "Deletion.rate.per.base", "Insertion.rate.per.base",
                          "X..of.reads.mapped.to.multiple.loci", 
                          "X..of.reads.mapped.to.too.many.loci", 
                          "X..of.reads.unmapped..too.many.mismatches", 
                          "X..of.reads.unmapped..too.short", 
                          "X..of.reads.unmapped..other", 
                          "X..of.chimeric.reads")

  character_columns <- c("libraryBatch", "sequencingBatch", "libraryPrep", "runType",
                        "notes", "Study", "projid", "Mapping.speed..Million.of.reads.per.hour",
                        "Uniquely.mapped.reads..", "Mismatch.rate.per.base...", 
                        "Deletion.rate.per.base", "Insertion.rate.per.base",
                        "X..of.reads.mapped.to.multiple.loci", 
                        "X..of.reads.mapped.to.too.many.loci", 
                        "X..of.reads.unmapped..too.many.mismatches", 
                        "X..of.reads.unmapped..too.short", 
                        "X..of.reads.unmapped..other", 
                        "X..of.chimeric.reads")

  # Apply conversions
  tech_covariates[percentage_columns] <- lapply(tech_covariates[percentage_columns], convert_percent_to_numeric)
  tech_covariates[character_columns] <- lapply(tech_covariates[character_columns], as.factor)

  for (tech_cov in tech_covar_list){
    if(is.numeric(tech_covariates[[tech_cov]])){
      tech_covariates[[tech_cov]] <- as.numeric(tech_covariates[[tech_cov]])}
    else{
      tech_covariates[[tech_cov]] <- factor(tech_covariates[[tech_cov]])
      print(tech_cov) 
    }
  }
  sapply(tech_covariates, class)


  rna_filtered_ordered <- rna_filtered[, rownames(tech_covariates)]

  #check if there is any NA in zscored count matrix and metadata
  dds <- DESeqDataSetFromMatrix(countData=rna_filtered_ordered, colData=metadata_filtered, design=~1)
  rld <- vst(dds, blind=FALSE)

  na_per_column <- colSums(is.na(tech_covariates))
  print(na_per_column)
  na_per_row <- rowSums(is.na(tech_covariates))
  print(na_per_row)

  tech_covariates_clean <- tech_covariates[, !colnames(tech_covariates) %in% colnames(tech_covariates)[colSums(is.na(tech_covariates)) > 10]] %>%
    filter(sequencingBatch != '0, 6, 7') 
  if (tissue_abbreviation == 'pcc'){
    tech_covariates_clean <- tech_covariates_clean %>%
      filter(rownames(tech_covariates) != 'RISK_204')}
  tech_covariates_clean <- na.omit(tech_covariates_clean)
  columns_to_keep <- colnames(assay(rld)) %in% rownames(tech_covariates_clean)
  cleaned_corrected_counts <- t(assay(rld)[, columns_to_keep])
  print("cleaned metadata dimensions:")
  dim(tech_covariates_clean)
  print("cleaned rna dimensions:")
  dim(cleaned_corrected_counts)

  #filtering low expression genes
  low_exp_genes <- read.csv(paste0(scratch_path, "/processed_data/genes_pass_in_all_tissues.csv"))
  cleaned_filtered_corrected_counts <- cleaned_corrected_counts[, colnames(cleaned_corrected_counts) %in% low_exp_genes$gene]

  matched_indices <- match(rownames(tech_covariates_clean),rownames(metadata_filtered))
  design_matrix <- metadata_filtered[matched_indices[!is.na(matched_indices)], c('msex','cogdx'), drop = TRUE]

  tech_covariates_clean_filtered <- tech_covariates_clean %>%
    dplyr::select(-sequencingBatch, -assay, -projid, -notes,-isStranded,-runType,-Study,-Number.of.splices..Non.canonical,-X..of.reads.unmapped..too.many.mismatches,-Number.of.chimeric.reads,-X..of.chimeric.reads) %>%
    dplyr::select(AlignmentSummaryMetrics__PCT_PF_READS_ALIGNED, RnaSeqMetrics__PCT_INTERGENIC_BASES, AlignmentSummaryMetrics__PF_READS_ALIGNED, Number.of.input.reads, RnaSeqMetrics__PF_ALIGNED_BASES, AlignmentSummaryMetrics__PF_HQ_ALIGNED_BASES, libraryPrep, RnaSeqMetrics__MEDIAN_5PRIME_BIAS, Mapping.speed..Million.of.reads.per.hour, RnaSeqMetrics__PCT_UTR_BASES, RIN, AlignmentSummaryMetrics__READS_ALIGNED_IN_PAIRS, AlignmentSummaryMetrics__PF_HQ_ERROR_RATE, RnaSeqMetrics__INTRONIC_BASES, Average.mapped.length, X..of.reads.unmapped..too.short, AlignmentSummaryMetrics__PCT_ADAPTER, AlignmentSummaryMetrics__PF_MISMATCH_RATE, AlignmentSummaryMetrics__PCT_READS_ALIGNED_IN_PAIRS, RnaSeqMetrics__CORRECT_STRAND_READS, RnaSeqMetrics__PCT_INTRONIC_BASES)
  tech_covariates_clean_filtered$libraryPrep <- as.integer(factor(tech_covariates_clean_filtered$libraryPrep))
  covars_numeric <- as.data.frame(lapply(tech_covariates_clean_filtered, function(x) as.numeric(as.character(x))))

  result <- limma::removeBatchEffect(t(cleaned_filtered_corrected_counts), batch=as.vector(tech_covariates_clean$sequencingBatch), batch2 = NULL, covariates=covars_numeric, design = design_matrix) 

  write.csv(result, paste0(scratch_path, "/processed_data/", tissue_abbreviation, "_RNAseq_counts_filtered_corrected.csv"))
}