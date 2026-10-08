library(gprofiler2)
library(tidyverse)
library(fgsea)
library(yaml)

env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch
code_path <- env$paths$code
#select an example parcel/feature to generate plots for
reference_parcel <- 'middletemporal'
reference_feature <- 'GrayVol_norm'

##pull list of markers from Seurat findallmarkers output
seaad_classlevel <- read.csv(paste0(code_path,"/data/seaad_class_markers.csv"))
seaad_classlevel_topn <- seaad_classlevel %>%
    rename('subclass' = 'cluster') %>%
    filter(p_val_adj == 0 & avg_log2FC > 3 & pct.1 > 0.5) %>% 
    group_by(gene) %>%
    slice_max(pct.1) %>%
    slice_max(avg_log2FC) %>%
    ungroup()
seaad_classlevel_topn_lists <- seaad_classlevel_topn %>%
    split(x = .$gene, f = .$subclass)
#table(seaad_classlevel_topn$subclass)
#print('no duplicates?')
#length(seaad_classlevel_topn$gene) == length(unique(seaad_classlevel_topn$gene))

##make list of gene names from ensembl ids
tt_example <- read.csv(paste0(scratch_path, "/processed_data/toptables/dkt/temporalcortex/superiortemporal.ThickAvg.csv"))
ensembl_ids <- tt_example$X
symbols <- gconvert(ensembl_ids)[c('name', 'input')]

##creating directory/data structure
tt_dir <- paste0(scratch_path, "/processed_data/toptables/dkt/")
#gsea_results_dir <- paste0(scratch_path, "/processed_data/celltypeGSEA/")
#dir.create(gsea_results_dir, recursive = TRUE, showWarnings = FALSE)

regions <- c("dlpfc", "tc", "fc", "hcn", "pcc")
mri_features <- c("CurvInd", "SurfArea_norm", "GrayVol_norm", "ThickAvg")
 
tissue_abbreviations <- c(
  dlpfc = "dorsolateralprefrontalcortex",
  tc = "temporalcortex",
  fc = "frontalcortex",
  pcc = "posteriorcingulatecortex",
  hcn = "Headofcaudatenucleus"
)

##make results dataframes for each cell type
for (cell in unique(seaad_classlevel_topn$subclass)){
    assign(paste0(cell, "_results"), data.frame(rna_region = NA, mri_feature = NA, dkt_parcel = NA, nes = NA, p = NA))
}

##loop through parcels and get NES for each cell type
for (rna_region in regions){
    toptables <- list.files(paste0(tt_dir,tissue_abbreviations[[rna_region]]))
    for (tt in toptables){
        #extract region and feature from file name
        reg <- strsplit(tt,'.', fixed = TRUE)[[1]][[1]]
        feat <- strsplit(tt,'.', fixed = TRUE)[[1]][[2]]
        #read tt and rank genes
        tt_data <- read.csv(paste0(tt_dir,tissue_abbreviations[[rna_region]],"/",tt)) %>%
            mutate(ranking = logFC*-log10(P.Value))
        rankingsdf <- tt_data %>%
            rename('input' = 'X') %>%
            left_join(symbols, by = 'input') %>%
            dplyr::select(c('name', 'ranking'))
        generankings <- rankingsdf[!duplicated(rankingsdf['name']),] %>%
            drop_na() %>%
            remove_rownames() %>%
            column_to_rownames('name') 
        grlist <- as.list(generankings)$ranking
        names(grlist) <- row.names(generankings)
        generankings <- grlist
        #checks in case p is +/- inf
        max_ranking <- max(generankings[is.finite(generankings)])
        min_ranking <- min(generankings[is.finite(generankings)])
        generankings <- replace(generankings, generankings > max_ranking, max_ranking * 10)
        generankings <- replace(generankings, generankings < min_ranking, min_ranking * 10)
        #sort genes by ranking
        generankings <- sort(generankings, decreasing = TRUE) 
        #gsea for all cell types
        GSEAresult <- fgsea(pathways = seaad_classlevel_topn_lists,
                 stats = generankings,
                 scoreType = 'std', 
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) 
        
        for (cell in unique(seaad_classlevel_topn$subclass)){
            cell_result <- filter(GSEAresult, pathway == cell)
            data_to_bind <- data.frame(rna_region = paste0(rna_region), 
                                       mri_feature = paste0(feat), 
                                       dkt_parcel = paste0(reg), 
                                       nes = cell_result$NES,
                                       p = cell_result$pval)
            all_cell_results <- get(paste0(cell, "_results"))
            updated_results <- rbind(all_cell_results, data_to_bind)
            assign(paste0(cell, "_results"), updated_results)
        }

        if (reg == reference_parcel && feat == reference_feature){
            plot <- plotGseaTable(seaad_classlevel_topn_lists, generankings, GSEAresult)
                ggsave(file = paste0(scratch_path, "/figures/supplemental/11_GSEAtable_",rna_region, ".",reg,".",feat,".svg"),
                    plot = plot, width = 4, height = 10, device = svg)
        }
        
    }
}

dir.create(paste0(scratch_path, "/processed_data/gsea_result/"), recursive = TRUE, showWarnings = FALSE)
for (cell in unique(seaad_classlevel_topn$subclass)){
    results <- get(paste0(cell, "_results"))%>% filter(!is.na(rna_region))
    write.csv(results, paste0(scratch_path, "/processed_data/gsea_result/", cell, "_results.csv"))
}