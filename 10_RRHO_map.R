library(RRHO2)
library(yaml)
library(tidyverse)
library(patchwork)
library(circlize)
library(ggseg)
library(ggsegFreeSurfer)
env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

refs <- "tissue, parcel, feature
dorsolateralprefrontalcortex, parahippocampal, CurvInd
frontalcortex, lateraloccipital, ThickAvg
temporalcortex, parsorbitalis, SurfArea_norm
Headofcaudatenucleus, pericalcarine, ThickAvg
"
refs <- read.csv(text = refs, strip.white = TRUE)

tt_dir <- paste0(scratch_path, "/processed_data/toptables/dkt/")
mri_features <- c("CurvInd", "SurfArea_norm", "GrayVol_norm", "ThickAvg")
regions <- c("dorsolateralprefrontalcortex","temporalcortex","frontalcortex","posteriorcingulatecortex","Headofcaudatenucleus")

compare_RRHO <- function(ref_toptable, comparison){
    tt2 <- read.csv(comparison) %>%
        mutate(rank = -log10(P.Value) * logFC) %>%
        dplyr::select(gene = X, rank) %>%
        arrange(desc(rank))

    rrho_output <- RRHO2_initialize(ref_toptable, tt2,
        log10.ind = TRUE, method = "hyper")
    return(data.frame(rrho_max = max(rrho_output$hypermat, na.rm = TRUE),
        rrho_mean = mean(rrho_output$hypermat, na.rm = TRUE)
    ))
}

run_ref <- function(ref_tissue, ref_parcel, ref_feature) {
    ref_fullname <- paste(ref_tissue, ref_parcel, ref_feature, sep = '_')
    print(ref_fullname)
    rrho_file <- paste0(scratch_path, "/processed_data/rrho/", ref_fullname, '.csv')
    dir.create(dirname(rrho_file), recursive = TRUE, showWarnings = FALSE)

    ref_tt <- read.csv(paste0(tt_dir, ref_tissue, '/', ref_parcel, '.', ref_feature, '.csv')) %>%
        mutate(rank = -log10(P.Value) * logFC) %>%
        dplyr::select(gene = X, rank) %>%
        arrange(desc(rank))

    output <- data.frame(tissue = NA, parcel = NA, feature = NA, rrho_max = NA, rrho_mean = NA)
    for (region in regions) {
        for (mri_feature in mri_features) {
            comparison_tts <- list.files(paste0(tt_dir, region, '/'), pattern = paste0(mri_feature, "\\.csv$"), full.names = TRUE)
            for (file in comparison_tts) {
                parcel <- str_remove(basename(file), paste0("\\.", mri_feature, "\\.csv$"))
                print(paste(region, parcel, mri_feature))
                rrho_result <- compare_RRHO(ref_tt, file)
                print(rrho_result)
                result <- data.frame(tissue = region, parcel = parcel, feature = mri_feature, rrho_result)
                output <- rbind(output, result)
            }
        }
    }

    output <- dplyr::filter(output, !is.na(tissue))
    write.csv(output, rrho_file, row.names = FALSE)
}

for (r in seq_len(nrow(refs))) {
    run_ref(refs$tissue[r], refs$parcel[r], refs$feature[r])
}