library(RRHO2)
library(yaml)
library(tidyverse)
library(parallel)

env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))
message("Using ", n_cores, " cores")

refs <- "tissue, parcel, feature
dorsolateralprefrontalcortex, parahippocampal, CurvInd
frontalcortex, lateraloccipital, ThickAvg
temporalcortex, parsorbitalis, SurfArea_norm
Headofcaudatenucleus, lateraloccipital, GrayVol_norm
"
refs <- read.csv(text = refs, strip.white = TRUE)
refs$fullname <- paste(refs$tissue, refs$parcel, refs$feature, sep = "_")

tt_dir <- paste0(scratch_path, "/processed_data/toptables/dkt/")
out_dir <- paste0(scratch_path, "/processed_data/rrho/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mri_features <- c("CurvInd", "SurfArea_norm", "GrayVol_norm", "ThickAvg")
regions <- c("dorsolateralprefrontalcortex", "temporalcortex", "frontalcortex",
             "posteriorcingulatecortex", "Headofcaudatenucleus")

read_ranked <- function(path) {
    read.csv(path) %>%
        mutate(rank = -log10(P.Value) * logFC) %>%
        dplyr::select(gene = X, rank) %>%
        arrange(desc(rank))
}

#load refs
ref_tts <- setNames(
    lapply(seq_len(nrow(refs)), function(r) {
        read_ranked(paste0(tt_dir, refs$tissue[r], "/", refs$parcel[r], ".", refs$feature[r], ".csv"))
    }),
    refs$fullname
)

# list of tasks , one per reference
tasks <- expand_grid(region = regions, feature = mri_features) %>%
    mutate(file = map2(region, feature, ~ list.files(
        paste0(tt_dir, .x, "/"), pattern = paste0(.y, "\\.csv$"), full.names = TRUE))) %>%
    unnest(file) %>%
    mutate(parcel = str_remove(basename(file), paste0("\\.", feature, "\\.csv$")))
message(nrow(tasks), " comparison files x ", nrow(refs), " references")

#compare every parcel to reference
run_task <- function(i) {
    t <- tasks[i, ]
    tt2 <- read_ranked(t$file)
    map_dfr(names(ref_tts), function(ref_name) {
        res <- tryCatch({
            rr <- RRHO2_initialize(ref_tts[[ref_name]], tt2, log10.ind = TRUE, method = "hyper")
            c(max(rr$hypermat, na.rm = TRUE), mean(rr$hypermat, na.rm = TRUE))
        }, error = function(e) c(NA_real_, NA_real_))
        data.frame(ref = ref_name, tissue = t$region, parcel = t$parcel, feature = t$feature,
                   rrho_max = res[1], rrho_mean = res[2])
    })
}

results <- mclapply(seq_len(nrow(tasks)), run_task,
                    mc.cores = n_cores, mc.preschedule = FALSE)

#failchecker
failed <- vapply(results, inherits, logical(1), "try-error")
if (any(failed)) warning(sum(failed), " tasks failed: ", paste(tasks$file[failed], collapse = ", "))
results <- bind_rows(results[!failed])

#make output csv for each reference parcel
for (ref_name in names(ref_tts)) {
    results %>%
        filter(ref == ref_name) %>%
        dplyr::select(-ref) %>%
        write.csv(paste0(out_dir, ref_name, ".csv"), row.names = FALSE)
}