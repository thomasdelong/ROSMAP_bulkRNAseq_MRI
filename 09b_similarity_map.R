library(yaml)
library(tidyverse)
library(ggseg)
library(ggsegFreeSurfer)
env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

rankmaps_dir <- paste0(scratch_path, "/processed_data/similarity_maps/")
dir.create(rankmaps_dir, recursive = TRUE, showWarnings = FALSE)
pi1maps_dir <- paste0(scratch_path, "/processed_data/pi1maps/")
tt_dir <- paste0(scratch_path, "/processed_data/toptables/dkt/")
tissue_abbreviations <- c(
  dlpfc = "dorsolateralprefrontalcortex",
  tc = "temporalcortex",
  fc = "frontalcortex",
  pcc = "posteriorcingulatecortex",
  hcn = "Headofcaudatenucleus"
)

random_shuffle <- function(x, y, n_perm){
    observed_pearson <- cor(x, y, method = "pearson")

    perm_pearson <- numeric(n_perm)
    for (i in seq_len(n_perm)){
        perm_pearson[i] <- cor(x, sample(y), method = 'pearson')
    }

    p_pearson  <- mean(abs(perm_pearson) >= abs(observed_pearson))

    list(observed_pearson = observed_pearson,
        p_pearson = p_pearson)
}

rank_toptable <- function(tt) {
  tt %>%
    filter(!is.na(P.Value), !is.na(logFC)) %>%
    mutate(rank = -log10(P.Value) * logFC) %>%
    dplyr::select(gene = X, rank)
}

compute_rank_similarity <- function(region_abbrev, mri_feature, n_perm = 1000) {
  dir_path <- paste0(tt_dir, tissue_abbreviations[[region_abbrev]])
  pattern  <- paste0("\\.", mri_feature, "\\.csv$")
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(invisible(NULL))

  # reference = region with max pi1 for this tissue/feature
  pi1_map <- read.csv(paste0(pi1maps_dir, region_abbrev, "_", mri_feature, ".csv"))
  ref_region <- pi1_map$dkt_region[which.max(pi1_map$pi1)]
  ref_file <- files[str_remove(basename(files), pattern) == ref_region]
  ref_ranked <- rank_toptable(read.csv(ref_file))

  results <- map_dfr(files, function(f) {
    dkt_region <- str_remove(basename(f), pattern)

    if (dkt_region == ref_region) {
      return(tibble(dkt_region = dkt_region, r_obs = NA_real_, p_value = NA_real_))
    }

    ranked <- rank_toptable(read.csv(f))

    merged <- inner_join(ref_ranked, ranked, by = "gene", suffix = c("_ref", "_other"))

    x <- merged$rank_ref
    y <- merged$rank_other

    perm_result <- random_shuffle(x, y, n_perm)

    tibble(dkt_region = dkt_region,
           r_obs = perm_result$observed_pearson,
           p_value = perm_result$p_pearson)
  })

  # reference region gets p = 0 and R= 1
  results$p_value[results$dkt_region == ref_region] <- 0
  results$r_obs[results$dkt_region == ref_region] <- 1

    print(results)

  results$label <- paste0("lh_", results$dkt_region)
  write.csv(results, paste0(rankmaps_dir, region_abbrev, "_", mri_feature, ".csv"),
            row.names = FALSE)
  results
}

plot_rank_similarity <- function(region_abbrev, mri_feature) {
  d <- compute_rank_similarity(region_abbrev, mri_feature)
  if (is.null(d)) return(invisible(NULL))
  #d <- d %>% mutate(r_obs = case_when(dkt_region == 'entorhinal' ~ 1, dkt_region != 'entorhinal' ~ 0))
  p_r <- ggplot() +
    geom_brain(data = d,
               mapping = aes(fill = .data$r_obs),
               hemi = "left",
               view = c("lateral", "medial"),
               atlas = dkt(),
               position = position_brain(hemi ~ view),
               color = "black") +
    scale_fill_gradient2(low = "white", mid = "white", high = "purple",
                          midpoint = 0.5, name = "shuffle\nR",
                          na.value = "grey") +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 10, 0, 0), units = "points"))
  ggsave(file = paste0(scratch_path, paste0("/figures/supplemental/09b_",region_abbrev, "_",mri_feature,"_R.svg")),
       plot = p_r, width = 10, height = 4, device = svg)
  assign(paste0(region_abbrev, "_", mri_feature, "_rank_plot"), p_r, envir = parent.frame())

  p_p <- ggplot() +
    geom_brain(data = d,
               mapping = aes(fill = .data$p_value),
               hemi = "left",
               view = c("lateral", "medial"),
               atlas = dkt(),
               position = position_brain(hemi ~ view),
               color = "black") +
    scale_fill_gradient2(low = "purple", mid = "white", high = "white",
                          midpoint = 0.05, name = "shuffle\np",
                          na.value = "grey") +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 10, 0, 0), units = "points"))
  ggsave(file = paste0(scratch_path, paste0("/figures/supplemental/09b_",region_abbrev, "_",mri_feature,"_p.svg")),
       plot = p_p, width = 10, height = 4, device = svg)
  assign(paste0(region_abbrev, "_", mri_feature, "_p_plot"), p_p, envir = parent.frame())
}

plot_rank_similarity('tc', 'ThickAvg')
plot_rank_similarity('pcc', 'CurvInd')