library(tidyverse)
library(yaml)
library(broom)
library(ggseg)
library(ggsegFreeSurfer)

env <- read_yaml(here::here("01_environment.yml"))
data_path <- env$paths$data
scratch_path <- env$paths$scratch

stats_dir <- file.path(scratch_path, "processed_data", "covariates")
fig_dir <- file.path(scratch_path, "figures", "main")
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

covar <- "ceradsc"
mri_features <- c("ThickAvg", "SurfArea_norm", "GrayVol_norm", "CurvInd")

feature_names <- c(
  ThickAvg      = "Thickness",
  SurfArea_norm = "Surface Area",
  GrayVol_norm  = "Volume",
  CurvInd       = "Curvature Index"
)

dkt <- read.csv(file.path(scratch_path, "processed_data", "DKT_stats_combat_averaged.csv"), check.names = FALSE) %>%
  rename(projid = Subject) %>%
  group_by(projid) %>%
  filter(Session == max(Session)) %>%
  ungroup() %>%
  dplyr::select(-Session)

rosmaster_df <- read.csv(file.path(data_path, "ROSMAP_RNA", "Metadata", "rosmaster.csv"), check.names = FALSE) %>%
  dplyr::select(-1) %>%
  filter(projid %in% dkt$projid, !is.na(.data[[covar]])) %>%
  dplyr::select(projid, all_of(covar))

print(paste0("donors with MRI features: ", n_distinct(dkt$projid)))
print(paste0("donors with MRI features and ", covar, ": ", nrow(rosmaster_df)))

fit_feature <- function(mri_feature) {
  feature_data <- dkt %>%
    dplyr::select(projid, ends_with(paste0(".", mri_feature))) %>%
    rename_with(~ str_remove(.x, fixed(paste0(".", mri_feature))), -projid)

  rosmaster_df %>%
    inner_join(feature_data, by = "projid") %>%
    pivot_longer(cols = -c(projid, all_of(covar)), names_to = "parcel", values_to = "value") %>%
    group_by(parcel) %>%
    group_modify(function(d, key) {
      model <- lm(as.formula(paste("value ~", covar)), data = d)
      coefs <- tidy(model)
      tibble(
        n         = nobs(model),
        beta      = coefs$estimate[coefs$term == covar],
        p_value   = coefs$p.value[coefs$term == covar],
        intercept = coefs$estimate[coefs$term == "(Intercept)"],
        r_squared = glance(model)$r.squared
      )
    }) %>%
    ungroup() %>%
    mutate(p_fdr = p.adjust(p_value, method = "BH"),
           feature = mri_feature,
           label = paste0("lh_", parcel))
}

region_stats <- map(set_names(mri_features), fit_feature)

for (mri_feature in mri_features) {
  write.csv(region_stats[[mri_feature]],
            file.path(stats_dir, paste0(covar, "_", mri_feature, ".csv")),
            row.names = FALSE)
  print(paste0(mri_feature, ": beta range ", paste(signif(range(region_stats[[mri_feature]]$beta), 3), collapse = " to "),
               ", ", sum(region_stats[[mri_feature]]$p_fdr < 0.05), " parcels FDR < 0.05"))
}

plot_cognition_map <- function(stats, mri_feature) {
  lim <- max(abs(stats$beta), na.rm = TRUE)

  ggplot() +
    geom_brain(data = stats,
               mapping = aes(fill = .data$beta),
               hemi = "left",
               view = c("lateral", "medial"),
               atlas = dkt(),
               position = position_brain(hemi ~ view),
               color = "black") +
    scale_fill_gradient2(low = '#0f12d0', mid = "white", high = "#b22222", midpoint = 0,
                         limits = c(-lim, lim),
                         name = "Beta",
                         na.value = "grey80",
                         guide = guide_colorbar(direction = "horizontal", barwidth = 10)) +
    labs(title = paste0(feature_names[[mri_feature]], " vs. ", covar)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "bottom")
}

for (mri_feature in mri_features) {
  p <- plot_cognition_map(region_stats[[mri_feature]], mri_feature)
  ggsave(file = file.path(fig_dir, paste0("07e_brainmap_", covar, "_", mri_feature, ".svg")),
         plot = p, width = 10, height = 4, device = svg)
}