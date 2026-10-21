library(tidyverse)
library(limma)
library(yaml)
library(patchwork)
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

regions <- c("dlpfc", "tc", "fc", "hcn", "pcc")
mri_features <- c("CurvInd", "SurfArea_norm", "GrayVol_norm", "ThickAvg")

tissue_abbreviations <- c(
  dlpfc = "dorsolateralprefrontalcortex",
  tc = "temporalcortex",
  fc = "frontalcortex",
  pcc = "posteriorcingulatecortex",
  hcn = "Headofcaudatenucleus"
)

plot_rrho <- function(rrho_data, region_abbrev, mri_feature) {
  region_data <- rrho_data %>%
    filter(tissue == tissue_abbreviations[[region_abbrev]], feature == mri_feature) %>%
    mutate(label = paste0("lh_", parcel))

  p <- ggplot() +
    geom_brain(data = region_data,
               mapping = aes(fill = .data$rrho_mean),
               hemi = "left",
               view = c("lateral", "medial"),
               atlas = dkt(),
               position = position_brain(hemi ~ view),
               color = "black") +
    scale_fill_gradient2(low = "#191970", high = "#b22222",
                          limits = c(0, max(rrho_data$rrho_mean)),
                          name = "Max RRHO -log10(p)",
                          na.value = "black") +
    theme_void() +
    theme(legend.position = "none", plot.margin = unit(c(0, 10, 0, 0), units = "points"))

  assign(paste0(region_abbrev, "_", mri_feature, "_plot"), p, envir = parent.frame())
}

title_panel <- function(label, right_margin = FALSE) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold") +
    theme_void()
  if (right_margin) p <- p + theme(plot.margin = unit(c(0, 10, 0, 0), units = "points"))
  p
}

blanktitle <- title_panel("")

row_labels <- c(
  ThickAvg      = "Average\nThickness",
  SurfArea_norm = "Normalized\nArea",
  GrayVol_norm  = "Normalized\nVolume",
  CurvInd       = "Curvature\nIndex"
)

run_ref <- function(ref_tissue, ref_parcel, ref_feature) {
  ref_fullname <- paste(ref_tissue, ref_parcel, ref_feature, sep = '_')
  print(ref_fullname)

  rrho_file <- paste0(scratch_path, "/processed_data/rrho/", ref_fullname, '.csv')
  rrho_data <- read.csv(rrho_file)

  for (region in regions) {
    for (mri_feature in mri_features) {
      plot_rrho(rrho_data, region, mri_feature)
      print(paste0(region, " ", mri_feature))
    }
  }

  col_means <- rrho_data %>%
    group_by(tissue) %>%
    summarise(m = mean(rrho_mean, na.rm = TRUE)) %>%
    deframe()
  row_means <- rrho_data %>%
    group_by(feature) %>%
    summarise(m = mean(rrho_mean, na.rm = TRUE)) %>%
    deframe()

  coltitle <- function(abbrev) {
    title_panel(sprintf("%s\n(%.2f)", toupper(abbrev),
                        col_means[[tissue_abbreviations[[abbrev]]]]),
                right_margin = TRUE)
  }
  rowtitle <- function(feature) {
    title_panel(sprintf("%s\n(%.2f)", row_labels[[feature]], row_means[[feature]]))
  }

  plot <- (blanktitle | coltitle("dlpfc") | coltitle("fc") | coltitle("tc") | coltitle("pcc") | coltitle("hcn")) /
    (rowtitle("ThickAvg") | dlpfc_ThickAvg_plot | fc_ThickAvg_plot | tc_ThickAvg_plot | pcc_ThickAvg_plot | hcn_ThickAvg_plot) /
    (rowtitle("SurfArea_norm") | dlpfc_SurfArea_norm_plot | fc_SurfArea_norm_plot | tc_SurfArea_norm_plot | pcc_SurfArea_norm_plot | hcn_SurfArea_norm_plot) /
    (rowtitle("GrayVol_norm") | dlpfc_GrayVol_norm_plot | fc_GrayVol_norm_plot | tc_GrayVol_norm_plot | pcc_GrayVol_norm_plot | hcn_GrayVol_norm_plot) /
    (rowtitle("CurvInd") | dlpfc_CurvInd_plot | fc_CurvInd_plot | tc_CurvInd_plot | pcc_CurvInd_plot | hcn_CurvInd_plot)

  ggsave(file = paste0(scratch_path, "/figures/supplemental/10_brainmaps_", ref_fullname, ".svg"),
         plot = plot, width = 10, height = 4, device = svg)
}

for (r in seq_len(nrow(refs))) {
  run_ref(refs$tissue[r], refs$parcel[r], refs$feature[r])
}