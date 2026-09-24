library(tidyverse)
library(limma)
library(yaml)
library(patchwork)
library(ggseg)
library(ggsegFreeSurfer)
env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

source(here::here("10_reference.R"))

rrho_file <- paste0(scratch_path, paste0("/processed_data/rrho/", ref_fullname, '.csv')) 
rrho_data <- read.csv(rrho_file) #%>%
    # mutate(rrho_max = if_else(tissue == ref_tissue & parcel == ref_parcel & feature == ref_feature,
    #                          NA, rrho_max))


regions <- c("dlpfc", "tc", "fc", "hcn", "pcc")
mri_features <- c("CurvInd", "SurfArea_norm", "GrayVol_norm", "ThickAvg")
 
tissue_abbreviations <- c(
  dlpfc = "dorsolateralprefrontalcortex",
  tc = "temporalcortex",
  fc = "frontalcortex",
  pcc = "posteriorcingulatecortex",
  hcn = "Headofcaudatenucleus"
)

plot_rrho <- function(region_abbrev, mri_feature) {
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
    scale_fill_gradient2(low = "white", high = "purple",
                          limits = c(0, max(rrho_data$rrho_mean)),
                          name = "Max RRHO -log10(p)",
                          na.value = "black") +
    theme_void() +
    theme(legend.position = "none", plot.margin = unit(c(0, 10, 0, 0), units = "points"))

  assign(paste0(region_abbrev, "_", mri_feature, "_plot"), p, envir = parent.frame())
}

for (region in regions) {
  for (mri_feature in mri_features) {
    plot_rrho(region, mri_feature)
    print(paste0(region, " ", mri_feature))
  }
}

blanktitle <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "", fontface = "bold") +
  theme_void()
coltitle_dlpfc <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "DLPFC", fontface = "bold") +
  theme_void() +
  theme(plot.margin = unit(c(0, 10, 0, 0), units = "points"))
coltitle_fc <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "FC", fontface = "bold") +
  theme_void() +
  theme(plot.margin = unit(c(0, 10, 0, 0), units = "points"))
coltitle_tc <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "TC", fontface = "bold") +
  theme_void() +
  theme(plot.margin = unit(c(0, 10, 0, 0), units = "points"))
coltitle_pcc <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "PCC", fontface = "bold") +
  theme_void() +
  theme(plot.margin = unit(c(0, 10, 0, 0), units = "points"))
coltitle_hcn <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "HCN", fontface = "bold") +
  theme_void() +
  theme(plot.margin = unit(c(0, 10, 0, 0), units = "points"))
 
rowtitle_thickness <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Average\nThickness", fontface = "bold") +
  theme_void()
rowtitle_area <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Normalized\nArea", fontface = "bold") +
  theme_void()
rowtitle_volume <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Normalized\nVolume", fontface = "bold") +
  theme_void()
rowtitle_curvature <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Curvature\nIndex", fontface = "bold") +
  theme_void()
 
 
plot <- (blanktitle | coltitle_dlpfc | coltitle_fc | coltitle_tc | coltitle_pcc | coltitle_hcn) /
  (rowtitle_thickness | dlpfc_ThickAvg_plot | fc_ThickAvg_plot | tc_ThickAvg_plot | pcc_ThickAvg_plot | hcn_ThickAvg_plot) /
  (rowtitle_area | dlpfc_SurfArea_norm_plot | fc_SurfArea_norm_plot | tc_SurfArea_norm_plot | pcc_SurfArea_norm_plot | hcn_SurfArea_norm_plot) /
  (rowtitle_volume | dlpfc_GrayVol_norm_plot | fc_GrayVol_norm_plot | tc_GrayVol_norm_plot | pcc_GrayVol_norm_plot | hcn_GrayVol_norm_plot) /
  (rowtitle_curvature | dlpfc_CurvInd_plot | fc_CurvInd_plot | tc_CurvInd_plot | pcc_CurvInd_plot | hcn_CurvInd_plot)

ggsave(file = paste0(scratch_path, "/figures/supplemental/10_brainmaps_",ref_fullname,".svg"),
       plot = plot, width = 10, height = 4, device = svg)
