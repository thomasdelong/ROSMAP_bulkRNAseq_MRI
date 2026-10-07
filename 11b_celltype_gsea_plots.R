library(tidyverse)
library(yaml)
library(patchwork)
library(ggseg)
library(ggsegFreeSurfer)
library(gprofiler2)
env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch
reference_parcel <- 'middletemporal'
reference_feature <- 'GrayVol_norm'

regions <- c("dlpfc", "tc", "fc", "hcn", "pcc")
mri_features <- c("CurvInd", "SurfArea_norm", "GrayVol_norm", "ThickAvg")
 
tissue_abbreviations <- c(
  dlpfc = "dorsolateralprefrontalcortex",
  tc = "temporalcortex",
  fc = "frontalcortex",
  pcc = "posteriorcingulatecortex",
  hcn = "Headofcaudatenucleus"
)

nes_colors <- c("#191970", "white", "white", "white", "#b22222")
nes_values <- scales::rescale(c(-3,-2, 0, 2, 3))

colorbar <- ggplot(data.frame(nes = seq(-3, 3, length.out = 200)),
                   aes(x = 1, y = nes, fill = nes)) +
  geom_raster() +
  scale_fill_gradientn(colors = nes_colors, values = nes_values,
                       limits = c(-3, 3)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(position = "right", expand = c(0, 0),
                     breaks = c(-3, -2, 0, 2, 3)) +
  labs(y = "NES") +
  theme_void() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 9, face = "bold", angle = 0,
                                    vjust = 0.5))

plot_nes <- function(region_abbrev, feature) {
  region_data <- results %>%
    filter(rna_region == region_abbrev, mri_feature == feature) %>%
    mutate(label = paste0("lh_", dkt_parcel))

  p <- ggplot() +
    geom_brain(data = region_data,
               mapping = aes(fill = .data$nes),
               hemi = "left",
               view = c("lateral", "medial"),
               atlas = dkt(),
               position = position_brain(hemi ~ view),
               color = "black") +
    # scale_fill_gradient2(low = "#191970", high = "#b22222",
    #                       limits = c(min(results$nes), max(results$nes)),
    #                       name = "NES",
    #                       na.value = "black") +
    scale_fill_gradientn(colors = nes_colors,
            values = nes_values)+
    theme_void() +
    theme(legend.position = "none", plot.margin = unit(c(0, 10, 0, 0), units = "points"))

  assign(paste0(region_abbrev, "_", feature, "_plot"), p, envir = parent.frame())
}
print(list.files(paste0(scratch_path, "/processed_data/gsea_result/")))
reference_parcel_results <- data.frame(cell_type = NA, rna_region = NA, nes = NA, p = NA)
for (result_file in list.files(paste0(scratch_path, "/processed_data/gsea_result/"))){
    results <- read.csv(paste0(scratch_path, "/processed_data/gsea_result/", result_file))
    
    for (region in regions) {
        for (mri_feature in mri_features) {
            plot_nes(region, mri_feature)
            #print(paste0(region, " ", mri_feature))
        }
    }

    cell <- gsub('_results.csv', '', paste0(result_file))

    celltitle <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste0(cell), fontface = "bold") +
    theme_void()

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
    
    
    plot <- (celltitle) /
    (blanktitle | coltitle_dlpfc | coltitle_fc | coltitle_tc | coltitle_pcc | coltitle_hcn) /
    (rowtitle_thickness | dlpfc_ThickAvg_plot | fc_ThickAvg_plot | tc_ThickAvg_plot | pcc_ThickAvg_plot | hcn_ThickAvg_plot) /
    (rowtitle_area | dlpfc_SurfArea_norm_plot | fc_SurfArea_norm_plot | tc_SurfArea_norm_plot | pcc_SurfArea_norm_plot | hcn_SurfArea_norm_plot) /
    (rowtitle_volume | dlpfc_GrayVol_norm_plot | fc_GrayVol_norm_plot | tc_GrayVol_norm_plot | pcc_GrayVol_norm_plot | hcn_GrayVol_norm_plot) /
    (rowtitle_curvature | dlpfc_CurvInd_plot | fc_CurvInd_plot | tc_CurvInd_plot | pcc_CurvInd_plot | hcn_CurvInd_plot) 

    plot <- plot | (plot_spacer() / colorbar / plot_spacer() +
                        plot_layout(heights = c(2, 3, 2)))
    plot <- plot + plot_layout(widths = c(25, 1))

    print(cell)
    ggsave(file = paste0(scratch_path, "/figures/supplemental/11_brainmaps_",cell,".svg"),
        plot = plot, width = 10, height = 4, device = svg)

  ref_res <- results %>%
    filter(mri_feature == reference_feature) %>%
    filter(dkt_parcel == reference_parcel) %>%
    select(-c('X','mri_feature', 'dkt_parcel')) %>%
    mutate(cell_type = cell)
  reference_parcel_results <- rbind(reference_parcel_results, ref_res)
}



reference_parcel_results <- reference_parcel_results %>%
  drop_na(cell_type) %>%
  mutate(significant = case_when(
      p < 0.0001 ~ "***",
      p < 0.001  ~ "**",
      p < 0.05   ~ "*",
      TRUE          ~ ""),
    region_upper = factor(toupper(rna_region), levels = c("HCN", "PCC", "TC", "FC", "DLPFC")),
    cell_type = factor(cell_type, levels = unique(cell_type))
  )
print(reference_parcel_results)

p <- ggplot(reference_parcel_results, aes(y = region_upper, x = cell_type, fill = nes)) +
  geom_tile() +
  geom_text(aes(label = significant), size = 10, vjust = 0.7) +
  scale_fill_gradient2(low = "blue", high = "orange") +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0.1),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  labs(y = "RNA Sample Region", x = "Cell Type Marker Enrichment")
p

ggsave(file = paste0(scratch_path, "/figures/supplemental/11_heatmap_",reference_parcel, '.', reference_feature,".svg"),
    plot = p, width = 10, height = 6, device = svg)

