library(tidyverse)
library(limma)
library(yaml)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
# if (!require("ggseg", quietly = TRUE))
#   install.packages("ggseg")
library(ggseg)
# if (!require("ggsegFreeSurfer", quietly = TRUE))
#   install.packages('ggsegFreeSurfer', repos = c('https://ggsegverse.r-universe.dev', 'https://cloud.r-project.org'))
library(ggsegFreeSurfer)
env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch
#print(dkt())

tt_dir <- paste0(scratch_path, "/processed_data/toptables/dkt/")
pi1maps_dir <- paste0(scratch_path, "/processed_data/pi1maps/")
dir.create(pi1maps_dir, recursive = TRUE, showWarnings = FALSE)
 
regions <- c("dlpfc", "tc", "fc", "hcn", "pcc")
mri_features <- c("CurvInd", "SurfArea_norm", "GrayVol_norm", "ThickAvg")
 
tissue_abbreviations <- c(
  dlpfc = "dorsolateralprefrontalcortex",
  tc = "temporalcortex",
  fc = "frontalcortex",
  pcc = "posteriorcingulatecortex",
  hcn = "Headofcaudatenucleus"
)

plot_pi1 <- function(region_abbrev, mri_feature) {
  dir_path <- paste0(tt_dir, tissue_abbreviations[[region_abbrev]], "/")
  pattern <- paste0("\\.", mri_feature, "\\.csv$")
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
 
  feature_truenull <- map_dfr(files, function(f) {
    dkt_region <- str_remove(basename(f), pattern)
    tt <- read.csv(f)
    p <- tt$P.Value[!is.na(tt$P.Value)]
    pi0 <- propTrueNull(p)
    tibble(dkt_region = dkt_region, pi0 = pi0, pi1 = 1 - pi0)
  })
 
  feature_truenull$label <- paste0("lh_", feature_truenull$dkt_region)
 
  assign("highest", filter(feature_truenull, pi1 == max(pi1, na.rm = TRUE)), envir = parent.frame())

  #original ggsegDKT code from before it got merged into ggsegFreeSurfer
  # plot_data <- brain_join(feature_truenull, dkt, "label") %>%
  #   filter(hemi == "left") %>%
  #   drop_na(region)
  # #print(head(plot_data))
  # print(dim(plot_data))
  # p <- ggplot(plot_data) +
  #   geom_brain(aes(fill = pi1), atlas = dkt, hemi = "left", view = "lateral") +
  #   scale_fill_gradient2(low = "white", high = "purple",
  #   limits = c(0, 0.415), name = "\u03c01") +
  #   theme_void() +
  #   theme(legend.position = "none", plot.margin = unit(c(0, 10, 0, 0), units = "points"))
  p <- ggplot() +
  geom_brain(data = feature_truenull,
             mapping = aes(fill = .data$pi1),
             hemi = "left",
             view = c("lateral", "medial"),
             atlas = dkt(),
             position = position_brain(hemi ~ view),
             color = "black", linewidth = 0.01) +
  scale_fill_gradient2(low = "white", high = "purple",
                        limits = c(0, 0.425), name = "\u03c01", #catch 22 with the scale here. need to rewrite so the pi1 calculations happen in a different loop. or just run it once then set it manually if anything is out of scale
                        na.value = "grey") +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(0, 10, 0, 0), units = "points"))
  assign(paste0(region_abbrev, "_", mri_feature, "_plot"), p, envir = parent.frame())
 
  write.csv(feature_truenull, paste0(pi1maps_dir, region_abbrev, "_", mri_feature, ".csv"), row.names = FALSE)
}

top_pi1 <- 0

for (region in regions) {
  for (mri_feature in mri_features) {
    plot_pi1(region, mri_feature)
    print(paste0(region, " ", mri_feature))
    print(highest)
    if (highest$pi1 > top_pi1){
      top_pi1 <- highest$pi1
    }
  }
}
 
print(paste('overall highest =', top_pi1))

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

ggsave(file = paste0(scratch_path, "/figures/main/09_brainmaps.svg"),
       plot = plot, width = 10, height = 4, device = svg)

n_perm = 1000
set.seed(1)

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

pi1_dataframes <- list()
i = 0
for (file in list.files(pi1maps_dir)){
    name <- str_replace(file, '.csv', '')
    data <- read.csv(paste0(pi1maps_dir, file)) %>%  
        select(c('label', 'pi1'))

    assign(name, data) 
    i = i + 1
    pi1_dataframes[[i]] <- paste0(name)
}

results <- data.frame(matrix(NA, nrow = length(pi1_dataframes), 
                ncol = length(pi1_dataframes),
                dimnames = list(pi1_dataframes, pi1_dataframes))
)
results_r <- data.frame(matrix(NA, nrow = length(pi1_dataframes), 
                ncol = length(pi1_dataframes),
                dimnames = list(pi1_dataframes, pi1_dataframes))
)

for (df in pi1_dataframes){
    for (df2 in pi1_dataframes){
        results[df, df2] <- random_shuffle(get(df)$pi1, get(df2)$pi1, 10000)$p_pearson
        results_r[df, df2] <- cor(get(df)$pi1, get(df2)$pi1, method = 'pearson')
    }
}

results_p_thickness <- results %>%
    select(contains('Thick')) %>%
    t() %>%
    as.data.frame() %>%
    select(contains('Thick'))


col_fun = colorRamp2(c(0, 0.05, 1), c("red", "white", "white"))
col_fun(seq(0, 1))
heatmap_thickness <- Heatmap(as.matrix(results_p_thickness), col = col_fun,
    heatmap_legend_param = list(at = c(0, 0.01, 0.05),
                                title = 'random \nshuffle \np-value'),
                                cluster_rows = FALSE, cluster_columns = FALSE)

results_p_volume <- results %>%
    select(contains('Vol')) %>%
    t() %>%
    as.data.frame() %>%
    select(contains('Vol'))
heatmap_volume <- Heatmap(as.matrix(results_p_volume), col = col_fun,
    heatmap_legend_param = list(at = c(0, 0.01, 0.05),
                                title = 'random \nshuffle \np-value'),
                                cluster_rows = FALSE, cluster_columns = FALSE)

results_p_surf <- results %>%
    select(contains('Surf')) %>%
    t() %>%
    as.data.frame() %>%
    select(contains('Surf'))
heatmap_surfacearea <- Heatmap(as.matrix(results_p_surf), col = col_fun,
    heatmap_legend_param = list(at = c(0, 0.01, 0.05),
                                title = 'random \nshuffle \np-value'),
                                cluster_rows = FALSE, cluster_columns = FALSE)

results_p_curv <- results %>%
    select(contains('Curv')) %>%
    t() %>%
    as.data.frame() %>%
    select(contains('Curv'))
heatmap_curvature <- Heatmap(as.matrix(results_p_curv), col = col_fun,
    heatmap_legend_param = list(at = c(0, 0.01, 0.05),
                                title = 'random \nshuffle \np-value'),
                                cluster_rows = FALSE, cluster_columns = FALSE)

save_heatmap <- function(heatmap, filepath, width = 7, height = 7){
  svg(filepath, width = width, height = height)
  draw(heatmap)
  dev.off()
}
save_heatmap(heatmap_thickness, paste0(scratch_path, "/figures/supplemental/09_thickness_shuffle_heatmap.svg"))
save_heatmap(heatmap_volume, paste0(scratch_path, "/figures/supplemental/09_volume_shuffle_heatmap.svg"))
save_heatmap(heatmap_surfacearea, paste0(scratch_path, "/figures/supplemental/09_surfacearea_shuffle_heatmap.svg"))
save_heatmap(heatmap_curvature, paste0(scratch_path, "/figures/supplemental/09_curvature_shuffle_heatmap.svg"))




# tccurvind <- read.csv(paste0(scratch_path, "/processed_data/toptables/dkt/temporalcortex/middletemporal.CurvInd.csv"))
# svg(filename=paste0(scratch_path, "/figures/supplemental/09_TCmiddletemporalcurvind.svg"), width=10, height=3)
# hist(tccurvind$P.Value, main = paste("TC Genes vs. Middle Temporal Curvature Index p-values"), breaks = 100, border = 'transparent', col = '#9c61ba', ylim = c(0,1500), xlab = 'p-value (100 bins)')
# dev.off()

# hcncurvind <- read.csv(paste0(scratch_path, "/processed_data/toptables/dkt/Headofcaudatenucleus/middletemporal.CurvInd.csv"))
# svg(filename=paste0(scratch_path, "/figures/supplemental/09_HCNmiddletemporalcurvind.svg"), width=10, height=3)
# hist(hcncurvind$P.Value, main = paste("HCN Genes vs. Middle Temporal Curvature Index p-values"), breaks = 100, border = 'transparent', col = '#9c61ba', ylim = c(0,1500), xlab = 'p-value (100 bins)')
# dev.off()

pcccurvind <- read.csv(paste0(scratch_path, "/processed_data/toptables/dkt/posteriorcingulatecortex/entorhinal.CurvInd.csv"))
svg(filename=paste0(scratch_path, "/figures/main/09_PCCentorhinalcurvind.svg"), width=10, height=3)
hist(pcccurvind$P.Value, main = paste("PCC Genes vs. Entorhinal Curvature Index p-values"), breaks = 100, border = 'transparent', col = '#9c61ba', ylim = c(0,1500), xlab = 'p-value (100 bins)')
dev.off()

hcncurvind <- read.csv(paste0(scratch_path, "/processed_data/toptables/dkt/Headofcaudatenucleus/entorhinal.CurvInd.csv"))
svg(filename=paste0(scratch_path, "/figures/main/09_HCNentorhinalcurvind.svg"), width=10, height=3)
hist(hcncurvind$P.Value, main = paste("HCN Genes vs. Entorhinal Curvature Index p-values"), breaks = 100, border = 'transparent', col = '#9c61ba', ylim = c(0,1500), xlab = 'p-value (100 bins)')
dev.off()