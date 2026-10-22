library(tidyverse)
library(limma)
library(yaml)
library(patchwork)
library(ggseg)
library(ggsegFreeSurfer)
env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

refs <- "label, tissue, parcel, feature
DLPFC - Parahippocampal Curvature, dorsolateralprefrontalcortex, parahippocampal, CurvInd
FC - Lateral Occipital Thickness, frontalcortex, lateraloccipital, ThickAvg
TC - Pars Orbitalis Surface Area, temporalcortex, parsorbitalis, SurfArea_norm
HCN - Lateral Occipital Volume, Headofcaudatenucleus, lateraloccipital, GrayVol_norm
"
refs <- read.csv(text = refs, strip.white = TRUE) %>%
  mutate(fullname = paste(tissue, parcel, feature, sep = "_"))

regions <- c("dlpfc", "tc", "fc", "hcn", "pcc")
mri_features <- c("CurvInd", "SurfArea_norm", "GrayVol_norm", "ThickAvg")

tissue_abbreviations <- c(
  dlpfc = "dorsolateralprefrontalcortex",
  tc = "temporalcortex",
  fc = "frontalcortex",
  pcc = "posteriorcingulatecortex",
  hcn = "Headofcaudatenucleus"
)
tissue_labels <- setNames(toupper(names(tissue_abbreviations)), tissue_abbreviations)
tissue_order <- c("DLPFC", "FC", "TC", "PCC", "HCN")

row_labels <- c(
  ThickAvg      = "Average\nThickness",
  SurfArea_norm = "Normalized\nArea",
  GrayVol_norm  = "Normalized\nVolume",
  CurvInd       = "Curvature\nIndex"
)
feature_labels <- str_replace(row_labels, "\n", " ") %>% setNames(names(row_labels))

all_rrho <- pmap(list(refs$label, refs$tissue, refs$parcel, refs$feature, refs$fullname),
                 \(lab, t, p, f, full) {
  read.csv(paste0(scratch_path, "/processed_data/rrho/", full, ".csv")) %>%
    #filter out the reference parcel since you get a super high value comparing it to itself
    filter(!(tissue == t & parcel == p & feature == f)) %>%
    mutate(ref = lab)
}) %>%
  bind_rows() %>%
  mutate(ref = factor(ref, levels = refs$label))

#function to compare the observed mean RRHO value across maps for a given tissue
#or feature and compare it to a shuffled null value. 
perm_summary <- function(df, group_col, n_perm = 1000) {
#observed mean
  obs <- df %>%
    group_by(ref, group = .data[[group_col]]) %>%
    summarise(obs = mean(rrho_mean, na.rm = TRUE), .groups = "drop")
#null 
  null <- map(seq_len(n_perm), \(i) {
    df %>%
      group_by(ref) %>%
      mutate(group = sample(.data[[group_col]])) %>%
      group_by(ref, group) %>%
      summarise(null = mean(rrho_mean, na.rm = TRUE), .groups = "drop")
  }) %>%
    bind_rows()
#make p by comparing observed to null distribution
  null %>%
    left_join(obs, by = c("ref", "group")) %>%
    group_by(ref, group) %>%
    summarise(obs = first(obs),
              ratio = first(obs) / mean(null),
              z = (first(obs) - mean(null)) / sd(null),
              p = (sum(null >= first(obs)) + 1) / (n() + 1),
              .groups = "drop") %>%
    mutate(p_adj = p.adjust(p, method = "BH"))
}

set.seed(1)
tissue_summary <- perm_summary(all_rrho, "tissue")
feature_summary <- perm_summary(all_rrho, "feature")

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
    scale_fill_gradient2(low = "white", high = "#b22222",
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

run_ref <- function(ref_label, ref_fullname) {
  print(ref_fullname)

  rrho_data <- filter(all_rrho, ref == ref_label)

  for (region in regions) {
    for (mri_feature in mri_features) {
      plot_rrho(rrho_data, region, mri_feature)
      print(paste0(region, " ", mri_feature))
    }
  }

  col_means <- tissue_summary %>% filter(ref == ref_label) %>% select(group, obs) %>% deframe()
  row_means <- feature_summary %>% filter(ref == ref_label) %>% select(group, obs) %>% deframe()

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
  run_ref(refs$label[r], refs$fullname[r])
}

z_lim <- max(abs(c(tissue_summary$z, feature_summary$z)), na.rm = TRUE)
obs_max <- max(c(tissue_summary$obs, feature_summary$obs))
obs_min <- min(c(tissue_summary$obs, feature_summary$obs))
ratio_max <- max(c(tissue_summary$ratio, feature_summary$ratio))
ratio_min <- min(c(tissue_summary$ratio, feature_summary$ratio))

# plot_summary <- function(s, xlab) {
#   ggplot(s, aes(group, ref, fill = z)) +
#     geom_tile(color = "white") +
#     geom_text(aes(label = paste0(sprintf("%.1f", obs), if_else(p_adj < 0.05, "*", ""))),
#               size = 3) +
#     scale_fill_gradient2(low = "white", mid = "white", high = "#b22222",
#                          limits = c(-z_lim, z_lim), name = "z vs\nshuffled") +
#     scale_y_discrete(labels = \(x) str_wrap(x, 20), limits = rev) +
#     scale_x_discrete(labels = \(x) str_wrap(x, 10)) +
#     labs(x = xlab, y = NULL) +
#     theme_minimal(base_size = 10) +
#     theme(panel.grid = element_blank())
# }

plot_summary <- function(s, xlab) {
  ggplot(s, aes(group, ref, fill = ratio)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(if_else(p_adj < 0.05, "*", ""))),
              size = 10) +
    scale_fill_gradient2(low = "white", mid = "white", high = "#b22222",
                        midpoint = 1, 
                        limits = c(ratio_min, ratio_max), name = "Observered/Shuffled\nRRHO -log10(p) ratio") +
    scale_y_discrete(labels = \(x) str_wrap(x, 20), limits = rev) +
    scale_x_discrete(labels = \(x) str_wrap(x, 10)) +
    labs(x = xlab, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank())
}

tissue_plot_data <- tissue_summary %>%
  mutate(group = factor(tissue_labels[group], levels = tissue_order))
feature_plot_data <- feature_summary %>%
  mutate(group = factor(feature_labels[group], levels = feature_labels[c("ThickAvg", "SurfArea_norm", "GrayVol_norm", "CurvInd")]))

summary_plot <- (plot_summary(tissue_plot_data, "Tissue") |
  (plot_summary(feature_plot_data, "Feature") + theme(axis.text.y = element_blank()))) +
  plot_layout(guides = "collect", widths = c(5, 4))

ggsave(file = paste0(scratch_path, "/figures/supplemental/10d_rrho_summary.svg"),
       plot = summary_plot, width = 9, height = 4, device = svg)