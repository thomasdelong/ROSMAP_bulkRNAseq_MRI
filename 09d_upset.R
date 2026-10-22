library(yaml)
library(tidyverse)
library(ComplexUpset)
env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

sets <- "label, tissue, parcel, feature
DLPFC - Parahippocampal Curvature, dorsolateralprefrontalcortex, parahippocampal, CurvInd
FC - Lateral Occipital Thickness, frontalcortex, lateraloccipital, ThickAvg
HCN - Lateral Occipital Volume, Headofcaudatenucleus, lateraloccipital, GrayVol_norm
TC - Pars Orbitalis Surface Area, temporalcortex, parsorbitalis, SurfArea_norm
"
sets <- read.csv(text = sets, strip.white = TRUE)

sig_sets <- pmap(list(sets$tissue, sets$parcel, sets$feature), \(t, p, f) {
  read.csv(paste0(scratch_path, "/processed_data/toptables/dkt/", t, "/", p, ".", f, ".csv")) %>%
    filter(adj.P.Val < 0.1) %>%
    pull(X)
}) %>%
  set_names(sets$label)

all_genes <- unique(unlist(sig_sets))
membership_df <- data.frame(gene = all_genes, check.names = FALSE)
for (s in names(sig_sets)) membership_df[[s]] <- all_genes %in% sig_sets[[s]]

p <- upset(
  membership_df,
  names(sig_sets),
  name = "Significant gene sets",
  min_size = 10,
  base_annotations = list(
    'Intersection size' = intersection_size(counts = FALSE)
  ),
  set_sizes = FALSE,
  height_ratio = 0.3,
  stripes = 'white',
  themes = theme_void()
)

ggsave(file = paste0(scratch_path, "/figures/main/09d_upset_sig_genes.svg"),
       plot = p, width = 8, height = 4, device = svg)