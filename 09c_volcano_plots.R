library(yaml)
library(tidyverse)
library(gprofiler2)
library(ggrepel)
library(ggrastr)

env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

out_dir <- file.path(scratch_path, "figures", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Highest-pi1 parcel for each RNA region
pi1maps_dir <- paste0(scratch_path, "/processed_data/pi1maps/")
all_pi1 <- list.files(pi1maps_dir, full.names = TRUE) %>%
  map_dfr(function(f) {
    filename <- str_remove(basename(f), "\\.csv$")
    parts <- str_split_fixed(filename, "_", 2)
    read.csv(f) %>% mutate(RNA_region = parts[1], mri_feature = parts[2])
  })

tissue_abbreviations <- data.frame(
  RNA_region = c("dlpfc", "tc", "fc", "pcc", "hcn"),
  tissue_name = c("dorsolateralprefrontalcortex", "temporalcortex", "frontalcortex",
                  "posteriorcingulatecortex", "Headofcaudatenucleus"))
feature_abbreviations <- data.frame(
  mri_feature = c("ThickAvg", "SurfArea_norm", "GrayVol_norm", "CurvInd"),
  mri_feature_name = c("Thickness", "Surface Area", "Volume", "Curvature Index"))

highest_by_RNA_region <- all_pi1 %>%
  group_by(RNA_region) %>%
  filter(pi1 == max(pi1, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(tissue_abbreviations, by = "RNA_region") %>%
  left_join(feature_abbreviations, by = "mri_feature") %>%
  mutate(RNA_region_upper = toupper(RNA_region))
print(highest_by_RNA_region)

example_for_names <- read.csv(paste0(scratch_path, "/processed_data/toptables/dkt/dorsolateralprefrontalcortex/caudalanteriorcingulate.CurvInd.csv"))
genenames <- gconvert(example_for_names$X, filter_na = FALSE) %>% dplyr::select(target, name)


region_colors <- c(dlpfc = "#e90000", fc = "#0041b2", tc = "#e19a00", pcc = "#04801b", hcn = "#D55E00")

top_tts <- highest_by_RNA_region %>%
  pmap(function(RNA_region, tissue_name, dkt_region, mri_feature, ...) {
    read.csv(file.path(scratch_path, "processed_data/toptables/dkt", tissue_name,
                       paste0(dkt_region, ".", mri_feature, ".csv"))) %>%
      dplyr::select(X, P.Value, adj.P.Val, logFC, B) %>%
      merge(genenames, by.x = "X", by.y = "target") %>%
      mutate(BH_1 = adj.P.Val < 0.1)
  }) %>%
  set_names(highest_by_RNA_region$RNA_region)

# highest gene for each
top_gene_per_region <- highest_by_RNA_region %>%
  mutate(
    gene_id = map_chr(RNA_region, ~ top_tts[[.x]] %>% arrange(P.Value) %>% slice(1) %>% pull(X)),
    gene    = map_chr(RNA_region, ~ top_tts[[.x]] %>% arrange(P.Value) %>% slice(1) %>% pull(name)),
    color   = unname(region_colors[RNA_region]),
    region_feature = paste0(dkt_region, ".", mri_feature)
  )
print(top_gene_per_region %>%
        dplyr::select(RNA_region, dkt_region, mri_feature, gene, gene_id, color))


genes_to_label <- top_gene_per_region %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(X = gene_id)

volcano_plots <- list()
for (tissue in top_gene_per_region$RNA_region) {
  parcel <- top_gene_per_region %>% filter(RNA_region == tissue)
  tt <- top_tts[[tissue]]
  plot_color <- parcel$color

  labeled <- tt %>% inner_join(genes_to_label, by = "X")

  p <- ggplot(tt, aes(x = logFC, y = -log10(P.Value))) +
    geom_point_rast(aes(color = BH_1), raster.dpi = 600) +
    labs(color = "q < 0.1") +
    xlab("β") +
    ylim(0,8) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = plot_color)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) +
    geom_text_repel(data = labeled, aes(label = name),
                    show.legend = FALSE, box.padding = 0.75, max.overlaps = Inf) +
    ggtitle(paste(parcel$RNA_region_upper, parcel$dkt_region, parcel$mri_feature_name))

  volcano_plots[[tissue]] <- p
  ggsave(file.path(out_dir, paste0("09c_volcano_", tissue, ".pdf")), p, width = 5, height = 4)
}

# data for scatter plots
combined_data <- read.csv(file.path(scratch_path, "processed_data/07_combined_wholebrain_metadata.csv"), check.names = FALSE)
rnaseq_subset_filtered <- read.csv(file.path(scratch_path, "processed_data/07_rnaseq_subset.csv"), row.names = 1, check.names = FALSE)
dkt_latest <- read.csv(file.path(scratch_path, "processed_data/DKT_stats_combat_averaged.csv")) %>%
  rename(projid = Subject) %>%
  group_by(projid) %>%
  filter(Session == max(Session)) %>%
  ungroup()

covariates <- c("age_diff", "age_death", "master_brain_side", "msex")

get_scatter_data <- function(rna_tissue, region_feature, gene_id) {
  data <- combined_data[, -1] %>%
    filter(gsub(" ", "", tissue) == rna_tissue) %>%
    left_join(dkt_latest[, c("projid", region_feature)], by = "projid") %>%
    column_to_rownames("specimenID")

  df <- data.frame(
    x = data[[region_feature]],
    y = as.numeric(unlist(rnaseq_subset_filtered[gene_id, rownames(data)])),
    data[, covariates]
  )
  df[complete.cases(df), ]
}

# scatter for each gene
scatter_plots <- list()
for (i in seq_len(nrow(top_gene_per_region))) {
  g <- top_gene_per_region[i, ]
  stats <- top_tts[[g$RNA_region]] %>% filter(X == g$gene_id)

  df <- get_scatter_data(rna_tissue = g$tissue_name,
                         region_feature = g$region_feature,
                         gene_id = g$gene_id)

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(color = g$color, alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, color = g$color, fill = g$color, alpha = 0.15) +
    labs(title = paste0(g$gene, " (", g$RNA_region_upper, ")"),
         x = paste(g$dkt_region, g$mri_feature_name),
         y = paste(g$gene, "expression")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())
  scatter_plots[[g$gene]] <- p
  ggsave(file.path(out_dir, paste0("09c_scatter_", g$RNA_region, "_", g$gene, ".pdf")), p, width = 5, height = 4)
}