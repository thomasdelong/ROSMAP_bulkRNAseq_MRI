library(RRHO2)
library(yaml)
library(tidyverse)
library(patchwork)

env <- read_yaml(here::here("01_environment.yml"))
scratch_path <- env$paths$scratch

sets <- "label, t, dkt, f
DLPFC - Parahippocampal Curvature, dorsolateralprefrontalcortex, parahippocampal, CurvInd
FC - Lateral Occipital Thickness, frontalcortex, lateraloccipital, ThickAvg
TC - Pars Orbitalis Surface Area, temporalcortex, parsorbitalis, SurfArea_norm
HCN - Lateral Occipital Volume, Headofcaudatenucleus, lateraloccipital, GrayVol_norm
"
sets <- read.csv(text = sets, strip.white = TRUE)

# gene ranking
tt_rank <- function(t, dkt, f) {
  file.path(scratch_path, "processed_data/toptables/dkt", t,
            paste0(dkt, ".", f, ".csv")) %>%
    read.csv() %>%
    mutate(rank = -log10(P.Value) * logFC) %>%
    dplyr::select(gene = X, rank) %>%
    arrange(desc(rank))
}

ranked <- pmap(list(sets$t, sets$dkt, sets$f), tt_rank) %>%
  set_names(sets$label)

# get unique pairs and compare with RRHO
pairs <- combn(nrow(sets), 2, simplify = FALSE)

rrho_results <- map(pairs, \(p) {
  RRHO2_initialize(ranked[[p[1]]], ranked[[p[2]]],
                   labels = sets$label[p],
                   log10.ind = TRUE, method = "hyper")
}) %>%
  set_names(map_chr(pairs, \(p) paste(sets$label[p], collapse = " vs ")))

# Mean hypermat per comparison (replaces the print() in the original)
rrho_means <- tibble(
  comparison = names(rrho_results),
  mean_hyper = map_dbl(rrho_results, \(r) mean(r$hypermat, na.rm = TRUE))
)
print(rrho_means)

# have to do the RRHO plots in ggplot so they're easy to combine
jet <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
"yellow", "#FF7F00", "red", "#7F0000")

fill_lims <- range(unlist(map(rrho_results, "hypermat")), na.rm = TRUE)

rrho_gg <- function(obj, lims, tag = NULL) {
  m  <- obj$hypermat
  df <- tibble(y = rep(seq_len(nrow(m)), ncol(m)),
               x = rep(seq_len(ncol(m)), each = nrow(m)),
               z = as.vector(m))
  ggplot(df, aes(x, y, fill = z)) +
    geom_raster() +
    scale_fill_gradientn(colours = jet, limits = lims,
                         name = expression("-"*log[10](p))) +
    coord_fixed(expand = FALSE) +
    labs(title = tag) +
    theme_minimal(base_size = 10) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10,
                                    margin = margin(b = 2)),
          plot.margin = margin(2, 2, 2, 2))
}

# Text-only panel for outer labels
label_panel <- function(txt, angle = 0, width = 22) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = str_wrap(txt, width),
             angle = angle, fontface = "bold", size = 3.5) +
    coord_cartesian(clip = "off") +
    theme_void()
}

# putting plots together
n <- nrow(sets)
grid <- rep(list(plot_spacer()), n^2)

# column and row labels, trims one so it's not too big (ie for four sets it should be three rows/cols)
for (j in 2:n) grid[[j]] <- label_panel(sets$label[j])

for (i in 1:(n - 1)) grid[[i * n + 1]] <- label_panel(sets$label[i], angle = 90)

# heatmaps
for (k in seq_along(pairs)) {
  i <- pairs[[k]][1]; j <- pairs[[k]][2]
  grid[[i * n + j]] <- rrho_gg(rrho_results[[k]], lims = fill_lims,
                               tag = sprintf("%.2f", rrho_means$mean_hyper[k]))
}

rrho_grid <- wrap_plots(grid, ncol = n, byrow = TRUE,
                        widths = c(0.3, rep(1, n - 1)),
                        heights = c(0.3, rep(1, n - 1))) +
  plot_layout(guides = "collect")
ggsave(file = paste0(scratch_path, "/figures/main/10e_rrho_grid.svg"),
       plot = rrho_grid, width = 8, height = 8, device = svg)