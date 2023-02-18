#' @title Figure 1B Cell Atlas
#' @author Kirk Gosik
#' @description
#'
#'


library(data.table)
library(tidyverse)

figure_path <- "results/figures/"
data_path <- "data/"
plot_df <- fread(paste0(data_path, "manuscript-plot-data.csv"), data.table = FALSE)

## subsetting to ILC populations
plot_df <- plot_df %>%
  filter(louvain_labels %in% (c(1,2,3,4,5,8))) %>%
  mutate(cell_types = factor(louvain_labels)) %>%
  mutate(cell_types = fct_recode(cell_types,
                                 "LTi-like/CCR6_low" = "1",
                                 "ILC3/RORgt_low" = "2",
                                 "ILC2" = "3",
                                 "ILC3/RORgt_high" = "4",
                                 "LTi-like/CCR6_high" = "5",
                                 "ILC1" = "8"))
  
label_df <- plot_df %>%
  dplyr::group_by(cell_types) %>%
  dplyr:: summarise(X_tsne1 = median(X_tsne1),
            X_tsne2 = median(X_tsne2),
            X_umap1 = median(X_umap1),
            X_umap2 = median(X_umap2))


## tSNE Atlas #####################

p_tsne <- plot_df %>%
  ggplot(aes(X_tsne1, X_tsne2, color = cell_types)) + 
    geom_point(shape = 46) + 
    geom_label(mapping = aes(X_tsne1, X_tsne2, label = cell_types), data = label_df) +
    theme_void() +
    theme(legend.position = "none") +
    labs(color = "Cell Type") + 
    guides(color = guide_legend(override.aes = list(size=5, shape=20)))


ggsave(filename = "./results/figures/Figure-1B-tsne-cell-types.pdf",
       plot = p_tsne,
       dpi = 330,
       width = 7,
       height = 7)


## UMAP Atlas #####################

p_umap <- plot_df %>%
  ggplot(aes(X_umap1, X_umap2, color = cell_types)) + 
  geom_point(shape = 46) + 
  geom_label(mapping = aes(X_umap1, X_umap2, label = cell_types), data = label_df) +
  theme_void() +
  theme(legend.position = "none") +
  labs(color = "Cell Type") + 
  guides(color = guide_legend(override.aes = list(size=5, shape=20)))

ggsave(filename = "./results/figures/Figure-1B-umap-cell-types.pdf",
       plot = p_umap,
       dpi = 330,
       width = 7,
       height = 7)


s