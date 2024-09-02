#' @title Figure 3B Topic Scores Embeddings
#' @author Kirk Gosik
#' @description
#'
#'


library(data.table)
library(ggplot2)
library(dplyr)

source("setup.R")
figure_dir <- paste0(results_dir, "figures/")
plot_df <- fread(paste0(data_dir, "manuscript-plot-data.csv.gz"), data.table = FALSE)


## tSNE ########################################

lapply(c(0,3,8,11,12,16), function(i) {
  
  ## topic plots
  p1 <- plot_df %>%
    filter(louvain_labels %in% (c(1,2,3,4,5,8))) %>%
    ggplot(aes_string("X_tsne1", "X_tsne2", color = paste0("topic", i))) + 
    geom_point(shape = 46) + 
    theme_void() +
    theme(legend.position = "right") +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    labs(title = paste0("Topic ", i),
         color = "")
  
  
  ggsave(filename = paste0(figure_path, "Figure-3B-tsne-topic", i, ".pdf"),
         plot = p1,
         dpi = 330,
         width = 7,
         height = 7)
  
})


## UMAP ########################################

## QTL Loci Genes
# topic 0: Bmp15
# topic 3: Sstr4, Thbd
# topic 8: 2610035D17Rik, Gm5845, 4933434M16Rik, Sox9(?)
# topic 11: Gm43781, Gm25513

lapply(c(0:19), function(i) {
  
  ## topic plots
  p1 <- plot_df %>%
    filter(louvain_labels %in% (c(1,2,3,4,5,8))) %>%
    ggplot(aes_string("X_umap1", "X_umap2", color = paste0("topic", i))) + 
    geom_point(shape = 46) + 
    theme_void() +
    theme(legend.position = "right") +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    labs(title = paste0("Topic ", i),
         color = "")
  
  
  ggsave(filename = paste0(figure_path, "Figure-3B-umap-topic", i, ".pdf"),
         plot = p1,
         dpi = 330,
         width = 7,
         height = 7)
  
})