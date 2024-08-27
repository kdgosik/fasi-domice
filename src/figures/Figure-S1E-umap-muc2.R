#' @title Figure S1E Gene UMAPs
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

## UMAP ########################################

## genes
lapply(c("Muc2"), function(i) {
  
  ## topic plots
  p1 <- plot_df %>%
    # filter(louvain_labels %in% (c(1,2,3,4,5,8))) %>%
    ggplot(aes_string("X_umap1", "X_umap2", color = i)) + 
    geom_point(shape = 46) + 
    theme_void() +
    theme(legend.position = "right") +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    labs(title = i,
         color = "")
  
  
  ggsave(filename = paste0(figure_dir, "gene-umaps/Figure-7D-umap-all-cells-gene-", i, ".pdf"),
         plot = p1,
         dpi = 330,
         width = 7,
         height = 7)
  
})


