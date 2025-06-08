library(data.table)
library(tibble)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

library(Seurat)

source('/workspace/fasi-domice/setup.R')
figure_dir <- "/workspace/fasi-domice/results/"

plot_df <- fread(paste0(data_dir, "manuscript-plot-data.csv.gz"), data.table=FALSE)


# obj <- readRDS(paste0(data_dir, "allchannels/ilc-seurat-object.rds"))
obj <- readRDS(paste0(data_dir, "allchannels/allchannels.rds"))

genes_to_check <- c("Cebpb","Ebp","Muc2", "Stub1", "Ift20")
gene_df <- data.frame(t(obj@assays$RNA$counts[genes_to_check,])) %>%
  # dplyr::rename(Cebpb = `obj.assays.RNA.counts..Cebpb....`) %>%
  tibble::rownames_to_column("index")
# dim(gene_df)
# colnames(gene_df)
# head(gene_df)

plot_subset_df <- plot_df %>%
  dplyr::select(index, called_cell_types_new, X_umap1, X_umap2) %>%
  left_join(gene_df)
# head(plot_subset_df)

write.csv(plot_subset_df, paste0(figure_dir, "additional-gene-umap.csv"))

plot_subset_df_summary <- plot_subset_df %>%
  # filter(louvain_labels %in% (c(1,2,3,4,5,8))) %>%
  group_by(called_cell_types_new) %>%
  summarise(across(everything(), ~ mean(.), .names = "mean_{.col}"),
            across(everything(), ~ mean(. != 0) * 100, .names = "pct_nonzero_{.col}"))



write.csv(plot_subset_df_summary, paste0(figure_dir, "additional-gene-summary.csv"))

# p1 <- ggplot(plot_subset_df, aes(X_umap1, X_umap2, color = Cebpb)) + geom_point(shape=46)
# ggsave('results/umap_cebpb.png', p1)


## genes
lapply(genes_to_check, function(i) {
  
  ## gene plots
  p1 <- plot_subset_df %>%
    # filter(louvain_labels %in% (c(1,2,3,4,5,8))) %>%
    ggplot(aes_string("X_umap1", "X_umap2", color = i)) + 
    geom_point(shape = 46) + 
    theme_void() +
    theme(legend.position = "right") +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    labs(title = i,
         color = "")
  
  ggsave(filename = paste0(figure_dir, "umap-all-cells-gene-", i, ".png"),
         plot = p1,
         dpi = 330,
         width = 7,
         height = 7)
  
  ggsave(filename = paste0(figure_dir, "umap-all-cells-gene-", i, ".pdf"),
         plot = p1,
         dpi = 330,
         width = 7,
         height = 7)
  
})
