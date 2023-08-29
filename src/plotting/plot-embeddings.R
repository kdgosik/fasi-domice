library(data.table)
library(ggplot2)
library(dplyr)


project_path <- "/workspace/fasi-domice/"


plotdf <- fread(paste0(data_dir, "manuscript-plot-data.csv.gz"))

p_mouse <- plotdf %>%
  filter(louvain_labels %in% c(1:5,8),
         BestCall != "AMB") %>%
  ggplot(aes(X_umap1, X_umap2, color = BestCall)) + 
    geom_point(shape = 46) +
    theme_void() + 
    theme(legend.position = "none")
ggsave(filename = "umap_mouseid.pdf",
       plot = p_mouse,
       dpi = 330)


p_batch <- plotdf %>%
  filter(louvain_labels %in% c(1:5,8),
         BestCall != "AMB") %>%
  mutate(batch = str_remove(BestCall, "_[0-9]$")) %>%
  ggplot(aes(X_umap1, X_umap2, color = batch)) + 
  geom_point(shape = 46) +
  theme_void() + 
  theme(legend.position = "none") +
  theme(legend.position = "none") + scale_color_manual(values = godsnot_64)
ggsave(filename = "umap_batchid.pdf",
       plot = p_batch,
       dpi = 330)