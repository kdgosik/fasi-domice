library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)


project_path <- "/workspace/fasi-domice/"


plotdf <- fread(paste0(data_dir, "manuscript-plot-data.csv.gz"))

p_mouse <- plotdf %>%
  filter(louvain_labels %in% c(1:5,8),
         BestCall != "AMB") %>%
  ggplot(aes(X_umap1, X_umap2, color = BestCall)) + 
    geom_point(shape = 46, alpha = 0.25) +
    theme_void() + 
    labs(color = "sample") +
    theme(legend.title = element_text(size = 4),
          # legend.position = c(1, .5),
          legend.spacing = unit(0,'cm'),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          panel.spacing = unit(0, "cm")) +
    guides(color = guide_legend(ncol = 8,
                                keywidth=0,
                                keyheight=0,
                                label.theme = element_text(size = 2))) +
  scale_color_manual(values = c(rep(godsnot_64,3)))
ggsave(filename = "umap_mouseid_w_legend.pdf",
       plot = p_mouse,
       dpi = 330)


p_batch <- plotdf %>%
  filter(louvain_labels %in% c(1:5,8),
         BestCall != "AMB") %>%
  mutate(batch = str_remove(BestCall, "_[0-9]$")) %>%
  ggplot(aes(X_umap1, X_umap2, color = batch)) + 
  geom_point(shape = 46) +
  theme_void() + 
  labs(color = "batch") +
  theme(legend.title = element_text(size = 6)) +
  guides(color = guide_legend(ncol = 3,
                              keywidth=0.1,
                              keyheight=0.1,
                              label.theme = element_text(size = 6))) +
  scale_color_manual(values = godsnot_64)
ggsave(filename = "umap_batchid_w_legend.pdf",
       plot = p_batch,
       dpi = 330)