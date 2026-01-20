library(Seurat)
library(data.table)
library(tibble)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)


source('/workspace/fasi-domice/setup.R')
figure_dir <- "/workspace/fasi-domice/results/"
cell_type_vec <- c("ILC1", "ILC2", "ILC3", "ILC3(LTi-like)","Enterocyte","NK","NKT","DCs")
genes_to_check <- c("Cebpb","Ebp","Muc2", "Stub1", "Ift20")


plot_df <- fread(paste0(data_dir, "manuscript-plot-data.csv.gz"), data.table=FALSE)


# obj <- readRDS(paste0(data_dir, "allchannels/ilc-seurat-object.rds"))
obj <- readRDS(paste0(data_dir, "allchannels/allchannels.rds"))
# obj_subset <-  subset(obj, subset = called_cell_types_new %in% cell_type_vec)

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
            across(everything(), ~ mean(. != 0) * 100, .names = "pct_nonzero_{.col}"),
            across(
              where(is.numeric), 
              ~ sum(. , na.rm = TRUE) / sum(. != 0, na.rm = TRUE), 
              .names = "mean_nonzero_{.col}"
            ))



write.csv(plot_subset_df_summary, paste0(figure_dir, "additional-gene-summary.csv"))

# p1 <- ggplot(plot_subset_df, aes(X_umap1, X_umap2, color = Cebpb)) + geom_point(shape=46)
# ggsave('results/umap_cebpb.png', p1)



# 1. Set the active identity to your cell type column
Idents(obj) <- "called_cell_types_new"

# 2. Subset using 'idents' argument (cleaner syntax)
# obj_subset <- subset(obj, idents = cell_type_vec)
obj_subset <-  subset(obj, subset = called_cell_types_new %in% cell_type_vec)

# 3. Now DotPlot will automatically use these idents
DotPlot(obj_subset, features = genes_to_check)
ggsave('/workspace/fasi-domice/dotplot-check.png')


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



genes_to_check <- c("Muc2", "Stub1", "Ift20")


# 1. Choose your genes and pull the data
plot_data <- plot_df %>%
  # dplyr::select(index, called_cell_types_new, X_umap1, X_umap2) %>%
  dplyr::select(index, called_cell_types_new) %>%
  dplyr::filter(called_cell_types_new %in% cell_type_vec) %>%
  dplyr::left_join(gene_df) %>%
  tidyr::pivot_longer(cols = all_of(genes_to_check), names_to = "gene", values_to = "expression")

# 2. Summarize: Calculate Percent Expressed and Average Expression
summary_data <- plot_data %>%
  group_by(called_cell_types_new, gene) %>%
  summarize(
    avg_exp = mean(expression),
    pct_exp = sum(expression > 0) / n() * 100,
    .groups = "drop"
  )


dot_plot <- ggplot(summary_data, aes(x = gene, y = called_cell_types_new)) +
  # Use geom_point to create the dots
  geom_point(aes(size = pct_exp, color = avg_exp)) +
  
  # Set the color gradient (similar to Seurat's default)
  scale_color_gradient(low = "lightgrey", high = "blue") +
  
  # Adjust the size scale to look like a DotPlot
  scale_size_continuous(range = c(1, 8)) +
  
  # Clean up the theme
  theme_minimal() +
  labs(
    x = "Features",
    y = "Cell Type",
    color = "Avg. Expression",
    size = "% Expressed"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )


ggsave(filename = '/workspace/fasi-domice/dotplot_test.png',
       dot_plot)
