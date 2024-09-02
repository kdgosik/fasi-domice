#' @title Figure 5B polygenic distribution
#' @author Kirk Gosik
#' @description 
#'
#'

library(data.table)
library(tidyverse)
library(ggpubr)

project_path <-"/home/rstudio/"
figure_path <- paste0(project_path, "results/figures/")
data_path <- paste0(project_path, "data/")
source(paste0(figure_path, "helpers.R"))


dir("data", "polygenic", full.names = TRUE)
## remake polygenic plots in ED5
df <- fread(paste0(data_path, "eqtl/polygenic/polygenic-scaled-compare.csv.gz"), data.table = FALSE) %>%
  dplyr::rename(ILC3 = NCR1ILC3,
                `LTi-like` = `Lti ILC3`)
dflong <- df %>% dplyr::select(-V1) %>% pivot_longer(cols = ILC1:ILC3)


# make_polygenic_plot <- function(gene_name) {
#   
#   out_plot <- dflong %>% 
#     filter(gene == gene_name) %>%
#     mutate(mouseid = as.numeric(factor(mouse_id))) %>%
#     ggplot(aes(x = value, fill = name)) + 
#       geom_dotplot(stackgroups = TRUE, method = "histodot", binwidth = 0.1)  +
#       # scale_y_continuous(NULL, breaks = NULL) +
#       facet_wrap(~name) +
#       theme(legend.position = "none") +
#       labs(title = paste0(gene_name, ": Polygenic Score Distribution by Mouse"),
#            x = "",
#            y = "density")
#   
#   
#   ggsave(filename = paste0("results/figures/ED5B-", gene_name, "-polygenic-by-celltype.pdf"),
#          plot = out_plot,
#          width = 7,
#          height = 5,
#          dpi = 330)
# 
# }




make_polygenic_plot <- function(gene_name, cell_type) {
  
  
  get_mice <- dflong$mouse_id[dflong$value > 1.5 & dflong$name == cell_type & dflong$gene == gene_name]
  
  
  out_plot <- dflong %>% 
    filter(gene == gene_name) %>%
    mutate(breaks = cut(value, 30)) %>% 
    group_by(name, breaks) %>% 
    mutate(count = n():1,
           value = mean(value)) %>% 
    mutate(mouse2 = ifelse(mouse_id %in% get_mice, 1, 0)) %>%
    ggplot(aes(x = value, y = count, color = mouse2)) + 
      geom_point(size = 2) +
      # scale_color_gradient(low = "white", high = "red") +
      facet_wrap(~name) +
      theme_void() +
      theme(legend.position = "none") +
      labs(title = paste0(gene_name, ": Polygenic Score Distribution by Mouse"),
           subtitle = "High Lineage Specific Gene Highlighted",
           x = "Polygenic Score",
           y = "Count")
  
  out_plot
  
}


ilc1_out_plot <- make_polygenic_plot(gene_name = "Tbx21", "ILC1")
ggsave(filename = paste0("results/figures/Figure-5C-", gene_name, "-polygenic-by-celltype-highlight.pdf"),
       plot = ilc1_out_plot,
       width = 7,
       height = 5,
       dpi = 330)


ilc2_kit_out_plot <- make_polygenic_plot(gene_name = "Kit", "ILC2")
ggsave(filename = paste0("results/figures/Figure-5C-", gene_name, "-polygenic-by-celltype-highlight.pdf"),
       plot = ilc2_kit_out_plot,
       width = 7,
       height = 5,
       dpi = 330)


ilc2_gata3_out_plot <- make_polygenic_plot(gene_name = "Gata3", "ILC2")
ggsave(filename = paste0("results/figures/ED5B-", gene_name, "-polygenic-by-celltype-highlight.pdf"),
       plot = ilc2_gata3_out_plot,
       width = 7,
       height = 5,
       dpi = 330)


ilc3_out_plot <- make_polygenic_plot(gene_name = "Ncr1", "ILC3")
ggsave(filename = paste0("results/figures/ED5B-", gene_name, "-polygenic-by-celltype-highlight.pdf"),
       plot = ilc3_out_plot,
       width = 7,
       height = 5,
       dpi = 330)

lti_out_plot <- make_polygenic_plot(gene_name = "Rorc", "LTi-like")
ggsave(filename = paste0("results/figures/ED5B-", gene_name, "-polygenic-by-celltype-highlight.pdf"),
       plot = lti_out_plot,
       width = 7,
       height = 5,
       dpi = 330)

