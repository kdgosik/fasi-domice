#' Figure 2C Proportion qtl
#' 
#' 

library(data.table)
library(dplyr)
library(ggplot2)
# library(ggpubr)


figure_path <- paste0("results/figures/")
data_path <- "data/"
source(paste0(figure_path, "helpers.R"))

## Read marker map
marker_map <- readRDS(paste0(data_path, "genotype/Regev_map_20171221.rds"))

## ILC Cell Proportions ###############
plot_files <- c("ILC3_stressed_vs_non_qtl.csv.gz", 
                "LTi_stressed_vs_non_qtl.csv.gz")
names(plot_files) <- c("RORgt_high ILC3 /RORgt_low ILC3",  "CCR6_high LTi /CCR6_low LTi")


i <- 1
plot_list <- Map(function(f, n) {
  
  plot_out <- fread(paste0("results/proportion/", f)) %>%
    dplyr::mutate(
      padj = p.adjust(p_values, method = "hochberg"),
      lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
    ) %>%
    gwas_plot(., marker_map = marker_map, plot_title = n, cutoff = 6) + 
    geom_hline(yintercept = 6, color = "red", linetype = "dashed")
  
  ggsave(filename = paste0(figure_path, "Figure-2C-", gsub(".csv.gz",".pdf", f)),
         plot = plot_out,
         height = 7,
         width = 7,
         dpi = 330)
  
  i <- i + 1
  
}, plot_files, names(plot_files))
