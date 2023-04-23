#' Figure 2B Proportion qtl
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
plot_files <- c("gwas-ilc1-ilc2-results.csv.gz", "gwas-ilc3-lti-results.csv.gz")
names(plot_files) <- c("ILC1 vs ILC2", "ILC3 vs LTi")

i <- 1
plot_list <- Map(function(f, n) {
  
  plot_out <- fread(paste0("results/proportions/", f)) %>%
    dplyr::mutate(
      padj = p.adjust(p_values, method = "hochberg"),
      lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
    ) %>%
    gwas_plot(., marker_map = marker_map, plot_title = n, cutoff = 6) + 
      geom_hline(yintercept = 6, 
                 color = "red", 
                 linetype = "dashed")
  
  ggsave(filename = paste0(figure_path, "Figure-2B-", gsub(".csv.gz",".pdf", f)),
         plot = plot_out,
         height = 7,
         width = 7,
         dpi = 330)
  
  i <- i + 1
  
}, plot_files, names(plot_files))



# FAM21 - Chromosome 6: 116,208,038-116,262,686 forward strand.
# IL-20RA - Chromosome 10: 19,712,570-19,760,053 forward strand.
f <- "gwas-ilc1-ilc2-results.csv.gz"
gwas_data <- fread(paste0(results_dir, "proportions/", f)) %>%
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  )

colordf <- gwas_data %>%
  filter(chr == 6, between(pos, 116208038, 116262686)) %>% # FAM21
  filter(chr == 10, between(pos, 19712570, 19760053)) %>% # IL20ra

  
gwas_plot(gwas_data = gwas_data, plot_title = n, cutoff = 6) + 
  geom_hline(yintercept = 6, 
             color = "red", 
             linetype = "dashed")

