#' @title Figure 1D eQTL map
#' @author Kirk Gosik
#' 
#' @description Cell type specific barplot of eQTL types
#'
#'


library(data.table)
# library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(ggplot2)



data_path <- "data/"
results_path <- "results/"
figure_path <- paste0("./results/figures/")

## read in genesets
# dir('~/Google Drive/My Drive/projects/domice/paper of QTL study/Revised materials of ILC-QTL paper cell Science format/data/results/eqtl', "qtl-loci-by-genes-", full.names = T) 
ilc1_eqtl_loci_by_gene <- read_csv(paste0(results_dir, "eqtl/qtl-loci-by-genes-ILC1.csv"))
ilc2_eqtl_loci_by_gene <- read_csv(paste0(results_dir, "eqtl/qtl-loci-by-genes-ILC2.csv"))
ilc3_eqtl_loci_by_gene <- read_csv(paste0(results_dir, "eqtl/qtl-loci-by-genes-ILC3.csv"))
lti_eqtl_loci_by_gene <- read_csv(paste0(results_dir, "eqtl/qtl-loci-by-genes-LTi.csv"))



bar_plot_df <- list("ILC1" = ilc1_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
                    "ILC2" = ilc2_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
                    "ILC3" = ilc3_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
                    "LTi" = lti_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect)) %>% 
  bind_rows(.id = "cell_type") %>%
  mutate(cis_effect = factor(cis_effect, levels = 0:1, labels = c("trans", "cis")))



ggplot(bar_plot_df, aes(cell_type, n, group = cis_effect, fill = cis_effect)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() + 
  theme(legend.position = "top") +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "",
       x = "Cell Type",
       y = "Number of loci",
       fill = "Loci Type")





