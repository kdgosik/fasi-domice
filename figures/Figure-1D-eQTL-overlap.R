#' @title Figure 1C eQTL Overlap
#' @author Kirk Gosik
#' @description
#'
#'


library(data.table)
library(tidyverse)
library(ggvenn)

figure_path <- "results/figures/"
data_path <- "data/"
source(paste0(figure_path, "helpers.R"))

## Genes in eQTL loci ######
vars <- fread(paste0(data_path, "vars.csv"))

p1 <- ggvenn(list("ILC1"=na.omit(vars$index[vars[["ilc1_eqtl_loci_genes_cv"]]==1]),
            "ILC2"=na.omit(vars$index[vars[["ilc2_eqtl_loci_genes_cv"]]==1]),
            "ILC3"=na.omit(vars$index[vars[["ilc3_eqtl_loci_genes_cv"]]==1]),
            "LTi"=na.omit(vars$index[vars[["lti_eqtl_loci_genes_cv"]]==1])),
       set_name_size = 10) +
  labs(title = "eQTL: 3224") +
  theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5, vjust = -5))



ggsave(filename = "./results/figures/Figure1C-eQTL-overlap.pdf",
       plot = p1,
       dpi = 330,
       width = 7,
       height = 7)



## eGenes ######

p2 <- ggvenn(list("ILC1"=na.omit(vars$index[vars[["ilc1_egenes_cv"]]==1]),
                  "ILC2"=na.omit(vars$index[vars[["ilc2_egenes_cv"]]==1]),
                  "ILC3"=na.omit(vars$index[vars[["ilc3_egenes_cv"]]==1]),
                  "LTi"=na.omit(vars$index[vars[["lti_egenes_cv"]]==1])),
             set_name_size = 10) +
  labs(title = "eGene: 833") +
  theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5, vjust = -5))



ggsave(filename = "./results/figures/Figure1C-eGenes-overlap.pdf",
       plot = p2,
       dpi = 330,
       width = 7,
       height = 7)
