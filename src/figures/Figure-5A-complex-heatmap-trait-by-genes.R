# remotes::install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(data.table)
# library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(ggvenn)
library(igraph)

project_path <- my_path <-"/workspace/fasi-domice/"
# project_path <- "~/Documents/projects/domice/"
figure_path <- paste0(project_path, "results/figures/")
data_path <- paste0(project_path, "data/")
# source(paste0(figure_path, "helpers.R"))

library(circlize)
col_fun = colorRamp2(c(0, 1), c("lightgray", "red"))
col_fun(seq(-3, 3))



vars <- fread(paste0(data_path, "allchannels/vars.csv"), data.table = FALSE)
ccre <- fread(paste0(data_path, "references/GM_SNPS_Consequence_cCRE.csv"), data.table = FALSE)

trait <- fread(paste0(results_dir, "trait_by_loci_10kb_window.csv" ))

trait <- trait %>%
  mutate(loci_chr_num = case_when(loci_chr == "X"~20, TRUE ~ as.numeric(loci_chr)),
         loci_start_pad = str_pad(loci_start, width = 9, side = "left", pad = "0")) %>%
  unite("loci_order", c(loci_chr_num, loci_start_pad), sep = "", remove=FALSE) %>%
  mutate(loci_order = as.numeric(loci_order)) %>%
  mutate(loci_order_rank = rank(loci_order, ties.method = "min"))

ggplot(trait, aes(loci_order_rank, trait, fill = count)) + geom_area()

# mat <- vars %>%
#   dplyr::select(213:218, 220:225, 227:232, 234:239, 260:283, 287:315) %>%
#   dplyr::mutate(topic2_qtl_genes = 0,
#                 topic5_qtl_genes = 0,
#                 topic7_qtl_genes = 0,
#                 topic9_qtl_genes = 0,
#                 topic10_qtl_genes = 0,
#                 topic12_qtl_genes = 0,
#                 topic16_qtl_genes = 0,
#                 topic17_qtl_genes =0) %>%
#   as.matrix
# mat <- mat[,-1]

mat <- vars %>%
  dplyr::select(217,224,231,238, 260:277, 287,288, 292:315) %>%
  # dplyr::mutate(topic2_qtl_genes = 0,
  #               topic5_qtl_genes = 0,
  #               topic7_qtl_genes = 0,
  #               topic9_qtl_genes = 0,
  #               topic10_qtl_genes = 0,
  #               topic12_qtl_genes = 0,
  #               topic16_qtl_genes = 0,
  #               topic17_qtl_genes =0) %>%
  dplyr::select(contains("eqtl"), 
         starts_with("topic"),
         "ilc3_activated_qtl_genes", "lti_activated_qtl_genes",
         "ilc1_ilc2_prop_qtl_genes", "ilc1_ilc3_prop_qtl_genes",
         "ilc2_lti_prop_qtl_genes", "ilc1_lti_prop_qtl_genes",
         "ilc2_ilc3_prop_qtl_genes",  "ilc3_lti_prop_qtl_genes", 
         ends_with("_steady")) %>%
  as.matrix

colnames(mat) <- str_replace(colnames(mat), "_genes_steady|_genes|_loci_genes_cv", "")
colnames(mat) <- str_replace(colnames(mat), "topic", "Topic")
colnames(mat) <- str_replace(colnames(mat), "ilc", "ILC")
colnames(mat) <- str_replace(colnames(mat), "lti", "LTi")
colnames(mat)[25:36] <- toupper(colnames(mat)[25:36])

rowmat <- vars %>%
  dplyr::select(ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed)


colmat <- data.frame(trait = c(rep("eqtl", 4), rep("topic", 12), rep("proportion", 8), 
                               rep("cytokine", 12)))
ca <- HeatmapAnnotation(df = colmat, which = "column")

ha <- Heatmap(mat)
ra <- HeatmapAnnotation(df = rowmat, which = "row")

rbar <- rowAnnotation(genes = row_anno_barplot(rowSums(mat)))
cbar <- HeatmapAnnotation(traits = anno_barplot(colSums(mat, na.rm=T)),
                          trait = colmat$trait)


# Heatmap(mat, name = "genes",
#         col = col_fun,
#         top_annotation = cbar, width = unit(2, "cm"), 
#         row_annotation = rbar,
#         row_order = NULL,
#         column_order = NULL)


pdf("Figure-5A-polygenic-heatmap-effects.pdf")
Heatmap(mat, 
        col = col_fun, 
        name = "genes", 
        top_annotation = cbar, 
        row_order = NULL, 
        column_order = NULL) + rbar
dev.off()
