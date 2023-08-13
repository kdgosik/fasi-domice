#' @title Waterfall of QTL Architecture
#' @author Kirk Gosik
#' @description 
#'
#'
#' @references Figure-3b-complex-heatmap-trait-by-genes.R


## setup #####################
# devtools::install_github("jokergoo/ComplexHeatmap")

library(ComplexHeatmap)
library(data.table)
# library(tidyverse)
library(dplyr)
library(tidyr)
# library(GenVisR)
library(igraph)
library(circlize)

figure_path <- "results/figures/"
data_path <- "data/"
source(paste0(figure_path, "helpers.R"))


col_fun = colorRamp2(c(0, 1), c("lightgray", "red"))
col_fun(seq(-3, 3))


## read data #########################
vars <- fread(paste0(data_path, "vars.csv"))
ccre <- fread(paste0(data_path, "references/GM_SNPS_Consequence_cCRE.csv"))


# ## Notes:
#   - Use gene, consequence columns
#   - instead of samples use trait/phenotype
#   - put in long format

plotdata <- ccre %>%
  dplyr::select(1,2,3,19, 29, 37:97) %>%
  pivot_longer(cols = c(-ensembl_gene, -conseq_clean)) %>%
  filter(value == 1) %>% 
  filter(!is.na(ensembl_gene), !is.na(conseq_clean)) %>%
  dplyr::rename(sample = name, gene = ensembl_gene, variant_class = conseq_clean)

# Create a vector to save mutation priority order for plotting
mutation_priority <- as.character(unique(plotdata$variant_class))

# Create an initial plot
pdf("Figure-7A-waterfall-plot.pdf")
waterfall(plotdata,
          main_geneLabSize = 0,
          fileType = "Custom", 
          variant_class_order = mutation_priority)
dev.off()



mat <- plotdata <- ccre %>%
  dplyr::select(1,2,3,19, 29, 37:97) %>% 
  dplyr::select(6:66) %>% 
  as.matrix()
Heatmap(mat, name = "mat", cluster_rows = FALSE, cluster_columns = FALSE) # turn off row clustering
## reformat #########################

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
  select(contains("eqtl"), 
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

## plot ###############################


## TODO: change from example

## example data
mutationData <- read.delim("http://genomedata.org/gen-viz-workshop/GenVisR/BKM120_Mutation_Data.tsv")


# Reformat the mutation data for waterfall()
mutationData <- mutationData[,c("patient", "gene.name", "trv.type", "amino.acid.change")]
colnames(mutationData) <- c("sample", "gene", "variant_class", "amino.acid.change")

# Create a vector to save mutation priority order for plotting
mutation_priority <- as.character(unique(mutationData$variant_class))

# Create an initial plot
waterfall(mutationData, fileType = "Custom", variant_class_order=mutation_priority)
