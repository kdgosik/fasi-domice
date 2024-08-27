#' @title Figure 4 A.QTL intersection among four phenotypic traits
#' @author Kirk Gosik
#' 
#' @description 
#'
#'
#'



## setup ###############################
library(UpSetR)
library(tidyverse)
library(qtl2)
library(ComplexUpset)


source("/home/rstudio/src/plotting-functions.R")

datadir <- "/home/rstudio/data/"
resultsdir <- "/home/rstudio/results/"

## marker map
marker_map <- readRDS(paste0(datadir, "genotype/Regev_map_20171221.rds"))


## reading data ########################
ccre <- readr::read_csv(paste0(datadir, "references/GM_SNPS_Consequence_cCRE.csv"))
lods <- readr::read_csv(paste0(resultsdir, "qtl-steady-cytokines-lods.csv.gz"))
vars <- fread(paste0(datadir, "vars.csv"), data.table = FALSE)

# plot_qtl(qtl_output = lods, marker_map = marker_map, varcolumn = 'IFNg')
# 
# mydata <- ccre %>%
#   # dplyr::select(starts_with("topic")) 
#   dplyr::select(starts_with("topic"), ends_with("steady"))
#   
# # ComplexUpset::upset(data = mydata[,1:5], colnames(mydata)[1:5])
# ComplexUpset::upset(data = mydata, colnames(mydata))

# topics - qtl - 260:271
# cell proportions - 272:277
# eqtl loci genes - 
# anitbodies - 278:283
# cytokine - steady - 292:303
# cytokine - allergy - 304:315



## Cytokine steady and eQTL
plotdata <- vars[, c(217, 224,231,238, 292:303)]
colnames(plotdata) <- gsub("_qtl_genes_steady", "", colnames(plotdata))
colnames(plotdata)[1:4] <- c("ILC1 eQTLs", "ILC2 eQTLs", "ILC3 eQTLs", "LTi eQTLs")
colnames(plotdata)[5:16] <- toupper(colnames(plotdata)[5:16])

pdf("results/figures/Figure-6_EDF-upset-cytokine-eqtls.pdf")
ComplexUpset::upset(
  plotdata,
  colnames(plotdata),
  min_degree = 2,
  min_size = 4,
  width_ratio=0.1,
  name = "intersections"
) + labs(title = "Intersection of eQTL and Cytokine QTL Loci")
dev.off()



## eQTL and topics
plotdata <- vars[, c(217, 224,231,238, 260:271)]
colnames(plotdata)[1:4] <- c("ILC1 eQTLs", "ILC2 eQTLs", "ILC3 eQTLs", "LTi eQTLs")
colnames(plotdata) <- gsub("_qtl_genes", "", colnames(plotdata))

pdf("results/figures/Figure-6D-upset-eQTL-topics.pdf")
ComplexUpset::upset(
  plotdata,
  colnames(plotdata),
  min_degree = 2,
  min_size = 6,
  width_ratio=0.1,
  name = "intersections"
) + labs(title = "Intersection of Topic and eQTL Loci")
dev.off()




## Cytokine cytokine and topics
plotdata <- vars[, c(217,224,231,238,292:303,260:271)]
colnames(plotdata) <- gsub("_steady", "", gsub("_qtl_genes", "", colnames(plotdata)))
metadata <- data.frame(column = colnames(plotdata),
                       type = sapply(strsplit(colnames(plotdata), " "), "[", 2)) %>%
  mutate(column = gsub("genes", "topics", column),
         type = gsub("genes", "topics", type))


pdf("results/figures/Figure-6E-upset-cytokine-topics.pdf")
ComplexUpset::upset(
  plotdata,
  colnames(plotdata),
  min_degree = 2,
  min_size = 3,
  width_ratio=0.1,
  name = "intersections"
) + labs(title = "Intersection of Topic and Cytokine QTL Loci")
dev.off()





## Cytokine steady and proportion qtls
plotdata <- vars[, c(292:303, 274, 275, 276, 277, 287, 288)]
colnames(plotdata) <- gsub("_steady", "", gsub("_qtl_genes", "", colnames(plotdata)))
colnames(plotdata)[13:18] <- c("ILC1 vs ILC2", "ILC1 vs ILC3", "ILC2 vs ILC3", "ILC3 vs LTi", "ILC2 vs LTi", "ILC1 vs LTi")
colnames(plotdata)[1:12] <- toupper(colnames(plotdata)[1:12])


pdf("results/figures/Figure-6C-upset-cytokine-proportion.pdf")
ComplexUpset::upset(
  plotdata,
  colnames(plotdata),
  min_degree = 2,
  min_size = 2,
  width_ratio=0.1,
  name = "intersections"
) + labs(title = "Intersection of Proportion and Cytokine QTL Loci")
dev.off()





## eQTL and proportion
plotdata <- vars[, c(217, 224, 231, 238, 274, 275, 276, 277, 287, 288)]
colnames(plotdata)[1:4] <- c("ILC1 eQTLs", "ILC2 eQTLs", "ILC3 eQTLs", "LTi eQTLs")
colnames(plotdata)[5:10] <- c("ILC1 vs ILC2", "ILC1 vs ILC3", "ILC2 vs ILC3", "ILC3 vs LTi", "ILC2 vs LTi", "ILC1 vs LTi")

pdf("results/figures/Figure-6_EDC-upset-eQTL-proportions.pdf")
ComplexUpset::upset(
  plotdata,
  colnames(plotdata),
  min_degree = 2,
  min_size = 4,
  width_ratio=0.1,
  name = "intersections"
) + labs(title = "Intersection of eQTL of Proportions Loci")
dev.off()





## topics and proportion
plotdata <- vars[, c(260:271, 274, 275, 276, 277, 287, 288)]
colnames(plotdata) <- gsub("_qtl_genes", "", colnames(plotdata))
colnames(plotdata)[13:18] <- c("ILC1 vs ILC2", "ILC1 vs ILC3", "ILC2 vs ILC3", "ILC3 vs LTi", "ILC2 vs LTi", "ILC1 vs LTi")

pdf("results/figures/Figure-6_EDD-upset-proportion-topics.pdf")
ComplexUpset::upset(
  plotdata,
  colnames(plotdata),
  min_degree = 2,
  min_size = 4,
  width_ratio=0.1,
  name = "intersections"
) + labs(title = "Intersection of Topic and proportion Loci")
dev.off()
