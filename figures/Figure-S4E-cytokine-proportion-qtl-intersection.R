#' @title Figure S4E.QTL intersection cytokine and proportion
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


## Cytokine steady and proportion qtls
plotdata <- vars[, c(292:303, 274, 275, 276, 277, 287, 288)]
colnames(plotdata) <- gsub("_steady", "", gsub("_qtl_genes", "", colnames(plotdata)))
colnames(plotdata)[13:18] <- c("ILC1 vs ILC2", "ILC1 vs ILC3", "ILC2 vs ILC3", "ILC3 vs LTi", "ILC2 vs LTi", "ILC1 vs LTi")
colnames(plotdata)[1:12] <- toupper(colnames(plotdata)[1:12])


pdf("results/figures/Figure-S4E-upset-cytokine-proportion.pdf")
ComplexUpset::upset(
  plotdata,
  colnames(plotdata),
  min_degree = 2,
  min_size = 2,
  width_ratio=0.1,
  name = "intersections"
) + labs(title = "Intersection of Proportion and Cytokine QTL Loci")
dev.off()



