#' @title Figure S4D A.QTL intersection among cytokine X eQTL
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


## Cytokine steady and eQTL
plotdata <- vars[, c(217, 224,231,238, 292:303)]
colnames(plotdata) <- gsub("_qtl_genes_steady", "", colnames(plotdata))
colnames(plotdata)[1:4] <- c("ILC1 eQTLs", "ILC2 eQTLs", "ILC3 eQTLs", "LTi eQTLs")
colnames(plotdata)[5:16] <- toupper(colnames(plotdata)[5:16])

pdf("results/figures/Figure-S4D-upset-cytokine-eqtls.pdf")
ComplexUpset::upset(
  plotdata,
  colnames(plotdata),
  min_degree = 2,
  min_size = 4,
  width_ratio=0.1,
  name = "intersections"
) + labs(title = "Intersection of eQTL and Cytokine QTL Loci")
dev.off()


