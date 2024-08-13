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


## eQTL and programs
plotdata <- vars[, c(217, 224,231,238, 260:271)]
colnames(plotdata)[1:4] <- c("ILC1 eQTLs", "ILC2 eQTLs", "ILC3 eQTLs", "LTi eQTLs")
colnames(plotdata) <- gsub("_qtl_genes", "", colnames(plotdata))

pdf("results/figures/Figure-4A-upset-eQTL-programs.pdf")
ComplexUpset::upset(
  plotdata,
  colnames(plotdata),
  min_degree = 2,
  min_size = 6,
  width_ratio=0.1,
  name = "intersections"
) + labs(title = "Intersection of Program and eQTL Loci")
dev.off()
