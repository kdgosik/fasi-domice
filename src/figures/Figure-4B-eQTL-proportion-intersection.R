#' @title Figure 4 B eQTL intersection among four phenotypic traits
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

## eQTL and proportion
plotdata <- vars[, c(217, 224, 231, 238, 274, 275, 276, 277, 287, 288)]
colnames(plotdata)[1:4] <- c("ILC1 eQTLs", "ILC2 eQTLs", "ILC3 eQTLs", "LTi eQTLs")
colnames(plotdata)[5:10] <- c("ILC1 vs ILC2", "ILC1 vs ILC3", "ILC2 vs ILC3", "ILC3 vs LTi", "ILC2 vs LTi", "ILC1 vs LTi")

pdf("results/figures/Figure-4B-upset-eQTL-proportions.pdf")
ComplexUpset::upset(
  plotdata,
  colnames(plotdata),
  min_degree = 2,
  min_size = 4,
  width_ratio=0.1,
  name = "intersections"
) + labs(title = "Intersection of eQTL of Proportions Loci")
dev.off()
