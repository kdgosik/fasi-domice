#' Checking to see if the loci containing Ccl17 affect ILC3 stressed, Cd48 and Frmd4b
#' and Nmu loci
#'
#'
#'
#'
#'

library(data.table)
# library(tidyverse)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(readr)
# library(ggvenn)
# library(igraph)
# library(ggforce)


source("../fasi-domice/setup.R")

# GM_Snps meta data
ccre <- fread(paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"), data.table = FALSE)

# gene level tracking
vars <- read_csv("/workspace/fasi-domice/data/allchannels/vars.csv")



### Ccl17 (ENSMUSG00000031780)
## Ccl17 SNPs - Chromosome 8: 95,537,081-95,538,664 forward strand.
start_irange <- 95000000
end_irange <- 96000000
chr_str <- "chr8"
chr_num <- "8"
gen <- "mm10"



vars %>% 
  filter(chr == chr_num, start > start_irange, end < end_irange)  %>% 
  dplyr::select(index, ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed)





## Nmu - Chr5:76333495..76363777
start_irange <- 76000000
end_irange <- 77000000
chr_str <- "chr5"
chr_num <- "5"
gen <- "mm10"


vars %>% 
  filter(chr == chr_num, start > start_irange, end < end_irange)  %>% 
  dplyr::select(index, ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed)


