#' @title Figure 4B eQTL proportion intersection plot
#' @author Kirk Gosik
#' @description

library(tidyverse)
library(UpSetR)
library(ComplexUpset)

library(readr)
library(dplyr)
library(tidyr)


project_dir <- "/Users/kirkgosik/Google Drive/My Drive/projects/"
domice_dir <- paste0(project_dir, "domice/paper of QTL study/Revised materials of ILC-QTL paper cell Science format/")
data_dir <- paste0(domice_dir, "data/")
results_dir <- paste0(data_dir, "results/")


vars <- read_csv(paste0(data_dir, 'allchannels/vars.csv'))

intrinsic <- vars %>%
  dplyr::select(gene_ids, ilc1_expressed, ilc2_expressed, 
                ilc3_expressed, lti_expressed) %>%
  pivot_longer(-gene_ids) %>%
  filter(value == 1) %>%
  distinct(gene_ids) %>%
  mutate(intrinsic = 1)


ccre <- read_csv(paste0(data_dir, 'references/GM_SNPS_Consequence_cCRE.csv'))
ccre_summary_with_ensembl <- ccre %>%
  mutate(pos_floor = as.character(1e6*floor(pos / 1e6))) %>% 
  unite("loci", c(chr, pos_floor)) %>%
  group_by(loci, ensembl_gene) %>%
  summarise(across(contains("qtl"), max, na.rm=T)) %>%
  left_join(dplyr::select(vars, 
                          ensembl_gene = gene_ids, 
                          ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed))


ccre_summary <- ccre %>%
  mutate(pos_floor = as.character(1e6*floor(pos / 1e6))) %>% 
  unite("loci", c(chr, pos_floor)) %>%
  group_by(loci) %>%
  summarise(across(contains("qtl"), max, na.rm=T),
            marker = first(marker),
            gene_ids = first(ensembl_gene)) %>%
  left_join(intrinsic) %>%
  replace_na(list(intrinsic = 0))



## Figure 4B: eQTL x proportion QTL  #########################

eqtl_proportion_upsetdf <- ccre_summary %>%
  ungroup() %>%
  dplyr::select(contains("prop"), contains('eqtl_loci_cv')) 

colnames(eqtl_proportion_upsetdf) <- c("ILC1 vs ILC2","ILC2 vs ILC3",
                                       "ILC3 vs LTi","ILC1 vs LTi",
                                       "ILC2 vs LTi","ILC1 vs ILC3",
                                       "ILC1 eQTL", "ILC2 eQTL", "ILC3 eQTL", "LTi-like eQTL")


eqtl_proportion_upset<- ComplexUpset::upset(data = eqtl_proportion_upsetdf, 
                                            intersect = colnames(eqtl_proportion_upsetdf),
                                            min_degree = 2,
                                            min_size = 6,
                                            set_size = FALSE) + 
  ggplot2::labs(title = "eQTL x proportion QTL")


## save plot
ggplot2::ggsave(filename =  "results/figures/Figure-4B-eQTL-proportion-upset.pdf",
                plot = eqtl_proportion_upset,
                width = 7,
                height = 5,
                dpi = 330)

