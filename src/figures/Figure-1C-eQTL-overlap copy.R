#' @title Figure 1C eQTL Overlap
#' @author Kirk Gosik
#' @description
#'
#'


library(data.table)
library(tidyverse)
library(ggvenn)

figure_path <- paste0("./results/figures/")
data_path <- "data/"
source(paste0(figure_path, "helpers.R"))

## Genes in eQTL loci ######
# vars <- fread(paste0(data_path, "vars.csv"))
## read in genesets
# dir('~/Google Drive/My Drive/projects/domice/paper of QTL study/Revised materials of ILC-QTL paper cell Science format/data/results/eqtl', "qtl-loci-by-genes-", full.names = T) 
ilc1_eqtl_loci_by_gene <- read_csv("results/eqtl/qtl-loci-by-genes-ILC1.csv")
ilc2_eqtl_loci_by_gene <- read_csv("results/eqtl/qtl-loci-by-genes-ILC2.csv")
ilc3_eqtl_loci_by_gene <- read_csv("results/eqtl/qtl-loci-by-genes-ILC3.csv")
lti_eqtl_loci_by_gene <- read_csv("results/eqtl/qtl-loci-by-genes-LTi.csv")

df <- bind_rows(list(ilc1 = ilc1_eqtl_loci_by_gene,
                     ilc2 = ilc2_eqtl_loci_by_gene,
                     ilc3 = ilc3_eqtl_loci_by_gene,
                     lti = lti_eqtl_loci_by_gene),.id='cell_type')


df_list_eqtls <- list('ILC1' = unique(ilc1_eqtl_loci_by_gene$loci), 
                'ILC2' = unique(ilc2_eqtl_loci_by_gene$loci), 
                'ILC3' = unique(ilc3_eqtl_loci_by_gene$loci),
                'LTi' = unique(lti_eqtl_loci_by_gene$loci))


p1 <- ggvenn(df_list_eqtls, set_name_size = 10) +
  labs(title = paste0("eQTL: ", n_distinct(df$loci))) +
  theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5, vjust = -5))



ggsave(filename = "./results/figures/Figure1C-eQTL-overlap.pdf",
       plot = p1,
       dpi = 330,
       width = 7,
       height = 7)



## eGenes ######
df_list_egenes <- list('ILC1' = unique(ilc1_eqtl_loci_by_gene$gene), 
                       'ILC2' = unique(ilc2_eqtl_loci_by_gene$gene),
                       'ILC3' = unique(ilc3_eqtl_loci_by_gene$gene),
                       'LTi' = unique(lti_eqtl_loci_by_gene$gene))

p2 <- ggvenn(df_list_egenes, set_name_size = 10) +
  labs(title = paste0("eGene: ", n_distinct(df$gene))) +
  theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5, vjust = -5))



ggsave(filename = "./results/figures/Figure1C-eGenes-overlap.pdf",
       plot = p2,
       dpi = 330,
       width = 7,
       height = 7)








## xxArchive #################################
## eQTL Associations #########################

## eQTL counts by gene
lod_cutoff <- 6


## ILC1
ilc1_eqtl <- fread("results/qtl-plot-lods-ILC1-cv.csv.gz", data.table = FALSE)


ilc1_counts <- ilc1_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  mutate(ILC1 = as.numeric(lod > lod_cutoff)) %>% 
  unite("loci", gene:xpos_floor) %>% 
  dplyr::select(-lod) %>%
  filter(ILC1 == 1)



## ILC2
ilc2_eqtl <- fread("results/qtl-plot-lods-ILC2-cv.csv.gz", data.table = FALSE)    



ilc2_counts <- ilc2_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  mutate(ILC2 = as.numeric(lod > lod_cutoff)) %>% 
  unite("loci", gene:xpos_floor) %>% 
  dplyr::select(-lod) %>%
  filter(ILC2 == 1)



## ILC3
ilc3_eqtl <- fread("results/qtl-plot-lods-Lti ILC3-cv.csv.gz", data.table = FALSE)


ilc3_counts <- ilc3_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  mutate(ILC3 = as.numeric(lod > lod_cutoff)) %>% 
  unite("loci", gene:xpos_floor) %>% 
  dplyr::select(-lod) %>%
  filter(ILC3 == 1)



## LTi
lti_eqtl <- fread("results/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz", data.table = FALSE)


lti_counts <- lti_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  mutate(LTi = as.numeric(lod > lod_cutoff)) %>% 
  unite("loci", gene:xpos_floor) %>% 
  dplyr::select(-lod) %>%
  filter(LTi == 1)



p1 <- ggvenn(list("ILC1" = ilc1_counts$loci, 
                  "ILC2" = ilc2_counts$loci, 
                  "ILC3" = ilc3_counts$loci, 
                  "LTi" = lti_counts$loci),
             set_name_size = 10) +
  labs(title = "eQTL") +
  theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5, vjust = -5))
