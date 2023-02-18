#' @title eQTL Associations
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


project_dir <- "/Users/kirkgosik/Google Drive/My Drive/projects/"
domice_dir <- paste0(project_dir, "domice/paper of QTL study/Revised materials of ILC-QTL paper cell Science format/")
data_dir <- paste0(domice_dir, "data/")
results_dir <- paste0(data_dir, "results/")


ensembl <- fread(paste0(data_path, "references/Mus_musculus.GRCm38.102.gtf"),
                 col.names = c("seqid", "type", "source", "start", "end", 
                               "V6", "strand", "V8", "annotation"),
                 data.table = FALSE) %>%
  filter(seqid %in% c(as.character(1:19), "X"), 
         str_detect(annotation, 'gene_biotype \"protein_coding\"'),
         source == "gene") %>%
  mutate(gene_id = str_extract(annotation, 'ENSMUSG[0-9]{1,12}'),
         gene_name = str_extract(annotation, '[A-Z][a-z]{1,4}[0-9]{1,4}'),
         chr = paste0("chr", seqid))


geneids <- ensembl %>% dplyr::select(marker_gene_symbol = gene_name, ensembl_gene = gene_id)

## eQTL Map #################################

## Plot file list
plot_files <- c(
  paste0("results/qtl-plot-lods-ILC1-cv.csv.gz"),
  paste0("results/qtl-plot-lods-ILC2-cv.csv.gz"),
  paste0("results/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz"),
  paste0("results/qtl-plot-lods-Lti ILC3-cv.csv.gz")
)


## Loading Marker Map
cat("Loading SNPs...", "\n")
# load(paste0(geno_path, "GM_snps.Rdata"))
GM_snps <- read_csv(paste0(data_path, "references/GM_SNPS_Consequence_cCRE.csv"))


## eQTL Associations #########################

## eQTL counts by gene
lod_cutoff <- 6


## ILC1 ##########
ilc1_eqtl <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-ILC1-cv.csv.gz"), 
                   data.table = FALSE) %>%
  dplyr::rename(eGene = gene, eGene_id = id) %>%
  left_join(
    {GM_snps %>% 
        mutate(marker_pos = pos) %>% 
        select(marker = marker, 
               marker_pos, 
               conseq_clean,
               # type,
               ensembl_gene)}
  ) %>%
  left_join(
    {ensembl %>% dplyr::select(eGene_start = start, 
                               eGene_end = end, 
                               eGene_id = gene_id)}
  ) %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (eGene_start - 1e6) & marker_pos < (eGene_end + 1e6)})) %>%
  mutate(cell_type = "ILC1",
         value_adj = ifelse(cis_effect==1, value + 6, value)) %>%
  left_join(geneids) %>%
  filter(value_adj > 10) %>%
  dplyr::select(-value, -value_adj)



write.csv(ilc1_eqtl, paste0(results_dir, "eqtl/ilc1_eqtl-genes_by_egenes.csv"))



## ILC2 ########
ilc2_eqtl <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-ILC2-cv.csv.gz"), 
                   data.table = FALSE)  %>%
  dplyr::rename(eGene = gene, eGene_id = id) %>%
  left_join(
    {GM_snps %>% 
        mutate(marker_pos = pos) %>% 
        select(marker = marker, 
               marker_pos, 
               conseq_clean,
               # type,
               ensembl_gene)}
  ) %>%
  left_join(
    {ensembl %>% dplyr::select(eGene_start = start, 
                               eGene_end = end, 
                               eGene_id = gene_id)}
  ) %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (eGene_start - 1e6) & marker_pos < (eGene_end + 1e6)})) %>%
  mutate(cell_type = "ILC2",
         value_adj = ifelse(cis_effect==1, value + 6, value)) %>%
  left_join(geneids) %>%
  filter(value_adj > 10) %>%
  dplyr::select(-value, -value_adj)

write.csv(ilc2_eqtl, paste0(results_dir, "eqtl/ilc2_eqtl-genes_by_egenes.csv"))

## ILC3 #### 
# "results/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz" 
ilc3_eqtl <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz"), 
                   data.table = FALSE) %>%
  dplyr::rename(eGene = gene, eGene_id = id) %>%
  left_join(
    {GM_snps %>% 
        mutate(marker_pos = pos) %>% 
        select(marker = marker, 
               marker_pos, 
               conseq_clean,
               # type,
               ensembl_gene)}
  ) %>%
  left_join(
    {ensembl %>% dplyr::select(eGene_start = start, 
                               eGene_end = end, 
                               eGene_id = gene_id)}
  ) %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (eGene_start - 1e6) & marker_pos < (eGene_end + 1e6)})) %>%
  mutate(cell_type = "ILC3",
         value_adj = ifelse(cis_effect==1, value + 6, value)) %>%
  left_join(geneids) %>%
  filter(value_adj > 10) %>%
  dplyr::select(-value, -value_adj)

write.csv(ilc3_eqtl, paste0(results_dir, "eqtl/ilc3_eqtl-genes_by_egenes.csv"))


## LTi #### 
# "results/qtl-plot-lods-Lti ILC3-cv.csv.gz"
lti_eqtl <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-Lti ILC3-cv.csv.gz"), 
                  data.table = FALSE) %>%
  dplyr::rename(eGene = gene, eGene_id = id) %>%
  left_join(
    {GM_snps %>% 
        mutate(marker_pos = pos) %>% 
        select(marker = marker, 
               marker_pos, 
               conseq_clean,
               # type,
               ensembl_gene)}
  ) %>%
  left_join(
    {ensembl %>% dplyr::select(eGene_start = start, 
                               eGene_end = end, 
                               eGene_id = gene_id)}
  ) %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (eGene_start - 1e6) & marker_pos < (eGene_end + 1e6)})) %>%
  mutate(cell_type = "LTi-like",
         value_adj = ifelse(cis_effect==1, value + 6, value)) %>%
  left_join(geneids) %>%
  filter(value_adj > 10) %>%
  dplyr::select(-value, -value_adj)


write.csv(lti_eqtl, paste0(results_dir, "eqtl/lti_eqtl-genes_by_egenes.csv"))
