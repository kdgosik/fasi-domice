#' @title Figure 7A Topic Scores Embeddings
#' @author Kirk Gosik
#' @description
#'
#'


library(data.table)
library(tidyverse)
library(ggvenn)
library(igraph)
library(ggforce)

project_path <- my_path <-"./"

figure_path <- paste0(project_path, "results/figures/")
data_path <- paste0(project_path, "data/")
source(paste0(figure_path, "helpers.R"))

# GM_Snps meta data
ccre <- fread(paste0(data_path, "references/GM_SNPS_Consequence_cCRE.csv"), data.table = FALSE)


## Gene Data
ensembl <- as.data.frame(readGFF(paste0(data_path, "references/Mus_musculus.GRCm38.102.gtf"))) %>%
  filter(seqid %in% c(as.character(1:19), "X"), gene_biotype == "protein_coding", type == "gene") %>%
  mutate(chr = paste0("chr", seqid))


## eQTL counts by gene ####
lod_cutoff <- 6


## ILC1 ####
ilc1_eqtl <- fread("results/eqtl/qtl-plot-lods-ILC1-cv.csv.gz", data.table = FALSE)


## creating loci by gene long table
ilc1_eqtl_loci_by_gene <- ilc1_eqtl %>%
  left_join(select(ccre, marker,pos)) %>%
  mutate(marker_pos_floor = floor(pos/1e6)) %>% 
  group_by(gene, marker_chr, marker_pos_floor) %>% 
  summarise(lod = max(value),
            marker_pos = median(pos)) %>%
  left_join(
    {ensembl %>% dplyr::select(gene_start = start, gene_end = end, gene = gene_name, gene_chr = seqid)}
  ) %>%
  ungroup() %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (gene_start - 1e6) & marker_pos < (gene_end + 1e6)})) %>%
  mutate(value_adj = ifelse(cis_effect==1, lod + 6, lod)) %>%
  filter(value_adj > 10) %>%
  unite("loci", c(marker_chr, marker_pos), remove = FALSE)


## ILC2
ilc2_eqtl <- fread("results/eqtl/qtl-plot-lods-ILC2-cv.csv.gz", data.table = FALSE)

## creating loci by gene long table
ilc2_eqtl_loci_by_gene <- ilc2_eqtl %>%
  left_join(select(ccre, marker,pos)) %>%
  mutate(marker_pos_floor = floor(pos/1e6)) %>% 
  group_by(gene, marker_chr, marker_pos_floor) %>% 
  summarise(lod = max(value),
            marker_pos = median(pos)) %>%
  left_join(
    {ensembl %>% dplyr::select(gene_start = start, gene_end = end, gene = gene_name, gene_chr = seqid)}
  ) %>%
  ungroup() %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (gene_start - 1e6) & marker_pos < (gene_end + 1e6)})) %>%
  mutate(value_adj = ifelse(cis_effect==1, lod + 6, lod)) %>%
  filter(value_adj > 10) %>%
  unite("loci", c(marker_chr, marker_pos), remove = FALSE)



## ILC3
ilc3_eqtl <- fread("results/eqtl/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz", data.table = FALSE)

## creating loci by gene long table
ilc2_eqtl_loci_by_gene <- ilc2_eqtl %>%
  left_join(select(ccre, marker,pos)) %>%
  mutate(marker_pos_floor = floor(pos/1e6)) %>% 
  group_by(gene, marker_chr, marker_pos_floor) %>% 
  summarise(lod = max(value),
            marker_pos = median(pos)) %>%
  left_join(
    {ensembl %>% dplyr::select(gene_start = start, gene_end = end, gene = gene_name, gene_chr = seqid)}
  ) %>%
  ungroup() %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (gene_start - 1e6) & marker_pos < (gene_end + 1e6)})) %>%
  mutate(value_adj = ifelse(cis_effect==1, lod + 6, lod)) %>%
  filter(value_adj > 10) %>%
  unite("loci", c(marker_chr, marker_pos), remove = FALSE)



## LTi
lti_eqtl <- fread("results/eqtl/qtl-plot-lods-Lti ILC3-cv.csv.gz", data.table = FALSE)

## creating loci by gene long table
lti_eqtl_loci_by_gene <- lti_eqtl %>%
  left_join(select(ccre, marker,pos)) %>%
  mutate(marker_pos_floor = floor(pos/1e6)) %>% 
  group_by(gene, marker_chr, marker_pos_floor) %>% 
  summarise(lod = max(value),
            marker_pos = median(pos)) %>%
  left_join(
    {ensembl %>% dplyr::select(gene_start = start, gene_end = end, gene = gene_name, gene_chr = seqid)}
  ) %>%
  ungroup() %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (gene_start - 1e6) & marker_pos < (gene_end + 1e6)})) %>%
  mutate(value_adj = ifelse(cis_effect==1, lod + 6, lod)) %>%
  filter(value_adj > 10) %>%
  unite("loci", c(marker_chr, marker_pos), remove = FALSE)




bar_plot_df <- list("ILC1" = ilc1_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
                    "ILC2" = ilc2_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
                    "ILC3" = ilc3_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
                    "LTi" = lti_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect)) %>% 
  bind_rows(.id = "cell_type") %>%
  mutate(cis_effect = factor(cis_effect, levels = 0:1, labels = c("trans", "cis")))



p1 <- ggplot(bar_plot_df, aes(cell_type, n, group = cis_effect, fill = cis_effect)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() + 
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "Proportion of trans-eQTL Loci",
       x = "Cell Type",
       y = "Number of Loci",
       fill = "Loci Type")



ggsave(filename = paste0(figure_path, "Figure-7A-proportion-trans-eQTL-loci.pdf"),
       plot = p1,
       dpi = 330,
       width = 8.5,
       height = 5)
