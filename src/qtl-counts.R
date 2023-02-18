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

write.csv(ilc1_eqtl_loci_by_gene,"results/eqtl/qtl-loci-by-gene-lods-ILC1.csv")
ilc1_eqtl_loci_by_gene %>% 
  dplyr::select(gene, gene_chr, gene_start, gene_end, 
                loci, loci_chr = marker_chr, loci_pos = marker_pos, cis_effect) %>%
  write.csv("results/eqtl/qtl-loci-by-genes-ILC1.csv")


# Counting number of QTLs per eGene
ilc1_eqtl_count <- ilc1_eqtl_loci_by_gene %>%
  group_by(gene) %>%
  summarise(qtl_count = sum(lod > lod_cutoff)) %>%
  mutate(cell_type ="ILC1",
         trait = paste0("ILC1 eQTL: ", gene)) %>%
  dplyr::rename(name = gene)


ilc1_counts <- ilc1_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  mutate(ILC1 = as.numeric(lod > lod_cutoff)) %>% 
  unite("loci", gene:xpos_floor) %>% 
  dplyr::select(-lod) %>%
  filter(ILC1 == 1)



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

write.csv(ilc2_eqtl_loci_by_gene,"results/eqtl/qtl-loci-by-gene-lods-ILC2.csv")
ilc2_eqtl_loci_by_gene %>% 
  dplyr::select(gene, gene_chr, gene_start, gene_end, 
                loci, loci_chr = marker_chr, loci_pos = marker_pos, cis_effect) %>%
  write.csv("results/eqtl/qtl-loci-by-genes-ILC2.csv")



# Counting number of QTLs per eGene
ilc2_eqtl_count <- ilc2_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  group_by(gene) %>%
  summarise(qtl_count = sum(lod > lod_cutoff)) %>%
  mutate(cell_type ="ILC2",
         trait = paste0("ILC2 eQTL: ", gene)) %>%
  dplyr::rename(name = gene)


ilc2_counts <- ilc2_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  mutate(ILC2 = as.numeric(lod > lod_cutoff)) %>% 
  unite("loci", gene:xpos_floor) %>% 
  dplyr::select(-lod) %>%
  filter(ILC2 == 1)



## ILC3
ilc3_eqtl <- fread("results/eqtl/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz", data.table = FALSE)

## creating loci by gene long table
ilc3_eqtl_loci_by_gene <- ilc3_eqtl %>%
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

write.csv(ilc3_eqtl_loci_by_gene,"results/eqtl/qtl-loci-by-gene-lods-ILC3.csv")
ilc3_eqtl_loci_by_gene %>% 
  dplyr::select(gene, gene_chr, gene_start, gene_end, 
                loci, loci_chr = marker_chr, loci_pos = marker_pos, cis_effect) %>%
  write.csv("results/eqtl/qtl-loci-by-genes-ILC3.csv")

# Counting number of QTLs per eGene
ilc3_eqtl_count <- ilc3_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  group_by(gene) %>%
  summarise(qtl_count = sum(lod > lod_cutoff)) %>%
  mutate(cell_type ="ILC3",
         trait = paste0("ILC3 eQTL: ", gene)) %>%
  dplyr::rename(name = gene)


ilc3_counts <- ilc3_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  mutate(ILC3 = as.numeric(lod > lod_cutoff)) %>% 
  unite("loci", gene:xpos_floor) %>% 
  dplyr::select(-lod) %>%
  filter(ILC3 == 1)



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

write.csv(lti_eqtl_loci_by_gene,"results/eqtl/qtl-loci-by-gene-lods-LTi.csv")
lti_eqtl_loci_by_gene %>% 
  dplyr::select(gene, gene_chr, gene_start, gene_end, 
                loci, loci_chr = marker_chr, loci_pos = marker_pos, cis_effect) %>%
  write.csv("results/eqtl/qtl-loci-by-genes-LTi.csv")
# Counting number of QTLs per eGene
lti_eqtl_count <- lti_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  group_by(gene) %>%
  summarise(qtl_count = sum(lod > lod_cutoff)) %>%
  mutate(cell_type ="LTi-like",
         trait = paste0("LTi-like eQTL: ", gene)) %>%
  dplyr::rename(name = gene)


lti_counts <- lti_eqtl %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  mutate(LTi = as.numeric(lod > lod_cutoff)) %>% 
  unite("loci", gene:xpos_floor) %>% 
  dplyr::select(-lod) %>%
  filter(LTi == 1)





bar_plot_df <- list("ILC1" = ilc1_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
                    "ILC2" = ilc2_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
                    "ILC3" = ilc3_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
                    "LTi" = lti_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect)) %>% 
  bind_rows(.id = "cell_type") %>%
  mutate(cis_effect = factor(cis_effect, levels = 0:1, labels = c("trans", "cis")))



ggplot(bar_plot_df, aes(cell_type, n, group = cis_effect, fill = cis_effect)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() + 
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "")
  


## proportion qtl counts ####

ilc1_ilc2 <- fread("results/proportion/gwas-ilc1-ilc2-results.csv.gz", data.table = FALSE)

# Counting number of QTL
ilc1_ilc2_count <- ilc1_ilc2 %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(chr, xpos_floor) %>% 
  summarise(lod = max(lods)) %>%
  ungroup() %>%
  summarise(qtl_count = sum(lod > 100)) %>%
  mutate(name = "ILC1 vs ILC2",
         trait = "ILC1 vs ILC2 QTL",
         cell_type = NA)



ilc1_ilc3 <- fread("results/proportion/gwas-ilc1-ilc3-results.csv.gz", data.table = FALSE) 

# Counting number of QTL
ilc1_ilc3_count <- ilc1_ilc3 %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(chr, xpos_floor) %>% 
  summarise(lod = max(lods)) %>%
  ungroup() %>%
  summarise(qtl_count = sum(lod > 100)) %>%
  mutate(name = "ILC1 vs ILC3",
         trait = "ILC1 vs ILC3 QTL",
         cell_type = NA)



ilc1_lti <- fread("results/proportion/gwas-ilc1-lti-results.csv.gz", data.table = FALSE)

# Counting number of QTL
ilc1_lti_count <- ilc1_lti %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(chr, xpos_floor) %>% 
  summarise(lod = max(lods)) %>%
  ungroup() %>%
  summarise(qtl_count = sum(lod > 100)) %>%
  mutate(name = "ILC1 vs LTi-like",
         trait = "ILC1 vs LTi-like QTL",
         cell_type = NA)



ilc2_ilc3 <- fread("results/proportion/gwas-ilc2-ilc3-results.csv.gz", data.table = FALSE)

# Counting number of QTL
ilc2_ilc3_count <- ilc2_ilc3 %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(chr, xpos_floor) %>% 
  summarise(lod = max(lods)) %>%
  ungroup() %>%
  summarise(qtl_count = sum(lod > 100)) %>%
  mutate(name = "ILC2 vs ILC3",
         trait = "ILC2 vs ILC3 QTL",
         cell_type = NA)



ilc2_lti <- fread("results/proportion/gwas-ilc2-lti-results.csv.gz", data.table = FALSE)

# Counting number of QTL
ilc2_lti_count <- ilc2_lti %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(chr, xpos_floor) %>% 
  summarise(lod = max(lods)) %>%
  ungroup() %>%
  summarise(qtl_count = sum(lod > 100)) %>%
  mutate(name= "ILC2 vs LTi",
         trait = "ILC2 vs LTi QTL",
         cell_type = NA)



ilc3_lti <- fread("results/proportion/gwas-ilc3-lti-results.csv.gz", data.table = FALSE)

# Counting number of QTL
ilc3_lti_count <- ilc3_lti %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(chr, xpos_floor) %>% 
  summarise(lod = max(lods)) %>%
  ungroup() %>%
  summarise(qtl_count = sum(lod > 100)) %>%
  mutate(name = "ILC3 vs LTi",
         trait = "ILC3 vs LTi QTL",
         cell_type = NA)


ilc3_stressed <- fread("results/proportion/ILC3_stressed_vs_non_qtl.csv.gz", data.table = FALSE)

# Counting number of QTL
ilc3_stressed_count <- ilc3_stressed %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(chr, xpos_floor) %>% 
  summarise(lod = max(lods)) %>%
  ungroup() %>%
  summarise(qtl_count = sum(lod > 200)) %>%
  mutate(name = "ILC3 vs LTi",
         trait = "ILC3 vs LTi QTL",
         cell_type = NA)



lti_stressed <- fread("results/proportion/LTi_stressed_vs_non_qtl.csv.gz", data.table = FALSE)

# Counting number of QTL
lti_stressed_count <- lti_stressed %>% 
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(chr, xpos_floor) %>% 
  summarise(lod = max(lods)) %>%
  ungroup() %>%
  summarise(qtl_count = sum(lod > 200)) %>%
  mutate(name = "ILC3 vs LTi",
         trait = "ILC3 vs LTi QTL",
         cell_type = NA)


ilc3_stressed_count + lti_stressed_count




## topic qtl counts ####

topics <- fread("results/topics/qtl-topic-lods.csv.gz", data.table = FALSE) %>%
  dplyr::select(-V1) %>%
  left_join(ccre) %>%
  dplyr::select(-ends_with("_qtl")) %>%
  mutate(pos_floor = floor(pos/1e6)) %>% 
  group_by(chr, pos_floor) %>% 
  summarise(across(starts_with("topic"), .fns = max)) %>%
  dplyr::select(chr, pos_floor, starts_with("topic")) %>%
  pivot_longer(-c(chr, pos_floor)) %>%
  group_by(name) %>%
  summarise(qtl_count = sum(value > 6)) %>%
  dplyr::mutate(trait = name,
                cell_type = NA)



## cytokine qtl counts

cytokines <- fread("results/cytokine/qtl-steady-cytokines-lods.csv.gz", data.table = FALSE) %>%
  dplyr::select(-V1) %>%
  left_join(ccre) %>%
  dplyr::select(1:15) %>%
  mutate(pos_floor = floor(pos/1e6)) %>% 
  group_by(chr, pos_floor) %>% 
  summarise(across(.cols = everything(), .fns = max)) %>%
  dplyr::select(-marker, -pos) %>%
  pivot_longer(-c(chr, pos_floor)) %>%
  group_by(name) %>%
  summarise(qtl_count = sum(value > 6)) %>%
  dplyr::mutate(trait = name,
                cell_type = NA)




outdf <- bind_rows(list(ilc1_eqtl_count, ilc1_ilc2_count, ilc1_ilc3_count,
                        ilc1_lti_count,  ilc2_eqtl_count, ilc2_ilc3_count, 
                        ilc2_lti_count, ilc3_eqtl_count,
                        ilc3_lti_count, lti_eqtl_count, topics, cytokines))
write.csv(outdf, "qtl-counts.csv", row.names = F)

