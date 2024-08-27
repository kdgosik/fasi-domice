library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

project_path <- my_path <-"./domice/"
project_path <-"~/Documents/projects/domice/"
figure_path <- paste0(project_path, "results/figures/")
data_path <- paste0(project_path, "data/")
source(paste0(figure_path, "helpers.R"))


vars <- fread(paste0(data_dir, "allchannels/vars.csv"), data.table = F)
ccre <- fread(paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"), data.table = F)

## add polygenic vs monogenic feature? ##################


## eQTLs #######

## ILC1 #########################

# ILC1 eGenes but loci gene is not expresssed in ILC1

# getting all egenes assocatied with ILC1
ilc1_egenes <- vars %>% filter(ilc1_egenes_cv == 1) %>% pull(index)

# reading in LOD scores
ilc1_lods <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-ILC1-cv.csv.gz"), data.table = FALSE) 

# Counting number of QTLs per eGene
ilc1_polygenic <- ilc1_lods %>% 
  filter(gene %in% ilc1_egenes) %>%
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  group_by(gene) %>%
  summarise(ilc1_polygenic_egene = as.numeric(sum(lod > 5)>1)) %>%
  dplyr::rename(index = gene)


## do i want eGenes or eqtl_loci_genes?

ilc1 <- vars %>% 
  left_join(ilc1_polygenic) %>%
  dplyr::count(ilc1_expressed, ilc1_egenes_cv, ilc1_polygenic_egene, ilc1_eqtl_loci_genes_cv, ilc1_ciseqtl_assoc_cv) %>% 
  # dplyr::filter(ilc1_eqtl_loci_genes_cv == 1) %>%
  dplyr::rename(intrinsic = ilc1_expressed,
                egene = ilc1_egenes_cv,
                loci_gene = ilc1_eqtl_loci_genes_cv,
                polygenic = ilc1_polygenic_egene,
                cis_effect = ilc1_ciseqtl_assoc_cv) %>%
  dplyr::mutate(cell_type = "ILC1",
                percent = round(n / sum(n), 4))
                




## ILC2 #########################

# ILC2 eGenes but loci gene is not expresssed in ILC2

# getting all egenes assocatied with ILC2
ilc2_egenes <- vars %>% filter(ilc2_egenes_cv == 1) %>% pull(index)

# reading in LOD scores
ilc2_lods <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-ILC2-cv.csv.gz"), data.table = FALSE) 

# Counting number of QTLs per eGene
ilc2_polygenic <- ilc2_lods %>% 
  filter(gene %in% ilc2_egenes) %>%
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  group_by(gene) %>%
  summarise(ilc2_polygenic_egene = as.numeric(sum(lod > 5)>1)) %>%
  dplyr::rename(index = gene)


## do i want eGenes or eqtl_loci_genes?

ilc2 <- vars %>% 
  left_join(ilc2_polygenic) %>%
  dplyr::count(ilc2_expressed, ilc2_egenes_cv, ilc2_polygenic_egene, ilc2_eqtl_loci_genes_cv, ilc2_ciseqtl_assoc_cv) %>% 
  # dplyr::filter(ilc2_eqtl_loci_genes_cv == 1) %>%
  dplyr::rename(intrinsic = ilc2_expressed,
                egene = ilc2_egenes_cv,
                loci_gene = ilc2_eqtl_loci_genes_cv,
                polygenic = ilc2_polygenic_egene,
                cis_effect = ilc2_ciseqtl_assoc_cv) %>%
  dplyr::mutate(cell_type = "ILC2",
                percent = round(n / sum(n), 4))





## ILC3 #########################

# ILC3 eGenes but loci gene is not expresssed in ilc3

# getting all egenes assocatied with ilc3
ilc3_egenes <- vars %>% filter(ilc3_egenes_cv == 1) %>% pull(index)

# reading in LOD scores
ilc3_lods <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-NCR1\\+\\ ILC3-cv.csv.gz"), data.table = FALSE) 

# Counting number of QTLs per eGene
ilc3_polygenic <- ilc3_lods %>% 
  filter(gene %in% ilc3_egenes) %>%
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  group_by(gene) %>%
  summarise(ilc3_polygenic_egene = as.numeric(sum(lod > 5)>1)) %>%
  dplyr::rename(index = gene)


## do i want eGenes or eqtl_loci_genes?

ilc3 <- vars %>% 
  left_join(ilc3_polygenic) %>%
  dplyr::count(ilc3_expressed, ilc3_egenes_cv, ilc3_polygenic_egene, ilc3_eqtl_loci_genes_cv, ilc3_ciseqtl_assoc_cv) %>% 
  # dplyr::filter(ilc3_eqtl_loci_genes_cv == 1) %>%
  dplyr::rename(intrinsic = ilc3_expressed,
                egene = ilc3_egenes_cv,
                loci_gene = ilc3_eqtl_loci_genes_cv,
                polygenic = ilc3_polygenic_egene,
                cis_effect = ilc3_ciseqtl_assoc_cv) %>%
  dplyr::mutate(cell_type = "ILC3",
                percent = round(n / sum(n), 4))






## LTi-like #########################

# LTi eGenes but loci gene is not expresssed in LTi

# getting all egenes assocatied with LTi
lti_egenes <- vars %>% filter(lti_egenes_cv == 1) %>% pull(index)

# reading in LOD scores
lti_lods <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-Lti ILC3-cv.csv.gz"), data.table = FALSE) 

# Counting number of QTLs per eGene
lti_polygenic <- lti_lods %>% 
  filter(gene %in% lti_egenes) %>%
  mutate(xpos_floor = floor(xpos)) %>% 
  group_by(gene, marker_chr, xpos_floor) %>% 
  summarise(lod = max(value)) %>%
  group_by(gene) %>%
  summarise(lti_polygenic_egene = as.numeric(sum(lod > 5)>1)) %>%
  dplyr::rename(index = gene)


## do i want eGenes or eqtl_loci_genes?

lti <- vars %>% 
  left_join(lti_polygenic) %>%
  dplyr::count(lti_expressed, lti_egenes_cv, lti_polygenic_egene, lti_eqtl_loci_genes_cv, lti_ciseqtl_assoc_cv) %>% 
  # dplyr::filter(lti_eqtl_loci_genes_cv == 1) %>%
  dplyr::rename(intrinsic = lti_expressed,
                egene = lti_egenes_cv,
                loci_gene = lti_eqtl_loci_genes_cv,
                polygenic = lti_polygenic_egene,
                cis_effect = lti_ciseqtl_assoc_cv) %>%
  dplyr::mutate(cell_type = "LTi",
                percent = round(n / sum(n), 4))




df_ilc_poly <- list(ilc1, ilc2, ilc3, lti) %>% bind_rows()
write.csv(df_ilc_poly, paste0(data_path, "intrinsic-polygenic-by-celltype.csv"), row.names = F)


## proportion QTL: polygenic #######


prop_poly <- vars %>% 
  dplyr::select(index,272:277,287,288) %>% 
  pivot_longer(-index) %>% 
  filter(value==1) %>% 
  group_by(name) %>% 
  summarise(qtl_count = sum(value),
            group = "proportion")



## using a family-wide error rate #####

# gwas_lods <- fread(paste0(project_path, "results/gwas-ilc1-ilc2-results.csv.gz"), data.table = F) %>%
#   mutate(fwer = quantile(p_values, 0.01),
#          xpos_floor = floor(xpos)) %>%
#   filter(p_values < fwer) %>%
#   summarise(qtl_count = n_distinct(xpos_floor))





## topic QTL: polygenic #######

topic_poly <- fread(paste0(domice_dir, "results/topics/topic-qtl-lods.csv.gz")) %>%
  left_join({
    ccre %>% 
      dplyr::select(marker, chr, pos)
      }) %>%
  mutate(pos_floor = round(pos/1e6)) %>%
  dplyr::select(-V1, -pos) %>%
  pivot_longer(-c(marker,chr,pos_floor)) %>%
  group_by(name, chr, pos_floor) %>%
  summarise(lods = max(value)) %>%
  group_by(name) %>%
  summarise(qtl_count = sum(lods > 6),
            group = "topic")




## cytokine qtl: polygenic ###########################

cytokine_poly <- fread(paste0(data_path, "qtl-steady-cytokines-lods.csv.gz"), data.table = F) %>%
  left_join({
    ccre %>% 
      dplyr::select(marker, chr, pos)
  }) %>%
  mutate(pos_floor = round(pos/1e6)) %>%
  dplyr::select(-V1, -pos) %>%
  pivot_longer(-c(marker,chr,pos_floor)) %>%
  group_by(name, chr, pos_floor) %>%
  summarise(lods = max(value)) %>%
  group_by(name) %>%
  summarise(qtl_count = sum(lods > 6),
            group = "cytokine")



df_poly <- list(prop_poly, topic_poly, cytokine_poly) %>% 
  bind_rows() %>%
  mutate(polygenic = ifelse(qtl_count > 1, "polygenic", "monogenic"))


write.csv(df_poly, paste0(data_path, "proportion-topic-cytokine-qtl-count.csv"), row.names = F)


df_poly %>% count(group, polygenic) %>%
  write.csv(., paste0(data_path, "proportion-topic-cytokine-polygenic-count.csv"), row.names = F)





cytokine <- fread(paste0(data_path, "qtl-steady-cytokines-lods.csv.gz"), data.table = F) %>%
  left_join({
    ccre %>% 
      dplyr::select(marker, chr, pos, ensembl_gene)
  }) %>%
  left_join({
    vars %>% dplyr::select(ensembl_gene = gene_ids, index)
  }) %>% 
  select(3:14,18) %>% 
  pivot_longer(-index)


cytokine %>% 
  filter(value > 6, !is.na(index)) %>% 
  count(index, name) %>%
  dplyr::rename(gene = index,
                trait = name) %>%
  write.csv(., paste0(data_path, "cytokine-genelist.csv"), row.names = F)


