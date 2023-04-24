library(data.table)
library(rtracklayer)
# library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(ggvenn)
library(igraph)
library(ggforce)

# BiocManager::install(c("GenomicAlignments", "GenomicRanges"))

# project_dir <- "/Users/kirkgosik/Google Drive/My Drive/projects/"
# domice_dir <- paste0(project_dir, "domice/paper of QTL study/Revised materials of ILC-QTL paper cell Science format/")
# data_dir <- paste0(domice_dir, "data/")
# results_dir <- paste0(data_dir, "results/")
# 
# figure_path <- paste0(project_path, "results/figures/")
# data_path <- paste0(project_path, "data/")
# source(paste0(figure_path, "helpers.R"))



# GM_Snps meta data
ccre <- fread(paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"), data.table = FALSE)


## Gene Data
ensembl <- as.data.frame(readGFF(paste0(data_dir, "references/Mus_musculus.GRCm38.102.gtf"))) %>%
  filter(seqid %in% c(as.character(1:19), "X"), gene_biotype == "protein_coding", type == "gene") %>%
  mutate(chr = paste0("chr", seqid))


## eQTL counts by gene ####
lod_cutoff <- 6


## ILC1 ####
ilc1_eqtl <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-ILC1-cv.csv.gz"), 
                   data.table = FALSE)

## creating loci by gene long table
ilc1_eqtl_loci_by_gene <- ilc1_eqtl %>%
  left_join(select(ccre, marker, pos)) %>%
  mutate(gene_chr = as.character(gene_chr),
         marker_pos = pos) %>%
  left_join(
    {ensembl %>% dplyr::select(gene_start = start, 
                               gene_end = end, 
                               gene = gene_name, 
                               gene_chr = seqid) %>%
        mutate(gene_chr = as.character(gene_chr))}
  ) %>%
  ungroup() %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (gene_start - 1e6) & marker_pos < (gene_end + 1e6)})) %>%
  mutate(value_adj = ifelse(cis_effect==1, value + 6, value),
         strand = "+") %>%
  filter(value_adj > 10) %>%
  dplyr::select(marker, gene, value, marker_chr, pos, strand,
                id, gene_chr, gene_start, gene_end, 
                cis_effect, value_adj)

ilc1eqtlgrl <- makeGRangesListFromDataFrame(ilc1_eqtl_loci_by_gene,
                                            split.field = "gene",
                                            names.field = "marker",
                                            seqnames = "marker_chr",
                                            start.field = "pos",
                                            end.field = "pos",
                                            strand = "strand",
                                            keep.extra.columns = TRUE)

## length of all genes in grl
ilc1_eqtl_list <- lapply(reduce(resize(sort(ilc1eqtlgrl),10000)), length)



# write.csv(ilc1_eqtl_loci_by_gene,"results/eqtl/qtl-loci-by-gene-lods-ILC1.csv")
# ilc1_eqtl_loci_by_gene %>% 
#   dplyr::select(gene, gene_chr, gene_start, gene_end, 
#                 loci, loci_chr = marker_chr, loci_pos = marker_pos, cis_effect) %>%
#   write.csv("results/eqtl/qtl-loci-by-genes-ILC1.csv")
# 
# 
# # Counting number of QTLs per eGene
# ilc1_eqtl_count <- ilc1_eqtl_loci_by_gene %>%
#   group_by(gene) %>%
#   summarise(qtl_count = sum(lod > lod_cutoff)) %>%
#   mutate(cell_type ="ILC1",
#          trait = paste0("ILC1 eQTL: ", gene)) %>%
#   dplyr::rename(name = gene)
# 
# 
# ilc1_counts <- ilc1_eqtl %>% 
#   mutate(xpos_floor = floor(xpos)) %>% 
#   group_by(gene, marker_chr, xpos_floor) %>% 
#   summarise(lod = max(value)) %>%
#   mutate(ILC1 = as.numeric(lod > lod_cutoff)) %>% 
#   unite("loci", gene:xpos_floor) %>% 
#   dplyr::select(-lod) %>%
#   filter(ILC1 == 1)



## ILC2 ####
ilc2_eqtl <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-ILC2-cv.csv.gz"), 
                   data.table = FALSE)

## creating loci by gene long table
ilc2_eqtl_loci_by_gene <- ilc2_eqtl %>%
  left_join(select(ccre, marker, pos)) %>%
  mutate(gene_chr = as.character(gene_chr),
         marker_pos = pos) %>%
  left_join(
    {ensembl %>% dplyr::select(gene_start = start, 
                               gene_end = end, 
                               gene = gene_name, 
                               gene_chr = seqid) %>%
        mutate(gene_chr = as.character(gene_chr))}
  ) %>%
  ungroup() %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (gene_start - 1e6) & marker_pos < (gene_end + 1e6)})) %>%
  mutate(value_adj = ifelse(cis_effect==1, value + 6, value)) %>%
  filter(value_adj > 10) %>%
  dplyr::select(marker, gene, value, marker_chr, pos, 
                id, gene_chr, gene_start, gene_end, cis_effect, value_adj)

ilc2eqtlgrl <- makeGRangesListFromDataFrame(ilc2_eqtl_loci_by_gene,
                                            split.field = "gene",
                                            names.field = "marker",
                                            seqnames="marker_chr",
                                            start.field="pos",
                                            end.field="pos",
                                            strand = "strand",
                                            keep.extra.columns = TRUE)

## length of all genes in grl
ilc2_eqtl_list <- lapply(reduce(resize(sort(ilc2eqtlgrl),10000)), length)

# write.csv(ilc2_eqtl_loci_by_gene,"results/eqtl/qtl-loci-by-gene-lods-ILC2.csv")
# ilc2_eqtl_loci_by_gene %>% 
#   dplyr::select(gene, gene_chr, gene_start, gene_end, 
#                 loci, loci_chr = marker_chr, loci_pos = marker_pos, cis_effect) %>%
#   write.csv("results/eqtl/qtl-loci-by-genes-ILC2.csv")
# 
# 
# 
# # Counting number of QTLs per eGene
# ilc2_eqtl_count <- ilc2_eqtl %>% 
#   mutate(xpos_floor = floor(xpos)) %>% 
#   group_by(gene, marker_chr, xpos_floor) %>% 
#   summarise(lod = max(value)) %>%
#   group_by(gene) %>%
#   summarise(qtl_count = sum(lod > lod_cutoff)) %>%
#   mutate(cell_type ="ILC2",
#          trait = paste0("ILC2 eQTL: ", gene)) %>%
#   dplyr::rename(name = gene)
# 
# 
# ilc2_counts <- ilc2_eqtl %>% 
#   mutate(xpos_floor = floor(xpos)) %>% 
#   group_by(gene, marker_chr, xpos_floor) %>% 
#   summarise(lod = max(value)) %>%
#   mutate(ILC2 = as.numeric(lod > lod_cutoff)) %>% 
#   unite("loci", gene:xpos_floor) %>% 
#   dplyr::select(-lod) %>%
#   filter(ILC2 == 1)



## ILC3 ####
ilc3_eqtl <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-NCR1\\+\\ ILC3-cv.csv.gz"), 
                   data.table = FALSE)

## creating loci by gene long table
ilc3_eqtl_loci_by_gene <- ilc3_eqtl %>%
  left_join(select(ccre, marker, pos)) %>%
  mutate(gene_chr = as.character(gene_chr),
         marker_pos = pos) %>%
  left_join(
    {ensembl %>% dplyr::select(gene_start = start, 
                               gene_end = end, 
                               gene = gene_name, 
                               gene_chr = seqid) %>%
        mutate(gene_chr = as.character(gene_chr))}
  ) %>%
  ungroup() %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (gene_start - 1e6) & marker_pos < (gene_end + 1e6)})) %>%
  mutate(value_adj = ifelse(cis_effect==1, value + 6, value)) %>%
  filter(value_adj > 10) %>%
  dplyr::select(marker, gene, value, marker_chr, pos, 
                id, gene_chr, gene_start, gene_end, cis_effect, value_adj)

ilc3eqtlgrl <- makeGRangesListFromDataFrame(ilc3_eqtl_loci_by_gene,
                                            split.field = "gene",
                                            names.field = "marker",
                                            seqnames="marker_chr",
                                            start.field="pos",
                                            end.field="pos",
                                            strand = "strand",
                                            keep.extra.columns = TRUE)

## length of all genes in grl
ilc3_eqtl_list <- lapply(reduce(resize(sort(ilc3eqtlgrl),10000)), length)

# write.csv(ilc3_eqtl_loci_by_gene,"results/eqtl/qtl-loci-by-gene-lods-ILC3.csv")
# ilc3_eqtl_loci_by_gene %>% 
#   dplyr::select(gene, gene_chr, gene_start, gene_end, 
#                 loci, loci_chr = marker_chr, loci_pos = marker_pos, cis_effect) %>%
#   write.csv("results/eqtl/qtl-loci-by-genes-ILC3.csv")
# 
# # Counting number of QTLs per eGene
# ilc3_eqtl_count <- ilc3_eqtl %>% 
#   mutate(xpos_floor = floor(xpos)) %>% 
#   group_by(gene, marker_chr, xpos_floor) %>% 
#   summarise(lod = max(value)) %>%
#   group_by(gene) %>%
#   summarise(qtl_count = sum(lod > lod_cutoff)) %>%
#   mutate(cell_type ="ILC3",
#          trait = paste0("ILC3 eQTL: ", gene)) %>%
#   dplyr::rename(name = gene)
# 
# 
# ilc3_counts <- ilc3_eqtl %>% 
#   mutate(xpos_floor = floor(xpos)) %>% 
#   group_by(gene, marker_chr, xpos_floor) %>% 
#   summarise(lod = max(value)) %>%
#   mutate(ILC3 = as.numeric(lod > lod_cutoff)) %>% 
#   unite("loci", gene:xpos_floor) %>% 
#   dplyr::select(-lod) %>%
#   filter(ILC3 == 1)



## LTi #####
lti_eqtl <- fread(paste0(results_dir, "eqtl/qtl-plot-lods-Lti ILC3-cv.csv.gz"), 
                         data.table = FALSE)

## creating loci by gene long table
lti_eqtl_loci_by_gene <- lti_eqtl %>%
  left_join(select(ccre, marker, pos)) %>%
  mutate(gene_chr = as.character(gene_chr),
         marker_pos = pos) %>%
  left_join(
    {ensembl %>% dplyr::select(gene_start = start, 
                               gene_end = end, 
                               gene = gene_name, 
                               gene_chr = seqid) %>%
        mutate(gene_chr = as.character(gene_chr))}
  ) %>%
  ungroup() %>%
  mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
                                   {marker_pos > (gene_start - 1e6) & marker_pos < (gene_end + 1e6)})) %>%
  mutate(value_adj = ifelse(cis_effect==1, value + 6, value)) %>%
  filter(value_adj > 10) %>%
  dplyr::select(marker, gene, value, marker_chr, pos, 
                id, gene_chr, gene_start, gene_end, cis_effect, value_adj)

ltieqtlgrl <- makeGRangesListFromDataFrame(lti_eqtl_loci_by_gene,
                                           split.field = "gene",
                                           names.field = "marker",
                                           seqnames="marker_chr",
                                           start.field="pos",
                                           end.field="pos",
                                           strand = "strand",
                                           keep.extra.columns = TRUE)

## length of all genes in grl
lti_eqtl_list <- lapply(reduce(resize(sort(ltieqtlgrl),10000)), length)


# write.csv(lti_eqtl_loci_by_gene,"results/eqtl/qtl-loci-by-gene-lods-LTi.csv")
# lti_eqtl_loci_by_gene %>% 
#   dplyr::select(gene, gene_chr, gene_start, gene_end, 
#                 loci, loci_chr = marker_chr, loci_pos = marker_pos, cis_effect) %>%
#   write.csv("results/eqtl/qtl-loci-by-genes-LTi.csv")
# # Counting number of QTLs per eGene
# lti_eqtl_count <- lti_eqtl %>% 
#   mutate(xpos_floor = floor(xpos)) %>% 
#   group_by(gene, marker_chr, xpos_floor) %>% 
#   summarise(lod = max(value)) %>%
#   group_by(gene) %>%
#   summarise(qtl_count = sum(lod > lod_cutoff)) %>%
#   mutate(cell_type ="LTi-like",
#          trait = paste0("LTi-like eQTL: ", gene)) %>%
#   dplyr::rename(name = gene)
# 
# 
# lti_counts <- lti_eqtl %>% 
#   mutate(xpos_floor = floor(xpos)) %>% 
#   group_by(gene, marker_chr, xpos_floor) %>% 
#   summarise(lod = max(value)) %>%
#   mutate(LTi = as.numeric(lod > lod_cutoff)) %>% 
#   unite("loci", gene:xpos_floor) %>% 
#   dplyr::select(-lod) %>%
#   filter(LTi == 1)





# bar_plot_df <- list("ILC1" = ilc1_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
#                     "ILC2" = ilc2_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
#                     "ILC3" = ilc3_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect),
#                     "LTi" = lti_eqtl_loci_by_gene %>% group_by(cis_effect) %>% distinct(loci) %>% count(cis_effect)) %>% 
#   bind_rows(.id = "cell_type") %>%
#   mutate(cis_effect = factor(cis_effect, levels = 0:1, labels = c("trans", "cis")))
# 
# 
# 
# ggplot(bar_plot_df, aes(cell_type, n, group = cis_effect, fill = cis_effect)) + 
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_minimal() + 
#   scale_fill_brewer(palette = "Dark2") +
#   labs(title = "")
  


## proportion qtl counts ####

ilc1_ilc2 <- fread(paste0(results_dir,"proportions/gwas-ilc1-ilc2-results.csv.gz"), 
                   data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.999)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand, xpos)

ilc1_ilc2gr <- makeGRangesFromDataFrame(ilc1_ilc2,
                                        seqnames = "chr",
                                        start.field = "pos",
                                        end.field = "pos",
                                        strand = "strand",
                                        keep.extra.columns = TRUE)

## length of all genes in grl
ilc1_ilc2_count <- length(reduce(resize(sort(ilc1_ilc2gr), 1000000)))



ilc1_ilc3 <- fread(paste0(results_dir,"proportions/gwas-ilc1-ilc3-results.csv.gz"), 
                   data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.999)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)

ilc1_ilc3gr <- makeGRangesFromDataFrame(ilc1_ilc3,
                                        seqnames="chr",
                                        start.field="pos",
                                        end.field="pos",
                                        strand = "strand",
                                        keep.extra.columns = TRUE,
                                        na.rm = TRUE)

## length of all genes in grl
ilc1_ilc3_count <- length(reduce(resize(sort(ilc1_ilc3gr), 1000000)))




ilc1_lti <- fread(paste0(results_dir,"proportions/gwas-ilc1-lti-results.csv.gz"), 
                   data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.999)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)

ilc1_ltigr <- makeGRangesFromDataFrame(ilc1_lti,
                                       seqnames="chr",
                                       start.field="pos",
                                       end.field="pos",
                                       strand = "strand",
                                       keep.extra.columns = TRUE,
                                       na.rm = TRUE)

## length of all genes in grl
ilc1_lti_count <- length(reduce(resize(sort(ilc1_ltigr), 1000000)))




ilc2_ilc3 <- fread(paste0(results_dir,"proportions/gwas-ilc2-ilc3-results.csv.gz"), 
                  data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.999)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)

ilc2_ilc3gr <- makeGRangesFromDataFrame(ilc2_ilc3,
                                       seqnames="chr",
                                       start.field="pos",
                                       end.field="pos",
                                       strand = "strand",
                                       keep.extra.columns = TRUE,
                                       na.rm = TRUE)

## length of all genes in grl
ilc2_ilc3_count <- length(reduce(resize(sort(ilc2_ilc3gr), 1000000)))



# quantile(ilc2_lti$lods, 0.999)
ilc2_lti <- fread(paste0(results_dir,"proportions/gwas-ilc2-lti-results.csv.gz"), 
                  data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.999)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)

ilc2_ltigr <- makeGRangesFromDataFrame(ilc2_lti,
                                       seqnames="chr",
                                       start.field="pos",
                                       end.field="pos",
                                       strand = "strand",
                                       keep.extra.columns = TRUE,
                                       na.rm = TRUE)

## length of all genes in grl
ilc2_lti_count <- length(reduce(resize(sort(ilc2_ltigr), 1000000)))



ilc3_lti <- fread(paste0(results_dir,"proportions/gwas-ilc3-lti-results.csv.gz"), 
                  data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.999)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)

ilc3_ltigr <- makeGRangesFromDataFrame(ilc3_lti,
                                       seqnames="chr",
                                       start.field="pos",
                                       end.field="pos",
                                       strand = "strand",
                                       keep.extra.columns = TRUE,
                                       na.rm = TRUE)

## length of all genes in grl
ilc3_lti_count <- length(reduce(resize(sort(ilc3_ltigr), 1000000)))





ilc3_stressed <- fread(paste0(results_dir,"proportions/ILC3_stressed_vs_non_qtl.csv.gz"), 
                       data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.999)) %>%
  dplyr::select(marker, chr, betas, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)

ilc3_stressedgr <- makeGRangesFromDataFrame(ilc3_stressed,
                                            seqnames="chr",
                                            start.field="pos",
                                            end.field="pos",
                                            strand = "strand",
                                            keep.extra.columns = TRUE,
                                            na.rm = TRUE)

# Counting number of QTL
## length of all genes in grl
ilc3_stressed_count <- length(reduce(resize(sort(ilc3_stressedgr), 1000000)))




lti_stressed <- fread(paste0(results_dir,"proportions/LTi_stressed_vs_non_qtl.csv.gz"), 
                      data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.999)) %>%
  dplyr::select(marker, chr, betas, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)

lti_stressedgr <- makeGRangesFromDataFrame(lti_stressed,
                                           seqnames="chr",
                                           start.field="pos",
                                           end.field="pos",
                                           strand = "strand",
                                           keep.extra.columns = TRUE,
                                           na.rm = TRUE)

# Counting number of QTL
## length of all genes in grl
lti_stressed_count <- length(reduce(resize(sort(lti_stressedgr), 1000000)))




## TODO: update beyond this point ####
## topic qtl counts ####

topics <- fread(paste0(results_dir, "topics/qtl-topic-lods.csv.gz"), 
                data.table = FALSE) %>%
  dplyr::select(-V1) %>%
  left_join(ccre) %>%
  dplyr::select(1:23) %>%
  pivot_longer(-c(marker, chr, pos)) %>%
  mutate(strand = "+") %>%
  filter(!is.na(pos), !is.na(chr), value > 6)



topicqtlgrl <- makeGRangesListFromDataFrame(topics,
                                            split.field = "name",
                                            names.field = "marker",
                                            seqnames = "chr",
                                            start.field = "pos",
                                            end.field = "pos",
                                            strand = "strand",
                                            na.rm = TRUE,
                                            keep.extra.columns = TRUE)

## length of all genes in grl
topicqtl_list <- lapply(reduce(resize(sort(topicqtlgrl),10000)), length)


## cytokine qtl counts

cytokines <- fread(paste0(results_dir,"cytokines/qtl-cytokines-steady-lods.csv.gz"), 
                   data.table = FALSE) %>%
  dplyr::select(-V1) %>%
  left_join(ccre) %>%
  dplyr::select(1:15) %>%
  pivot_longer(-c(marker, chr, pos)) %>%
  mutate(strand = "+") %>%
  filter(!is.na(pos), !is.na(chr), value > 6)



cytokinesqtlgrl <- makeGRangesListFromDataFrame(cytokines,
                                            split.field = "name",
                                            names.field = "marker",
                                            seqnames = "chr",
                                            start.field = "pos",
                                            end.field = "pos",
                                            strand = "strand",
                                            na.rm = TRUE,
                                            keep.extra.columns = TRUE)

## length of all genes in grl
cytokineqtl_list <- lapply(reduce(resize(sort(cytokinesqtlgrl),50000)), length)




# outdf <- bind_rows(list(ilc1_eqtl_count, ilc1_ilc2_count, ilc1_ilc3_count,
#                         ilc1_lti_count,  ilc2_eqtl_count, ilc2_ilc3_count, 
#                         ilc2_lti_count, ilc3_eqtl_count,
#                         ilc3_lti_count, lti_eqtl_count, topics, cytokines))
# write.csv(outdf, "qtl-counts.csv", row.names = F)


## Counting ##########

## ILC1
ilc1_eqtl_list <- lapply(reduce(resize(sort(ilc1eqtlgrl),10000)), length)

## ILC2
ilc2_eqtl_list <- lapply(reduce(resize(sort(ilc2eqtlgrl),10000)), length)

## ILC3
ilc3_eqtl_list <- lapply(reduce(resize(sort(ilc3eqtlgrl),10000)), length)

## LTi
lti_eqtl_list <- lapply(reduce(resize(sort(ltieqtlgrl),10000)), length)

## ILC1 vs ILC2
ilc1_ilc2_count <- length(reduce(resize(sort(ilc1_ilc2gr), 10000)))

## ILC1 vs ILC3
ilc1_ilc3_count <- length(reduce(resize(sort(ilc1_ilc3gr), 10000)))

## ILC1 vs LTi
ilc1_ilc3_count <- length(reduce(resize(sort(ilc1_ilc3gr), 10000)))

## ILC2 vs ILC3
ilc2_ilc3_count <- length(reduce(resize(sort(ilc2_ilc3gr), 10000)))

## ILC2 vs LTi
ilc2_lti_count <- length(reduce(resize(sort(ilc2_ltigr), 10000)))

## ILC3 vs LTi
ilc3_lti_count <- length(reduce(resize(sort(ilc3_ltigr), 10000)))

## ILC3 Stressed
ilc3_stressed_count <- length(reduce(resize(sort(ilc3_stressedgr), 10000)))

## LTi Stressed
lti_stressed_count <- length(reduce(resize(sort(lti_stressedgr), 10000)))

## topic
topicqtl_list <- lapply(reduce(resize(sort(topicqtlgrl),10000)), length)

## cytokine
cytokineqtl_list <- lapply(reduce(resize(sort(cytokinesqtlgrl),10000)), length)




