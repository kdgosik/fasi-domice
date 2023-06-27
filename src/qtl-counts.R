library(data.table)
library(rtracklayer)
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

## Gene Data
# ensembl <- as.data.frame(readGFF(paste0(data_dir, "references/Mus_musculus.GRCm38.102.gtf"))) %>%
#   filter(seqid %in% c(as.character(1:19), "X"), gene_biotype == "protein_coding", type == "gene") %>%
#   mutate(chr = paste0("chr", seqid))

ensembl <- readGFF(paste0(data_dir, "references/Mus_musculus.GRCm38.102.gtf")) %>%
  filter(seqid %in% c(as.character(1:19), "X"), 
         gene_biotype == "protein_coding", 
         type == "gene") %>%
  dplyr::select(seqid, start, end, strand, gene_name)


mm10 <- makeGRangesFromDataFrame(ensembl, keep.extra.columns = T)

ccregr <- ccre %>%
  dplyr::select(chr, pos, marker) %>%
  filter(!is.na(chr)) %>%
  makeGRangesFromDataFrame(.,
                           seqnames.field = "chr",
                           start.field = "pos",
                           end.field = "pos",
                           na.rm=TRUE, 
                           keep.extra.columns = T)



## eQTL counts by gene ####
lod_cutoff <- 6


## ILC1 ####

# Step 1:  
# Step 2:
# Step 3:
# Step 4:
# Step 5:
# Step 6:


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
                cis_effect, value_adj) %>%
  mutate(gene = paste0("ILC1_", gene))

ilc1eqtlgrl <- makeGRangesListFromDataFrame(ilc1_eqtl_loci_by_gene,
                                            split.field = "gene",
                                            names.field = "marker",
                                            seqnames = "marker_chr",
                                            start.field = "pos",
                                            end.field = "pos",
                                            strand = "strand",
                                            keep.extra.columns = TRUE)

## length of all genes in grl
ilc1_eqtl_list <- lapply(reduce(resize(sort(ilc1eqtlgrl), 500000)), length)

## count monogenic vs polygenic
ilc1_polygenic_counts <- lapply(seq_along(ilc1eqtlgrl), function(i) { 
  sum(GenomicRanges::countOverlaps(reduce(resize(sort(ilc1eqtlgrl), 500000)),
                                   reduce(resize(sort(ilc1eqtlgrl), 500000))[[i]])>0) 
})

## e.g.
# subsetByOverlaps(reduce(resize(sort(ilc1eqtlgrl), 500000)),
#                  +                  reduce(resize(sort(ilc1eqtlgrl), 500000))[[203]])

ilc1_polygenic_df <- lapply(seq_along(ilc1eqtlgrl), function(i) { 
  subsetByOverlaps(reduce(resize(sort(ilc1eqtlgrl), 500000)),
                                   reduce(resize(sort(ilc1eqtlgrl), 500000))[[i]]) %>%
    as.data.frame() %>%
    distinct(group_name, .keep_all = TRUE)
}) %>%
  bind_rows(.,.id = "loci_id") %>%
  mutate(loci_id = paste0('loci_', str_pad(loci_id, width=3, side = 'left', pad='0'))) %>%
  dplyr::rename(eGene = group_name,
                loci_chr = seqnames,
                loci_start = start,
                loci_end = end) %>%
  dplyr::select(-width) %>%
  group_by(loci_id) %>%
  mutate(eGene_count = n_distinct(eGene)) %>%
  ungroup()

## making GRange List of loci_id to find overlaps
ilc1_polygenic_dfgrl <- makeGRangesListFromDataFrame(ilc1_polygenic_df,
                                               split.field = "loci_id",
                                               names.field = "eGene",
                                               seqnames = "loci_chr",
                                               start.field = "loci_start",
                                               end.field = "loci_end",
                                               strand = "strand",
                                               keep.extra.columns = TRUE)

## creating loci gene annotations
ilc1_outdf <- sapply(names(ilc1_polygenic_dfgrl), function(i) { 
  subsetByOverlaps(mm10, ilc1_polygenic_dfgrl[[i]]) %>% 
    as.data.frame() 
  }, USE.NAMES = TRUE, simplify = FALSE ) %>% 
  bind_rows(., .id='loci_id') %>%
  dplyr::select(loci_id,
                eQTL_loci_gene_chr = seqnames,
                eQTL_loci_gene_start = start,
                eQTL_loci_gene_end = end,
                eQTL_loci_gene = gene_name)


ilc1_eqtl_loci_by_gene_outdf <- full_join(ilc1_polygenic_df, ilc1_outdf) %>% 
  dplyr::select(-group) %>%
  mutate(cell_type = "ILC1")
write.csv(ilc1_eqtl_loci_by_gene_outdf, paste0(results_dir, "eqtl/qtl-loci-by-gene-ILC1.csv"))


ilc1_eqtl_polygeneicdf <- ilc1_eqtl_loci_by_gene_outdf %>%
  filter(eGene_count > 1)
write.csv(ilc1_eqtl_polygeneicdf, paste0(results_dir, "eqtl/qtl-loci-by-gene-ILC1_polygenic_only.csv"))

qtl_loci_by_gene_ILC1 <- read_csv("/workspace/fasi-domice/results/eqtl/qtl-loci-by-gene-ILC1.csv")

loci_gene_expressed_ilc1 <- qtl_loci_by_gene_ILC1 %>% 
  filter(eQTL_loci_gene %in% vars$index[vars$ilc1_expressed==1]) %>% 
  pull(loci_id) %>% 
  n_distinct()

loci_total_ilc1 <- qtl_loci_by_gene_ILC1 %>% 
  # filter(eQTL_loci_gene %in% vars$index[vars$ilc1_expressed==1]) %>% 
  pull(loci_id) %>% 
  n_distinct()


loci_gene_expressed_ilc1 / loci_total_ilc1


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
  filter(value_adj > 10,
         strand = "+") %>%
  dplyr::select(marker, gene, value, marker_chr, pos, 
                id, gene_chr, gene_start, gene_end, cis_effect, value_adj) %>%
  mutate(gene = paste0("ILC2_", gene))

ilc2eqtlgrl <- makeGRangesListFromDataFrame(ilc2_eqtl_loci_by_gene,
                                            split.field = "gene",
                                            names.field = "marker",
                                            seqnames="marker_chr",
                                            start.field="pos",
                                            end.field="pos",
                                            strand = "strand",
                                            keep.extra.columns = TRUE)

## length of all genes in grl
ilc2_eqtl_list <- lapply(reduce(resize(sort(ilc2eqtlgrl),500000)), length)


## count monogenic vs polygenic
ilc2_polygenic_counts <- lapply(seq_along(ilc2eqtlgrl), function(i) { 
  sum(GenomicRanges::countOverlaps(reduce(resize(sort(ilc2eqtlgrl), 500000)),
                                   reduce(resize(sort(ilc2eqtlgrl), 500000))[[i]])>0) 
})



ilc2_polygenic_df <- lapply(seq_along(ilc2eqtlgrl), function(i) { 
  subsetByOverlaps(reduce(resize(sort(ilc2eqtlgrl), 500000)),
                   reduce(resize(sort(ilc2eqtlgrl), 500000))[[i]]) %>%
    as.data.frame() %>%
    distinct(group_name, .keep_all = TRUE)
}) %>%
  bind_rows(.,.id = "loci_id") %>%
  mutate(loci_id = paste0('loci_', str_pad(loci_id, width=3, side = 'left', pad='0'))) %>%
  dplyr::rename(eGene = group_name,
                loci_chr = seqnames,
                loci_start = start,
                loci_end = end) %>%
  dplyr::select(-width) %>%
  group_by(loci_id) %>%
  mutate(eGene_count = n_distinct(eGene)) %>%
  ungroup()

## making GRange List of loci_id to find overlaps
ilc2_polygenic_dfgrl <- makeGRangesListFromDataFrame(ilc2_polygenic_df,
                                                split.field = "loci_id",
                                                names.field = "eGene",
                                                seqnames = "loci_chr",
                                                start.field = "loci_start",
                                                end.field = "loci_end",
                                                strand = "strand",
                                                keep.extra.columns = TRUE)

## creating loci gene annotations
ilc2_outdf <- sapply(names(ilc2_polygenic_dfgrl), function(i) { 
  subsetByOverlaps(mm10, ilc2_polygenic_dfgrl[[i]]) %>% 
    as.data.frame() 
}, USE.NAMES = TRUE, simplify = FALSE ) %>% 
  bind_rows(., .id='loci_id') %>%
  dplyr::select(loci_id,
                eQTL_loci_gene_chr = seqnames,
                eQTL_loci_gene_start = start,
                eQTL_loci_gene_end = end,
                eQTL_loci_gene = gene_name)


ilc2_eqtl_loci_by_gene_outdf <- full_join(ilc2_polygenic_df, ilc2_outdf) %>% 
  dplyr::select(-group) %>%
  mutate(cell_type = "ILC2")
write.csv(ilc2_eqtl_loci_by_gene_outdf, paste0(results_dir, "eqtl/qtl-loci-by-gene-ILC2.csv"))


ilc2_eqtl_polygeneicdf <- ilc2_eqtl_loci_by_gene_outdf %>%
  filter(eGene_count > 1)
write.csv(ilc2_eqtl_polygeneicdf, paste0(results_dir, "eqtl/qtl-loci-by-gene-ILC2_polygenic_only.csv"))


## count loci genes expressed in cell type
qtl_loci_by_gene_ILC2 <- read_csv("/workspace/fasi-domice/results/eqtl/qtl-loci-by-gene-ILC2.csv")

loci_gene_expressed_ilc2 <- qtl_loci_by_gene_ILC2 %>% 
  filter(eQTL_loci_gene %in% vars$index[vars$ilc2_expressed==1]) %>% 
  pull(loci_id) %>% 
  n_distinct()

loci_total_ilc2 <- qtl_loci_by_gene_ILC2 %>% 
  # filter(eQTL_loci_gene %in% vars$index[vars$ilc1_expressed==1]) %>% 
  pull(loci_id) %>% 
  n_distinct()
loci_gene_expressed_ilc2 / loci_total_ilc2



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
  filter(value_adj > 10,
         strand = "+") %>%
  dplyr::select(marker, gene, value, marker_chr, pos, 
                id, gene_chr, gene_start, gene_end, cis_effect, value_adj) %>%
  mutate(gene = paste0("ILC3_", gene))

ilc3eqtlgrl <- makeGRangesListFromDataFrame(ilc3_eqtl_loci_by_gene,
                                            split.field = "gene",
                                            names.field = "marker",
                                            seqnames="marker_chr",
                                            start.field="pos",
                                            end.field="pos",
                                            strand = "strand",
                                            keep.extra.columns = TRUE)

## length of all genes in grl
ilc3_eqtl_list <- lapply(reduce(resize(sort(ilc3eqtlgrl), 500000)), length)


## count monogenic vs polygenic
ilc3_polygenic_counts <- lapply(seq_along(ilc3eqtlgrl), function(i) { 
  sum(GenomicRanges::countOverlaps(reduce(resize(sort(ilc3eqtlgrl), 500000)),
                                   reduce(resize(sort(ilc3eqtlgrl), 500000))[[i]])>0) 
})



ilc3_polygenic_df <- lapply(seq_along(ilc3eqtlgrl), function(i) { 
  subsetByOverlaps(reduce(resize(sort(ilc3eqtlgrl), 500000)),
                   reduce(resize(sort(ilc3eqtlgrl), 500000))[[i]]) %>%
    as.data.frame() %>%
    distinct(group_name, .keep_all = TRUE)
}) %>%
  bind_rows(.,.id = "loci_id") %>%
  mutate(loci_id = paste0('loci_', str_pad(loci_id, width=3, side = 'left', pad='0'))) %>%
  dplyr::rename(eGene = group_name,
                loci_chr = seqnames,
                loci_start = start,
                loci_end = end) %>%
  dplyr::select(-width) %>%
  group_by(loci_id) %>%
  mutate(eGene_count = n_distinct(eGene)) %>%
  ungroup()

## making GRange List of loci_id to find overlaps
ilc3_polygenic_dfgrl <- makeGRangesListFromDataFrame(ilc3_polygenic_df,
                                                     split.field = "loci_id",
                                                     names.field = "eGene",
                                                     seqnames = "loci_chr",
                                                     start.field = "loci_start",
                                                     end.field = "loci_end",
                                                     strand = "strand",
                                                     keep.extra.columns = TRUE)

## creating loci gene annotations
ilc3_outdf <- sapply(names(ilc3_polygenic_dfgrl), function(i) { 
  subsetByOverlaps(mm10, ilc3_polygenic_dfgrl[[i]]) %>% 
    as.data.frame() 
}, USE.NAMES = TRUE, simplify = FALSE ) %>% 
  bind_rows(., .id='loci_id') %>%
  dplyr::select(loci_id,
                eQTL_loci_gene_chr = seqnames,
                eQTL_loci_gene_start = start,
                eQTL_loci_gene_end = end,
                eQTL_loci_gene = gene_name)


ilc3_eqtl_loci_by_gene_outdf <- full_join(ilc3_polygenic_df, ilc3_outdf) %>% 
  dplyr::select(-group) %>%
  mutate(cell_type = "ILC3")
write.csv(ilc3_eqtl_loci_by_gene_outdf, paste0(results_dir, "eqtl/qtl-loci-by-gene-ILC3.csv"))


ilc3_eqtl_polygeneicdf <- ilc3_eqtl_loci_by_gene_outdf %>%
  filter(eGene_count > 1)
write.csv(ilc3_eqtl_polygeneicdf, paste0(results_dir, "eqtl/qtl-loci-by-gene-ILC3_polygenic_only.csv"))

## count loci genes expressed in cell type
qtl_loci_by_gene_ILC3 <- read_csv("/workspace/fasi-domice/results/eqtl/qtl-loci-by-gene-ILC3.csv")

loci_gene_expressed_ilc3 <- qtl_loci_by_gene_ILC3 %>% 
  filter(eQTL_loci_gene %in% vars$index[vars$ilc3_expressed==1]) %>% 
  pull(loci_id) %>% 
  n_distinct()

loci_total_ilc3 <- qtl_loci_by_gene_ILC3 %>% 
  # filter(eQTL_loci_gene %in% vars$index[vars$ilc1_expressed==1]) %>% 
  pull(loci_id) %>% 
  n_distinct()
loci_gene_expressed_ilc3 / loci_total_ilc3



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
  filter(value_adj > 10,
         strand = "+") %>%
  dplyr::select(marker, gene, value, marker_chr, pos, 
                id, gene_chr, gene_start, gene_end, cis_effect, value_adj) %>%
  mutate(gene = paste0("LTi_", gene))

ltieqtlgrl <- makeGRangesListFromDataFrame(lti_eqtl_loci_by_gene,
                                           split.field = "gene",
                                           names.field = "marker",
                                           seqnames="marker_chr",
                                           start.field="pos",
                                           end.field="pos",
                                           strand = "strand",
                                           keep.extra.columns = TRUE)

## length of all genes in grl
lti_eqtl_list <- lapply(reduce(resize(sort(ltieqtlgrl), 500000)), length)


## count monogenic vs polygenic
lti_polygenic_counts <- lapply(seq_along(ltieqtlgrl), function(i) { 
  sum(GenomicRanges::countOverlaps(reduce(resize(sort(ltieqtlgrl), 500000)),
                                   reduce(resize(sort(ltieqtlgrl), 500000))[[i]])>0) 
})



lti_polygenic_df <- lapply(seq_along(ltieqtlgrl), function(i) { 
  subsetByOverlaps(reduce(resize(sort(ltieqtlgrl), 500000)),
                   reduce(resize(sort(ltieqtlgrl), 500000))[[i]]) %>%
    as.data.frame() %>%
    distinct(group_name, .keep_all = TRUE)
}) %>%
  bind_rows(.,.id = "loci_id") %>%
  mutate(loci_id = paste0('loci_', str_pad(loci_id, width=3, side = 'left', pad='0'))) %>%
  dplyr::rename(eGene = group_name,
                loci_chr = seqnames,
                loci_start = start,
                loci_end = end) %>%
  dplyr::select(-width) %>%
  group_by(loci_id) %>%
  mutate(eGene_count = n_distinct(eGene)) %>%
  ungroup()

## making GRange List of loci_id to find overlaps
lti_polygenic_dfgrl <- makeGRangesListFromDataFrame(lti_polygenic_df,
                                                     split.field = "loci_id",
                                                     names.field = "eGene",
                                                     seqnames = "loci_chr",
                                                     start.field = "loci_start",
                                                     end.field = "loci_end",
                                                     strand = "strand",
                                                     keep.extra.columns = TRUE)

## creating loci gene annotations
lti_outdf <- sapply(names(lti_polygenic_dfgrl), function(i) { 
  subsetByOverlaps(mm10, lti_polygenic_dfgrl[[i]]) %>% 
    as.data.frame() 
}, USE.NAMES = TRUE, simplify = FALSE ) %>% 
  bind_rows(., .id='loci_id') %>%
  dplyr::select(loci_id,
                eQTL_loci_gene_chr = seqnames,
                eQTL_loci_gene_start = start,
                eQTL_loci_gene_end = end,
                eQTL_loci_gene = gene_name)


lti_eqtl_loci_by_gene_outdf <- full_join(lti_polygenic_df, lti_outdf) %>% 
  dplyr::select(-group) %>%
  mutate(cell_type = "LTi-like")
write.csv(lti_eqtl_loci_by_gene_outdf, paste0(results_dir, "eqtl/qtl-loci-by-gene-LTi.csv"))


lti_eqtl_polygeneicdf <- lti_eqtl_loci_by_gene_outdf %>%
  filter(eGene_count > 1)
write.csv(lti_eqtl_polygeneicdf, paste0(results_dir, "eqtl/qtl-loci-by-gene-LTi_polygenic_only.csv"))

## count loci genes expressed in cell type
qtl_loci_by_gene_LTi <- read_csv("/workspace/fasi-domice/results/eqtl/qtl-loci-by-gene-LTi.csv")

loci_gene_expressed_lti <- qtl_loci_by_gene_LTi %>% 
  filter(eQTL_loci_gene %in% vars$index[vars$lti_expressed==1]) %>% 
  pull(loci_id) %>% 
  n_distinct()

loci_total_lti <- qtl_loci_by_gene_LTi %>% 
  # filter(eQTL_loci_gene %in% vars$index[vars$ilc1_expressed==1]) %>% 
  pull(loci_id) %>% 
  n_distinct()
loci_gene_expressed_lti / loci_total_lti


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
ilc1_ilc2_count <- length(reduce(resize(sort(ilc1_ilc2gr), 500000)))



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
ilc1_ilc3_count <- length(reduce(resize(sort(ilc1_ilc3gr), 500000)))




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



## within proportion ####
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

prop_count_list <- list("ILC1_ILC2"=ilc1_ilc2, 
                        "ILC1_ILC3"=ilc1_ilc3,
                        "ILC1_LTi"=ilc1_lti,
                        "ILC2_ILC3"=ilc2_ilc3,
                        "ILC2_LTi"=ilc2_lti,
                        "ILC3_LTi"=ilc3_lti,
                        "ILC3_activated"=ilc3_stressed,
                        "LTi_activated"=lti_stressed) %>%
  bind_rows(.,.id = "trait")


prop_count_list <- list("ILC1_ILC2"=ilc1_ilc2_count, 
                        "ILC1_ILC3"=ilc1_ilc3_count,
                        "ILC1_LTi"=ilc1_lti_count,
                        "ILC2_ILC3"=ilc2_ilc3_count,
                        "ILC2_LTi"=ilc2_lti_count,
                        "ILC3_LTi"=ilc3_lti_count,
                        "ILC3_activated"=ilc3_stressed_count,
                        "LTi_activated"=lti_stressed_count)

prop_gr_list <- list("ILC1_ILC2"=ilc1_ilc2gr, 
                  "ILC1_ILC3"=ilc1_ilc3gr,
                  "ILC1_LTi"=ilc1_ltigr,
                  "ILC2_ILC3"=ilc2_ilc3gr,
                  "ILC2_LTi"=ilc2_ltigr,
                  "ILC3_LTi"=ilc3_ltigr,
                  "ILC3_activated"=ilc3_stressedgr,
                  "LTi_activated"=lti_stressedgr)

prop_grl <- as(prop_list, "GRangesList")


## length of all genes in grl
prop_grl_list <- lapply(reduce(resize(sort(prop_grl), 500000)), length)


## count monogenic vs polygenic
prop_polygenic_counts <- lapply(seq_along(prop_grl), function(i) { 
  sum(GenomicRanges::countOverlaps(reduce(resize(sort(prop_grl), 500000)),
                                   reduce(resize(sort(prop_grl), 500000))[[i]])>0) 
})



prop_polygenic_df <- lapply(seq_along(prop_grl), function(i) { 
  subsetByOverlaps(reduce(resize(sort(prop_grl), 500000)),
                   reduce(resize(sort(prop_grl), 500000))[[i]]) %>%
    as.data.frame() %>%
    distinct(group_name, .keep_all = TRUE)
}) %>%
  bind_rows(.,.id = "loci_id") %>%
  mutate(loci_id = paste0('loci_', str_pad(loci_id, width=3, side = 'left', pad='0'))) %>%
  dplyr::rename(proportion = group_name,
                loci_chr = seqnames,
                loci_start = start,
                loci_end = end) %>%
  dplyr::select(-width) %>%
  group_by(loci_id) %>%
  mutate(proportion_count = n_distinct(proportion)) %>%
  ungroup()

## making GRange List of loci_id to find overlaps
prop_polygenic_dfgrl <- makeGRangesListFromDataFrame(prop_polygenic_df,
                                                     split.field = "loci_id",
                                                     names.field = "proportion",
                                                     seqnames = "loci_chr",
                                                     start.field = "loci_start",
                                                     end.field = "loci_end",
                                                     strand = "strand",
                                                     keep.extra.columns = TRUE)

## creating loci gene annotations
prop_outdf <- sapply(names(prop_polygenic_dfgrl), function(i) { 
  subsetByOverlaps(mm10, prop_polygenic_dfgrl[[i]]) %>% 
    as.data.frame() 
}, USE.NAMES = TRUE, simplify = FALSE ) %>% 
  bind_rows(., .id='loci_id') %>%
  dplyr::select(loci_id,
                propQTL_loci_gene_chr = seqnames,
                propQTL_loci_gene_start = start,
                propQTL_loci_gene_end = end,
                propQTL_loci_gene = gene_name)


prop_loci_by_gene_outdf <- full_join(prop_polygenic_df, prop_outdf) %>% 
  dplyr::select(-group)


write.csv(prop_loci_by_gene_outdf, paste0(results_dir, "proportions/qtl-loci-by-gene-proportions.csv"))


prop_polygeneicdf <- prop_loci_by_gene_outdf %>%
  filter(proportion_count > 1)
write.csv(prop_polygeneicdf, paste0(results_dir, "proportions/qtl-loci-by-gene-proportions_polygenic_only.csv"))


prop_qtl_count <- prop_loci_by_gene_outdf %>%
  filter(propQTL_loci_gene %in% vars$index[vars$ilc1_expressed==1] |
         propQTL_loci_gene %in% vars$index[vars$ilc2_expressed==1] |
         propQTL_loci_gene %in% vars$index[vars$ilc3_expressed==1] |
         propQTL_loci_gene %in% vars$index[vars$lti_expressed==1]) %>%
  pull(loci_id) %>%
  n_distinct()

prop_qtl_total_count <- prop_loci_by_gene_outdf %>%
  # filter(propQTL_loci_gene %in% vars$index[vars$ilc1_expressed==1] |
  #          propQTL_loci_gene %in% vars$index[vars$ilc2_expressed==1] |
  #          propQTL_loci_gene %in% vars$index[vars$ilc3_expressed==1] |
  #          propQTL_loci_gene %in% vars$index[vars$lti_expressed==1]) %>%
  pull(loci_id) %>%
  n_distinct()




## topic qtl counts ####

topics <- fread(paste0(results_dir, "topics/qtl-topic-lods.csv.gz"), 
                data.table = FALSE) %>%
  dplyr::select(-V1) %>%
  left_join(ccre) %>%
  dplyr::select(1:23) %>%
  pivot_longer(-c(marker, chr, pos)) %>%
  mutate(strand = "+") %>%
  filter(!is.na(pos), !is.na(chr), value > 6) %>%
  dplyr::rename(trait = name)



topicqtlgrl <- makeGRangesListFromDataFrame(topics,
                                            split.field = "trait",
                                            names.field = "marker",
                                            seqnames = "chr",
                                            start.field = "pos",
                                            end.field = "pos",
                                            strand = "strand",
                                            na.rm = TRUE,
                                            keep.extra.columns = TRUE)

## length of all genes in grl
topicqtl_list <- lapply(reduce(resize(sort(topicqtlgrl),500000)), length)


## count monogenic vs polygenic
topic_polygenic_counts <- lapply(seq_along(topicqtlgrl), function(i) { 
  sum(GenomicRanges::countOverlaps(reduce(resize(sort(topicqtlgrl), 500000)),
                                   reduce(resize(sort(topicqtlgrl), 500000))[[i]])>0) 
})



topic_polygenic_df <- lapply(seq_along(topicqtlgrl), function(i) { 
  
  fovdf <- findOverlaps(reduce(resize(sort(topicqtlgrl), 500000)),
                        reduce(resize(sort(topicqtlgrl), 500000))[[i]]) %>%
    as.data.frame()
  
  keep_ids <- sort(unique(c(fovdf$queryHits, fovdf$subjectHits)))
  
  ovdf <- subsetByOverlaps(reduce(resize(sort(topicqtlgrl), 500000)),
                           reduce(resize(sort(topicqtlgrl), 500000))[[i]]) %>%
    as.data.frame()
  
  outdf <- ovdf[keep_ids, ] %>% mutate(loci = fovdf$subjectHits)
  outdf
  
}) %>%
  bind_rows(.,.id = "loci_id") %>%
  mutate(loci_id = paste0('loci_', str_pad(loci_id, width=3, side = 'left', pad='0'), loci)) %>%
  dplyr::rename(topic = group_name,
                loci_chr = seqnames,
                loci_start = start,
                loci_end = end) %>%
  dplyr::select(-width) %>%
  group_by(loci_id) %>%
  mutate(topic_count = n_distinct(topic)) %>%
  ungroup()

## making GRange List of loci_id to find overlaps
topic_polygenic_dfgrl <- makeGRangesListFromDataFrame(topic_polygenic_df,
                                                     split.field = "loci_id",
                                                     names.field = "topic",
                                                     seqnames = "loci_chr",
                                                     start.field = "loci_start",
                                                     end.field = "loci_end",
                                                     strand = "strand",
                                                     keep.extra.columns = TRUE)

## creating loci gene annotations
topic_outdf <- sapply(names(topic_polygenic_dfgrl), function(i) { 
  subsetByOverlaps(mm10, topic_polygenic_dfgrl[[i]]) %>% 
    as.data.frame() 
}, USE.NAMES = TRUE, simplify = FALSE ) %>% 
  bind_rows(., .id='loci_id') %>%
  dplyr::select(loci_id,
                topicQTL_loci_gene_chr = seqnames,
                topicQTL_loci_gene_start = start,
                topicQTL_loci_gene_end = end,
                topicQTL_loci_gene = gene_name)


topic_loci_by_gene_outdf <- full_join(topic_polygenic_df, topic_outdf) %>% 
  dplyr::select(-group)


write.csv(topic_loci_by_gene_outdf, paste0(results_dir, "topics/qtl-loci-by-gene-topics.csv"))


topic_polygeneicdf <- topic_loci_by_gene_outdf %>%
  filter(topic_count > 1)
write.csv(topic_polygeneicdf, paste0(results_dir, "topics/qtl-loci-by-gene-topics_polygenic_only.csv"))



## cytokine qtl counts

cytokines <- fread(paste0(results_dir,"cytokines/qtl-cytokines-steady-lods.csv.gz"), 
                   data.table = FALSE) %>%
  dplyr::select(-V1) %>%
  left_join(ccre) %>%
  dplyr::select(1:15) %>%
  pivot_longer(-c(marker, chr, pos)) %>%
  mutate(strand = "+") %>%
  filter(!is.na(pos), !is.na(chr), value > 6) %>%
  dplyr::rename(trait = name)



cytokineqtlgrl <- makeGRangesListFromDataFrame(cytokines,
                                            split.field = "name",
                                            names.field = "marker",
                                            seqnames = "chr",
                                            start.field = "pos",
                                            end.field = "pos",
                                            strand = "strand",
                                            na.rm = TRUE,
                                            keep.extra.columns = TRUE)

## length of all genes in grl
cytokineqtl_list <- lapply(reduce(resize(sort(cytokineqtlgrl), 5000000)), length)


## count monogenic vs polygenic
cytokine_polygenic_counts <- lapply(seq_along(cytokineqtlgrl), function(i) { 
  sum(GenomicRanges::countOverlaps(reduce(resize(sort(cytokineqtlgrl), 500000)),
                                   reduce(resize(sort(cytokineqtlgrl), 500000))[[i]])>0) 
})



cytokine_polygenic_df <- lapply(seq_along(cytokineqtlgrl), function(i) { 
  
  fovdf <- GenomicRanges::findOverlaps(reduce(resize(sort(cytokineqtlgrl), 500000)),
                        reduce(resize(sort(cytokineqtlgrl), 500000))[[i]]) %>%
    as.data.frame()
  
  keep_ids <- sort(unique(c(fovdf$queryHits, fovdf$subjectHits)))
  
  ovdf <- subsetByOverlaps(reduce(resize(sort(cytokineqtlgrl), 500000)),
                           reduce(resize(sort(cytokineqtlgrl), 500000))[[i]]) %>%
    as.data.frame()
  
  outdf <- ovdf[keep_ids, ] %>% mutate(loci = fovdf$subjectHits)
  outdf
  
}) %>%
  bind_rows(.,.id = "loci_id") %>%
  mutate(loci_id = paste0('loci_', str_pad(loci_id, width=3, side = 'left', pad='0'), loci)) %>%
  # mutate(loci_id = paste0('loci_', str_pad(loci_id, width=3, side = 'left', pad='0'))) %>%
  dplyr::rename(cytokine = group_name,
                loci_chr = seqnames,
                loci_start = start,
                loci_end = end) %>%
  dplyr::select(-width) %>%
  group_by(loci_id) %>%
  mutate(cytokine_count = n_distinct(cytokine)) %>%
  ungroup()

## making GRange List of loci_id to find overlaps
cytokine_polygenic_dfgrl <- makeGRangesListFromDataFrame(cytokine_polygenic_df,
                                                      split.field = "loci_id",
                                                      names.field = "cytokine",
                                                      seqnames = "loci_chr",
                                                      start.field = "loci_start",
                                                      end.field = "loci_end",
                                                      strand = "strand",
                                                      keep.extra.columns = TRUE)

## creating loci gene annotations
cytokine_outdf <- sapply(names(cytokine_polygenic_dfgrl), function(i) { 
  subsetByOverlaps(mm10, cytokine_polygenic_dfgrl[[i]]) %>% 
    as.data.frame() 
}, USE.NAMES = TRUE, simplify = FALSE ) %>% 
  bind_rows(., .id='loci_id') %>%
  dplyr::select(loci_id,
                cytokineQTL_loci_gene_chr = seqnames,
                cytokineQTL_loci_gene_start = start,
                cytokineQTL_loci_gene_end = end,
                cytokineQTL_loci_gene = gene_name)


cytokine_loci_by_gene_outdf <- full_join(cytokine_polygenic_df, cytokine_outdf) %>% 
  dplyr::select(-group)


write.csv(cytokine_loci_by_gene_outdf, paste0(results_dir, "cytokines/qtl-loci-by-gene-cytokines.csv"))


cytokine_polygeneicdf <- cytokine_loci_by_gene_outdf %>%
  filter(cytokine_count > 1)
write.csv(cytokine_polygeneicdf, paste0(results_dir, "cytokines/qtl-loci-by-gene-cytokines_polygenic_only.csv"))


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



## Combining to Dataset ##########

list(ilc1_eqtl = ilc1_eqtl_loci_by_gene %>% 
       dplyr::select(trait = gene,
                     marker,
                     chr = marker_chr,
                     strand,
                     value),
     ilc2_eqtl = ilc2_eqtl_loci_by_gene %>% 
       dplyr::select(trait = gene,
                     marker,
                     chr = marker_chr,
                     strand,
                     value),
     ilc3_eqtl = ilc3_eqtl_loci_by_gene %>% 
       dplyr::select(trait = gene,
                     marker,
                     chr = marker_chr,
                     strand,
                     value),
     lti_eqtl = lti_eqtl_loci_by_gene %>% 
       dplyr::select(trait = gene,
                     marker,
                     chr = marker_chr,
                     strand,
                     value))



