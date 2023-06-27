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
  dplyr::select(seqid, start, end, strand, gene_name, gene_id)


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
  mutate(value_adj = ifelse(cis_effect==1, value + 6, value),
         strand = "+") %>%
  filter(value_adj > 10) %>%
  dplyr::select(marker, gene, value, marker_chr, pos, strand,
                id, gene_chr, gene_start, gene_end, cis_effect, value_adj) %>%
  mutate(gene = paste0("ILC2_", gene))


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
  mutate(value_adj = ifelse(cis_effect==1, value + 6, value),
         strand = "+") %>%
  filter(value_adj > 10) %>%
  dplyr::select(marker, gene, value, marker_chr, pos, strand,
                id, gene_chr, gene_start, gene_end, cis_effect, value_adj) %>%
  mutate(gene = paste0("ILC3_", gene))



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
  mutate(value_adj = ifelse(cis_effect==1, value + 6, value),
         strand = "+") %>%
  filter(value_adj > 10) %>%
  dplyr::select(marker, gene, value, marker_chr, pos, strand,
                id, gene_chr, gene_start, gene_end, cis_effect, value_adj) %>%
  mutate(gene = paste0("LTi_", gene))

## proportion qtl counts ####

ilc1_ilc2 <- fread(paste0(results_dir,"proportions/gwas-ilc1-ilc2-results.csv.gz"), 
                   data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.995)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand, xpos)


ilc1_ilc3 <- fread(paste0(results_dir,"proportions/gwas-ilc1-ilc3-results.csv.gz"), 
                   data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.995)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)


ilc1_lti <- fread(paste0(results_dir,"proportions/gwas-ilc1-lti-results.csv.gz"), 
                   data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.995)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)



ilc2_ilc3 <- fread(paste0(results_dir,"proportions/gwas-ilc2-ilc3-results.csv.gz"), 
                  data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.995)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)

# quantile(ilc2_lti$lods, 0.999)
ilc2_lti <- fread(paste0(results_dir,"proportions/gwas-ilc2-lti-results.csv.gz"), 
                  data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.995)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)

ilc3_lti <- fread(paste0(results_dir,"proportions/gwas-ilc3-lti-results.csv.gz"), 
                  data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.995)) %>%
  dplyr::select(marker, chr, betas, se, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)


## within proportion ####
ilc3_stressed <- fread(paste0(results_dir,"proportions/ILC3_stressed_vs_non_qtl.csv.gz"), 
                       data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.995)) %>%
  dplyr::select(marker, chr, betas, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)


lti_stressed <- fread(paste0(results_dir,"proportions/LTi_stressed_vs_non_qtl.csv.gz"), 
                      data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  ) %>%
  dplyr::filter(lods > quantile(lods, 0.995)) %>%
  dplyr::select(marker, chr, betas, p_values, 
                lods, chr, pos, cM, ensembl_gene, 
                consequence, strand)

## combining
prop_qtl_df <- list("ILC1_ILC2"=ilc1_ilc2, 
                        "ILC1_ILC3"=ilc1_ilc3,
                        "ILC1_LTi"=ilc1_lti,
                        "ILC2_ILC3"=ilc2_ilc3,
                        "ILC2_LTi"=ilc2_lti,
                        "ILC3_LTi"=ilc3_lti,
                        "ILC3_activated"=ilc3_stressed,
                        "LTi_activated"=lti_stressed) %>%
  bind_rows(.,.id = "trait")


## topic qtl counts ####

topics <- fread(paste0(results_dir, "topics/qtl-topic-lods.csv.gz"), 
                data.table = FALSE) %>%
  dplyr::select(-V1) %>%
  left_join(ccre) %>%
  dplyr::select(1:23) %>%
  pivot_longer(-c(marker, chr, pos)) %>%
  mutate(strand = "+",
         chr = as.character(chr)) %>%
  filter(!is.na(pos), !is.na(chr), value > 6) %>%
  dplyr::rename(trait = name)




## cytokine qtl counts

cytokines <- fread(paste0(results_dir,"cytokines/qtl-cytokines-steady-lods.csv.gz"), 
                   data.table = FALSE) %>%
  dplyr::select(-V1) %>%
  left_join(ccre) %>%
  dplyr::select(1:15) %>%
  pivot_longer(-c(marker, chr, pos)) %>%
  mutate(strand = "+",
         chr = as.character(chr)) %>%
  filter(!is.na(pos), !is.na(chr), value > 6) %>%
  dplyr::rename(trait = name)



## Combining to Dataset ##########

combined_qtl_df <- list(ilc1_eqtl = ilc1_eqtl_loci_by_gene %>% 
       dplyr::select(trait = gene,
                     marker,
                     chr = marker_chr,
                     pos,
                     strand) %>%
       mutate(chr = as.character(chr)),
     ilc2_eqtl = ilc2_eqtl_loci_by_gene %>% 
       dplyr::select(trait = gene,
                     marker,
                     pos,
                     chr = marker_chr,
                     strand) %>%
       mutate(chr = as.character(chr)),
     ilc3_eqtl = ilc3_eqtl_loci_by_gene %>% 
       dplyr::select(trait = gene,
                     marker,
                     pos,
                     chr = marker_chr,
                     strand) %>%
       mutate(chr = as.character(chr)),
     lti_eqtl = lti_eqtl_loci_by_gene %>% 
       dplyr::select(trait = gene,
                     marker,
                     pos,
                     chr = marker_chr,
                     strand) %>%
       mutate(chr = as.character(chr)),
     prop_qtl_df %>% 
       dplyr::select(trait,
                     marker,
                     pos,
                     chr,
                     strand),
     topics %>% 
       dplyr::select(trait,
                     marker,
                     chr,
                     pos,
                     strand),
     cytokines %>% 
       dplyr::select(trait,
                     marker,
                     chr,
                     pos,
                     strand)) %>% 
  bind_rows() %>%
  filter(!is.na(pos)) 


combined_qtl_df %>%
  left_join({
    ccre %>%
      mutate(chr = as.character(chr)) %>% 
      dplyr::select(marker, chr, pos, ensembl_gene)
    }) %>%
  left_join(dplyr::select(ensembl, ensembl_gene = gene_id, eQTL_loci_gene_name = gene_name)) %>%
  left_join(dplyr::select(vars, eQTL_loci_gene_name = index, 
                          ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed)) %>%
write.csv(., "trait_by_loci_no_window.csv", row.names = FALSE)


combined_qtl_df %>%
  group_by(marker) %>%
  mutate(count = n_distinct(trait)) %>%
  ungroup() %>%
  filter(count > 1) %>%
  left_join({
    ccre %>%
      mutate(chr = as.character(chr)) %>% 
      dplyr::select(marker, chr, pos, ensembl_gene)
  }) %>%
  left_join(dplyr::select(ensembl, ensembl_gene = gene_id, eQTL_loci_gene_name = gene_name)) %>%
  left_join(dplyr::select(vars, eQTL_loci_gene_name = index, 
                          ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed)) %>%
  write.csv(., "trait_by_loci_no_window_polygenic_only.csv", row.names = FALSE)


combined_qtl_gr <- makeGRangesFromDataFrame(combined_qtl_df,
                         seqnames="chr",
                         start.field="pos",
                         end.field="pos",
                         strand = "strand",
                         keep.extra.columns = TRUE,
                         na.rm = TRUE)

combined_qtl_grl <- makeGRangesListFromDataFrame(combined_qtl_df,
                                            split.field = "trait",
                                            names.field = "marker",
                                            seqnames = "chr",
                                            start.field = "pos",
                                            end.field = "pos",
                                            strand = "strand",
                                            keep.extra.columns = TRUE)



qtl_count_df <- combined_qtl_df %>%
  group_by(marker) %>%
  summarise(count = n_distinct(trait)) %>%
  left_join(dplyr::select(ccre, marker, chr, pos, ensembl_gene)) %>%
  left_join(dplyr::select(ensembl, ensembl_gene = gene_id, qtl_loci_gene = gene_name)) %>%
  left_join(dplyr::select(vars, qtl_loci_gene = index, 
                          ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed))




## start window search
reduce(resize(sort(combined_qtl_grl), 500000))


test_grl <- reduce(resize(sort(combined_qtl_grl), 10000))
outdf <- lapply(seq_along(test_grl), function(i) {
  
  grl_subset <- test_grl[[i]]
  
  lapply(seq_along(grl_subset), function(j) {
    
    # find overlapping pairs between test loci and other loci
    p <- findOverlapPairs(test_grl, test_grl[[i]][j])
    # filter just to hits
    as.data.frame(pintersect(p)) %>%
      filter(hit) %>%
      mutate(loci = paste0(names(test_grl)[[i]], "_loci", 
                           str_pad(j, pad = 0, side = "left", width = 3)),
             i = i,
             j = j)
    
  }) %>% bind_rows()
  
}) %>% bind_rows(.,.id = "trait") %>%
  ## pivot wider to group by position of commmon overlapping loci
  pivot_wider(id_cols = c(seqnames, start, end,i,j), 
              names_from = loci,
              values_from = hit,
              values_fn = n_distinct) %>% 
  ## pivot longer and filter NAs to get only overlapping loci for each trait
  pivot_longer(-c(seqnames, start, end,i,j)) %>% 
  filter(!is.na(value)) %>% 
  group_by(seqnames, start,end) %>%
  ## count the loci
  mutate(count = n()) %>% 
  ## filter to above 1 to show multiple
  # filter(count > 1) %>%
  mutate(loci = paste0("loci_", str_pad(cur_group_id(), pad = 0, side = "left", width = 3))) %>%
  dplyr::rename(trait = name)


outgr <- outdf %>% distinct(seqnames, start, end) %>% as(.,"GRanges")


annotdf <- lapply(seq_along(outgr), function(i) {
  
  tryCatch({
    tmpdf <- as.data.frame(subsetByOverlaps(as(ensembl,"GRanges"), outgr[i]))
    out <- data.frame(cbind(as.data.frame(outgr[i]), tmpdf[,c("gene_name", "gene_id")]))
    out %>%
      dplyr::select(seqnames, start, end, gene_name, gene_id)
  }, error = function(e) NULL)

}) %>% bind_rows()

fulldf <- full_join(outdf, annotdf, by = c("seqnames","start","end")) %>%
  dplyr::rename(loci_chr = seqnames,
                loci_start = start,
                loci_end = end,
                eQTL_loci_gene_name = gene_name,
                eQTL_loci_gene_id = gene_id) %>%
  mutate(trait = str_remove(trait, "_loci[0-9]{3}")) %>%
  dplyr::select(trait, loci, loci_chr, loci_start, loci_end, 
                eQTL_loci_gene_name, eQTL_loci_gene_id,
                count) %>%
  left_join(dplyr::select(vars, eQTL_loci_gene_name = index, 
                          ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed))


write.csv(fulldf, "trait_by_loci_10kb_window.csv", row.names = FALSE)

fulldf %>%
  dplyr::filter(count > 1) %>%
  write.csv(., "trait_by_loci_10kb_window_polygenic_only.csv", row.names = FALSE)







## check Markers #####

# UNCHS005636
combined_qtl_df %>% 
  filter(marker == "UNCHS005636") %>% 
  count(trait) %>% 
  separate(trait, c("cell_type", "gene")) %>% 
  pivot_wider(id_cols = gene, names_from = cell_type, values_from = n) %>%
  left_join(dplyr::select(vars, gene = index, 
                          ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed))

combined_qtl_df %>% filter(marker %in% c("ICR1041", "ICR5041", "UNCHS005636")) %>%  View



## run through ENRICHR
combined_qtl_df %>% 
  filter(marker == "JAX00511533") %>% 
  count(trait) %>% 
  separate(trait, c("cell_type", "gene")) %>% 
  filter(cell_type == "LTi") %>% 
  pull(gene) %>% 
  str_c(., collapse = "\n") %>%
  cat()


combined_qtl_df %>% 
  filter(marker == "JAX00511533") %>% 
  count(trait) %>% 
  separate(trait, c("cell_type", "gene")) %>% 
  pivot_wider(id_cols = gene, names_from = cell_type, values_from = n) %>%
  left_join(dplyr::select(vars, gene = index, 
                          ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed))
