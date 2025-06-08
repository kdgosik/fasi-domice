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
source("../fasi-domice/src/qtl-database.R")




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




# ## start window search
# reduce(resize(sort(combined_qtl_grl), 500000))


test_grl <- reduce(resize(sort(combined_qtl_grl), 1000))
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


write.csv(fulldf, "trait_by_loci_500kb_window.csv", row.names = FALSE)

fulldf %>%
  dplyr::filter(count > 1) %>%
  write.csv(., "trait_by_loci_500kb_window_polygenic_only.csv", row.names = FALSE)




