#' Figure 1 E Muc2 Gene Model
#' 
#' 


library(data.table)
library(tidyverse)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
# library(BSgenome.Mmusculus.UCSC.mm10)

my_path <- "/home/rstudio/"
### Muc2 (ENSMUSG00000025515)
## Chromosome 7: 141,690,340-141,754,693 forward strand.
start_irange <- 141000000
end_irange <- 142000000
chr_str <- "chr7"
chr_num <- "7"
gen <- "mm10"

ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
## check results
# ilc1 <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-ILC1-cv.csv.gz"), data.table = FALSE)
# ilc1 %>% left_join(ccre) %>% filter(marker_chr == chr_num, between(pos, start_irange, end_irange))

# ilc2 <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-ILC2-cv.csv.gz"), data.table = FALSE)
# ilc2 %>% left_join(ccre) %>% filter(marker_chr == chr_num, between(pos, start_irange, end_irange))

# ilc3 <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz"), data.table = FALSE)
# ilc3 %>% left_join(ccre) %>% filter(marker_chr == chr_num, between(pos, start_irange, end_irange))

# lti <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-Lti ILC3-cv.csv.gz"), data.table = FALSE)
# lti %>% left_join(ccre) %>% filter(marker_chr == chr_num, between(pos, start_irange, end_irange))


## Stub1, Ift20
ilc1 <- fread(file = paste0(my_path, "data/eqtl/qtl-lods-ILC1-cv.csv.gz"),
              select = c("marker", "Stub1", "Ift20"),
              data.table = FALSE)

ccre <- ccre %>% left_join(ilc1)



txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
keepStandardChromosomes(txdb, pruning.mode = "coarse")

## make tracks ###############################3

# strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)



pos_snps <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)
names_snps <-  ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(marker)

snps_atrack <- AnnotationTrack(start = pos_snps,
                               width = rep(1, length(pos_snps)),
                               chromosome = chr_str,
                               group = names_snps,
                               genome = gen,
                               name = "SNPs")

start_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% pull(start)
end_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% pull(end)
names_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% pull(type.y)

ccre_atrack <- AnnotationTrack(start = start_ccre,
                               width = end_ccre - start_ccre,
                               chromosome = chr_str,
                               group = names_ccre,
                               genome = gen,
                               name = "cCREs")


grtrack <- GeneRegionTrack(txdb, 
                           genome = gen,
                           chromosome = chr_str, 
                           name = "Gene Model",
                           geneAnnotation = "symbol")




create_eGene_track <- function(gene_name, ccre, chr_num, gen) {
  
  start_pos <- ccre %>%
    arrange(pos) %>%
    dplyr::rename(gene_lod = matches(gene_name)) %>%
    filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(gene_lod)) %>% pull(pos)


  data_lod <- ccre %>%
    arrange(pos) %>%
    dplyr::rename(gene_lod = matches(gene_name)) %>%
    filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(gene_lod)) %>% pull(gene_lod)
  
  out_dtrack <- DataTrack(data = data_lod, 
                           start = start_pos,
                           end = start_pos+1, 
                           chromosome = chr_num, 
                           genome = gen,
                           name = paste0("ILC1-", gene_name))
  
  out_dtrack

}


## Stub1, Ift20
Stub1_dtrack <- create_eGene_track(gene_name = "Stub1", ccre = ccre, chr_num = chr_num, gen = gen)
Ift20_dtrack <- create_eGene_track(gene_name = "Ift20", ccre = ccre, chr_num = chr_num, gen = gen)

pdf(paste0(my_path, "results/figures/Figure-6-gviz-genemodel-Muc2.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, grtrack, ccre_atrack,
                Stub1_dtrack, Ift20_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()



## UMAP Muc2 ##########################

## setup ###############################
source("/home/rstudio/src/plotting-functions.R")

data_path <- "data/"
results_path <- "results/"

plotdf <- readr::read_csv(paste0(data_path, "manuscript-plot-data.csv"))
Muc2_expression <- fread(paste0(data_path, "expression/Gene-list-M-expression.csv"),
                         select = c("V1", "Muc2"))


plotdf <- plotdf %>%
  dplyr::select("index","X_umap1", "X_umap2") %>%
  left_join(dplyr::rename(Muc2_expression, index = V1))



p1 <- plotdf %>% 
  ggplot(aes_string("X_umap1", "X_umap2", color = "Muc2")) + 
  geom_point(shape = 46) +
  scale_color_gradient(low = "lightgrey", high="darkblue") +
  theme_void() +
  labs(title = "Muc2")

ggsave(filename = paste0(results_path, "figures/Figure-1E-umap-all-cells-Muc2.pdf"),
       plot = p1,
       width = 7,
       height = 5,
       dpi = 330)
