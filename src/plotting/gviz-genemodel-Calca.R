#' Figure ED 6 Calca Gene Model
#' 
#' 


library(data.table)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
my_path <- "/home/rstudio/"
### Calca (ENSMUSG00000030669)
## Chromosome 7: 114,631,478-114,636,357 reverse strand.
start_irange <- 114500000
end_irange <- 115000000
chr_str <- "chr7"
chr_num <- "7"
gen <- "mm10"

ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
## check results
# ilc2 <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-ILC2-cv.csv.gz"), data.table = FALSE)
# ilc2 %>% left_join(ccre) %>% filter(marker_chr == chr_num, between(pos, start_irange, end_irange))


egene_list <- c("Ebp", "Eif4a2", "Actr10", "Cd81")
ilc2 <- fread(paste0(my_path, "data/eqtl/qtl-lods-ILC2-cv.csv.gz"),
              select = c("marker", "Ebp", "Eif4a2", "Actr10", "Cd81"),
              data.table = FALSE)


ccre <- ccre %>% left_join(ilc2)

# ilc3_gwas <- fread(paste0(my_path, "results/proportion/ILC3_stressed_vs_non_qtl.csv.gz"))


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
keepStandardChromosomes(txdb, pruning.mode = "coarse")

## make tracks ###############################3

strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
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






# start_ilc3_gwas <- end_ilc3_gwas <- ilc3_gwas %>%
#   arrange(pos) %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)
# 
# data_ilc3_gwas <- ilc3_gwas %>%
#   arrange(pos) %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(lods)/12
# 
# ilc3_dtrack <- DataTrack(data = data_ilc3_gwas,
#                         start = start_ilc3_gwas,
#                         end = end_ilc3_gwas+1,
#                         chromosome = chr_num,
#                         genome = gen,
#                         name = "ILC3 Activated")




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
                           name = gene_name)
  
  out_dtrack

}

egene1_dtrack <- create_eGene_track(gene_name = egene_list[1], ccre = ccre, chr_num = chr_num, gen = gen)
egene2_dtrack <- create_eGene_track(gene_name = egene_list[2], ccre = ccre, chr_num = chr_num, gen = gen)
egene3_dtrack <- create_eGene_track(gene_name = egene_list[3], ccre = ccre, chr_num = chr_num, gen = gen)
egene4_dtrack <- create_eGene_track(gene_name = egene_list[4], ccre = ccre, chr_num = chr_num, gen = gen)

pdf(paste0(my_path, "results/figures/Figure-6-gviz-genemodel-calca.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, grtrack, ccre_atrack,
                egene1_dtrack, egene2_dtrack, egene3_dtrack, egene4_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()


## check that 3 eGenes lods are highly correlated in the desired region
ccre  %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange)) %>% 
  dplyr::select(matches(egene_list)) %>% 
  cor(.,use = "pairwise.complete.obs")

## UMAP

## setup ###############################
source("/home/rstudio/src/plotting-functions.R")

datadir <- "/home/rstudio/data/"
resultsdir <- "/home/rstudio/results/"

plotdf <- readr::read_csv(paste0(datadir, "manuscript-plot-data.csv"))
egene_expression <- fread(paste0(datadir, "expression/Gene-list-E-expression.csv"),
                         select = c("V1", "Ebp", "Eif4a2"))
cd81_expression <- fread(paste0(datadir, "expression/Gene-list-C-expression.csv"),
                         select = c("V1", "Cd81"))
actr10_expression <- fread(paste0(datadir, "expression/Gene-list-A-expression.csv"),
                         select = c("V1", "Actr10"))

plotdf <- plotdf %>%
  dplyr::select("index","X_umap1", "X_umap2") %>%
  left_join(dplyr::rename(egene_expression, index = V1)) %>%
  left_join(dplyr::rename(cd81_expression, index = V1))



p1 <- plotdf %>% 
  ggplot(aes_string("X_umap1", "X_umap2", color = "Cd81")) + 
  geom_point(shape = 46) +
  scale_color_gradient(low = "lightgrey", high="darkblue") +
  theme_void() +
  labs(title = "Cd81")

ggsave(filename = paste0(resultsdir, "figures/Figure-6-umap-all-cells-Cd81.pdf"),
       plot = p1,
       width = 7,
       height = 5,
       dpi = 330)


p2 <- plotdf %>% 
  ggplot(aes_string("X_umap1", "X_umap2", color = "Ebp")) + 
  geom_point(shape = 46) +
  scale_color_gradient(low = "lightgrey", high="darkblue") +
  theme_void() +
  labs(title = "Ebp")

ggsave(filename = paste0(resultsdir, "figures/Figure-6-umap-all-cells-Ebp.pdf"),
       plot = p2,
       width = 7,
       height = 5,
       dpi = 330)
