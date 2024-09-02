#' @title Figure S6F Vipr2 Gene Model
#' @author Kirk Gosik
#' @description
#'

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
my_path <- "/home/rstudio/"
my_path <- "/workspace/fasi-domice/"
data_dir <- paste0(my_path, "data/")
results_dir <- paste0(my_path, "results/")


### Vipr2 (ENSMUSG00000011171)
## Chromosome 12: 116,077,726-116,146,261 forward strand.
start_irange <- 116000000
end_irange <- 116500000
chr_str <- "chr12"
chr_num <- "12"
gen <- "mm10"

ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
## check results
# lti <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-Lti ILC3-cv.csv.gz"), data.table = FALSE)
# lti %>% left_join(ccre) %>% filter(marker_chr == chr_num, between(pos, start_irange, end_irange))


## Slx1b, Rpl35, Ablim1
lti <- fread(paste0(my_path, "data/eqtl/qtl-lods-Lti ILC3-cv.csv.gz"),
              select = c("marker", "Slx1b", "Rpl35", "Ablim1"),
              data.table = FALSE)

ccre <- ccre %>% left_join(lti)


ensembl <- readGFF(paste0(data_dir, "references/Mus_musculus.GRCm38.102.gtf")) %>%
  filter(seqid %in% c(as.character(1:19), "X"), 
         gene_biotype == "protein_coding",
         str_detect(transcript_name, "-201"),
         type %in% c("gene", "exon")) %>%
  dplyr::select(chromosome = seqid, start, end, strand, feature = type,
                gene = gene_id, 
                exon = exon_id, 
                transcript = transcript_id, 
                symbol = gene_name) %>%
  mutate(width = abs(start - end))

grtrack <- GeneRegionTrack(ensembl,
                           chromosome = chr_num, 
                           genome = "mm10", 
                           transcriptAnnotation = "symbol")


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

# start_ccre <- ccre %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% pull(start)
# end_ccre <- ccre %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% pull(end)
# names_ccre <- ccre %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% pull(type.y)
# 
# ccre_atrack <- AnnotationTrack(start = start_ccre,
#                                width = end_ccre - start_ccre,
#                                chromosome = chr_str,
#                                group = names_ccre,
#                                genome = gen,
#                                name = "cCREs")


lti_gwas <- fread(paste0(results_dir,"proportions/LTi_stressed_vs_non_qtl.csv.gz"), 
                      data.table = FALSE) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  )


start_lti_gwas <- end_lti_gwas <- lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_lti_gwas <- lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(lods)/12

lti_dtrack <- DataTrack(data = data_lti_gwas, 
                        start = start_lti_gwas,
                        end = end_lti_gwas+1, 
                        chromosome = chr_num, 
                        genome = gen,
                        name = "ILC3 Activated")


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

## Slx1b, Rpl35, Ablim1
Slx1b_dtrack <- create_eGene_track(gene_name = "Slx1b", ccre = ccre, chr_num = chr_num, gen = gen)
Rpl35_dtrack <- create_eGene_track(gene_name = "Rpl35", ccre = ccre, chr_num = chr_num, gen = gen)
Ablim1_dtrack <- create_eGene_track(gene_name = "Ablim1", ccre = ccre, chr_num = chr_num, gen = gen)


ht <- HighlightTrack(trackList = list(grtrack, Slx1b_dtrack, Rpl35_dtrack),
                     start = 116180000-10000, end = 116180000+10000,
                     chromosome = chr_num)



pdf(paste0(my_path, "results/figures/gviz-genemodel-Vipr2.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()



## check that 3 eGenes lods are highly correlated in the desired region
## Slx1b, Rpl35, Ablim1
ccre  %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange)) %>% 
  dplyr::select(Slx1b, Rpl35, Ablim1) %>% 
  cor(.,use = "pairwise.complete.obs")

## UMAP Slx1b, Rpl35, Ablim1

## setup ###############################
source("/home/rstudio/src/plotting-functions.R")

datadir <- "/home/rstudio/data/"
resultsdir <- "/home/rstudio/results/"

plotdf <- readr::read_csv(paste0(datadir, "manuscript-plot-data.csv"))
Slx1b_expression <- fread(paste0(datadir, "expression/Gene-list-S-expression.csv"),
                         select = c("V1", "Slx1b"))
Rpl35_expression <- fread(paste0(datadir, "expression/Gene-list-R-expression.csv"),
                         select = c("V1", "Rpl35"))
Ablim1_expression <- fread(paste0(datadir, "expression/Gene-list-A-expression.csv"),
                           select = c("V1", "Ablim1"))


plotdf <- plotdf %>%
  dplyr::select("index","X_umap1", "X_umap2") %>%
  left_join(dplyr::rename(Slx1b_expression, index = V1)) %>%
  left_join(dplyr::rename(Rpl35_expression, index = V1))



p1 <- plotdf %>% 
  ggplot(aes_string("X_umap1", "X_umap2", color = "Slx1b")) + 
  geom_point(shape = 46) +
  scale_color_gradient(low = "lightgrey", high="darkblue") +
  theme_void() +
  labs(title = "Slx1b")

ggsave(filename = paste0(resultsdir, "figures/Figure-6-umap-all-cells-Slx1b.pdf"),
       plot = p1,
       width = 7,
       height = 5,
       dpi = 330)


p2 <- plotdf %>% 
  ggplot(aes_string("X_umap1", "X_umap2", color = "Rpl35")) + 
  geom_point(shape = 46) +
  scale_color_gradient(low = "lightgrey", high="darkblue") +
  theme_void() +
  labs(title = "Rpl35")

ggsave(filename = paste0(resultsdir, "figures/Figure-6-umap-all-cells-Rpl35.pdf"),
       plot = p2,
       width = 7,
       height = 5,
       dpi = 330)
