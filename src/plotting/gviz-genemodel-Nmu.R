#' Figure 4D Nmu locus Gene Model
#' 
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

## setup ####################################
data_path <- "data/"
results_path <- "results/"
figure_path <- "results/figures/"
source(paste0(figure_path, "helpers.R"))


## read data #######3
ccre <- fread(paste0(data_path, "references/GM_SNPS_Consequence_cCRE.csv"))
ilc1 <- fread(paste0(results_path, "eqtl/qtl-plot-lods-ILC1-cv.csv.gz"))

## Nmu - Chr5:76333495..76363777
start_irange <- 76000000
end_irange <- 77000000
chr_str <- "chr5"
chr_num <- "5"
gen <- "mm10"

marker_list <- ccre %>%   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(marker)
## Sequence and Ideogram
strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
keepStandardChromosomes(txdb, pruning.mode="coarse")
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)

## Genome Track
grtrack <- GeneRegionTrack(txdb, 
                           genome = gen,
                           chromosome = chr_str, 
                           name = "Gene Model",
                           geneAnnotation = "symbol")

pos_snps <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
names_snps <-  ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker

snps_atrack <- AnnotationTrack(start = pos_snps,
                               width = rep(1, length(pos_snps)),
                               chromosome = chr_str,
                               group = names_snps,
                               genome = gen,
                               name = "SNPs")

start_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
end_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
names_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y

ccre_atrack <- AnnotationTrack(start = start_ccre,
                               width = end_ccre - start_ccre,
                               chromosome = chr_str,
                               group = names_ccre,
                               genome = gen,
                               name = "cCREs")



data_ilc1 <- ilc1 %>%
  dplyr::arrange(xpos) %>%
  dplyr::filter(marker %in% marker_list, gene == "Lman2") %>%
  dplyr::select(marker, value)

lman2_markers <- data_ilc1 %>% pull(marker)
data_ilc1 <- data_ilc1 %>% pull(value)

## Lman2 eGene (Chr13:55491646-55510596 bp, - strand)
# http://www.informatics.jax.org/marker/MGI:1914140
start_ilc1 <- end_ilc1 <- ccre %>%
  arrange(pos) %>%
  filter(marker %in% lman2_markers) %>% pull(pos)

ilc1_dtrack <- DataTrack(data = data_ilc1, 
                           start = start_ilc1,
                           end = end_ilc1+1, 
                           chromosome = chr_num, 
                           genome = gen,
                           name = "ILC1- Lman2")


pdf("Figure-6-genemodel-Nmu.pdf")
plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, ccre_atrack, ilc1_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()
