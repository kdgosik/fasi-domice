#' Figure 4A Sox9 Gene Model
#' 
#' 


library(data.table)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)

## setup ####################################3
data_path <- "data/"
results_path <- "results/"
figure_path <- "results/figures/"
source(paste0(figure_path, "helpers.R"))


## read data #######3
ccre <- fread(paste0(data_path, "references/GM_SNPS_Consequence_cCRE.csv"))
## topic8
topics <- fread(paste0(results_path, "topics/qtl-topic-lods.csv.gz"))

## topic 8
# topics <- fread(paste0(my_path, "results/topic-qtl-lods.csv"))
## Sox9 - Chr11:112782224-112787760
start_irange <- 112500000
end_irange <- 113250000
chr_str <- "chr11"
chr_num <- "11"
gen <- "mm10"


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

## topic8
start_topic8 <- end_topic8 <- topics %>%
  left_join({dplyr::select(ccre, marker, pos, chr)}) %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_topic8 <- topics %>%
  left_join({dplyr::select(ccre, marker, pos, chr)}) %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(topic8)

topic8_dtrack <- DataTrack(data = data_topic8, 
                           start = start_topic8,
                           end = end_topic8+1, 
                           chromosome = chr_num, 
                           genome = gen,
                           name = "Topic 8")


pdf("Figure-4A-genemodel-Sox9.pdf")
plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, ccre_atrack, topic8_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()
