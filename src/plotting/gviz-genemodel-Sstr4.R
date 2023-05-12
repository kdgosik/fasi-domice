library(data.table)
library(readxl)
library(dplyr)
library(qtl2)
library(Gviz)

## setup.R ###########
source("../fasi-domice/setup.R")

# GM_Snps meta data
ccre <- fread(paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"), data.table = FALSE)



## topic 3
# topics <- fread(paste0(my_path, "results/topic-qtl-lods.csv"))
## Sstr4 - Chr2:148395377-148396764
start_irange <- 148000000
end_irange <- 149000000
chr_str <- "chr2"
chr_num <- "2"
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

## topic1, topic3, topic6
topics <- fread(paste0(my_path, "results/topic-qtl-lods.csv"))

start_topic3 <- end_topic3 <- topics %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos

data_topic3 <- topics %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>% 
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>%
  dplyr::select(topic1, topic3) 
# dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% topic3

topic3_dtrack <- DataTrack(data = data_topic3,
                           start = start_topic3,
                           end = end_topic3,
                           chromosome = chr_num,
                           groups = c("Topic 1", "Topic 3"),
                           genome = gen,
                           name = "Topic")


plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, ccre_atrack, 
                topic3_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")