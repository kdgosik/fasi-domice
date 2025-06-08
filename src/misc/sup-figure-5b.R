#' Figure 5B Supplementary
#' 
#' 

project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"
figure_path <- paste0(project_path, "kirk/results/figures/")
source(paste0(figure_path, "helpers.R"))


vars <- fread(paste0(my_path, "data/scanpy/allchannels/vars.csv"), data.table = FALSE)
ccre <- fread(paste0(my_path, "data/GM_SNPS_Consequence_cCRE.csv"))
ccre[, recomb_rate := round(cM / (pos / 1e6), 4)]


library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
## IL-5 and topic19 Gene Model Plot #########################

## cis-eQTL for all major lineages
## Bcl2 Chr1:106538178-106714274 ENSMUSG00000057329
start_irange <- 106400000
end_irange <- 107000000
chr_str <- "chr1"
chr_num <- "1"
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
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)
names_snps <-  ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(marker)

snps_atrack <- AnnotationTrack(start = pos_snps,
                               width = rep(1, length(pos_snps)),
                               chromosome = chr_str,
                               group = names_snps,
                               genome = gen,
                               name = "SNPs")

start_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% 
  pull(start)
end_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% 
  pull(end)
names_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% 
  pull(type.y)

ccre_atrack <- AnnotationTrack(start = start_ccre,
                               width = end_ccre - start_ccre,
                               chromosome = chr_str,
                               group = names_ccre,
                               genome = gen,
                               name = "cCREs")

## topic19
topics <- fread(paste0(my_path, "results/topic-qtl-lods.csv"), 
                select=c("marker", "topic19"))

start_topic <- end_topic <- topics %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)

data_topic19 <- topics %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(topic19)

topic_dtrack <- DataTrack(data = data_topic19,
                          start = start_topic,
                          end = end_topic,
                          chromosome = chr_num,
                          genome = gen,
                          name = "Topic 19")



## IL-5 steady state
steady <- fread(paste0(my_path, "results/qtl-cytokines-steady-lods.csv"),
                select = c("marker", "IL5"))
start_steady <- end_steady <- steady %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)

data_steady <- steady %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(IL5)

steady_dtrack <- DataTrack(data = data_steady,
                           start = start_steady,
                           end = end_steady,
                           chromosome = chr_num,
                           genome = gen,
                           name = "IL-5")

## IL-5 allergy state
allergy <- fread(paste0(my_path, "results/qtl-cytokines-allergy-lods.csv"),
                 select = c("marker", "IL5"))
start_allergy <- end_allergy <- allergy %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)

data_allergy <- allergy %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(IL5)

allergy_dtrack <- DataTrack(data = data_allergy,
                            start = start_allergy,
                            end = end_allergy,
                            chromosome = chr_num,
                            genome = gen,
                            name = "IL-5 (Allergy)")



pdf(paste0(figure_path, "figure-5-extended-genemodel-Bcl2.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, ccre_atrack, 
                topic_dtrack, steady_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b",
           background.panel = "grey95", background.title = "darkblue")
dev.off()


