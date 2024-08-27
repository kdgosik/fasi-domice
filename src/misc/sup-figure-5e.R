#' Figure 5E Supplementary
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

## Itgav IL-4 IL-17A ####################################
## cis-eQTL for all major lineages
## Itgav Chr2:83724397-83806916  ENSMUSG00000027087
start_irange <- 83500000
end_irange <- 84000000
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


# ## topic19
# topics <- fread(paste0(my_path, "results/topic-qtl-lods.csv"), 
#                 select=c("marker", "topic19"))
# 
# start_topic <- end_topic <- topics %>%
#   dplyr::left_join(ccre, by = "marker") %>%
#   dplyr::arrange(pos) %>%
#   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
#   pull(pos)
# 
# data_topic19 <- topics %>%
#   dplyr::left_join(ccre, by = "marker") %>%
#   dplyr::arrange(pos) %>%
#   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
#   pull(topic19)
# 
# topic_dtrack <- DataTrack(data = data_topic19,
#                           start = start_topic,
#                           end = end_topic,
#                           chromosome = chr_num,
#                           genome = gen,
#                           name = "Topic 19")



## IL-4 allergy state
allergy <- fread(paste0(my_path, "results/qtl-cytokines-allergy-lods.csv"),
                 select = c("marker", "IL4", "IL17A"))
start_allergy <- end_allergy <- allergy %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)

data_il4 <- allergy %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(IL4)

il4_dtrack <- DataTrack(data = data_il4,
                        start = start_allergy,
                        end = end_allergy,
                        chromosome = chr_num,
                        genome = gen,
                        name = "IL-4")

## IL-17a allergy state
data_il17a <- allergy %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(IL17A)

il17a_dtrack <- DataTrack(data = data_il17a,
                          start = start_allergy,
                          end = end_allergy,
                          chromosome = chr_num,
                          genome = gen,
                          name = "IL-17A")


pdf(paste0(figure_path, "figure-5-extended-genemodel-Itgav.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, ccre_atrack, 
                il4_dtrack, il17a_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b",
           background.panel = "grey95", background.title = "darkblue")
dev.off()



