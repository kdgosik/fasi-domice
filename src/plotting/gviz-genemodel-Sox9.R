#' Sox9 Gene Model
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


## read data #######
ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
## topic8
topics <- fread(paste0(my_path, "results/topics/qtl-topic-lods.csv.gz"))

## topic 8
# topics <- fread(paste0(my_path, "results/topic-qtl-lods.csv"))
## Sox9 - Chr11:112782224-112787760
start_irange <- 112700000
end_irange <- 113200000
chr_str <- "chr11"
chr_num <- "11"
gen <- "mm10"


ensembl <- readGFF(paste0(my_path, "data/references/Mus_musculus.GRCm38.102.gtf")) %>%
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
                           chromosome = 18, 
                           genome = "mm10", 
                           transcriptAnnotation = "symbol")

# plotTracks(grtrack, start_irange, end_irange)

mm10 <- makeGRangesFromDataFrame(ensembl, keep.extra.columns = T)


## make tracks ###############################

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


## topic8
start_topic8 <- end_topic8 <- topics %>%
  left_join({dplyr::select(ccre, marker, pos, chr)}) %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)

data_topic8 <- topics %>%
  left_join({dplyr::select(ccre, marker, pos, chr)}) %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(topic8)

topic8_dtrack <- DataTrack(data = data_topic8, 
                           start = start_topic8,
                           end = end_topic8, 
                           chromosome = chr_num, 
                           genome = gen,
                           name = "Topic 8")


ht <- HighlightTrack(trackList = list(grtrack,
                                      topic8_dtrack),
                     start = 112940000-10000, end = 112940000+10000,
                     chromosome = chr_num)




pdf(paste0(my_path, "results/figures/gviz-genemodel-Sox9.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()
