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


## topic 3
# topics <- fread(paste0(my_path, "results/topic-qtl-lods.csv"))
## Sstr4 - Chr2:148395377-148396764
start_irange <- 148000000
end_irange <- 149000000
chr_str <- "chr2"
chr_num <- "2"
gen <- "mm10"


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

## topic1, topic3, topic6
topics <- fread(paste0(my_path, "results/topics/qtl-topic-lods.csv.gz"))

start_topic3 <- end_topic3 <- topics %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

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

ht <- HighlightTrack(trackList = list(grtrack, ccre_atrack,
                                      topic3_dtrack),
                     start = 148395377-10000, end = 148395377+10000,
                     chromosome = chr_num)

pdf(paste0(my_path, "/results/figures/genemodel-sstr4.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()

