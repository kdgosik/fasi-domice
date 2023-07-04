#' Figure 4E
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
project_path <- "/workspace/fasi-domice/"


ccre <- fread(paste0(project_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
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
topics <- fread(paste0(project_path, "results/topics/qtl-topic-lods.csv.gz"), 
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
steady <- fread(paste0(project_path, "results/cytokines/qtl-cytokines-steady-lods.csv.gz"),
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


ht <- HighlightTrack(trackList = list(grtrack, ccre_atrack,
                                      topic_dtrack, steady_dtrack),
                     start = 106735000-20000, end = 106735000+20000,
                     chromosome = 1)



pdf(paste0(project_path, "results/figures/gviz-genemodel-Bcl2.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()


