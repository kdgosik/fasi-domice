#' @title Figure 4F Gene genemodel Smad2
#' @author Kirk Gosik
#' @description
#'

#' Smad2 Gene Model
#' Chr18:76374651-76444034 bp, + strand
#' https://www.informatics.jax.org/marker/MGI:108051
#' 
#' Needs updated to be accurate, Calca version copy


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
data_dir <- paste0(my_path, "data")

### Smad2 (ENSMUSG00000024563)
## Chromosome 18: 76,241,580-76,310,963 forward strand.
start_irange <- 76200000
end_irange <- 76500000
chr_str <- "chr18"
chr_num <- "18"
gen <- "mm10"

ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
## check results

# IL17F IL9 cytokines
cytokines <- fread(paste0(my_path, "results/cytokines/qtl-cytokines-steady-lods.csv.gz"))

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





## IL-9 steady state
steady <- fread(paste0(my_path, "results/cytokines/qtl-cytokines-steady-lods.csv.gz"),
                select = c("marker", "IL9", "IL17F"))
start_steady <- end_steady <- steady %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)

data_il9_steady <- steady %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(IL9)

steady_il9_dtrack <- DataTrack(data = data_il9_steady,
                           start = start_steady,
                           end = end_steady,
                           chromosome = chr_num,
                           genome = gen,
                           name = "IL-9")


## IL-17F steady state
data_il17f_steady <- steady %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(IL17F)

steady_il17f_dtrack <- DataTrack(data = data_il17f_steady,
                               start = start_steady,
                               end = end_steady,
                               chromosome = chr_num,
                               genome = gen,
                               name = "IL-17F")



ht <- HighlightTrack(trackList = list(grtrack,
                                      steady_il9_dtrack,
                                      steady_il17f_dtrack),
                     start = 76350000-10000, end = 76350000+10000,
                     chromosome = 18)




pdf(paste0(my_path, "results/figures/gviz-genemodel-smad2.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()

