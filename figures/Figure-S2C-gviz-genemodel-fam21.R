#' @title Figure 2C Fam21 Gene Model
#' @author Kirk Gosik
#' @description
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

## cis-eQTL for all major lineages
## Washc2 (Fam21) Chr6:116208038-116262686  ENSMUSG00000024104 
start_irange <- 116000000
end_irange <- 117000000
chr_str <- "chr6"
chr_num <- "6"
gen <- "mm10"


## read data #######3
ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
ilc1_ilc2_gwas <- fread(paste0(my_path, 'results/proportions/gwas-ilc1-ilc2-results.csv.gz'))


# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# keepStandardChromosomes(txdb, pruning.mode = "coarse")

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

# plotTracks(grtrack, start_irange, end_irange)

mm10 <- makeGRangesFromDataFrame(ensembl, keep.extra.columns = T)


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



start_ilc1_ilc2_gwas <- end_ilc1_ilc2_gwas <- ilc1_ilc2_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_ilc1_ilc2_gwas <- ilc1_ilc2_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(lods)

ilc1_ilc2_dtrack <- DataTrack(data = data_ilc1_ilc2_gwas, 
                              start = start_ilc1_ilc2_gwas,
                              end = end_ilc1_ilc2_gwas, 
                              chromosome = chr_num, 
                              genome = gen,
                              name = "ILC1 vs ILC2")



ht <- HighlightTrack(trackList = list(grtrack, ccre_atrack,
                                      ilc1_ilc2_dtrack),
                     start = 116208038-10000, end = 116208038+10000,
                     chromosome = chr_num)


pdf("results/figures/gviz-genemodel-Fam21.pdf")
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()




