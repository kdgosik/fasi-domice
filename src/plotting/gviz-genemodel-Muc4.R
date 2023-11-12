#' Figure ED 6 Muc2 Gene Model
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
# library(BSgenome.Mmusculus.UCSC.mm10)

my_path <- "/home/rstudio/"
my_path <- "/workspace/fasi-domice/"
## Muc20 - ENSMUSG00000035638
### Muc4 (ENSMUSG00000118672)
## Chromosome 16: 32,735,886-32,782,391 forward strand..
start_irange <- 32600000
end_irange <- 32880000
chr_str <- "chr16"
chr_num <- "16"
gen <- "mm10"

ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
lti_gwas <- fread(paste0(my_path, "results/proportions/LTi_stressed_vs_non_qtl.csv.gz")) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  )


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
                           chromosome = 7, 
                           genome = "mm10", 
                           transcriptAnnotation = "symbol")
mm10 <- makeGRangesFromDataFrame(ensembl, keep.extra.columns = T)


## make tracks ###############################

# strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
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



start_lti_gwas <- end_lti_gwas <- lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_lti_gwas <- lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(lods)

lti_dtrack <- DataTrack(data = data_lti_gwas, 
                        start = start_lti_gwas,
                        end = end_lti_gwas+1, 
                        chromosome = chr_num, 
                        genome = gen,
                        name = "LTi Activated")


pdf(paste0(my_path, "results/figures/gviz-genemodel-Muc4.pdf"))
ht <- HighlightTrack(trackList = list(grtrack, lti_dtrack),
                     start = 32770000-10000, end = 32770000+10000,
                     chromosome = 16)
plotTracks(list(gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()

