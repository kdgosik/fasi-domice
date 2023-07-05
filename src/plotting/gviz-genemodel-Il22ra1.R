#' Il22ra1 Gene Model
#' Chr4: 135,728,172-135,752,140  bp, + strand
#' https://www.informatics.jax.org/marker/MGI:2663588
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

### Il22ra1 (ENSMUSG00000037157)
## Chromosome 4: 135,728,172-135,752,140 
start_irange <- 135600000
end_irange <- 135900000
chr_str <- "chr4"
chr_num <- "4"
gen <- "mm10"

ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))


ilc2_ilc3 <- fread(paste0(results_dir,"proportions/gwas-ilc2-ilc3-results.csv.gz"), 
                   data.table = FALSE) %>% 
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





start_ilc2_ilc3_gwas <- end_ilc2_ilc3_gwas <- ilc2_ilc3 %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_ilc2_ilc3_gwas <- ilc2_ilc3 %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(lods)

ilc2_ilc3_dtrack <- DataTrack(data = data_ilc2_ilc3_gwas,
                        start = start_ilc2_ilc3_gwas,
                        end = end_ilc2_ilc3_gwas+1,
                        chromosome = chr_num,
                        genome = gen,
                        name = "ILC2 vs ILC3")






ht <- HighlightTrack(trackList = list(grtrack, ccre_atrack,
                                      ilc2_ilc3_dtrack),
                     start = 135750000-10000, end = 135750000+10000,
                     chromosome = 4)




pdf(paste0(my_path, "results/figures/gviz-genemodel-Il22ra1.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()



