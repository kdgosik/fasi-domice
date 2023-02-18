#' Figure 3E Socs2 Gene Model
#' 
#' 


library(data.table)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)

## setup ####################################
data_path <- "data/"
results_path <- "results/"
figure_path <- "results/figures/"
source(paste0(figure_path, "helpers.R"))


# 95383903..95417321
# Chr10:95383903-95417321 bp, - strand
# From Ensembl annotation of GRCm38
start_irange <- 95000000
end_irange <- 96500000
chr_str <- "chr10"
chr_num <- "10"
gen <- "mm10"

## read data #######3
ccre <- fread(paste0(data_path, "references/GM_SNPS_Consequence_cCRE.csv"))
## topic6
topics <- fread(paste0(results_path, "topics/qtl-topic-lods.csv.gz"))


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
keepStandardChromosomes(txdb, pruning.mode = "coarse")

# gff_file <- paste0(data_path, "references/Mus_musculus.GRCm39.104.gff3.gz")
# gff <- readGFF(gff_file)
# socs2 <- subset(gff, Parent=="gene:ENSMUSG00000020027")
# gff_gr <- readGFFAsGRanges(gff_file)
# txdb <- makeTxDbFromGFF(gff_file, format="gff3")
# byexons <- exonsBy(txdb, by="gene")


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


grtrack <- GeneRegionTrack(txdb, 
                           genome = gen,
                           chromosome = chr_str, 
                           name = "Gene Model",
                           geneAnnotation = "symbol")




start_topic6 <- end_topic6 <- topics %>%
  left_join({dplyr::select(ccre, marker, pos, chr)}) %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_topic6 <- topics %>%
  left_join({dplyr::select(ccre, marker, pos, chr)}) %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(topic6)

topic6_dtrack <- DataTrack(data = data_topic6, 
                        start = start_topic6,
                        end = end_topic6+1, 
                        chromosome = chr_num, 
                        genome = gen,
                        name = "Topic 6")



pdf("Figure-3E-genemodel-Socs2.pdf")
plotTracks(list(itrack, gtrack, snps_atrack, grtrack, ccre_atrack, topic6_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()