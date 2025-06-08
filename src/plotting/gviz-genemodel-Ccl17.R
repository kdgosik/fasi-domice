#' Figure ED 6 Ccl17 Gene Model
#' 
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


### Ccl17 (ENSMUSG00000031780)
## Ccl17 SNPs - Chromosome 8: 95,537,081-95,538,664 forward strand.
start_irange <- 95000000
end_irange <- 96000000
chr_str <- "chr8"
chr_num <- "8"
gen <- "mm10"

## read data #######
ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"), data.table = FALSE)

## check results
# ilc3 <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz"), data.table = FALSE)

## Dap3, Myl12b, Tmsb4x
ilc2 <- fread(paste0(data_dir, "eqtl/qtl-lods-ILC2-cv.csv.gz"),
              select = c("marker", "Dap3", "Myl12b", "Tmsb4x"),
              data.table = FALSE)

## Frmd4b, Cd48
# ilc3 <- fread(paste0(data_dir, "eqtl/qtl-lods-NCR1+\ ILC3-cv.csv.gz"),
#               select = c("marker", "Frmd4b", "Cd48"),
#               data.table = FALSE)

ccre <- ccre %>% left_join(ilc2)

# ilc3_gwas <- fread(paste0(my_path, "results/proportions/ILC3_stressed_vs_non_qtl.csv.gz"))


# ilc2_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-ilc3-results.csv.gz'))
# ilc2_lti_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-lti-results.csv.gz'))

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



# start_ilc3_gwas <- end_ilc3_gwas <- ilc3_gwas %>%
#   arrange(pos) %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)
# 
# data_ilc3_gwas <- ilc3_gwas %>%
#   arrange(pos) %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(lods)/12
# 
# ilc3_dtrack <- DataTrack(data = data_ilc3_gwas, 
#                         start = start_ilc3_gwas,
#                         end = end_ilc3_gwas+1, 
#                         chromosome = chr_num, 
#                         genome = gen,
#                         name = "ILC3 Activated")



# start_Frmd4b <- end_Frmd4b <- ccre %>%
#   arrange(pos) %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(Frmd4b)) %>% pull(pos)
# 
# data_Frmd4b <- ccre %>%
#   arrange(pos) %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(Frmd4b)) %>% pull(Frmd4b)
# 
# Frmd4b_dtrack <- DataTrack(data = data_Frmd4b, 
#                            start = start_Frmd4b,
#                            end = end_Frmd4b+1, 
#                            chromosome = chr_num, 
#                            genome = gen,
#                            name = "Frmd4b")
# 
# 
# start_Cd48 <- end_Cd48 <- ccre %>%
#   arrange(pos) %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(Cd48)) %>% pull(pos)
# 
# data_Cd48 <- ccre %>%
#   arrange(pos) %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(Cd48)) %>% pull(Cd48)
# 
# Cd48_dtrack <- DataTrack(data = data_Cd48, 
#                            start = start_Cd48,
#                            end = end_Cd48+1, 
#                            chromosome = chr_num, 
#                            genome = gen,
#                            name = "Cd48")


start_Dap3 <- end_Dap3 <- ccre %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(Dap3)) %>% pull(pos)

data_Dap3 <- ccre %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(Dap3)) %>% pull(Dap3)

Dap3_dtrack <- DataTrack(data = data_Dap3, 
                         start = start_Dap3,
                         end = end_Dap3+1, 
                         chromosome = chr_num, 
                         genome = gen,
                         name = "Dap3")

# Myl12b
start_Myl12b <- end_Myl12b <- ccre %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(Myl12b)) %>% pull(pos)

data_Myl12b <- ccre %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(Myl12b)) %>% pull(Myl12b)

Myl12b_dtrack <- DataTrack(data = data_Myl12b, 
                         start = start_Myl12b,
                         end = end_Myl12b+1, 
                         chromosome = chr_num, 
                         genome = gen,
                         name = "Myl12b")


# Tmsb4x
start_Tmsb4x <- end_Tmsb4x <- ccre %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(Tmsb4x)) %>% pull(pos)

data_Tmsb4x <- ccre %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(Tmsb4x)) %>% pull(Tmsb4x)

Tmsb4x_dtrack <- DataTrack(data = data_Tmsb4x, 
                           start = start_Tmsb4x,
                           end = end_Myl12b+1, 
                           chromosome = chr_num, 
                           genome = gen,
                           name = "Tmsb4x")


# 94811601
ht <- HighlightTrack(trackList = list(grtrack, ccre_atrack,
                                      Dap3_dtrack, Myl12b_dtrack, Tmsb4x_dtrack),
                     start = 94819451-10000, end = 94819451+10000,
                     chromosome = chr_num)




pdf(paste0(my_path, "results/figures/gviz-genemodel-ccl17.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()




## UMAP Ccl17, Frmd4b and Cd48

## setup ###############################
source("/home/rstudio/src/plotting-functions.R")

datadir <- "/home/rstudio/data/"
resultsdir <- "/home/rstudio/results/"

plotdf <- readr::read_csv(paste0(datadir, "manuscript-plot-data.csv"))
cd48_expression <- fread(paste0(datadir, "expression/Gene-list-C-expression.csv"),
                         select = c("V1", "Ccl17", "Cd48"))
Frmd4b_expression <- fread(paste0(datadir, "expression/Gene-list-F-expression.csv"),
                         select = c("V1", "Frmd4b"))


plotdf <- plotdf %>%
  dplyr::select("index","X_umap1", "X_umap2") %>%
  left_join(dplyr::rename(cd48_expression, index = V1)) %>%
  left_join(dplyr::rename(Frmd4b_expression, index = V1))



p1 <- plotdf %>% 
  ggplot(aes_string("X_umap1", "X_umap2", color = "Cd48")) + 
  geom_point(shape = 46) +
  scale_color_gradient(low = "lightgrey", high="darkblue") +
  theme_void() +
  labs(title = "Cd48")

ggsave(filename = paste0(resultsdir, "figures/Figure-6-umap-all-cells-Cd48.pdf"),
       plot = p1,
       width = 7,
       height = 5,
       dpi = 330)


p2 <- plotdf %>% 
  ggplot(aes_string("X_umap1", "X_umap2", color = "Frmd4b")) + 
  geom_point(shape = 46) +
  scale_color_gradient(low = "lightgrey", high="darkblue") +
  theme_void() +
  labs(title = "Frmd4b")

ggsave(filename = paste0(resultsdir, "figures/Figure-6-umap-all-cells-Frmd4b.pdf"),
       plot = p2,
       width = 7,
       height = 5,
       dpi = 330)
