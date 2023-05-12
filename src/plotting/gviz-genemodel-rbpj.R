#' Figure 2D Rbpj Gene Model
#' 
#' 


library(data.table)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
my_path <- "/home/rstudio/domice/"
### Rbpj
## Rbpj SNPs - chr5 53553396 - 53661165
start_irange <- 53000000
end_irange <- 54000000
chr_str <- "chr5"
chr_num <- "5"
gen <- "mm10"

## read data ########
ccre <- fread(paste0(my_path, "data/GM_SNPS_Consequence_cCRE.csv"))
lti_gwas <- fread(paste0(my_path, "results/LTi_stressed_vs_non_qtl.csv.gz"))
ilc2_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-ilc3-results.csv.gz'))
ilc2_lti_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-lti-results.csv.gz'))


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
keepStandardChromosomes(txdb, pruning.mode = "coarse")

gff_file <- paste0(my_path, "data/Mus_musculus.GRCm39.104.gff3.gz")
gff <- readGFF(gff_file)
rbpj <- subset(gff, Parent=="gene:ENSMUSG00000039191")
gff_gr <- readGFFAsGRanges(gff_file)
txdb <- makeTxDbFromGFF(gff_file, format="gff3")
byexons <- exonsBy(txdb, by="gene")


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


grtrack <- GeneRegionTrack(txdb, 
                           genome = gen,
                           chromosome = chr_str, 
                           name = "Gene Model",
                           geneAnnotation = "symbol")






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


start_ilc2_ilc3_gwas <- end_ilc2_ilc3_gwas <- ilc2_ilc3_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_ilc2_ilc3_gwas <- ilc2_ilc3_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(lods)

ilc2_ilc3_dtrack <- DataTrack(data = data_ilc2_ilc3_gwas, 
                              start = start_ilc2_ilc3_gwas,
                              end = end_ilc2_ilc3_gwas, 
                              chromosome = chr_num, 
                              genome = gen,
                              name = "ILC2 vs ILC3")



start_ilc2_lti_gwas <- end_ilc2_lti_gwas <- ilc2_lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_ilc2_lti_gwas <- ilc2_lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(lods)

ilc2_lti_dtrack <- DataTrack(data = data_ilc2_lti_gwas, 
                             start = start_ilc2_lti_gwas,
                             end = end_ilc2_lti_gwas, 
                             chromosome = chr_num, 
                             genome = gen,
                             name = "ILC2 vs LTi-Like")

plotTracks(list(itrack, gtrack, snps_atrack, grtrack, ccre_atrack,
                lti_dtrack, ilc2_ilc3_dtrack, ilc2_lti_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
