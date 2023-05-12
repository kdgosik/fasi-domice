list.of.packages <- c("data.table", "tidyverse", "Gviz", "TxDb.Mmusculus.UCSC.mm10.knownGene", "rtracklayer",
                      "GenomicRanges", "GenomicFeatures", "BSgenome.Mmusculus.UCSC.mm10", "R.utils")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(c("DBI", "data.table", "tidyverse", "R.utils"))
  BiocManager::install(c("Gviz", "TxDb.Mmusculus.UCSC.mm10.knownGene", "rtracklayer", "GenomicRanges",
                         "GenomicFeatures", "BSgenome.Mmusculus.UCSC.mm10"))
}
library(data.table)
library(tidyverse)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DBI)

my_path <- "/home/rstudio/domice/"
### Ptpn5 http://www.informatics.jax.org/marker/MGI:97807
## PTPN5 Chr7:46727543-46783432 bp, - strand
start_irange <- 45000000
end_irange <- 49000000
chr_str <- "chr7"
chr_num <- "7"
gen <- "mm10"


## read data #######
ccre <- fread(paste0(my_path, "data/GM_SNPS_Consequence_cCRE.csv"))
marker_list <- ccre %>% filter(chr == chr_num, start > start_irange, end < end_irange) %>% pull(marker)

# lti_gwas <- fread(paste0(my_path, "results/LTi_stressed_vs_non_qtl.csv.gz"))
# ilc2_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-ilc3-results.csv.gz'))
# ilc2_lti_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-lti-results.csv.gz'))


vars <- fread(paste0(my_path, "data/vars.csv"))
cytokines <- fread(paste0(my_path, "data/qtl-steady-cytokines-lods.csv.gz"))
ilc1 <- fread(paste0(my_path, "data/qtl-lods-ILC1-cv.csv.gz"), 
              select = c("marker", "Il18r1", "Hspa14"), 
              data.table = FALSE) %>% 
  filter(marker %in% marker_list)
  

# ilc1_lods <- tbl(mydb, "ILC1-cv") %>%
#   filter(marker %in% marker_list, name %in% c("Il18r1","Hspa14")) %>%
#   collect()


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
keepStandardChromosomes(txdb, pruning.mode = "coarse")

# gff_file <- paste0(my_path, "data/Mus_musculus.GRCm39.104.gff3.gz")
# gff <- readGFF(gff_file)
# ptpn5 <- subset(gff, Parent=="gene:ENSMUSG00000030854") # ENSMUSG00000030854
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






start_ilc1 <- end_ilc1 <- ccre %>%
  filter(marker %in% ilc1$marker) %>%
  arrange(pos) %>% pull(pos)

# startdata_lti_gwas <- lti_gwas %>%
#   arrange(pos) %>%
#   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods

ilc1_il18a_dtrack <- DataTrack(data = ilc1$Il18r1, 
                        start = start_ilc1,
                        end = end_ilc1+1, 
                        chromosome = chr_num, 
                        genome = gen,
                        name = "ILC1 eGenes")


ilc1_Hspa14_dtrack <- DataTrack(data = ilc1$Hspa14, 
                         start = start_ilc1,
                         end = end_ilc1+1, 
                         chromosome = chr_num, 
                         genome = gen,
                         name = "ILC1 eGenes")


start_ifng<- end_ifng <- ccre %>%
  filter(marker %in% marker_list) %>%
  arrange(pos) %>% pull(pos)

data_ifng <- cytokines %>%
  filter(marker %in% marker_list) %>%
   pull(IFNg)

ifng_dtrack <- DataTrack(data = data_ifng, 
                         start = start_ifng[-2],
                         end = end_ifng[-2], 
                         chromosome = chr_num, 
                         genome = gen,
                         name = "IFNg")


pdf("genemodel-ptpn5.pdf")
plotTracks(list(itrack, gtrack, snps_atrack, grtrack, ccre_atrack,
                ilc1_il18a_dtrack, ilc1_Hspa14_dtrack, ifng_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()