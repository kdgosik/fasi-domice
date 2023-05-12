library(data.table)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
my_path <- "/home/rstudio/domice/"

## cis-eQTL for all major lineages
## Washc2 (Fam21) Chr6:116208038-116262686  ENSMUSG00000024104 
start_irange <- 116000000
end_irange <- 117000000
chr_str <- "chr6"
chr_num <- "6"
gen <- "mm10"


## read data #######3
ccre <- fread(paste0(my_path, "data/GM_SNPS_Consequence_cCRE.csv"))
ilc1_ilc2_gwas <- fread(paste0(my_path, 'results/gwas-ilc1-ilc2-results.csv.gz'))


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
keepStandardChromosomes(txdb, pruning.mode = "coarse")

# gff_file <- paste0(my_path, "data/Mus_musculus.GRCm39.104.gff3.gz")
# gff <- readGFF(gff_file)
# rbpj <- subset(gff, Parent=="gene:ENSMUSG00000039191")
# gff_gr <- readGFFAsGRanges(gff_file)
# txdb <- makeTxDbFromGFF(gff_file, format="gff3")
byexons <- exonsBy(txdb, by="gene")


## make tracks ###############################3

strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)


pos_snps <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
names_snps <-  ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker

snps_atrack <- AnnotationTrack(start = pos_snps,
                               width = rep(1, length(pos_snps)),
                               chromosome = chr_str,
                               group = names_snps,
                               genome = gen,
                               name = "SNPs")

start_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
end_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
names_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y

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




start_ilc1_ilc2_gwas <- end_ilc1_ilc2_gwas <- ilc1_ilc2_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos

data_ilc1_ilc2_gwas <- ilc1_ilc2_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods

ilc1_ilc2_dtrack <- DataTrack(data = data_ilc1_ilc2_gwas, 
                              start = start_ilc1_ilc2_gwas,
                              end = end_ilc1_ilc2_gwas, 
                              chromosome = chr_num, 
                              genome = gen,
                              name = "ILC1 vs ILC2")




pdf("domice/results/figures/genemodel-Fam21.pdf")
plotTracks(list(itrack, gtrack, snps_atrack, grtrack, ccre_atrack,
                ilc1_ilc2_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()




