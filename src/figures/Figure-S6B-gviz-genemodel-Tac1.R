#' Figure S6B Tac1 Gene Model
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

### Tac1 (ENSMUSG00000061762)
## Chromosome 6: 7,554,879-7,565,834 forward strand.
start_irange <- 7000000
end_irange <- 8000000
chr_str <- "chr6"
chr_num <- "6"
gen <- "mm10"

ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))

## check results
# ilc1 <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-ILC1-cv.csv.gz"), data.table = FALSE)
# ilc1 %>% left_join(ccre) %>% filter(marker_chr == chr_num, between(pos, start_irange, end_irange))

# ilc2 <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-ILC2-cv.csv.gz"), data.table = FALSE)
# ilc2 %>% left_join(ccre) %>% filter(marker_chr == chr_num, between(pos, start_irange, end_irange))

# ilc3 <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz"), data.table = FALSE)
# ilc3 %>% left_join(ccre) %>% filter(marker_chr == chr_num, between(pos, start_irange, end_irange))

# lti <- fread(paste0(my_path, "results/eqtl/qtl-plot-lods-Lti ILC3-cv.csv.gz"), data.table = FALSE)
# lti %>% left_join(ccre) %>% filter(marker_chr == chr_num, between(pos, start_irange, end_irange))


## Rpn2, Sec61b
lti <- fread(paste0(my_path, "data/eqtl/qtl-lods-Lti ILC3-cv.csv.gz"),
              select = c("marker", "Rpn2", "Sec61b"),
              data.table = FALSE)

ccre <- ccre %>% left_join(lti)

lti_gwas <- fread(paste0(my_path, "results/proportions/LTi_stressed_vs_non_qtl.csv.gz"))


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


## LTi GWAS
start_lti_gwas <- end_lti_gwas <- lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_lti_gwas <- lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(lods)/12

lti_dtrack <- DataTrack(data = data_lti_gwas, 
                        start = start_lti_gwas,
                        end = end_lti_gwas+1, 
                        chromosome = chr_num, 
                        genome = gen,
                        name = "LTi-like Activated")




create_eGene_track <- function(gene_name, ccre, chr_num, gen) {
  
  start_pos <- ccre %>%
    arrange(pos) %>%
    dplyr::rename(gene_lod = matches(gene_name)) %>%
    filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(gene_lod)) %>% pull(pos)


  data_lod <- ccre %>%
    arrange(pos) %>%
    dplyr::rename(gene_lod = matches(gene_name)) %>%
    filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(gene_lod)) %>% pull(gene_lod)
  
  out_dtrack <- DataTrack(data = data_lod, 
                           start = start_pos,
                           end = start_pos+1, 
                           chromosome = chr_num, 
                           genome = gen,
                           name = gene_name)
  
  out_dtrack

}


## Rpn2, Sec61b
Rpn2_dtrack <- create_eGene_track(gene_name = "Rpn2", ccre = ccre, chr_num = chr_num, gen = gen)
Sec61b_dtrack <- create_eGene_track(gene_name = "Sec61b", ccre = ccre, chr_num = chr_num, gen = gen)



ht <- HighlightTrack(trackList = list(grtrack,
                                      lti_dtrack,
                                      Rpn2_dtrack,
                                      Sec61b_dtrack),
                     start = 7520000-10000, end = 7520000+10000,
                     chromosome = 6)




pdf(paste0(my_path, "results/figures/gviz-genemodel-Tac1.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()

