#' @title Figure S4F PTPN5 Gene Model
#' @author Kirk Gosik
#' @description
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

### Ptpn5 http://www.informatics.jax.org/marker/MGI:97807
## PTPN5 Chromosome 7: 47,077,795-47,133,684 - strand
start_irange <- 46000000
end_irange <- 49000000
chr_str <- "chr7"
chr_num <- "7"
gen <- "mm10"


## read data #######
ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
marker_list <- ccre %>% filter(chr == chr_num, start > start_irange, end < end_irange) %>% pull(marker)

# lti_gwas <- fread(paste0(my_path, "results/LTi_stressed_vs_non_qtl.csv.gz"))
# ilc2_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-ilc3-results.csv.gz'))
# ilc2_lti_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-lti-results.csv.gz'))


vars <- fread(paste0(my_path, "data/allchannels/vars.csv"))

cytokines <- fread(paste0(my_path, "results/cytokines/qtl-cytokines-steady-lods.csv.gz"))
ilc1 <- fread(paste0(my_path, "data/eqtl/qtl-lods-ILC1-cv.csv.gz"), 
              select = c("marker", "Il18r1", "Hspa14"), 
              data.table = FALSE)
  
ccre <- ccre %>% left_join(ilc1)



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


start_ilc1 <- end_ilc1 <- ccre %>%
  filter(marker %in% ilc1$marker) %>%
  arrange(pos) %>% pull(pos)

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

ilc1_il18r1_dtrack <- create_eGene_track(gene_name = "Il18r1", ccre = ccre, chr_num = chr_num, gen = gen)
ilc1_Hspa14_dtrack <- create_eGene_track(gene_name = "Hspa14", ccre = ccre, chr_num = chr_num, gen = gen)



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


ht <- HighlightTrack(trackList = list(grtrack, ccre_atrack,
                                      ilc1_il18r1_dtrack, ilc1_Hspa14_dtrack, ifng_dtrack),
                     start = 47077795-10000, end = 47077795+10000,
                     chromosome = chr_num)


pdf(paste0("gviz-genemodel-ptpn5.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()