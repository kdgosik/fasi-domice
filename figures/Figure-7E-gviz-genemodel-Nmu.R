#' @title Figure 7E Nmu locus Gene Model
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

## setup ####################################
data_path <- paste0(my_path, "data/")
results_path <- paste0(my_path, "results/")
figure_path <- paste0(my_path, "results/figures/")
# source(paste0(figure_path, "helpers.R"))


## read data #######
ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
topics <- fread(paste0(my_path, "results/topics/qtl-topic-lods.csv.gz"), data.table = FALSE)

ilc1 <- fread(paste0(data_dir, "eqtl/qtl-lods-ILC1-cv.csv.gz"),
              select = c("marker", "Lman2"),
              data.table = FALSE)

ccre <- ccre %>% left_join(ilc1)

ilc1_lti_gwas <- fread(paste0(my_path, 'results/proportions/gwas-ilc1-lti-results.csv.gz')) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  )


ilc2_lti_gwas <- fread(paste0(my_path, 'results/proportions/gwas-ilc2-lti-results.csv.gz')) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  )

## Nmu - Chr5:76333495..76363777
start_irange <- 76000000
end_irange <- 77000000
chr_str <- "chr5"
chr_num <- "5"
gen <- "mm10"



marker_list <- ccre %>% 
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(marker)


## Sequence and Ideogram
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)

## Genome Track
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

pos_snps <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)
names_snps <-  ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(marker)

snps_atrack <- AnnotationTrack(start = pos_snps,
                               width = rep(1, length(pos_snps)),
                               chromosome = chr_str,
                               group = names_snps,
                               genome = gen,
                               name = "SNPs")

start_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% 
  pull(start)
end_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% 
  pull(end)
names_ccre <- ccre %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %>% 
  pull(type.y)

ccre_atrack <- AnnotationTrack(start = start_ccre,
                               width = end_ccre - start_ccre,
                               chromosome = chr_str,
                               group = names_ccre,
                               genome = gen,
                               name = "cCREs")


start_topic8 <- end_topic8 <- topics %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>%
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_topic8 <- topics %>%
  dplyr::left_join(ccre, by = "marker") %>%
  dplyr::arrange(pos) %>% 
  dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>%
  dplyr::select(topic8) 
# dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% topic3

topic8_dtrack <- DataTrack(data = data_topic8,
                           start = start_topic8,
                           end = end_topic8,
                           chromosome = chr_num,
                           genome = gen,
                           name = "Topic 8")


start_ilc1_lti_gwas <- end_ilc1_lti_gwas <- ilc1_lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(pos)

data_ilc1_lti_gwas <- ilc1_lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% pull(lods)

ilc1_lti_dtrack <- DataTrack(data = data_ilc1_lti_gwas, 
                             start = start_ilc1_lti_gwas,
                             end = end_ilc1_lti_gwas, 
                             chromosome = chr_num, 
                             genome = gen,
                             name = "ILC1 vs LTi-Like")


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

egene1_dtrack <- create_eGene_track(gene_name = "Lman2", ccre = ccre, chr_num = chr_num, gen = gen)


ht <- HighlightTrack(trackList = list(grtrack, ccre_atrack,
                                      egene1_dtrack,
                                      topic8_dtrack,
                                      ilc1_lti_dtrack,
                                      ilc2_lti_dtrack),
                     start = 76333495-20000, end = 76333495+10000,
                     chromosome = 5)




pdf("gviz-genemodel-Nmu.pdf")
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()
