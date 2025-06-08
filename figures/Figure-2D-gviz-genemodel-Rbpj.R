#' @title Figure 2D Rbpj Gene Model
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

### Rbpj
## Rbpj SNPs - chr5 53553396 - 53661165
start_irange <- 53000000
end_irange <- 54000000
chr_str <- "chr5"
chr_num <- "5"
gen <- "mm10"

## read data #######
ccre <- fread(paste0(my_path, "data/references/GM_SNPS_Consequence_cCRE.csv"))
lti_gwas <- fread(paste0(my_path, "results/proportions/LTi_stressed_vs_non_qtl.csv.gz")) %>% 
  mutate(strand = "+") %>% 
  dplyr::mutate(
    padj = p.adjust(p_values, method = "hochberg"),
    lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
  )
ilc2_ilc3_gwas <- fread(paste0(my_path, 'results/proportions/gwas-ilc2-ilc3-results.csv.gz')) %>% 
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

ht <- HighlightTrack(trackList = list(grtrack, ccre_atrack,
                                      lti_dtrack, ilc2_ilc3_dtrack, ilc2_lti_dtrack),
                     start = 53670000-10000, end = 53670000+10000,
                     chromosome = chr_num)

pdf(paste0(my_path, "/results/figures/gviz-genemodel-rbpj.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, ht),
           from = start_irange, to = end_irange, cex = 0.8, type = "b")
dev.off()



