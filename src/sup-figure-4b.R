#' Figure 4B Supplementary
#' 
#' 

project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"
figure_path <- paste0(project_path, "kirk/results/figures/")
source(paste0(figure_path, "helpers.R"))

plot_df <- fread(paste0(my_path, "data/manuscript-plot-data.csv"), data.table = FALSE)

## Gviz section #############

library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
ccre <- fread(paste0(my_path, "data/GM_SNPS_Consequence_cCRE.csv.gz"))
ccre[, recomb_rate := round(cM / (pos / 1e6), 4)]


## Rbpj SNPs - chr5 53553396 - 53661165
start_irange <- 53000000
end_irange <- 54000000
chr_str <- "chr5"
chr_num <- "5"
gen <- "mm10"


gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)

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


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
grtrack <- GeneRegionTrack(txdb,
                           genome = gen,
                           chromosome = chr_str,
                           name = "Gene Model",
                           geneAnnotation = "symbol")


strack <- SequenceTrack(Mmusculus, chromosome = chr_str)


lti_gwas <- fread(paste0(my_path, "results/LTi_stressed_vs_non_qtl.csv.gz"))
ilc2_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-ilc3-results.csv.gz'))
ilc2_lti_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-lti-results.csv.gz'))


start_lti_gwas <- end_lti_gwas <- lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)

data_lti_gwas <- lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(lods)

lti_dtrack <- DataTrack(data = data_lti_gwas,
                        start = start_lti_gwas,
                        end = end_lti_gwas+1,
                        chromosome = chr_num,
                        genome = gen,
                        name = "LTi Activated")


start_ilc2_ilc3_gwas <- end_ilc2_ilc3_gwas <- ilc2_ilc3_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)

data_ilc2_ilc3_gwas <- ilc2_ilc3_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(lods)

ilc2_ilc3_dtrack <- DataTrack(data = data_ilc2_ilc3_gwas,
                              start = start_ilc2_ilc3_gwas,
                              end = end_ilc2_ilc3_gwas,
                              chromosome = chr_num,
                              genome = gen,
                              name = "ILC2 vs ILC3")



start_ilc2_lti_gwas <- end_ilc2_lti_gwas <- ilc2_lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(pos)

data_ilc2_lti_gwas <- ilc2_lti_gwas %>%
  arrange(pos) %>%
  filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %>% 
  pull(lods)

ilc2_lti_dtrack <- DataTrack(data = data_ilc2_lti_gwas,
                             start = start_ilc2_lti_gwas,
                             end = end_ilc2_lti_gwas,
                             chromosome = chr_num,
                             genome = gen,
                             name = "ILC2 vs LTi-Like")


pdf(paste0(figure_path, "Fig.4C-genemodel-Rbpj.pdf"))
plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, ccre_atrack,
                lti_dtrack, ilc2_ilc3_dtrack, ilc2_lti_dtrack),
           from = start_irange, to = end_irange, cex = 0.8, type = "b",
           main = "Rbpj Gene Model", background.title = "darkblue")
dev.off()
