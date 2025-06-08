#' @title Figure 3B 
#' @author Kirk Gosik
#' @description
#'
#'

project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"
figure_path <- paste0(project_path, "kirk/results/figures/")
data_path <- paste0(project_path, "kirk/data/")
source(paste0(figure_path, "helpers.R"))

marker_map <- readRDS(paste0(geno_path, "Regev_map_20171221.rds"))
topics <- fread(paste0(my_path, "results/topic-qtl-lods.csv.gz"),
                data.table = FALSE)

i <- 1
for( col in paste0("topic", c(1,2,4,5,6,7,9,10,13,14,15,17,18,19)) ) {
  
  # col <- colnames(cytokines)[i]
  p1 <- plot_qtl(qtl_output = topics, marker_map = marker_map, varcolumn = col) + 
    geom_hline(yintercept = 6, col = "red")
  ggsave(filename = paste0(figure_path, "Supp.Fig.3B,", i, "-qtl-", col, ".pdf"),
         plot = p1,
         height = 7,
         width = 7,
         dpi = 330)
  i <- i + 1
  
}







## OLD CODE ##################
# library(Gviz)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# library(BSgenome.Mmusculus.UCSC.mm10)
#' ## Load Gene Data
#' vars <- fread(paste0(my_path, "data/scanpy/allchannels/vars.csv"), data.table = FALSE)
#' 
#' ## Load Marker Data
#' ccre <- fread(paste0(my_path, "data/GM_SNPS_Consequence_cCRE.csv"))
#' ccre[, recomb_rate := round(cM / (pos / 1e6), 4)]
#' 
#' 
#' ## Rbpj ################################
#' 
#' ## Rbpj SNPs - chr5 53553396 - 53661165
#' start_irange <- 53000000
#' end_irange <- 54000000
#' chr_str <- "chr5"
#' chr_num <- "5"
#' gen <- "mm10"
#' 
#' lti_gwas <- fread(paste0(my_path, "results/LTi_stressed_vs_non_qtl.csv")) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' ilc2_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-ilc3-results.csv')) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' ilc2_lti_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-lti-results.csv')) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' 
#' ## Sequence and Ideogram
#' strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' keepStandardChromosomes(txdb, pruning.mode="coarse")
#' gtrack <- GenomeAxisTrack()
#' itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)
#' 
#' ## Genome Track
#' grtrack <- GeneRegionTrack(txdb, 
#'                            genome = gen,
#'                            chromosome = chr_str, 
#'                            name = "Gene Model",
#'                            geneAnnotation = "symbol")
#' 
#' pos_snps <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' names_snps <-  ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker
#' 
#' snps_atrack <- AnnotationTrack(start = pos_snps,
#'                                width = rep(1, length(pos_snps)),
#'                                chromosome = chr_str,
#'                                group = names_snps,
#'                                genome = gen,
#'                                name = "SNPs")
#' 
#' start_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
#' end_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
#' names_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y
#' 
#' ccre_atrack <- AnnotationTrack(start = start_ccre,
#'                                width = end_ccre - start_ccre,
#'                                chromosome = chr_str,
#'                                group = names_ccre,
#'                                genome = gen,
#'                                name = "cCREs")
#' 
#' 
#' start_lti_gwas <- end_lti_gwas <- lti_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_lti_gwas <- lti_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' lti_dtrack <- DataTrack(data = data_lti_gwas, 
#'                         start = start_lti_gwas,
#'                         end = end_lti_gwas+1, 
#'                         chromosome = chr_num, 
#'                         genome = gen,
#'                         name = "LTi Activated")
#' 
#' 
#' start_ilc2_ilc3_gwas <- end_ilc2_ilc3_gwas <- ilc2_ilc3_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_ilc2_ilc3_gwas <- ilc2_ilc3_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' ilc2_ilc3_dtrack <- DataTrack(data = data_ilc2_ilc3_gwas, 
#'                               start = start_ilc2_ilc3_gwas,
#'                               end = end_ilc2_ilc3_gwas, 
#'                               chromosome = chr_num, 
#'                               genome = gen,
#'                               name = "ILC2 vs ILC3")
#' 
#' 
#' 
#' start_ilc2_lti_gwas <- end_ilc2_lti_gwas <- ilc2_lti_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_ilc2_lti_gwas <- ilc2_lti_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' ilc2_lti_dtrack <- DataTrack(data = data_ilc2_lti_gwas, 
#'                              start = start_ilc2_lti_gwas,
#'                              end = end_ilc2_lti_gwas, 
#'                              chromosome = chr_num, 
#'                              genome = gen,
#'                              name = "ILC2 vs LTi-Like")
#' 
#' plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, ccre_atrack,
#'                 lti_dtrack, ilc2_ilc3_dtrack, ilc2_lti_dtrack),
#'            from = start_irange, to = end_irange, cex = 0.8, type = "b")
#' 
#' 
#' ## STIM2 ##################
#' #' @title LTi Activated
#' #' @description Stromal interaction molecules STIM2 are endoplasmic reticulum (ER) 
#' #' Ca2+ sensors that initiate store-operated Ca 2+ entry (SOCE)(Chen et al., 2019) 
#' #' which is important for regulation of calcium-influx. In ILC-neuron axis study, 
#' #' it has been implicated that neuropeptide CGRP limits group 2 innate lymphoid 
#' #' cell responses and constrains type 2 inflammation. NMU controls ILC2s 
#' #' downstream of extracellular signal-regulated kinase and calcium-influx-dependent 
#' #' activation of both calcineurin and nuclear factor of activated T cells 
#' #' (NFAT)(Desvignes et al., 2015). In the immune system, STIM2 participates in 
#' #' T cell activation-induced production of interleukin2 (IL-2) and interferon gamma (IFNγ), 
#' #' probably by stabilization of NFAT residence in the nucleus, as well as in 
#' #' differentiation of naive T cells into Th17 lymphocytes, which presumably are 
#' #' important in early phases of autoimmune diseases(Oh-hora et al., 2008). 
#' #' STIM2 and its homologue STIM1 control the maintenance of CD8 memory and 
#' #' generation of LCMV-specific antibodies by regulating CD40L expression on 
#' #' CD4+ T cells in anti-viral immunity and Stim2-/- dysregulated CD40L 
#' #' expression (Shaw et al., 2014). ILC3 can express CD40 Ligand and 
#' #' MHCII (Robinette and Colonna, 2016)(Komlósi et al., 2018), which might be regulated 
#' #' by ER stress and Stim2 and conduct anti-viral immunity. Here we found the eQTL 
#' #' corresponding gene STIM2 implicated to be responsible for ILC3 plasticity.
#' 
#' 
#' ## Stim2 SNPs - chr5 53553396 - 53661165
#' start_irange <- 53500000
#' end_irange <- 54500000
#' chr_str <- "chr5"
#' chr_num <- "5"
#' gen <- "mm10"
#' 
#' 
#' lti_gwas <- fread(paste0(my_path, "results/LTi_stressed_vs_non_qtl.csv")) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' ilc2_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-ilc3-results.csv')) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' ilc2_lti_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-lti-results.csv')) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' 
#' ## Sequence and Ideogram
#' strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' keepStandardChromosomes(txdb, pruning.mode="coarse")
#' gtrack <- GenomeAxisTrack()
#' itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)
#' 
#' ## Genome Track
#' grtrack <- GeneRegionTrack(txdb, 
#'                            genome = gen,
#'                            chromosome = chr_str, 
#'                            name = "Gene Model",
#'                            geneAnnotation = "symbol")
#' 
#' pos_snps <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' names_snps <-  ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker
#' 
#' snps_atrack <- AnnotationTrack(start = pos_snps,
#'                                width = rep(1, length(pos_snps)),
#'                                chromosome = chr_str,
#'                                group = names_snps,
#'                                genome = gen,
#'                                name = "SNPs")
#' 
#' start_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
#' end_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
#' names_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y
#' 
#' ccre_atrack <- AnnotationTrack(start = start_ccre,
#'                                width = end_ccre - start_ccre,
#'                                chromosome = chr_str,
#'                                group = names_ccre,
#'                                genome = gen,
#'                                name = "cCREs")
#' 
#' 
#' start_lti_gwas <- end_lti_gwas <- lti_gwas %>%
#'   dplyr::arrange(pos) %>%
#'   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_lti_gwas <- lti_gwas %>%
#'   dplyr::arrange(pos) %>%
#'   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' lti_dtrack <- DataTrack(data = data_lti_gwas, 
#'                         start = start_lti_gwas,
#'                         end = end_lti_gwas+1, 
#'                         chromosome = chr_num, 
#'                         genome = gen,
#'                         name = "LTi Activated")
#' 
#' 
#' start_ilc2_ilc3_gwas <- end_ilc2_ilc3_gwas <- ilc2_ilc3_gwas %>%
#'   dplyr::arrange(pos) %>%
#'   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_ilc2_ilc3_gwas <- ilc2_ilc3_gwas %>%
#'   dplyr::arrange(pos) %>%
#'   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' ilc2_ilc3_dtrack <- DataTrack(data = data_ilc2_ilc3_gwas, 
#'                               start = start_ilc2_ilc3_gwas,
#'                               end = end_ilc2_ilc3_gwas, 
#'                               chromosome = chr_num, 
#'                               genome = gen,
#'                               name = "ILC2 vs ILC3")
#' 
#' 
#' 
#' start_ilc2_lti_gwas <- end_ilc2_lti_gwas <- ilc2_lti_gwas %>%
#'   dplyr::arrange(pos) %>%
#'   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_ilc2_lti_gwas <- ilc2_lti_gwas %>%
#'   dplyr::arrange(pos) %>%
#'   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' ilc2_lti_dtrack <- DataTrack(data = data_ilc2_lti_gwas, 
#'                              start = start_ilc2_lti_gwas,
#'                              end = end_ilc2_lti_gwas, 
#'                              chromosome = chr_num, 
#'                              genome = gen,
#'                              name = "ILC2 vs LTi-Like")
#' 
#' plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, ccre_atrack,
#'                 lti_dtrack, ilc2_ilc3_dtrack, ilc2_lti_dtrack),
#'            from = start_irange, to = end_irange, cex = 0.8, type = "b")
#' 
#' 
#' 
#' ## TMEM132c #################
#' #' @title  **ILC1 vs ILC3 proportion**
#' #' @description 5 proteins of the human TMEM132 family (TMEM132A, B, C, D and E). 
#' #' These are genes in which variants are enriched for individuals with hearing loss, 
#' #' panic disorder or cancer. Common variants within the TMEM132E gene are associated 
#' #' with insomnia symptoms(Li et al., 2014) common and rare variants near 
#' #' TMEM132D gene are robustly associated with panic disorder(Haaker et al., 2014); 
#' #' and variants near TMEM132B are associated with excessive daytime 
#' #' sleepiness (Lane et al., 2017). In healthy individuals, some of the TMEM132D 
#' #' non-coding variants exhibit higher anxiety scores and larger volumetric 
#' #' estimates of the amygdala and hippocampus, key neural structures associated 
#' #' with fear and anxiety(Haaker et al., 2014). In the mouse, anterior cingulate cortex 
#' #' TMEM132D expression correlates with anxiety-related behavior (Erhardt et al., 2010). 
#' #' Finally, mutations in TMEM132D are unusually frequent in small-cell lung 
#' #' cancer(Iwakawa et al., 2015)and in pancreatic cancer(Forbes et al., 2014). 
#' #' Disease mutations near, or functions of, TMEM132A or C have yet to be identified. 
#' #' Recent study has shown TMEM132c is highly expressed in human colorectal cancer. 
#' #' Here we found TEME132c is a SNP-gene with eQTL effect associated with 
#' #' ILC1/ILC3 relative composition in mouse and its expression correlates with…
#' #' Suggest the potential involvement of TEME protein in mediating the disease 
#' #' susceptibility to host defense and inflammation.
#' 
#' ## Tmem132c Chr5:127241826-127565790  ENSMUSG00000034324 
#' start_irange <- 127000000
#' end_irange <- 128000000
#' chr_str <- "chr5"
#' chr_num <- "5"
#' gen <- "mm10"
#' 
#' ilc1_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc1-ilc3-results.csv')) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' ## Sequence and Ideogram
#' strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' keepStandardChromosomes(txdb, pruning.mode="coarse")
#' gtrack <- GenomeAxisTrack()
#' itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)
#' 
#' ## Genome Track
#' grtrack <- GeneRegionTrack(txdb, 
#'                            genome = gen,
#'                            chromosome = chr_str, 
#'                            name = "Gene Model",
#'                            geneAnnotation = "symbol")
#' 
#' pos_snps <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' names_snps <-  ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker
#' 
#' snps_atrack <- AnnotationTrack(start = pos_snps,
#'                                width = rep(1, length(pos_snps)),
#'                                chromosome = chr_str,
#'                                group = names_snps,
#'                                genome = gen,
#'                                name = "SNPs")
#' 
#' start_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
#' end_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
#' names_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y
#' 
#' ccre_atrack <- AnnotationTrack(start = start_ccre,
#'                                width = end_ccre - start_ccre,
#'                                chromosome = chr_str,
#'                                group = names_ccre,
#'                                genome = gen,
#'                                name = "cCREs")
#' 
#' start_ilc1_ilc3_gwas <- end_ilc1_ilc3_gwas <- ilc1_ilc3_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_ilc1_ilc3_gwas <- ilc1_ilc3_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' ilc1_ilc3_dtrack <- DataTrack(data = data_ilc1_ilc3_gwas, 
#'                               start = start_ilc1_ilc3_gwas,
#'                               end = end_ilc1_ilc3_gwas, 
#'                               chromosome = chr_num, 
#'                               genome = gen,
#'                               name = "ILC1 vs ILC3")
#' 
#' 
#' plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, 
#'                 ccre_atrack, ilc1_ilc3_dtrack),
#'            from = start_irange, to = end_irange, cex = 0.8, type = "b")
#' 
#' 
#' 
#' 
#' 
#' ## Lgr6 ########################
#' #' @title ILC1 vs ILC3 proportion
#' #' @description the 7-transmembrane receptor, Lgr5 has recently gained prominence 
#' #' as a marker of Wnt-regulated adult stem cell populations in the hair-follicle, 
#' #' intestine and stomach. A closely-related protein, Lgr6 marks adult stem cells 
#' #' responsible for fueling the renewal of the sebaceous gland and skin(Leushacke 
#' #' and Barker, 2011). Lgr6 as a marker of adult stem cells in the skin(Barker et al., 2013). 
#' #' Lgr6 is expressed by airway smooth muscle cells (aSMCs) lining the bronchiolar 
#' #' epithelium. Secretion of Fgf10 by Lgr6 + aSMCs promotes differentiation of airway 
#' #' epithelial cells(Leung et al., 2018).While Igr6 is well implicated in colon cancer 
#' #' as an intestinal stem cell marker. We found the eQTL for this gene implicated here 
#' #' with the relative proportion of ILC1/ILC3 which further emphasize the importance of 
#' #' epithelium cell for the development of ILCs.
#' 
#' 
#' 
#' ## cis-eQTL for all major lineages
#' ## Lgr6 Chr1:134983301-135105276  ENSMUSG00000042793 
#' start_irange <- 134500000
#' end_irange <- 135500000
#' chr_str <- "chr1"
#' chr_num <- "1"
#' gen <- "mm10"
#' 
#' 
#' ilc1_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc1-ilc3-results.csv')) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' ## Sequence and Ideogram
#' strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' keepStandardChromosomes(txdb, pruning.mode="coarse")
#' gtrack <- GenomeAxisTrack()
#' itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)
#' 
#' ## Genome Track
#' grtrack <- GeneRegionTrack(txdb, 
#'                            genome = gen,
#'                            chromosome = chr_str, 
#'                            name = "Gene Model",
#'                            geneAnnotation = "symbol")
#' 
#' pos_snps <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' names_snps <-  ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker
#' 
#' snps_atrack <- AnnotationTrack(start = pos_snps,
#'                                width = rep(1, length(pos_snps)),
#'                                chromosome = chr_str,
#'                                group = names_snps,
#'                                genome = gen,
#'                                name = "SNPs")
#' 
#' start_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
#' end_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
#' names_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y
#' 
#' ccre_atrack <- AnnotationTrack(start = start_ccre,
#'                                width = end_ccre - start_ccre,
#'                                chromosome = chr_str,
#'                                group = names_ccre,
#'                                genome = gen,
#'                                name = "cCREs")
#' 
#' start_ilc1_ilc3_gwas <- end_ilc1_ilc3_gwas <- ilc1_ilc3_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_ilc1_ilc3_gwas <- ilc1_ilc3_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' ilc1_ilc3_dtrack <- DataTrack(data = data_ilc1_ilc3_gwas, 
#'                               start = start_ilc1_ilc3_gwas,
#'                               end = end_ilc1_ilc3_gwas, 
#'                               chromosome = chr_num, 
#'                               genome = gen,
#'                               name = "ILC1 vs ILC3")
#' 
#' 
#' plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, 
#'                 ccre_atrack, ilc1_ilc3_dtrack),
#'            from = start_irange, to = end_irange, cex = 0.8, type = "b")
#' 
#' 
#' 
#' 
#' ## Rrbp1 #####################
#' #' @title ILC1 vs ILC3 proportion
#' #' @description Rrbp1 is highly expressed in epithelium cells in small intestine. 
#' #' Ribosome-binding protein 1 (RRBP1) has been implicated in the regulation of 
#' #' unfolded protein response, which is involved in almost every aspect of cancer 
#' #' development. Further, Endoplasmic reticulum ribosome-binding protein 1, 
#' #' RRBP1, promotes progression of colorectal cancer and predicts an unfavourable 
#' #' prognosis(Pan et al., 2015). over the last decade have shown that the UPR 
#' #' plays a critical role in shaping immunity and inflammation, resulting in the 
#' #' recognition of the UPR as a key player in pathological processes including 
#' #' complex inflammatory, autoimmune and neoplastic diseases. Intestinal epithelium, 
#' #' with its many highly secretory cells, forms an important barrier and messenger 
#' #' between the luminal environment and the host immune system. It is not surprising, 
#' #' that numerous studies have associated ER stress and the UPR with intestinal diseases such as inflammatory bowel disease (IBD) and colorectal cancer (CRC). The UPR in immune cells including B cell, T cells, macrophage and epithelium cells, all of which has great crosstalk with innate lymphoid cells(Coleman and Haller, 2019). We detected the eQTLs for Rrbp1 gene which accounts for ILC1/ILC3 composition. Our observation again emphasized the significance of UPR in regulating innate immune response in terms of ILC1 and ILC3 composition in small intestine.
#' 
#' 
#' ## Rrbp1 Chr2:143947395-144011263  ENSMUSG00000027422 
#' start_irange <- 143500000
#' end_irange <- 144500000
#' chr_str <- "chr2"
#' chr_num <- "2"
#' gen <- "mm10"
#' 
#' ilc1_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc1-ilc3-results.csv')) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' ## Sequence and Ideogram
#' strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' keepStandardChromosomes(txdb, pruning.mode="coarse")
#' gtrack <- GenomeAxisTrack()
#' itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)
#' 
#' ## Genome Track
#' grtrack <- GeneRegionTrack(txdb, 
#'                            genome = gen,
#'                            chromosome = chr_str, 
#'                            name = "Gene Model",
#'                            geneAnnotation = "symbol")
#' 
#' pos_snps <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' names_snps <-  ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker
#' 
#' snps_atrack <- AnnotationTrack(start = pos_snps,
#'                                width = rep(1, length(pos_snps)),
#'                                chromosome = chr_str,
#'                                group = names_snps,
#'                                genome = gen,
#'                                name = "SNPs")
#' 
#' start_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
#' end_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
#' names_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y
#' 
#' ccre_atrack <- AnnotationTrack(start = start_ccre,
#'                                width = end_ccre - start_ccre,
#'                                chromosome = chr_str,
#'                                group = names_ccre,
#'                                genome = gen,
#'                                name = "cCREs")
#' 
#' 
#' start_ilc1_ilc3_gwas <- end_ilc1_ilc3_gwas <- ilc1_ilc3_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_ilc1_ilc3_gwas <- ilc1_ilc3_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' ilc1_ilc3_dtrack <- DataTrack(data = data_ilc1_ilc3_gwas, 
#'                               start = start_ilc1_ilc3_gwas,
#'                               end = end_ilc1_ilc3_gwas, 
#'                               chromosome = chr_num, 
#'                               genome = gen,
#'                               name = "ILC1 vs ILC3")
#' 
#' plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, 
#'                 ccre_atrack, ilc1_ilc3_dtrack),
#'            from = start_irange, to = end_irange, cex = 0.8, type = "b")
#' 
#' 
#' 
#' ## IL-20RA ######################
#' #' @title ILC3 vs LTi-like proportion
#' #' @description The interleukin-20-receptor I complex (IL-20-RI) is composed of 
#' #' two chains, IL20RA and IL20RB. Its ligands are the three members of the IL19 
#' #' subfamily of cytokines, IL-19, IL-20 and IL-24. These cytokines are important 
#' #' in the manifestation of psoriatic lesions and, recently, an association of 
#' #' polymorphisms of IL20 with psoriasis has been described(Kingo et al., 2008). 
#' #' Capture Hi-C Identifies a Novel Causal Gene, IL20RA, in the Pan-Autoimmune 
#' #' Genetic Susceptibility Region 6q23. Interleukin 20 receptor, alpha subunit is a 
#' #' subunit for the interleukin-20 receptor. The risk allele of the most likely 
#' #' causal SNP, rs6927172, is correlated with both a higher frequency of interactions 
#' #' and increased expression of IL20RA, along with a stronger binding of both the NFκB 
#' #' transcription factor and chromatin marks characteristic of active enhancers in 
#' #' T-cells(McGovern et al., 2016). Here, we found the eQTL for IL-20RA is responsible 
#' #' for distinguish between ILC3 and LTi-like ILC3.
#' 
#' 
#' ## Il20RA Chr10:19712570-19760053  ENSMUSG00000020007 
#' start_irange <- 19000000
#' end_irange <- 20000000
#' chr_str <- "chr10"
#' chr_num <- "10"
#' gen <- "mm10"
#' 
#' ilc3_lti_gwas <- fread(paste0(my_path, 'results/gwas-ilc3-lti-results.csv')) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' ## Sequence and Ideogram
#' strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' keepStandardChromosomes(txdb, pruning.mode="coarse")
#' gtrack <- GenomeAxisTrack()
#' itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)
#' 
#' ## Genome Track
#' grtrack <- GeneRegionTrack(txdb, 
#'                            genome = gen,
#'                            chromosome = chr_str, 
#'                            name = "Gene Model",
#'                            geneAnnotation = "symbol")
#' 
#' pos_snps <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' names_snps <-  ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker
#' 
#' snps_atrack <- AnnotationTrack(start = pos_snps,
#'                                width = rep(1, length(pos_snps)),
#'                                chromosome = chr_str,
#'                                group = names_snps,
#'                                genome = gen,
#'                                name = "SNPs")
#' 
#' start_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
#' end_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
#' names_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y
#' 
#' ccre_atrack <- AnnotationTrack(start = start_ccre,
#'                                width = end_ccre - start_ccre,
#'                                chromosome = chr_str,
#'                                group = names_ccre,
#'                                genome = gen,
#'                                name = "cCREs")
#' 
#' 
#' start_ilc3_lti_gwas <- end_ilc3_lti_gwas <- ilc3_lti_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_ilc3_lti_gwas <- ilc3_lti_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' ilc3_lti_dtrack <- DataTrack(data = data_ilc3_lti_gwas, 
#'                               start = start_ilc3_lti_gwas,
#'                               end = end_ilc3_lti_gwas, 
#'                               chromosome = chr_num, 
#'                               genome = gen,
#'                               name = "ILC3 vs LTi-like")
#' 
#' plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, 
#'                 ccre_atrack, ilc3_lti_dtrack),
#'            from = start_irange, to = end_irange, cex = 0.8, type = "b")
#' 
#' 
#' 
#' 
#' ## Fam21 ################
#' #' @title ILC1 vs ILC2 proportion
#' #' @description WASH exists in a multiprotein complex containing FAM21, 
#' #' which links WASH to endosomes and is required for WASH-dependent 
#' #' retromer-mediated sorting(Gomez and Billadeau, 2009; Jia et al., 2012). 
#' #' WASH deletion impairs the cell pool of NKp46+ILC3s. In NKp46+ ILC3s, 
#' #' WASH recruits Arid1a to the Ahr promoter thus activating AHR expression, 
#' #' indicating WASH-mediated AHR expression has a critical function in the maintenance 
#' #' of NKp46+ILC3s(Xia et al., 2017). While NKp46+ is also expressed in ILC1 and 
#' #' here in our data we found eQTLs for FAM21 tips the balance of relative proportion 
#' #' of IL1/ILC2 indicating its role in regulating ILC development.
#' 
#' 
#' ## Washc2 (Fam21) Chr6:116208038-116262686  ENSMUSG00000024104 
#' start_irange <- 116000000
#' end_irange <- 117000000
#' chr_str <- "chr6"
#' chr_num <- "6"
#' gen <- "mm10"
#' 
#' ilc1_ilc2_gwas <- fread(paste0(my_path, 'results/gwas-ilc1-ilc2-results.csv')) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' ## Sequence and Ideogram
#' strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' keepStandardChromosomes(txdb, pruning.mode="coarse")
#' gtrack <- GenomeAxisTrack()
#' itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)
#' 
#' ## Genome Track
#' grtrack <- GeneRegionTrack(txdb, 
#'                            genome = gen,
#'                            chromosome = chr_str, 
#'                            name = "Gene Model",
#'                            geneAnnotation = "symbol")
#' 
#' pos_snps <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' names_snps <-  ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker
#' 
#' snps_atrack <- AnnotationTrack(start = pos_snps,
#'                                width = rep(1, length(pos_snps)),
#'                                chromosome = chr_str,
#'                                group = names_snps,
#'                                genome = gen,
#'                                name = "SNPs")
#' 
#' start_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
#' end_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
#' names_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y
#' 
#' ccre_atrack <- AnnotationTrack(start = start_ccre,
#'                                width = end_ccre - start_ccre,
#'                                chromosome = chr_str,
#'                                group = names_ccre,
#'                                genome = gen,
#'                                name = "cCREs")
#' 
#' start_ilc1_ilc2_gwas <- end_ilc1_ilc2_gwas <- ilc1_ilc2_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_ilc1_ilc2_gwas <- lic1_ilc2_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' ilc1_ilc2_dtrack <- DataTrack(data = data_ilc1_ilc2_gwas, 
#'                              start = start_ilc1_ilc2_gwas,
#'                              end = end_ilc1_ilc2_gwas, 
#'                              chromosome = chr_num, 
#'                              genome = gen,
#'                              name = "ILC1 vs ILC2")
#' 
#' plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, 
#'                 ccre_atrack, ilc1_ilc2_dtrack),
#'            from = start_irange, to = end_irange, cex = 0.8, type = "b")
#' 
#' 
#' 
#' 
#' ## BANF2 ###########
#' #' @title ILC2 vs ILC3
#' #' @description BANF2 is termed for Barrier to autointegration factor 2, 
#' #' which has been associated with the abnormality of immune system physiology 
#' #' phenotype in GWAS datasets from the GWASdb SNP-Phenotype Associations dataset. 
#' #' This immune disease implicated SNP-gene is here to identified to be a eQTL-gene 
#' #' which influence the ILC development by tipping the balance of ILC2/ILC3 relative 
#' #' proportion.
#' 
#' 
#' ## Banf2 Chr2:144033102-144073979  ENSMUSG00000037307 
#' start_irange <- 144000000
#' end_irange <- 145000000
#' chr_str <- "chr2"
#' chr_num <- "2"
#' gen <- "mm10"
#' 
#' 
#' ilc2_ilc3_gwas <- fread(paste0(my_path, 'results/gwas-ilc2-ilc3-results.csv')) %>%
#'   dplyr::mutate(
#'     padj = p.adjust(p_values, method = "hochberg"),
#'     lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
#'   )
#' 
#' ## Sequence and Ideogram
#' strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' keepStandardChromosomes(txdb, pruning.mode="coarse")
#' gtrack <- GenomeAxisTrack()
#' itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)
#' 
#' ## Genome Track
#' grtrack <- GeneRegionTrack(txdb, 
#'                            genome = gen,
#'                            chromosome = chr_str, 
#'                            name = "Gene Model",
#'                            geneAnnotation = "symbol")
#' 
#' pos_snps <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' names_snps <-  ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker
#' 
#' snps_atrack <- AnnotationTrack(start = pos_snps,
#'                                width = rep(1, length(pos_snps)),
#'                                chromosome = chr_str,
#'                                group = names_snps,
#'                                genome = gen,
#'                                name = "SNPs")
#' 
#' start_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
#' end_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
#' names_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y
#' 
#' ccre_atrack <- AnnotationTrack(start = start_ccre,
#'                                width = end_ccre - start_ccre,
#'                                chromosome = chr_str,
#'                                group = names_ccre,
#'                                genome = gen,
#'                                name = "cCREs")
#' 
#' start_ilc2_ilc3_gwas <- end_ilc2_ilc3_gwas <- ilc2_ilc3_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_ilc2_ilc3_gwas <- lic1_ilc2_gwas %>%
#'   arrange(pos) %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% lods
#' 
#' ilc2_ilc3_dtrack <- DataTrack(data = data_ilc2_ilc3_gwas, 
#'                               start = start_ilc2_ilc3_gwas,
#'                               end = end_ilc2_ilc3_gwas, 
#'                               chromosome = chr_num, 
#'                               genome = gen,
#'                               name = "ILC2 vs ILC3")
#' 
#' plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack, 
#'                 ccre_atrack, ilc2_ilc3_dtrack),
#'            from = start_irange, to = end_irange, cex = 0.8, type = "b")
#' 
#' 
#' 
#' ## Gm43781 ###########
#' #' @title Topic 11 (ILC2)
#' #' @description Gm43781 predicted gene from MGI located in a significant loci for topic 11 
#' #' QTL
#' 
#' 
#' ## Gm43781 Chr5:54807334-54810014  ENSMUSG00000106126 
#' start_irange <- 54550000
#' end_irange <- 54950000
#' chr_str <- "chr5"
#' chr_num <- "5"
#' gen <- "mm10"
#' 
#' 
#' ## topic11
#' topics <- fread(paste0(my_path, "results/topic-qtl-lods.csv"))
#' 
#' start_topic <- end_topic <- topics %>%
#'   dplyr::left_join(ccre, by = "marker") %>%
#'   dplyr::arrange(pos) %>%
#'   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_topic <- topics %>%
#'   dplyr::left_join(ccre, by = "marker") %>%
#'   dplyr::arrange(pos) %>%
#'   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% topic11
#' 
#' topic_dtrack <- DataTrack(data = data_topic,
#'                     start = start_topic,
#'                     end = end_topic,
#'                     chromosome = chr_num,
#'                     genome = gen,
#'                     name = "Topic 11")
#' 
#' ## Sequence and Ideogram
#' strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' keepStandardChromosomes(txdb, pruning.mode="coarse")
#' gtrack <- GenomeAxisTrack()
#' itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)
#' 
#' ## Genome Track
#' # grtrack <- GeneRegionTrack(txdb, 
#' #                            genome = gen,
#' #                            chromosome = chr_str, 
#' #                            name = "Gene Model",
#' #                            geneAnnotation = "symbol")
#' 
#' ## Chr5:54807334-54810014
#' grtrack <- AnnotationTrack(
#'   range=GRanges(seqnames = chr_str, 
#'                 ranges = IRanges(start = 54807334,  end = 54810014),
#'                 names = c("Gm43781"),
#'                 strand = c("+")), 
#'   genome = gen, 
#'   chromosome = chr_str, name = "Gene Model")
#' 
#' 
#' pos_snps <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' names_snps <-  ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker
#' 
#' snps_atrack <- AnnotationTrack(start = pos_snps,
#'                                width = rep(1, length(pos_snps)),
#'                                chromosome = chr_str,
#'                                group = names_snps,
#'                                genome = gen,
#'                                name = "SNPs")
#' 
#' start_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
#' end_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
#' names_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y
#' 
#' ccre_atrack <- AnnotationTrack(start = start_ccre,
#'                                width = end_ccre - start_ccre,
#'                                chromosome = chr_str,
#'                                group = names_ccre,
#'                                genome = gen,
#'                                name = "cCREs")
#' 
#' 
#' 
#' plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack,
#'                 ccre_atrack, topic_dtrack),
#'            from = start_irange, to = end_irange, cex = 0.8, type = "b")
#' 
#' 
#' 
#' 
#' ## Topic 0
#' # chr2:171146966..171646966
#' start_irange <- 171000000
#' end_irange <- 172000000
#' chr_str <- "chr2"
#' chr_num <- "2"
#' gen <- "mm10"
#' 
#' 
#' ## topic0
#' topics <- fread(paste0(my_path, "results/topic-qtl-lods.csv"))
#' 
#' start_topic <- end_topic <- topics %>%
#'   dplyr::left_join(ccre, by = "marker") %>%
#'   dplyr::arrange(pos) %>%
#'   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' 
#' data_topic <- topics %>%
#'   dplyr::left_join(ccre, by = "marker") %>%
#'   dplyr::arrange(pos) %>%
#'   dplyr::filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% topic0
#' 
#' topic_dtrack <- DataTrack(data = data_topic,
#'                           start = start_topic,
#'                           end = end_topic,
#'                           chromosome = chr_num,
#'                           genome = gen,
#'                           name = "Topic 0")
#' 
#' ## Sequence and Ideogram
#' strack <- SequenceTrack(Mmusculus, chromosome = chr_str)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' keepStandardChromosomes(txdb, pruning.mode="coarse")
#' gtrack <- GenomeAxisTrack()
#' itrack <- IdeogramTrack(genome = gen, chromosome = chr_str)
#' 
#' ## Genome Track
#' grtrack <- GeneRegionTrack(txdb,
#'                            genome = gen,
#'                            chromosome = chr_str,
#'                            name = "Gene Model",
#'                            geneAnnotation = "symbol")
#' 
#' # ## Chr5:54807334-54810014
#' # grtrack <- AnnotationTrack(
#' #   range=GRanges(seqnames = chr_str, 
#' #                 ranges = IRanges(start = 54807334,  end = 54810014),
#' #                 names = c("Gm43781"),
#' #                 strand = c("+")), 
#' #   genome = gen, 
#' #   chromosome = chr_str, name = "Gene Model")
#' 
#' 
#' pos_snps <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% pos
#' names_snps <-  ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(pos)) %$% marker
#' 
#' snps_atrack <- AnnotationTrack(start = pos_snps,
#'                                width = rep(1, length(pos_snps)),
#'                                chromosome = chr_str,
#'                                group = names_snps,
#'                                genome = gen,
#'                                name = "SNPs")
#' 
#' start_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% start
#' end_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% end
#' names_ccre <- ccre %>%
#'   filter(chr == chr_num, between(pos, start_irange, end_irange), !is.na(type.y)) %$% type.y
#' 
#' ccre_atrack <- AnnotationTrack(start = start_ccre,
#'                                width = end_ccre - start_ccre,
#'                                chromosome = chr_str,
#'                                group = names_ccre,
#'                                genome = gen,
#'                                name = "cCREs")
#' 
#' 
#' 
#' plotTracks(list(itrack, gtrack, snps_atrack, grtrack, strack,
#'                 ccre_atrack, topic_dtrack),
#'            from = start_irange, to = end_irange, cex = 0.8, type = "b")
#' 
#' 
#' 
#' ## USCS Genome Example ########
#' from <- 54500000
#' to <- 55000000
#' chr_str <- "chr5"
#' chr_num <- "5"
#' gen <- "mm9"
#' # from <- 65921878
#' # to <- 65980988
#' knownGenes <- UcscTrack(genome = gen, chromosome = chr_str, 
#'                         track = "knownGene", from = from, to = to,
#'                         trackType = "GeneRegionTrack", 
#'                         rstarts = "exonStarts", rends = "exonEnds", 
#'                         gene = "name", symbol = "name", 
#'                         transcript = "name", strand = "strand", 
#'                         fill = "#8282d2", name = "UCSC Genes")
#' 
#' refGenes <- UcscTrack(genome = gen, chromosome = chr_str,
#'                       track = "xenoRefGene", from = from, to = to,
#'                       trackType = "GeneRegionTrack", 
#'                       rstarts = "exonStarts", rends = "exonEnds", 
#'                       gene = "name",  symbol = "name2", 
#'                       transcript = "name", strand = "strand",
#'                       fill = "#8282d2", stacking = "dense", 
#'                       name = "Other RefSeq")
#' 
#' ensGenes <- UcscTrack(genome = gen, chromosome = chr_str,
#'                       track = "ensGene", from = from, to = to,
#'                       trackType = "GeneRegionTrack", 
#'                       rstarts = "exonStarts", rends = "exonEnds",
#'                       gene = "name", symbol = "name2", 
#'                       transcript = "name", strand = "strand", 
#'                       fill = "#960000", name = "Ensembl Genes")
#' 
#' 
#' gencodeGenes <- UcscTrack(genome = gen, chromosome = chr_str,
#'                       track = "All GENCODE VM25", from = from, to = to,
#'                       trackType = "GeneRegionTrack", 
#'                       rstarts = "exonStarts", rends = "exonEnds",
#'                       gene = "name", symbol = "name2", 
#'                       transcript = "name", strand = "strand", 
#'                       fill = "#960000", name = "Ensembl Genes")
#' 
#' cpgIslands <- UcscTrack(genome = gen, chromosome = chr_str, 
#'                         track = "cpgIslandExt", from = from, to = to,
#'                         trackType = "AnnotationTrack", 
#'                         start = "chromStart", end = "chromEnd", 
#'                         id = "name", shape = "box", fill = "#006400", 
#'                         name = "CpG Islands")
#' 
#' snpLocations <-  UcscTrack(genome = gen, chromosome = chr_str,
#'                            track = "snp128", from = from, to = to,
#'                            trackType = "AnnotationTrack", 
#'                            start = "chromStart", end = "chromEnd", 
#'                            id = "name", feature = "func", 
#'                            strand = "strand", shape = "box", 
#'                            stacking = "dense", fill = "black",
#'                            name = "SNPs")
#' 
#' encodeLocations <- UcscTrack(genome = gen, chromosome = chr_str,
#'                              track = "ENCODE cCREs", from = from, to = to,
#'                              trackType = "AnnotationTrack", 
#'                              start = "chromStart", end = "chromEnd", 
#'                              id = "name", feature = "func", 
#'                              strand = "strand", shape = "box", 
#'                              stacking = "dense", fill = "black",
#'                              name = "cCREs")
#' 
#' # conservation <- UcscTrack(genome = gen, chromosome = chr_str,
#' #                           track = "Conservation", 
#' #                           table = "phyloP30wayPlacental",
#' #                           from = from, to = to, trackType = "DataTrack", 
#' #                           start = "start", end = "end", data = "score",
#' #                           type = "hist", window = "auto", 
#' #                           col.histogram = "darkblue", 
#' #                           fill.histogram = "darkblue",
#' #                           ylim = c(-3.7, 4), name = "Conservation")
#' 
#' # gcContent <- UcscTrack(genome = gen, chromosome = chr_str,
#' #                        track = "GC Percent", table = "gc5Base",
#' #                        from = from, to = to, trackType = "DataTrack", 
#' #                        start = "start", end = "end", data = "score",
#' #                        type = "hist", window = -1, windowSize = 1500, 
#' #                        fill.histogram = "black", col.histogram = "black",
#' #                        ylim = c(30, 70), name = "GC Percent")
#' 
#' axTrack <- GenomeAxisTrack()
#' idxTrack <- IdeogramTrack(genome=gen, chromosome=chr_str)
#' 
#' plotTracks(list(idxTrack, axTrack, knownGenes, refGenes, ensGenes, 
#'                 cpgIslands, gcContent, conservation, snpLocations), 
#'            from = from, to = to, showTitle = FALSE)
#' 
#' 
#' plotTracks(list(idxTrack, axTrack, knownGenes, refGenes, gencodeGenes, 
#'                 cpgIslands, gcContent, snpLocations), 
#'            from = from, to = to, showTitle = TRUE)
#' 
#' plotTracks(list(idxTrack, axTrack, knownGenes, refGenes, 
#'                 cpgIslands, gcContent, snpLocations), 
#'            from = from, to = to, showTitle = FALSE)
