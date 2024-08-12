#' @title Figure 1D eQTL map
#' @author Kirk Gosik
#' 
#' @description Cell type specific plots of eQTL by eGene genomic location
#'
#'


library(data.table)
library(tidyverse)

figure_path <- "results/figures/"
data_path <- "data/"
geno_path <- paste0(data_path, "genotype/")
source(paste0(figure_path, "helpers.R"))


## Plot file list
plot_files <- c(
  paste0("results/qtl-plot-lods-ILC1-cv.csv.gz"),
  paste0("results/qtl-plot-lods-ILC2-cv.csv.gz"),
  paste0("results/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz"),
  paste0("results/qtl-plot-lods-Lti ILC3-cv.csv.gz")
)


## Loading Marker Map
cat("Loading SNPs...", "\n")
# load(paste0(geno_path, "GM_snps.Rdata"))
GM_snps <- read_csv(paste0(data_path, "references/GM_SNPS_Consequence_cCRE.csv"))
marker_map <- readRDS(paste0(geno_path, "Regev_map_20171221.rds"))
xpos <- qtl2:::map_to_xpos(marker_map, gap = 25)
xchrbound <- qtl2:::map_to_boundaries(marker_map, gap = 25)

cat("Loading Genes...", "\n")
# ensembl2 <- fread(paste0(my_path, "data/reference/ensembl.Mus_musculus.GRCm38.93.csv"))
# ensembl2 <- ensembl2[chr %in% c(as.character(1:19), "X"), ]
# geneids <- fread(paste0(my_path, "data/reference/GeneID.csv"), col.names = c("symbol", "gene"))

ensembl <- as.data.frame(readGFF(paste0(data_path, "references/Mus_musculus.GRCm38.102.gtf"))) %>%
  filter(seqid %in% c(as.character(1:19), "X"), gene_biotype == "protein_coding") %>%
  mutate(chr = paste0("chr", seqid))

geneids <- ensembl %>% dplyr::select(symbol = gene_name, gene = gene_id)


## Creating Gene Map
gene_map <- split(ensembl, f = ensembl$chr)
for( nm in names(gene_map) ) {
  vec <- gene_map[[nm]]$start
  names(vec) <- gene_map[[nm]]$gene
  gene_map[[nm]] <- vec
}
gene_map <- gene_map[c(as.character(1:19), "X")]
ypos <- qtl2:::map_to_xpos(gene_map, gap = 25)
ychrbound <- qtl2:::map_to_boundaries(gene_map, gap = 25)

i <- 1
out <- list()
out_lods <- list()
df <- NULL
plot_title <- c("ILC1", "ILC2", "LTi", "ILC3")
## CoefVar: c(36:38,40)
## Mean: c(65:67,69)
for( f in plot_files ) {
  
  cat('Loading file: ', f, '\n')
  lods <- fread(f)
  
  cat('Identifying Cis-effects \n')
  lods <- lods %>%
    left_join(
      {GM_snps %>% 
          mutate(marker_pos = pos) %>% 
          select(marker = marker, 
                 marker_pos, 
                 conseq_clean,
                 # type,
                 ensembl_gene)}
    ) %>%
    left_join(
      {ensembl %>% dplyr::select(gene_start = start, gene_end = end, gene = gene_id)}
    ) %>%
    mutate(cis_effect = as.numeric(marker_chr == gene_chr & 
    {marker_pos > (gene_start - 1e6) & marker_pos < (gene_end + 1e6)})) %>%
    ## treating cis and trans effects separately
    mutate(value_cis = ifelse(cis_effect==1, value + 6, value), 
           cell_type = plot_title[i]) %>%
    left_join(geneids) %>%
    filter(value_cis > 10)
  
  cat('Plot File: ', plot_title[i], '\n')
  out_lods[[i]] <- lods
  i <- i + 1
  
}

out_df <- do.call(rbind, out_lods)
table(out_df$symbol[out_df$cis_effect==1], out_df$cell_type[out_df$cis_effect==1])


cat("Plotting LODs...", "\n")
p <- ggplot(out_lods[[1]], aes(x = xpos, y = ypos, color = cis_effect)) +
  # labs(x = "Chromosome of SNP", y = "Chromosome of Gene") +
  labs(x = "Chromosome of SNP", y = "Chromosome of Gene", title = plot_title[i]) + # title = paste0(plot_title, " (", data_subset[i], ")")) +
  geom_point(shape = 20) +
  scale_x_continuous(breaks = apply(xchrbound, 2, median), label = c(as.character(1:19), "X")) + 
  scale_y_continuous(breaks = apply(ychrbound, 2, median), label = c(as.character(1:19), "X")) + 
  scale_color_continuous(low="lightgrey", high="black") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")



# # conseq_clean, 
# p <- ggplot(lods, aes(x = xpos, y = ypos, color=ctcf)) +
#   labs(x = "Chromosome of SNP", y = "Chromosome of Gene", title = paste0(plot_title, " (", data_subset, ") ", geneset_name)) +
#   # geom_point(shape = c(46, 8)[(lods$geneset=="Yes") + 1]) +
#   geom_point(shape = 20) +
#   scale_x_continuous(breaks = apply(xchrbound, 2, median), label = c(as.character(1:19), "X")) +
#   scale_y_continuous(breaks = apply(ychrbound, 2, median), label = c(as.character(1:19), "X")) +
#   theme_bw() +
#   theme(panel.border = element_blank(), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.line = element_line(colour = "black"))




