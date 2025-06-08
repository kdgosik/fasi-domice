#' Figure 3d Supplementary
#' 
#' 


project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"
figure_path <- paste0(project_path, "kirk/results/figures/")
source(paste0(figure_path, "helpers.R"))

## load data
topic_lods <- fread(paste0(my_path, "results/topic-qtl-lods.csv"), data.table = FALSE)
ccre <- fread(paste0(my_path, "data/GM_SNPS_Consequence_cCRE.csv"))
ccre[, recomb_rate := round(cM / (pos / 1e6), 4)]


## topic_qtl #############
file_list <- dir(paste0(my_path, "results/founder-coefs"), pattern = "founder-coef.csv", full.names = TRUE)
# my_colors <- qtl2::CCcolors
# names(my_colors) <- LETTERS[1:8]
founderdf <- data.frame(founder_letter = LETTERS[1:8], 
                        founder_color = qtl2::CCcolors, 
                        founder = names(qtl2::CCcolors),
                        best_clean = c("A_J", "C57BL_6NJ", "129S1_SvImJ", "NOD_ShiLtJ",
                                       "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ"),
                        stringsAsFactors = FALSE)


marker_map <- readRDS(paste0(geno_path, "Regev_map_20171221.rds"))
xpos <- qtl2:::map_to_xpos(marker_map, gap = 25)
xchrbound <- qtl2:::map_to_boundaries(marker_map, gap = 25)

## Make Marker Reference Dataset
MarkerRef <- data.frame(marker = names(xpos), xpos = xpos, stringsAsFactors = FALSE) %>%
  left_join(
    {ccre %>% dplyr::select(marker, chr)}
  )


for( f in file_list ) {
  
  topic <- gsub("topic(.*?)-chr(.*?)-founder-coef.csv","\\1", basename(f))
  chr <- gsub("topic(.*?)-chr(.*?)-founder-coef.csv","\\2", basename(f))
  var_col <- paste0("topic", topic)
  
  coef <- read.csv(paste0(my_path, "results/founder-coefs/", var_col, "-chr", chr, "-founder-coef.csv"))
  plot_lods <- topic_lods %>% left_join(ccre) %>% left_join(MarkerRef)
  chr_markers <- ccre$marker[ccre$chr == chr]
  
  topic_qtl <- ggplot() +
    geom_line(data = plot_lods[plot_lods$marker %in% chr_markers,], 
              mapping = aes_string(x = "xpos", y = var_col)) +
    geom_hline(yintercept = 6, color = "red", linetype = "dashed") + 
    scale_x_continuous(breaks = apply(xchrbound, 2, median), label = c(as.character(1:19), "X")) + 
    labs(x = "", y = "LOD", title = "") +
    theme(
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    theme_minimal()
  
  
  
  p1 <- coef %>%
    dplyr::rename(marker = X) %>%
    left_join(data.frame(marker = names(xpos), xpos = xpos, stringsAsFactors = FALSE)) %>%
    dplyr::select(-intercept) %>%
    pivot_longer(cols = A:H, names_to = "founder_letter") %>%
    left_join(founderdf, by = "founder_letter") %>%
    ggplot(., aes(x = xpos, y = value, color = founder)) +
    geom_line() + 
    scale_color_manual(values = qtl2::CCcolors) +
    scale_x_continuous(breaks = apply(xchrbound, 2, median), label = c(as.character(1:19), "X")) +
    labs(x = "Chromosome of SNP", y = "LOD", title = paste0("Topic ", topic)) +
    theme(
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    theme_minimal()
  
  pout <- p1 / topic_qtl
  
  ggsave(filename = paste0(figure_path, "Supp.Fig.3D,",i,"-qtl-founder-allele-", var_col,"-chr",chr,".pdf"),
         plot = pout,
         width = 8.5,
         height = 5,
         dpi = 330)
  
}
