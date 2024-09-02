#' @title Figure 3H Founder program QTLs
#' @author Kirk Gosik
#' @description
#'


## Set project directory
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
  
  ggsave(filename = paste0(figure_path, "Figure-3H-,",i,"-qtl-founder-allele-", var_col,"-chr",chr,".pdf"),
         plot = pout,
         width = 8.5,
         height = 5,
         dpi = 330)
  
}










## OLD CODE ######################################
# cell_prop_geno_df <- read.csv(paste0(my_path, "results/cell-proportions-qtl-genotypes.csv"))
# plot_df <- fread(paste0(my_path, "data/manuscript-plot-data.csv"), data.table = FALSE)
# 
# 
# ## Rbpj significant marker - UNCHS014402
# ## Stim2 significant marker - UNCHS014413
# ## Fam21 (Washc2) significant marker - UNCHS018751
# ## Rrbp1 significant marker - UNCHS007079
# ## Tmem132c - significant marker - UNCHS015849
# ## Il20ra - significant marker - UNC17535571
# 
# cell_prop_order <- plot_df %>%
#   dplyr::left_join({ 
#     cell_prop_geno_df %>% dplyr::select(BestCall, founder_genotype = UNCHS014402)
#   }) %>%
#   filter(louvain_labels %in% c(1,2,3,4,5,8), !is.na(founder_genotype)) %>% 
#   group_by(founder_genotype, called_cell_types_new) %>% 
#   summarise(n = n()) %>% 
#   mutate(total = sum(n), prop = n / sum(n)) %>%
#   filter(called_cell_types_new == "ILC3(LTi-like)") %>%
#   arrange(prop) %$%
#   # filter(called_cell_types_new == "ILC1") %>%
#   # arrange(desc(prop)) %$%
#   founder_genotype %>%
#   as.character()
# 
# celltype_order <- c("NK", "ILC1", "ILC2", "ILC3", "ILC3(LTi-like)")
# plot_df <- plot_df %>%
#   dplyr::left_join({ 
#     cell_prop_geno_df %>% dplyr::select(BestCall, founder_genotype = UNCHS014402)
#   }) %>%
#   dplyr::mutate(founder_genotype = factor(founder_genotype, levels = cell_prop_order),
#                 called_cell_types_new = factor(called_cell_types_new, levels = celltype_order))
# 
# 
# pout <- plot_df %>% 
#   filter(louvain_labels %in% c(1,2,3,4,5,8), !is.na(founder_genotype)) %>% 
#   ggplot(., aes(founder_genotype, 
#                 color = called_cell_types_new, 
#                 fill = called_cell_types_new)) + 
#   geom_bar(position = "fill") + 
#   theme(axis.text.x = element_text(angle = 90)) +
#   labs(title = "Cell Proportion at Marker UNCHS014402 (Rbpj)")
# 
# ggsave(filename = paste0(figure_path, "figure-3d-barplot-UNCHS014402-proportion.pdf"),
#        plot = pout,
#        dpi = 330,
#        width = 7,
#        height = 5)
# 
