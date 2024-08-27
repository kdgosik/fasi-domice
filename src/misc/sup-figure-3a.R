#' Figure 3A Supplementary
#' 
#' 

project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"

figure_path <- paste0(project_path, "kirk/results/figures/")
source(paste0(figure_path, "helpers.R"))

plot_df <- fread(paste0(my_path, "data/manuscript-plot-data.csv.gz"), data.table = FALSE)


keep_topics <- c(1,2,4,5,6,7,9,10,13,14,15,17,18,19)
i <- 1
for( t in paste0("topic", keep_topics) ) {
  
  plot_out <- plot_df %>%
    dplyr::filter(louvain_labels %in% c(1, 2, 3, 4, 5, 8)) %>%
    ggplot(., aes_string("X_tsne1", "X_tsne2", color = t)) + 
    geom_point(shape = 46) + 
    scale_color_gradient(low = "lightgrey", high = "darkblue") +
    my_theme
  
  ggsave(filename = paste0(figure_path, "Supp.Fig.3A", i, "-tsne-", t, ".pdf"), 
         plot = plot_out,
         height = 7,
         width = 7,
         dpi = 330)
  
  i <- i + 1
  
}

















## OLD CODE ################
# topic_geno_df <- read.csv(paste0(my_path, "results/topic-qtl-genotypes.csv"))
# 
# plot_df <- plot_df %>%
#   left_join({
#     topic_geno_df
#   })
# 
# 
# 
# for( t in paste0("topic", 0:19) ) {
#   
#   topic_markers <- grep(t, colnames(topic_geno_df), value = TRUE)
#   
#   if( !is.null(topic_markers) ) {
#     
#     for( m in topic_markers ) {
#       plot_out <- plot_df %>%
#         dplyr::filter(!is.na(.data[[ m ]]),
#                       !grepl("Y", .data[[ m ]]),
#                       louvain_labels %in% c(1,2,3,4,5,8)) %>%
#         ggplot(., aes_string("X_tsne1", "X_tsne2", color = t)) + 
#         geom_point(shape = 46) + 
#         scale_color_gradient(low = "lightgrey", high = "darkblue") +
#         my_theme +
#         facet_wrap( m )
#       
#       ggsave(filename = paste0(figure_path, "figure-3a-tsne-ilcs-bygenotype-", m, ".pdf"), 
#              plot = plot_out,
#              height = 7,
#              width = 7,
#              dpi = 330)
#     }
#     
#   }
#   
# }
# 
# # # reorder is close to order, but is made to change the order of the factor levels.
# 
# for( t in paste0("topic", 0:19) ) {
#   
#   topic_markers <- grep(t, colnames(topic_geno_df), value = TRUE)
#   
#   if( !is.null(topic_markers) ) {
#     
#     for( m in topic_markers ) {
#       
#       genotype_order <- plot_df %>%
#         dplyr::filter(!is.na(.data[[ m ]]),
#                       !grepl("Y", .data[[ m ]])) %>%
#         dplyr::group_by( BestCall, .data[[ m ]] ) %>%
#         dplyr::summarise( topic = mean(.data[[ t ]] )) %>% 
#         dplyr::group_by( .data[[ m ]] ) %>% 
#         dplyr::summarise( topic = max(topic) ) %>% 
#         dplyr::arrange(topic) %>% 
#         .[[ m ]] %>%
#         as.character
#       
#       plot_df <- plot_df %>%
#         dplyr::mutate(founder_genotype = factor(.data[[ m ]], levels = genotype_order))
#       
#       plot_out <- plot_df %>%
#         dplyr::filter(!is.na(.data[[ m ]]),
#                       !grepl("Y", .data[[ m ]])) %>%
#         dplyr::group_by(BestCall, founder_genotype) %>% 
#         dplyr::summarise(topic_score = mean(.data[[ t ]])) %>%
#         ggplot(., aes(founder_genotype, topic_score)) + 
#         geom_point() +
#         theme(
#           panel.grid.minor=element_blank()
#         ) +
#         labs(title = paste(toupper(t), "at marker", strsplit(m, "_")[[1]][2]))
#       
#       ggsave(filename = paste0(figure_path, "figure-3a-scatter-topicscore-bygenotype-", m, ".pdf"), 
#              plot = plot_out, 
#              height = 7,
#              width = 7,
#              dpi = 330)
#       
#     }
#     
#   }
#   
# }













# 
# plot_df <- fread(paste0(my_path, "data/manuscript-plot-data.csv"), data.table = FALSE)
# 
# ccre <- fread(paste0(my_path, "data/GM_SNPS_Consequence_cCRE.csv"))
# ccre[, recomb_rate := round(cM / (pos / 1e6), 4)]
# 
# ## look at cCRE and cell_prop_geno_df to make plots for different proportions
# 
# # ccre
# # plot_df
# 
# ## get significant markers here
# # prop_markers <- unique(c(ccre$marker[ccre$ilc1_ilc2_prop_qtl == 1],
# #                   ccre$marker[ccre$ilc1_ilc3_prop_qtl == 1],
# #                   ccre$marker[ccre$ilc1_lti_prop_qtl == 1],
# #                   ccre$marker[ccre$ilc2_ilc3_prop_qtl == 1],
# #                   ccre$marker[ccre$ilc2_lti_prop_qtl == 1],
# #                   ccre$marker[ccre$ilc3_lti_prop_qtl == 1]))
# 
# cell_prop_geno_df <- read.csv(paste0(my_path, "results/cell-proportions-qtl-genotypes.csv"))
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
#     cell_prop_geno_df %>% dplyr::select(BestCall, founder_genotype = UNC17535571)
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
#     cell_prop_geno_df %>% dplyr::select(BestCall, founder_genotype = UNC17535571)
#   }) %>%
#   dplyr::mutate(founder_genotype = factor(founder_genotype, levels = cell_prop_order),
#                 called_cell_types_new = factor(called_cell_types_new, levels = celltype_order))
# 
# 
# plot_df %>% 
#   filter(louvain_labels %in% c(1,2,3,4,5,8), !is.na(founder_genotype)) %>% 
#   ggplot(., aes(founder_genotype, 
#                 color = called_cell_types_new, 
#                 fill = called_cell_types_new)) + 
#   geom_bar(position = "fill") + 
#   theme(axis.text.x = element_text(angle = 90)) +
#   labs(title = "Cell Proportion at Marker UNC17535571 (IL-20ra)")