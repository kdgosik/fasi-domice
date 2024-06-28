#' Figure Supp 3C
#' 
#' 

project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"

source(paste0(project_path, "kirk/results/figures/helpers.R"))


## load data
plot_df <- fread(paste0(my_path, "data/manuscript-plot-data.csv"), data.table = FALSE)

## tSNE ###############################
## Gata3
gata3 <- plot_df %>% 
  dplyr::filter(louvain_labels %in% c(1,2,3,4,5,8)) %>% 
  my_plot(.,"Gata3") +
  my_theme_nolegend + 
  labs(title = "Gata3", x = "", y = "") + 
  theme(plot.title = element_text(size = 24))


ggsave(filename = paste0(figure_path, "figure-2-tsne-gata3.pdf"),
       plot = gata3,
       width = 7,
       height = 7,
       dpi = 330)



## Gata3 BIOCARTA
gata3_biocarta <- plot_df %>% 
  dplyr::filter(louvain_labels %in% c(1,2,3,4,5,8)) %>% 
  my_plot(.,"BIOCARTA_GATA3_PATHWAY") + # , color = "orange") +
  my_theme_nolegend + 
  labs(title = "BIOCARTA GATA3 PATHWAY", x = "", y = "") + 
  theme(plot.title = element_text(size = 24))


ggsave(filename = paste0(figure_path, "figure-2-tsne-gata3-biocarta.pdf"),
       plot = gata3_biocarta,
       width = 7,
       height = 7,
       dpi = 330)


## Topic 11
topic11 <- plot_df %>%
  dplyr::filter(louvain_labels %in% c(1,2,3,4,5,8)) %>%
  my_plot(.,"topic11") + # , color = "orange") +
  my_theme_nolegend +
  labs(title = "Topic 11", x = "", y = "") +
  theme(plot.title = element_text(size = 24))


ggsave(filename = paste0(figure_path, "figure-2-tsne-topic11.pdf"),
       plot = topic11,
       width = 7,
       height = 7,
       dpi = 330)




## original plots #####################################################



UNCHS014430_aprobs <- fread(paste0(my_path, "results/UNCHS014430-aprobs.csv"),
                            col.names = c("BestCall", LETTERS[1:8]))

UNCHS014430_genotypes <-  fread(paste0(my_path, "results/UNCHS014430-genotypes.csv"),
                                col.names = c("BestCall", "UNCHS014430"))

left_join(UNCHS014430_aprobs) %>%
  left_join(UNCHS014430_genotypes) %>%
  
  
  topic11 <- ggplot(df, aes(x = X_tsne1, y = X_tsne2, color = topic11)) + 
  geom_point(shape = 46) + 
  my_theme_nolegened + 
  labs(title = 'topic11')


biocarta_gata3 <-  ggplot(df, aes(x = X_tsne1, y = X_tsne2, color = BIOCARTA_GATA3_PATHWAY)) + 
  geom_point(shape = 46) + 
  my_theme_nolegened + 
  labs(title = 'BIOCARTA_GATA3_PATHWAY')

gata3 <- ggplot(df, aes(x = X_tsne1, y = X_tsne2, color = Gata3)) + 
  geom_point(shape = 46) + 
  my_theme_nolegened + 
  labs(title = 'Gata3')


biocarta_gata3_quant <- plot_quantile(df, "BIOCARTA_GATA3_PATHWAY", "UNCHS014430")
topic11_quant <- plot_quantile(df, "topic11", "UNCHS014430")
gata3_quant <- plot_quantile(df, "Gata3", "UNCHS014430")
il4_quant <- plot_quantile(df, "Il4", "UNCHS014430")
il5_quant <- plot_quantile(df, "Il5", "UNCHS014430")
stim2_quant <- plot_quantile(df, "Stim2", "UNCHS014430")

plot_grid(topic11, gata3, biocarta_gata3,
          topic11_quant, gata3_quant, biocarta_gata3_quant,
          ncol = 3)