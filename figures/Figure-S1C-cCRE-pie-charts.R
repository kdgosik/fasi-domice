#' @title Supplementary Figure 1C cCRE Pie Charts
#' @author Kirk Gosik
#' @description
#'
#'

library(ggpubr)
library(tidyverse)
library(patchwork)


## setup ##################
figure_path <- paste0("./results/figures/")
data_path <- "data/"
source(paste0(figure_path, "helpers.R"))
setwd(project_path)

ccre <- fread(paste0(data_path, "GM_SNPS_Consequence_cCRE.csv"), data.table = FALSE) %>%
  mutate(ccre_type = ifelse(type.y == "CTCF-only,CTCF-bound", "CTCF-bound",
                            ifelse(str_detect(type.y, "ELS"), "Enhancer",
                                              ifelse(str_detect(type.y, "PLS"), "Promoter", "DNase-H3K4me3"))))

ccre <- fread(paste0(data_path, "GM_SNPS_Consequence_cCRE.csv"), data.table = FALSE) %>%
  mutate(ccre_type = ifelse(str_detect(type.y, "ELS"), "Enhancer",
                                   ifelse(str_detect(type.y, "PLS"), "Promoter", NA)))





## Pie ggplot function ##########
#Create a custom color scale
# library(RColorBrewer)
myColors <- brewer.pal(3,"Set2")
names(myColors) <- unique(na.omit(ccre$ccre_type))
# colScale <- scale_colour_manual(name = "grp",values = myColors)
pie_ggplot <- function( data ) {
  
  ggplot(data, aes(x="", y = value, fill = group)) + 
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_colour_manual(name = "group", values = myColors, aesthetics = "fill") +
    geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), 
                  label = paste0(prop, "%")), size=5)
    #geom_text(aes(y = ypos, label = group), color = "white", size = 6)
    #scale_fill_brewer(palette = "Set2", values = myColors[1:2])
  
}





# # Create Data ILC1 ccre
ilc1_eqtl_ccre <- table(ccre$ccre_type[ccre$ilc1_eqtl_loci_mean==1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "ILC1",
                region = "all")

ilc1_ccre_plot <- pie_ggplot(ilc1_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-ilc1-eqtl.pdf"), 
       plot = ilc1_ccre_plot, 
       dpi = 330)


## cis-eQTL
ilc1_cis_eqtl_ccre <- table(ccre$ccre_type[ccre$ilc1_ciseqtl_assoc_mean == 1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "ILC1",
                region = "cis" )

ilc1_cis_ccre_plot <- pie_ggplot(ilc1_cis_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-ilc1-cis-eqtl.pdf"), 
       plot = ilc1_cis_ccre_plot,
       dpi = 330)


## trans-eQTL
ilc1_trans_eqtl_ccre <- table(ccre$ccre_type[ccre$ilc1_ciseqtl_assoc_mean == 0 & ccre$ilc1_eqtl_loci_mean == 1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "ILC1",
                region = "trans" )

ilc1_trans_ccre_plot <- pie_ggplot(ilc1_trans_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-ilc1-trans-eqtl.pdf"),
       plot = ilc1_trans_ccre_plot,
       dpi = 320)



# # Create Data ILC2 ccre
ilc2_eqtl_ccre <- table(ccre$ccre_type[ccre$ilc2_eqtl_loci_mean==1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "ILC2",
                region = "all" )

ilc2_ccre_plot <- pie_ggplot(ilc2_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-ilc2-eqtl.pdf"), 
       plot = ilc2_ccre_plot,
       dpi = 320)


## cis-eQTL
ilc2_cis_eqtl_ccre <- table(ccre$ccre_type[ccre$ilc2_ciseqtl_assoc_mean == 1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "ILC2",
                region = "cis" )

ilc2_cis_ccre_plot <- pie_ggplot(ilc2_cis_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-ilc2-cis-eqtl.pdf"), 
       plot = ilc2_cis_ccre_plot,
       dpi = 320)


## trans-eQTL
ilc2_trans_eqtl_ccre <- table(ccre$ccre_type[ccre$ilc2_ciseqtl_assoc_mean == 0 & ccre$ilc2_eqtl_loci_mean == 1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "ILC2",
                region = "trans" )

ilc2_trans_ccre_plot <- pie_ggplot(ilc2_trans_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-ilc2-trans-eqtl.pdf"), 
       plot = ilc2_trans_ccre_plot,
       dpi = 320)


# # Create Data ILC3 ccre
ilc3_eqtl_ccre <- table(ccre$ccre_type[ccre$ilc3_eqtl_loci_mean==1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "ILC3",
                region = "all" )

ilc3_ccre_plot <- pie_ggplot(ilc3_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-ilc3-eqtl.pdf"), 
       plot = ilc3_ccre_plot,
       dpi = 320)


## cis-eQTL
ilc3_cis_eqtl_ccre <- table(ccre$ccre_type[ccre$ilc3_ciseqtl_assoc_mean == 1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "ILC3",
                region = "cis" )

ilc3_cis_ccre_plot <- pie_ggplot(ilc3_cis_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-ilc3-cis-eqtl.pdf"), 
       plot = ilc3_cis_ccre_plot,
       dpi = 320)


## trans-eQTL
ilc3_trans_eqtl_ccre <- table(ccre$ccre_type[ccre$ilc3_ciseqtl_assoc_mean == 0 & ccre$ilc3_eqtl_loci_mean == 1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "ILC3",
                region = "trans" )

ilc3_trans_ccre_plot <- pie_ggplot(ilc3_trans_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-ilc3-trans-eqtl.pdf"), 
       plot = ilc3_trans_ccre_plot,
       dpi = 320)


# # Create Data LTi ccre
lti_eqtl_ccre <- table(ccre$ccre_type[ccre$lti_eqtl_loci_mean==1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "LTi",
                region = "all" )

lti_ccre_plot <- pie_ggplot(lti_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-lti-eqtl.pdf"), 
       plot = lti_ccre_plot,
       dpi = 320)


## cis-eQTL
lti_cis_eqtl_ccre <- table(ccre$ccre_type[ccre$lti_ciseqtl_assoc_mean == 1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "LTi",
                region = "cis" )

lti_cis_ccre_plot <- pie_ggplot(lti_cis_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-lti-cis-eqtl.pdf"), 
       plot = lti_cis_ccre_plot,
       dpi = 320)


## trans-eQTL
lti_trans_eqtl_ccre <- table(ccre$ccre_type[ccre$lti_ciseqtl_assoc_mean == 0 & ccre$lti_eqtl_loci_mean == 1]) %>%
  as.data.frame() %>%
  dplyr::rename(value = Freq, group = Var1) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(prop = round(value / sum(value) * 100, 1)) %>%
  dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop,
                cell_type = "LTi",
                region = "trans" )

lti_trans_ccre_plot <- pie_ggplot(lti_trans_eqtl_ccre)
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-lti-trans-eqtl.pdf"),
       plot = lti_trans_ccre_plot,
       dpi = 320)




cis_row <- ilc1_cis_ccre_plot | ilc2_cis_ccre_plot | ilc3_cis_ccre_plot | lti_cis_ccre_plot
trans_row <- ilc1_trans_ccre_plot | ilc2_trans_ccre_plot | ilc3_trans_ccre_plot | lti_trans_ccre_plot


library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

p1 <-   gridExtra::grid.arrange(ilc1_cis_ccre_plot + theme(legend.position = "none") + labs(title = "ILC1"), 
                          ilc2_cis_ccre_plot + theme(legend.position = "none") + labs(title = "ILC2"), 
                          ilc3_cis_ccre_plot + theme(legend.position = "none") + labs(title = "ILC3"), 
                          lti_cis_ccre_plot + theme(legend.position = "none") + labs(title = "LTi"), 
                          nrow = 1, left = "cis")
p2 <- gridExtra::grid.arrange(ilc1_trans_ccre_plot + theme(legend.position = "none"), 
                          ilc2_trans_ccre_plot + theme(legend.position = "none"), 
                          ilc3_trans_ccre_plot + theme(legend.position = "none"), 
                          lti_trans_ccre_plot + theme(legend.position = "none"), 
                          nrow = 1, left = "trans")



legend <- get_legend(ilc1_cis_ccre_plot)
# 4. Arrange ggplot2 graphs with a specific width
grid.arrange(p1, p2, legend, ncol =3, widths=c(2.3, 2.3, 0.8))


pout <- grid.arrange(p1, p2, legend,
             widths = c(2, 1),
             layout_matrix = rbind(c(1, 3),
                                   c(2, NA)))
ggsave(filename = paste0(my_path, "results/figures/ED1F-pie-ccre-cis-trans-eqtl.pdf"),
       plot = pout,
       dpi = 320)



plotdf <- bind_rows(list(ilc1_cis_eqtl_ccre, ilc1_trans_eqtl_ccre,
                         ilc2_cis_eqtl_ccre, ilc2_trans_eqtl_ccre,
                         ilc3_cis_eqtl_ccre, ilc3_trans_eqtl_ccre,
                         lti_cis_eqtl_ccre, lti_trans_eqtl_ccre))

ggplot(plotdf, aes(x="ILC1", y = ypos, fill = group)) + 
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_colour_manual(name = "group", values = myColors, aesthetics = "fill") +
  geom_text(aes(y = prop/2 + c(0, cumsum(prop)[-length(value)]), 
                label = paste0(prop, "%")), size=5) +
  facet_grid(region ~ cell_type)

