library(ggpubr)

project_path <-"./domice/"
figure_path <- paste0("./results/figures/")
data_path <- paste0(project_path, "data/")
source(paste0(figure_path, "helpers.R"))
setwd(project_path)

lods <- fread("./data/qtl-steady-cytokines-lods.csv.gz")
load("./data/GM_snps.Rdata")
marker_map <- sapply(unique(GM_snps$chr)[1:20], function(x) {
  structure(GM_snps$pos[GM_snps$chr==x], names = GM_snps$marker[GM_snps$chr==x])
}, USE.NAMES = TRUE, simplify = FALSE)


plot_qtl(lods, marker_map, "IL5")


plot_qtl <- function(qtl_output, marker_map, varcolumn ) {
  
  xpos <- qtl2:::map_to_xpos(marker_map, gap = 25)
  xchrbound <- qtl2:::map_to_boundaries(marker_map, gap = 25)
  chr_breaks <- apply(xchrbound, 2, median)

  rectangles <- data.frame(
    xmin = apply(xchrbound, 2, min),
    xmax = apply(xchrbound, 2, max),
    ymin = 0,
    ymax = max(qtl_output[[varcolumn]])
    )[seq(2,20,2), ]


  plot_df <- as.data.frame(qtl_output) %>%
    left_join({
      data.frame(marker = names(xpos),
                 xpos = xpos)
    })

  plot_out <- ggplot() + 
    geom_rect(data = rectangles, 
              mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
              fill='gray80', alpha=0.8) +
    geom_line(data = plot_df, 
            mapping = aes_string(x = "xpos", y = varcolumn), shape = 46) +
    geom_hline(yintercept = 6, color = "red", linetype = "dashed") + 
  # geom_label(data = markers_df, 
  #            aes(x = xpos, y = plot_lods, label = marker), 
  #            size = 2) +
  # geom_label_repel(data = snp_lables, 
  #                  aes(label=as.factor(marker), alpha = 0.7), 
  #                  size = 2, force = 1.3)
    scale_x_continuous(breaks = chr_breaks, label = c(as.character(1:19), "X")) + 
    labs(x = "Chromosome of SNP", 
         y = "LOD", 
         title = varcolumn) +
    theme(
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  plot_out

}


i <- 1
for( a in colnames(lods)[3:14]) {
  cat("Using cytokine: ", a, "\n")
  p <- plot_qtl(lods, marker_map = marker_map, varcolumn = a) + 
    geom_hline(y_intercept = 6, color = "red", linetype = "dashed")
  
  ggsave(filename = paste0("./results/figures/ED4A-", i, "-", a, "-cytokine-qtl.pdf"),
         plot = p,
         dpi = 330,
         height = 5,
         width = 7)
  
  i <- i + 1
  
}
