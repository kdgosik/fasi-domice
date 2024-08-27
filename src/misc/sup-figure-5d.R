#' Figure 5D Supplementary
#' 
#' 

project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"
figure_path <- paste0(project_path, "kirk/results/figures/")
source(paste0(figure_path, "helpers.R"))

## Read marker map
marker_map <- readRDS(paste0(geno_path, "Regev_map_20171221.rds"))

cytokines <- fread(paste0(my_path, "results/qtl-allergy-cytokines-lods.csv.gz"),
                   data.table = FALSE)

i <- 1
plot_cytokines <- c("TNFa", "IL10", "IFNg","IL6", "IL2","IL9")
# plot_qtl(qtl_output = qtl, marker_map = marker_map, varcolumn = col)
for( col in plot_cytokines ) {
  
  # col <- colnames(cytokines)[i]
  p1 <- plot_qtl(qtl_output = cytokines, marker_map = marker_map, varcolumn = col) + 
    geom_hline(yintercept = 6, col = "red")
  ggsave(filename = paste0(figure_path, "Supp.Fig.5D", i, "-qtl-cytokine-allergy-", col, "-gwa.pdf"),
         plot = p1,
         height = 7,
         width = 7,
         dpi = 330)
  i <- i + 1
  
}


