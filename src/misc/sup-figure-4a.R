#' Figure 4A Supplementary
#' 
#' 

project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"
figure_path <- paste0(project_path, "kirk/results/figures/")
source(paste0(figure_path, "helpers.R"))


## Read marker map
marker_map <- readRDS(paste0(geno_path, "Regev_map_20171221.rds"))

## ILC Cell Proportions ###############
plot_files <- c("ILC3_stressed_vs_non_qtl.csv.gz", 
                "LTi_stressed_vs_non_qtl.csv.gz")
names(plot_files) <- c("RORgt+ILC3 /RORgt-ILC3",  "CCR6+LTi /CCR6-LTi")

i <- 1
plot_list <- Map(function(f, n) {
  
  plot_out <- fread(paste0(my_path, "results/", f)) %>%
    dplyr::mutate(
      padj = p.adjust(p_values, method = "hochberg"),
      lods = -log10(p.adjust(p_values, method = "hochberg")) / 10
    ) %>%
    gwas_plot(., marker_map = marker_map, plot_title = n) + 
    geom_hline(yintercept = 6, color = "red", linetype = "dashed")
  
  ggsave(filename = paste0(figure_path, "Supp.Fig.4A-", gsub(".csv.gz",".pdf", f)),
         plot = plot_out,
         height = 7,
         width = 7,
         dpi = 330)
  
  i <- i + 1
  
}, plot_files, names(plot_files))


