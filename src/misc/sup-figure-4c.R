#' Figure 4C Supplementary
#' 
#' 

project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"
figure_path <- paste0(project_path, "kirk/results/figures/")
source(paste0(figure_path, "helpers.R"))



library(UpSetR)

vars <- fread(paste0(my_path, "data/scanpy/allchannels/vars.csv"), data.table = FALSE)


upset(data = vars, 
      sets = c("ilc1_ilc2_prop_qtl_genes", "ilc1_ilc3_prop_qtl_genes", 
               "ilc1_lti_prop_qtl_genes", "ilc2_ilc3_prop_qtl_genes", 
               "ilc2_lti_prop_qtl_genes", "ilc3_lti_prop_qtl_genes"), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = "on")

