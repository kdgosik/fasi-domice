library(data.table)
library(dplyr)

source("/workspace/fasi-domice/src/helpers.R")
source("/workspace/fasi-domice/setup.R")


# GM_Snps meta data
ccre <- fread(paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"), data.table = FALSE)


## cytokines
cytokines <- fread(paste0(results_dir, "cytokines/qtl-cytokines-steady-lods.csv.gz"),
                   data.table = FALSE) %>%
  dplyr::select(-V1) %>%
  left_join(ccre) %>% 
  mutate(xpos=rank(chr*1e9+pos, ties.method = "first"))


plot_qtl(qtl_output = cytokines, varcolumn = "IL4", plot_title = "IL4")

for(i in colnames(cytokines)[2:13]) {
  
  p_out <- plot_qtl(qtl_output = cytokines, 
                    varcolumn = i, 
                    plot_title = i)
  
  ggsave(filename = paste0(results_dir, "figures/qtl_plot_cytokine_steady_", i, ".pdf"),
         plot = p_out,
         width = 7,
         height =5,
         dpi = 330)
  
}