library(data.table)
library(dplyr)

source("/workspace/fasi-domice/src/helpers.R")
source("/workspace/fasi-domice/setup.R")


# GM_Snps meta data
ccre <- fread(paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"), data.table = FALSE)


## topics
topics <- fread(paste0(results_dir, "topics/qtl-topic-lods.csv.gz"), 
                data.table = FALSE) %>%
  dplyr::select(-V1) %>%
  left_join(ccre) %>% 
  mutate(xpos=rank(chr*1e9+pos, ties.method = "first"))


plot_qtl(qtl_output = topics, varcolumn = "topic8", plot_title = "Topic 8")

for(i in 0:19) {
  
  p_out <- plot_qtl(qtl_output = topics, 
           varcolumn = paste0("topic", i), 
           plot_title = paste0("Topic ", i))
  
  ggsave(filename = paste0("qtl_plot_topic", i, ".pdf"),
         plot = p_out,
         width = 7,
         height =5,
         dpi = 330)
  
}
