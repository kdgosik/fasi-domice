library(tidyverse)
library(UpSetR)
library(ComplexUpset)



project_dir <- "/Users/kirkgosik/Google Drive/My Drive/projects/"
domice_dir <- paste0(project_dir, "domice/paper of QTL study/Revised materials of ILC-QTL paper cell Science format/")
data_dir <- paste0(domice_dir, "data/")
results_dir <- paste0(data_dir, "results/")


vars <- read_csv(paste0(data_dir, 'allchannels/vars.csv'))

intrinsic <- vars %>%
  dplyr::select(gene_ids, ilc1_expressed, ilc2_expressed, 
                ilc3_expressed, lti_expressed) %>%
  pivot_longer(-gene_ids) %>%
  filter(value == 1) %>%
  distinct(gene_ids) %>%
  mutate(intrinsic = 1)


ccre <- read_csv(paste0(data_dir, 'references/GM_SNPS_Consequence_cCRE.csv'))
ccre_summary_with_ensembl <- ccre %>%
  mutate(pos_floor = as.character(1e6*floor(pos / 1e6))) %>% 
  unite("loci", c(chr, pos_floor)) %>%
  group_by(loci, ensembl_gene) %>%
  summarise(across(contains("qtl"), max, na.rm=T)) %>%
  left_join(dplyr::select(vars, 
                          ensembl_gene = gene_ids, 
                          ilc1_expressed, ilc2_expressed, ilc3_expressed, lti_expressed))


ccre_summary <- ccre %>%
  mutate(pos_floor = as.character(1e6*floor(pos / 1e6))) %>% 
  unite("loci", c(chr, pos_floor)) %>%
  group_by(loci) %>%
  summarise(across(contains("qtl"), max, na.rm=T),
            marker = first(marker),
            gene_ids = first(ensembl_gene)) %>%
  left_join(intrinsic) %>%
  replace_na(list(intrinsic = 0))


## Figure 4A: eQTL x Topic QTL #########################

eqtl_topic_upsetdf <- ccre_summary %>%
  ungroup() %>%
  dplyr::select(starts_with("topic"), contains('eqtl_loci_cv'))

colnames(eqtl_topic_upsetdf) <- c(paste("Topic", 0:19), 
                                  "ILC1 eQTL", "ILC2 eQTL", "ILC3 eQTL", "LTi-like eQTL")


eqtl_topic_upset<- ComplexUpset::upset(data = eqtl_topic_upsetdf, 
                                       intersect = colnames(eqtl_topic_upsetdf),
                                       min_degree = 2,
                                       set_size = FALSE) + 
  ggplot2::labs(title = "eQTL x topic QTL")

## save plot
ggplot2::ggsave(filename =  "results/figures/Figure-4A-eQTL-topic-upset.pdf",
       plot = eqtl_topic_upset,
       width = 7,
       height = 5,
       dpi = 330)

## Figure 4B: eQTL x proportion QTL  #########################

eqtl_proportion_upsetdf <- ccre_summary %>%
  ungroup() %>%
  dplyr::select(contains("prop"), contains('eqtl_loci_cv')) 

colnames(eqtl_proportion_upsetdf) <- c("ILC1 vs ILC2","ILC2 vs ILC3",
                                       "ILC3 vs LTi","ILC1 vs LTi",
                                       "ILC2 vs LTi","ILC1 vs ILC3",
                                       "ILC1 eQTL", "ILC2 eQTL", "ILC3 eQTL", "LTi-like eQTL")


eqtl_proportion_upset<- ComplexUpset::upset(data = eqtl_proportion_upsetdf, 
                                            intersect = colnames(eqtl_proportion_upsetdf),
                                            min_degree = 2,
                                            min_size = 6,
                                            set_size = FALSE) + 
  ggplot2::labs(title = "eQTL x proportion QTL")


## save plot
ggplot2::ggsave(filename =  "results/figures/Figure-4B-eQTL-proportion-upset.pdf",
                plot = eqtl_proportion_upset,
                width = 7,
                height = 5,
                dpi = 330)


## Figure 4C: proportion QTL x Topic QTL #########################

topic_proportion_upsetdf <- ccre_summary %>%
  ungroup() %>%
  dplyr::select(contains("prop"), starts_with('topic')) 

colnames(topic_proportion_upsetdf) <- c("ILC1 vs ILC2","ILC2 vs ILC3",
                                       "ILC3 vs LTi","ILC1 vs LTi",
                                       "ILC2 vs LTi","ILC1 vs ILC3",
                                       paste("Topic", 0:19))


topic_proportion_upset<- ComplexUpset::upset(data = topic_proportion_upsetdf, 
                                            intersect = colnames(topic_proportion_upsetdf),
                                            min_degree = 2,
                                            min_size = 2,
                                            set_size = FALSE) + 
  ggplot2::labs(title = "topic QTL x proportion QTL")


## save plot
ggplot2::ggsave(filename =  "results/figures/Figure-4C-topic-proportion-upset.pdf",
                plot = topic_proportion_upset,
                width = 7,
                height = 5,
                dpi = 330)



## Figure 4H: trait type #########################

all_traits_upsetdf <- ccre_summary %>%
  ungroup() %>%
  rowwise() %>%
  mutate(Topic = as.numeric(sum(c_across(starts_with("topic")))>0),
         Proportion = as.numeric(sum(c_across(contains("prop")))>0),
         eQTL = as.numeric(sum(c_across(ends_with('eqtl_loci_cv')))>0),
         Cytokine = as.numeric(sum(c_across(ends_with('_steady')))>0)) %>%
  dplyr::select(Topic, Proportion, eQTL, Cytokine)

colnames(all_traits_upsetdf) <- c("Topic","Proportion", "eQTL","Cytokine")


all_traits_upset <- ComplexUpset::upset(data = all_traits_upsetdf, 
                                             intersect = colnames(all_traits_upsetdf),
                                             min_degree = 2,
                                             min_size = 2,
                                             set_size = FALSE) + 
  ggplot2::labs(title = "QTL Overlap")

## save plot
ggplot2::ggsave(filename =  "results/figures/Figure-4H-all-traits-upset.pdf",
                plot = all_traits_upsetdf,
                width = 7,
                height = 5,
                dpi = 330)

## Figure 4I: topic QTL x cytokine QTLs #########################

topic_cytokine_upsetdf <- ccre_summary %>%
  ungroup() %>%
  dplyr::select(contains("steady"), starts_with('topic'))

colnames(topic_cytokine_upsetdf) <- c("INFg","IL-5","TNF","IL-2",
                                      "IL-6","IL-4", "IL-10", "IL-9",
                                      "IL-17A", "IL-17F", "IL-22", "IL-13",
                                      paste("Topic", 0:19))


topic_cytokine_upset<- ComplexUpset::upset(data = topic_cytokine_upsetdf, 
                                           intersect = colnames(topic_cytokine_upsetdf),
                                           min_degree = 2,
                                           min_size = 1,
                                           set_size = FALSE) + 
  ggplot2::labs(title = "topic QTL x cytokine QTL")


## save plot
ggplot2::ggsave(filename =  "results/figures/Figure-4I-topic-cytokine-upset.pdf",
                plot = topic_cytokine_upset,
                width = 7,
                height = 5,
                dpi = 330)


## Figure ED-4C: eQTL x cytokine QTL #########################

eqtl_cytokine_upsetdf <- ccre_summary %>%
  ungroup() %>%
  dplyr::select(contains("steady"), contains('eqtl_loci_cv'))

colnames(eqtl_cytokine_upsetdf) <- c("INFg","IL-5","TNF","IL-2",
                                      "IL-6","IL-4", "IL-10", "IL-9",
                                      "IL-17A", "IL-17F", "IL-22", "IL-13",
                                      "ILC1 eQTL", "ILC2 eQTL", "ILC3 eQTL", "LTi-like eQTL")


eqtl_cytokine_upset <- ComplexUpset::upset(data = eqtl_cytokine_upsetdf, 
                                           intersect = colnames(eqtl_cytokine_upsetdf),
                                           min_degree = 2,
                                           min_size = 1,
                                           set_size = FALSE) + 
  ggplot2::labs(title = "eQTL x cytokine QTL")

## save plot
ggplot2::ggsave(filename =  "results/figures/Figure-ED-4C-eQTL-cytokine-upset.pdf",
                plot = eqtl_cytokine_upset,
                width = 7,
                height = 5,
                dpi = 330)


## Figure ED-4D: proportion QTL x cytokine QTL #########################

proportion_cytokine_upsetdf <- ccre_summary %>%
  ungroup() %>%
  dplyr::select(contains("steady"), contains("prop"))

colnames(proportion_cytokine_upsetdf) <- c("INFg","IL-5","TNF","IL-2",
                                           "IL-6","IL-4", "IL-10", "IL-9",
                                           "IL-17A", "IL-17F", "IL-22", "IL-13",
                                           "ILC1 vs ILC2","ILC2 vs ILC3",
                                           "ILC3 vs LTi","ILC1 vs LTi",
                                           "ILC2 vs LTi","ILC1 vs ILC3")


proportion_cytokine_upset<- ComplexUpset::upset(data = proportion_cytokine_upsetdf, 
                                           intersect = colnames(proportion_cytokine_upsetdf),
                                           min_degree = 2,
                                           min_size = 1,
                                           set_size = FALSE)

ggplot2::ggsave(filename =  "results/figures/Figure-ED-4D-proportion-cytokine-upset.pdf",
                plot = proportion_cytokine_upset,
                width = 7,
                height = 5,
                dpi = 330)



## intrinsic #############


## eQTL
eQTL_intrinsic <- ccre %>% 
  ungroup() %>% 
  rowwise() %>%
  mutate(eQTL = as.numeric(sum(c_across(ends_with('eqtl_loci_cv')))>0)) %>% 
  dplyr::select(ensembl_gene, eQTL) %>% 
  group_by(gene_ids = ensembl_gene) %>% 
  summarise(eQTL = max(eQTL)) %>% 
  left_join(intrinsic) %>% 
  count(eQTL, intrinsic)

## total proportion-QTL
sum(all_traits_upsetdf$eQTL)

# mutate(Topic = as.numeric(sum(c_across(starts_with("topic")))>0),
#        Proportion = as.numeric(sum(c_across(contains("prop")))>0),
#        eQTL = as.numeric(sum(c_across(ends_with('eqtl_loci_cv')))>0),

# Proportion QTL
proportion_intrinsic <- ccre %>% 
  ungroup() %>% 
  rowwise() %>%
  mutate(Proportion = as.numeric(sum(c_across(contains("prop")))>0)) %>% 
  dplyr::select(ensembl_gene, Proportion) %>% 
  group_by(gene_ids = ensembl_gene) %>% 
  summarise(Proportion = max(Proportion)) %>% 
  left_join(intrinsic) %>% 
  count(Proportion, intrinsic)

## total proportion-QTL
sum(all_traits_upsetdf$Proportion)


## topic
topic_intrinsic <- ccre %>% 
  ungroup() %>% 
  rowwise() %>%
  mutate(Topic = as.numeric(sum(c_across(starts_with("topic")))>0)) %>% 
  dplyr::select(ensembl_gene, Topic) %>% 
  group_by(gene_ids = ensembl_gene) %>% 
  summarise(Topic = max(Topic)) %>% 
  left_join(intrinsic) %>% 
  count(Topic, intrinsic)
 
## total topic
sum(rle(all_traits_upsetdf$Topic)$values == 1)
# 18

## xxArchive #############

ilc1_eqtl_loci_by_gene <- read_csv(paste0(results_dir, "eqtl/qtl-loci-by-genes-ILC1.csv"))
ilc1_qtl <- read_csv(paste0(results_dir, "eqtl/qtl-plot-lods-ILC1-cv.csv.gz"))
ilc1_qtl %>% 
  filter(value > 10) %>% 
  left_join(dplyr::select(ccre, marker, chr, pos)) %>% 
  mutate(pos_floor = as.character(1e6*floor(pos / 1e6))) %>%
  unite("loci", c(chr, pos_floor)) %>% 
  group_by(gene, loci) %>% 
  summarise(value = max(value)) %>% 
  count(gene) %>% 
  ungroup() %>% 
  summarise(mean(n))

ilc2_eqtl_loci_by_gene <- read_csv(paste0(results_dir, "eqtl/qtl-loci-by-genes-ILC2.csv"))
ilc2_qtl <- read_csv(paste0(results_dir, "eqtl/qtl-plot-lods-ILC2-cv.csv.gz"))
ilc2_qtl %>% 
  filter(value > 6) %>% 
  left_join(dplyr::select(ccre, marker, chr, pos)) %>% 
  mutate(pos_floor = as.character(1e6*floor(pos / 1e6))) %>%
  unite("loci", c(chr, pos_floor)) %>% 
  group_by(gene, loci) %>% 
  summarise(value = max(value)) %>% 
  count(gene) %>% 
  ungroup() %>% 
  summarise(mean(n))


ilc3_eqtl_loci_by_gene <- read_csv(paste0(results_dir, "eqtl/qtl-loci-by-genes-ILC3.csv"))
ilc3_qtl <- read_csv(paste0(results_dir, "eqtl/qtl-plot-lods-NCR1+ ILC3-cv.csv.gz"))
ilc3_qtl %>% 
  filter(value > 6) %>% 
  left_join(dplyr::select(ccre, marker, chr, pos)) %>% 
  mutate(pos_floor = as.character(1e6*floor(pos / 1e6))) %>%
  unite("loci", c(chr, pos_floor)) %>% 
  group_by(gene, loci) %>% 
  summarise(value = max(value)) %>% 
  count(gene) %>% 
  ungroup() %>% 
  summarise(mean(n))


lti_eqtl_loci_by_gene <- read_csv(paste0(results_dir, "eqtl/qtl-loci-by-genes-LTi.csv"))
lti_qtl <- read_csv(paste0(results_dir, "eqtl/qtl-plot-lods-Lti ILC3-cv.csv.gz"))
lti_qtl %>% 
  filter(value > 6) %>% 
  left_join(dplyr::select(ccre, marker, chr, pos)) %>% 
  mutate(pos_floor = as.character(1e6*floor(pos / 1e6))) %>%
  unite("loci", c(chr, pos_floor)) %>% 
  group_by(gene, loci) %>% 
  summarise(value = max(value)) %>% 
  count(gene) %>% 
  ungroup() %>% 
  summarise(mean(n))





ilc1_ilc2 <- read_csv(paste0(results_dir, "proportion/gwas-ilc1-ilc2-results.csv.gz"))
ilc1_ilc3 <- read_csv(paste0(results_dir, "proportion/gwas-ilc1-ilc3-results.csv.gz"))
ilc1_lti <- read_csv(paste0(results_dir, "proportion/gwas-ilc1-lti-results.csv.gz"))
ilc2_ilc3 <- read_csv(paste0(results_dir, "proportion/gwas-ilc2-ilc3-results.csv.gz"))
ilc2_lti <- read_csv(paste0(results_dir, "proportion/gwas-ilc2-lti-results.csv.gz"))
ilc3_lti <- read_csv(paste0(results_dir, "proportion/gwas-ilc3-lti-results.csv.gz"))

proportion_qtl_loci_genes_by_trait <- read_csv(paste0(results_dir, "proportion/proportion-qtl-loci-genes-by-trait.csv"))
proportion_qtl_loci_genes_by_trait %>% mutate(pos_floor = 1e6*floor(pos / 1e6)) %>% unite("loci", c(chr,pos_floor)) %>% pull(loci) %>% n_distinct


topics <- read_csv(paste0(results_dir, "topics/topic-qtl-lods.csv.gz"))

cytokines <- read_csv(paste0(results_dir, "cytokine/qtl-cytokines-steady-lods.csv.gz"))
cytokines_long <- cytokines %>% 
  dplyr::select(-`...1`) %>% 
  pivot_longer(-marker) %>% 
  filter(value > 6) %>% 
  left_join(dplyr::select(ccre, marker, chr, pos)) %>%
  unite('loci', c(chr, pos))



ccre_summary <- ccre_summary %>%
  rowwise() %>%
  mutate(topic_qtl = as.numeric(sum(c_across(starts_with("topic")))>0),
         prop_qtl = as.numeric(sum(c_across(contains("prop")))>0),
         eqtl = as.numeric(sum(c_across(ends_with('eqtl_loci_cv')))>0),
         cytokine_qtl = as.numeric(sum(c_across(ends_with('_steady')))>0))


upsetdf <- ccre_summary %>%
  dplyr::select(topic_qtl, prop_qtl, eqtl, cytokine_qtl)

upset(upsetdf)




## page 16 polygenic QTLs

df1 %>% filter(str_detect(trait,"allergy")) %>% summarise(mean(V1))
df1 %>% filter(str_detect(trait,"prop|activate"))  %>% summarise(mean(V1))
df1 %>% filter(str_detect(trait,"topic")) %>% summarise(mean(V1))
df1 %>% filter(str_detect(trait,"eqtl_loci_cv")) %>% summarise(mean(V1))

