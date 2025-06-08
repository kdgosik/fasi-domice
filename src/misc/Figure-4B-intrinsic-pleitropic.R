library(data.table)
library(tidyverse)
# library(ggvenn)
# library(igraph)

project_path <- my_path <-"./domice/"
project_path <-"~/Documents/projects/domice/"
figure_path <- paste0(project_path, "results/figures/")
data_path <- paste0(project_path, "data/")
source(paste0(figure_path, "helpers.R"))


vars <- fread(paste0(data_path, "vars.csv"), data.table = F)
ccre <- fread(paste0(data_path, "GM_SNPS_Consequence_cCRE.csv"), data.table = F)
colnames(ccre)[76] <- "ilc1_ilc3_prop_qtl"

## add polygenic vs monogenic feature? ##################


## eQTLs #######

## ILC1 #########################

# ILC1 eGenes but loci gene is not expresssed in ILC1

# # getting all egenes assocatied with ILC1
# ilc1_egenes <- vars %>% filter(ilc1_egenes_cv == 1) %>% pull(index)
# 
# # reading in LOD scores
# ilc1_lods <- fread(paste0(project_path, "results/qtl-plot-lods-ILC1-cv.csv.gz"), data.table = FALSE) 
# 
# # Counting number of QTLs per eGene
# ilc1_polygenic <- ilc1_lods %>% 
#   filter(gene %in% ilc1_egenes) %>%
#   mutate(xpos_floor = floor(xpos)) %>% 
#   group_by(gene, marker_chr, xpos_floor) %>% 
#   summarise(lod = max(value)) %>%
#   group_by(gene) %>%
#   summarise(ilc1_polygenic_egene = as.numeric(sum(lod > 5)>1)) %>%
#   dplyr::rename(index = gene)


ilc1_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(ilc1_pleitropic = sum(ilc1_eqtl_loci_cv)) %>% 
  filter(ilc1_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(ilc1_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


ilc1 <- vars %>% 
  left_join(ilc1_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(ilc1_expressed==1, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

ilc1_summary <- ilc1 %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>% 
  mutate(cell_type = "ILC1")


ilc1_plot <- ggplot(ilc1_summary, aes(intrinsic, n, fill = value)) + 
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "ILC1 Intrinsic Comparison of Pletropic QTLs",
         x = "",
         y = "percent",
         fill = "Pleitropic")

# ilc1 <- vars %>% 
#   left_join(ilc1_polygenic) %>%
#   dplyr::count(ilc1_expressed, ilc1_egenes_cv, ilc1_polygenic_egene, ilc1_eqtl_loci_genes_cv, ilc1_ciseqtl_assoc_cv) %>% 
#   # dplyr::filter(ilc1_eqtl_loci_genes_cv == 1) %>%
#   dplyr::rename(intrinsic = ilc1_expressed,
#                 egene = ilc1_egenes_cv,
#                 loci_gene = ilc1_eqtl_loci_genes_cv,
#                 polygenic = ilc1_polygenic_egene,
#                 cis_effect = ilc1_ciseqtl_assoc_cv) %>%
#   dplyr::mutate(cell_type = "ILC1",
#                 percent = round(n / sum(n), 4))
                




## ILC2 #########################

# ILC2 eGenes but loci gene is not expresssed in ILC2

ilc2_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(ilc2_pleitropic = sum(ilc2_eqtl_loci_cv)) %>% 
  filter(ilc2_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(ilc2_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


ilc2 <- vars %>% 
  left_join(ilc2_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(ilc2_expressed==1, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

ilc2_summary <- ilc2 %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>% 
  mutate(cell_type = "ILC2")


ilc2_plot <- ggplot(ilc2_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "ILC2 Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")





## ILC3 #########################

ilc3_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(ilc3_pleitropic = sum(ilc3_eqtl_loci_cv)) %>% 
  filter(ilc3_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(ilc3_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


ilc3 <- vars %>% 
  left_join(ilc3_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(ilc3_expressed==1, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

ilc3_summary <- ilc3 %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>%
  mutate(cell_type = "ILC3")


ilc3_plot <- ggplot(ilc3_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "ILC3 Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")







## LTi-like #########################


lti_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(lti_pleitropic = sum(lti_eqtl_loci_cv)) %>% 
  filter(lti_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(lti_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


lti <- vars %>% 
  left_join(lti_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(lti_expressed==1, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

lti_summary <- lti %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>%
  mutate(cell_type = "LTi-like")

lti_plot <- ggplot(lti_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "LTi-like Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")



df_ilc_poly <- list(ilc1_summary, ilc2_summary, ilc3_summary, lti_summary) %>% 
  bind_rows() %>%
  mutate(intrinsic = factor(intrinsic, levels = c("intrinsic", "extrinsic")))


write.csv(df_ilc_poly, paste0(data_path, "intrinsic-polygenic-by-celltype.csv"), row.names = F)


pout <- ggplot(df_ilc_poly, aes(value, n, fill = intrinsic)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~cell_type) +
  labs(title = "Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "") + 
  theme(axis.text.x = element_text(angle = 45))


ggsave(filename = paste0("results/figures/figure-4B-eqtl-intrinsic-pleitropic.pdf"),
       plot = pout,
       width = 7,
       height = 5,
       dpi = 5)

## proportion QTL #######

### ILC1 vs ILC2 #########

ilc1_ilc2_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(ilc1_ilc2_pleitropic = sum(ilc1_ilc2_prop_qtl)) %>% 
  filter(ilc1_ilc2_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(ilc1_ilc2_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


ilc1_ilc2 <- vars %>% 
  left_join(ilc1_ilc2_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(ilc1_expressed==1 | ilc2_expressed, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

ilc1_ilc2_summary <- ilc1_ilc2 %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>%
  mutate(comparison = "ILC1 vs ILC2")

ilc1_ilc3_plot <- ggplot(ilc1_ilc2_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "ILC1 vs ILC2 Proportion QTL Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")


### ILC1 vs ILC3 #########

ilc1_ilc3_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(ilc1_ilc3_pleitropic = sum(ilc1_ilc3_prop_qtl)) %>% 
  filter(ilc1_ilc3_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(ilc1_ilc3_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


ilc1_ilc3 <- vars %>% 
  left_join(ilc1_ilc3_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(ilc1_expressed==1 | ilc3_expressed, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

ilc1_ilc3_summary <- ilc1_ilc3 %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>%
  mutate(comparison = "ILC1 vs ILC3")

ilc1_ilc3_plot <- ggplot(ilc1_ilc3_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "ILC1 vs ILC3 Proportion QTL Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")


### ILC1 vs LTi #########

ilc1_lti_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(ilc1_lti_pleitropic = sum(ilc1_ilc3_prop_qtl)) %>% 
  filter(ilc1_lti_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(ilc1_lti_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


ilc1_lti <- vars %>% 
  left_join(ilc1_lti_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(ilc1_expressed==1, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

ilc1_lti_summary <- ilc1_lti %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>%
  mutate(comparison = "ILC1 vs LTi")

ilc1_lti_plot <- ggplot(ilc1_lti_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "ILC1 vs LTi Proportion QTL Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")


### ILC2 vs ILC3 #########

ilc2_ilc3_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(ilc2_ilc3_pleitropic = sum(ilc2_ilc3_prop_qtl)) %>% 
  filter(ilc2_ilc3_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(ilc2_ilc3_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


ilc2_ilc3 <- vars %>% 
  left_join(ilc2_ilc3_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(ilc2_expressed==1 | ilc3_expressed, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

ilc2_ilc3_summary <- ilc2_ilc3 %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>%
  mutate(comparison = "ILC2 vs ILC3")

ilc2_ilc3_plot <- ggplot(ilc2_ilc3_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "ILC2 vs ILC3 Proportion QTL Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")


### ILC2 vs LTi #########

ilc2_lti_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(ilc2_lti_pleitropic = sum(ilc2_ilc3_prop_qtl)) %>% 
  filter(ilc2_lti_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(ilc2_lti_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


ilc2_lti <- vars %>% 
  left_join(ilc2_lti_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(ilc2_expressed==1, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

ilc2_lti_summary <- ilc2_lti %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>%
  mutate(comparison = "ILC2 vs LTi")

ilc2_lti_plot <- ggplot(ilc2_lti_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "ILC2 vs LTi Proportion QTL Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")


### ILC3 vs LTi #########

ilc3_lti_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(ilc3_lti_pleitropic = sum(ilc3_lti_prop_qtl)) %>% 
  filter(ilc3_lti_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(ilc3_lti_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


ilc3_lti <- vars %>% 
  left_join(ilc3_lti_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(ilc3_expressed==1 | lti_expressed==1, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

ilc3_lti_summary <- ilc3_lti %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>%
  mutate(comparison = "ILC3 vs LTi")

ilc3_lti_plot <- ggplot(ilc3_lti_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "ILC3 vs LTi Proportion QTL Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")


### ilc3_activated_qtl #########

ilc3_activated_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(ilc3_activated_pleitropic = sum(ilc3_activated_qtl)) %>% 
  filter(ilc3_activated_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(ilc3_activated_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


ilc3_activated <- vars %>% 
  left_join(ilc3_activated_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(ilc3_expressed==1, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

ilc3_activated_summary <- ilc3_activated %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>%
  mutate(comparison = "ILC3 Rorgt+/-")

ilc3_activated_plot <- ggplot(ilc3_activated_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "ILC3 Rorgt +/-Proportion QTL Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")


### lti_activated_qtl #########

lti_activated_pleitropic <- ccre %>% 
  group_by(ensembl_gene) %>% 
  summarise(lti_activated_pleitropic = sum(lti_activated_qtl)) %>% 
  filter(lti_activated_pleitropic > 0) %>%
  mutate(pleitropic = ifelse(lti_activated_pleitropic == 1, "monotropic", "pleitropic")) %>%
  dplyr::rename(gene_ids = ensembl_gene)


lti_activated <- vars %>% 
  left_join(lti_activated_pleitropic) %>%
  dplyr::mutate(intrinsic = ifelse(lti_expressed==1, "intrinsic", "extrinsic")) %>%
  dplyr::select(gene_ids, pleitropic, intrinsic)

lti_activated_summary <- lti_activated %>% 
  filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
  count(intrinsic, pleitropic) %>% 
  pivot_longer(pleitropic) %>%
  mutate(comparison = "LTi-like Ccr6 +/-")

lti_activated_plot <- ggplot(lti_activated_summary, aes(intrinsic, n, fill = value)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "LTi-like Ccr6 +/- Proportion QTL Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "Pleitropic")

## using a family-wide error rate #####

# gwas_lods <- fread(paste0(project_path, "results/gwas-ilc1-ilc3-results.csv.gz"), data.table = F) %>%
#   mutate(fwer = quantile(p_values, 0.01),
#          xpos_floor = floor(xpos)) %>%
#   filter(p_values < fwer) %>%
#   summarise(qtl_count = n_distinct(xpos_floor))



df_poly_props <- list(ilc1_ilc2_summary, ilc1_ilc3_summary, ilc1_lti_summary,
     ilc2_ilc3_summary, ilc2_lti_summary, ilc3_lti_summary) %>%
  bind_rows() %>%
  mutate(intrinsic = factor(intrinsic, levels = c("intrinsic", "extrinsic")),
         comparison = str_replace(comparison, "_QTL", ""))

write.csv(df_poly_props, paste0(data_path, "intrinsic-polygenic-by-celltype.csv"), row.names = F)


pout <- ggplot(df_poly_props, aes(value, n, fill = intrinsic)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~comparison) +
  labs(title = "Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "") +
  theme(axis.text.x = element_text(angle = 45))


ggsave(filename = paste0("results/figures/figure-4B-proportion-qtls-intrinsic-pleitropic.pdf"),
       plot = pout,
       width = 7,
       height = 5,
       dpi = 5)


## topic QTL: polygenic #######

# topic_poly <- fread(paste0(project_path, "results/topic-qtl-lods.csv.gz")) %>%
#   left_join({
#     ccre %>% 
#       dplyr::select(marker, chr, pos)
#       }) %>%
#   mutate(pos_floor = round(pos/1e6)) %>%
#   dplyr::select(-V1, -pos) %>%
#   pivot_longer(-c(marker,chr,pos_floor)) %>%
#   group_by(name, chr, pos_floor) %>%
#   summarise(lods = max(value)) %>%
#   group_by(name) %>%
#   summarise(qtl_count = sum(lods > 6),
#             group = "topic")


create_pleitropic_compare <- function(varcol) {
  
  pleitropic_df <- ccre %>% 
    group_by(ensembl_gene) %>% 
    summarise(pleitropic_total = sum(.data[[varcol]])) %>% 
    filter(pleitropic_total > 0) %>%
    mutate(pleitropic = ifelse(pleitropic_total == 1, "monotropic", "pleitropic")) %>%
    dplyr::rename(gene_ids = ensembl_gene)
  
  intrinsic_df <- vars %>% 
    left_join(pleitropic_df) %>%
    dplyr::mutate(intrinsic = ifelse(ilc1_expressed == 1 | ilc2_expressed == 1 | ilc3_expressed == 1 | lti_expressed==1, "intrinsic", "extrinsic")) %>%
    dplyr::select(gene_ids, pleitropic, intrinsic)


  summary_df <- intrinsic_df %>% 
    filter(!is.na(intrinsic), !is.na(pleitropic)) %>%
    count(intrinsic, pleitropic) %>% 
    pivot_longer(pleitropic) %>%
    mutate(comparison = toupper(varcol))


  # plot_out <- ggplot(summary_df, aes(intrinsic, n, fill = value)) + 
  #   geom_bar(stat = "identity", position = "fill") +
  #   labs(title = "LTi-like Ccr6 +/- Proportion QTL Intrinsic Comparison of Pletropic QTLs",
  #        x = "",
  #        y = "percent",
  #        fill = "Pleitropic")
  
  summary_df

}


topic_pleitropic_df <- lapply(paste0("topic", 0:19, "_qtl"), create_pleitropic_compare) %>% 
  Filter(function(x) NROW(x) > 0,.) %>% 
  bind_rows() %>%
  mutate(intrinsic = factor(intrinsic, levels = c("intrinsic", "extrinsic")),
         comparison = str_replace(comparison, "_QTL", "")) %>%
  mutate(comparison = paste("Topic", str_extract(comparison, "[0-9]{1,2}"))) %>%
  mutate(comparison = factor(comparison, levels = paste("Topic", 0:19)))

write.csv(topic_pleitropic_df, paste0(data_path, "intrinsic-pleitropic-topics.csv"), row.names = F)


pout <- ggplot(topic_pleitropic_df, aes(value, n, fill = intrinsic)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~comparison) +
  labs(title = "Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "") + 
  theme(axis.text.x = element_text(angle=45))


ggsave(filename = paste0("results/figures/figure-4B-topic-qtls-intrinsic-pleitropic.pdf"),
       plot = pout,
       width = 7,
       height = 5,
       dpi = 5)




## cytokine qtl: polygenic ###########################


cytokine_pleitropic_df <- lapply(colnames(ccre)[85:96], create_pleitropic_compare) %>% 
  Filter(function(x) NROW(x) > 0,.) %>% 
  bind_rows() %>%
  mutate(intrinsic = factor(intrinsic, levels = c("intrinsic", "extrinsic")),
         comparison = str_replace(comparison, "_QTL_STEADY", ""))

write.csv(cytokine_pleitropic_df, paste0(data_path, "intrinsic-pleitropic-cytokines.csv"), row.names = F)


pout <- ggplot(cytokine_pleitropic_df, aes(value, n, fill = intrinsic)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~comparison) +
  labs(title = "Intrinsic Comparison of Pletropic QTLs",
       x = "",
       y = "percent",
       fill = "") + 
  theme(axis.text.x = element_text(angle=45))


ggsave(filename = paste0("results/figures/figure-4B-cytokines-qtls-intrinsic-pleitropic.pdf"),
       plot = pout,
       width = 7,
       height = 5,
       dpi = 5)




