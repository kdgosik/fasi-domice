#' Figure 6B Intrinsic by Trait bar plots
#'
#'
#'
#'
#'

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)


trait_by_loci_no_window <- fread(paste0(results_dir, "trait_by_loci_no_window.csv"))
trait_by_loci_10kb_window <- fread(paste0(results_dir, "trait_by_loci_10kb_window.csv"))
trait_by_loci_50kb_window <- fread(paste0(results_dir, "trait_by_loci_50kb_window.csv"))
trait_by_loci_500kb_window <- fread(paste0(results_dir, "trait_by_loci_500kb_window.csv"))




trait_by_loci_no_window %>% 
  replace_na(list(ilc1_expressed=0, ilc2_expressed=0, ilc3_expressed=0, lti_expressed=0)) %>% 
  mutate(ilc_expressed = as.numeric(ilc1_expressed+ilc2_expressed+ilc3_expressed+lti_expressed > 0))  %>%
  separate("trait", c("cell_type", "gene"), remove=FALSE) %>% 
  group_by(cell_type) %>% summarise(prop = mean(ilc_expressed))


trait_by_loci_50kb_window %>% 
  replace_na(list(ilc1_expressed=0, ilc2_expressed=0, ilc3_expressed=0, lti_expressed=0)) %>% 
  mutate(ilc_expressed = as.numeric(ilc1_expressed+ilc2_expressed+ilc3_expressed+lti_expressed > 0))  %>%
  separate("trait", c("cell_type", "gene"), remove=FALSE) %>% 
  group_by(cell_type) %>% 
  summarise(prop = mean(ilc_expressed)) %>%
  arrange(desc(prop))



trait_by_loci_50kb_window %>% 
  replace_na(list(ilc1_expressed=0, ilc2_expressed=0, ilc3_expressed=0, lti_expressed=0)) %>% 
  mutate(ilc_expressed = as.numeric(ilc1_expressed+ilc2_expressed+ilc3_expressed+lti_expressed > 0))  %>%
  separate("trait", c("cell_type", "gene"), remove=FALSE) %>% 
  group_by(cell_type, loci) %>% 
  summarise(ilc_expressed = max(ilc_expressed)) %>%
  group_by(cell_type) %>% 
  summarise(prop = mean(ilc_expressed)) %>%
  arrange(desc(prop))



outdf <- trait_by_loci_50kb_window %>% 
  replace_na(list(ilc1_expressed=0, ilc2_expressed=0, ilc3_expressed=0, lti_expressed=0)) %>% 
  mutate(ilc_expressed = as.numeric(ilc1_expressed+ilc2_expressed+ilc3_expressed+lti_expressed > 0))  %>%
  separate("trait", c("cell_type", "gene"), remove=FALSE) %>% 
  mutate(cell_type = case_when(str_detect(trait, "ILC[1-3]_ILC[1-3]|ILC[1-3]_LTi") ~ trait,
                               TRUE ~ cell_type)) %>%
  mutate(trait_group = case_when(str_detect(trait, "ILC[1-3]_ILC[1-3]|ILC[1-3]_LTi") ~ "proportion",
                                 str_detect(trait, "topic") ~ "topic",
                                 str_detect(trait, "ILC[1-3]") ~ "expression",
                                 cell_type == "LTi" ~ "expression",
                                 TRUE ~ "cytokine")) %>%
  group_by(trait_group, cell_type, loci) %>% 
  summarise(ilc_expressed = max(ilc_expressed)) %>%
  group_by(trait_group, cell_type) %>% 
  summarise(prop = mean(ilc_expressed)) %>%
  arrange(desc(prop))




plot_all_trait_groups <- outdf %>%
  ungroup() %>%
  mutate(extrinsic = 1-prop) %>%
  dplyr::select(trait_group, cell_type, intrinsic = prop, extrinsic) %>%
  pivot_longer(-c(trait_group, cell_type)) %>%
  mutate(trait_group = factor(trait_group, levels = c("expression", "proportion", "topic", "cytokine")),
         cell_type = str_pad(str_remove(cell_type, "topic"), pad = "0", width = 2, side = "left")) %>%
  ggplot(aes(cell_type, value, fill = name)) + 
    geom_bar(stat = "identity") +
    facet_wrap(~trait_group, scales = "free_x") +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = "bottom") +
    labs(title = "Intrinsic QTL Loci Genes",
         x = "trait",
         y = "proportion",
         fill = "location")


ggsave(file = "qtl_trait_groups_intrinsic_proportion.pdf",
       plot = plot_all_trait_groups,
       width = 7,
       height = 5)


plot_expression_trait_groups <- outdf %>%
  ungroup() %>%
  mutate(extrinsic = 1-prop) %>%
  dplyr::select(trait_group, cell_type, intrinsic = prop, extrinsic) %>%
  pivot_longer(-c(trait_group, cell_type)) %>%
  mutate(trait_group = factor(trait_group, levels = c("expression", "proportion", "topic", "cytokine")),
         cell_type = str_pad(str_remove(cell_type, "topic"), pad = "0", width = 2, side = "left")) %>%
  filter(str_detect(trait_group, "expression")) %>%
  ggplot(aes(cell_type, value, fill = name)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~trait_group, scales = "free_x") +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "bottom") +
  labs(title = "Intrinsic QTL Loci Genes",
       x = "trait",
       y = "proportion",
       fill = "location")


ggsave(file = "qtl_expression_trait_groups_intrinsic_proportion.pdf",
       plot = plot_expression_trait_groups,
       width = 7,
       height = 5)


plot_proportion_trait_groups <- outdf %>%
  ungroup() %>%
  mutate(extrinsic = 1-prop) %>%
  dplyr::select(trait_group, cell_type, intrinsic = prop, extrinsic) %>%
  pivot_longer(-c(trait_group, cell_type)) %>%
  mutate(trait_group = factor(trait_group, levels = c("expression", "proportion", "topic", "cytokine")),
         cell_type = str_pad(str_remove(cell_type, "topic"), pad = "0", width = 2, side = "left")) %>%
  filter(str_detect(trait_group, "proportion")) %>%
  ggplot(aes(cell_type, value, fill = name)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~trait_group, scales = "free_x") +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "bottom") +
  labs(title = "Intrinsic QTL Loci Genes",
       x = "trait",
       y = "proportion",
       fill = "location")


ggsave(file = "qtl_proportion_trait_groups_intrinsic_proportion.pdf",
       plot = plot_proportion_trait_groups,
       width = 7,
       height = 5)




plot_topic_trait_groups <- outdf %>%
  ungroup() %>%
  mutate(extrinsic = 1-prop) %>%
  dplyr::select(trait_group, cell_type, intrinsic = prop, extrinsic) %>%
  pivot_longer(-c(trait_group, cell_type)) %>%
  mutate(trait_group = factor(trait_group, levels = c("expression", "proportion", "topic", "cytokine")),
         cell_type = str_pad(str_remove(cell_type, "topic"), pad = "0", width = 2, side = "left")) %>%
  filter(str_detect(trait_group, "topic")) %>%
  ggplot(aes(cell_type, value, fill = name)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~trait_group, scales = "free_x") +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "bottom") +
  labs(title = "Intrinsic QTL Loci Genes",
       x = "trait",
       y = "proportion",
       fill = "location")


ggsave(file = "qtl_topic_trait_groups_intrinsic_proportion.pdf",
       plot = plot_topic_trait_groups,
       width = 7,
       height = 5)




poly_vs_mono <- trait_by_loci_500kb_window %>% 
  replace_na(list(ilc1_expressed=0, ilc2_expressed=0, ilc3_expressed=0, lti_expressed=0)) %>% 
  mutate(ilc_expressed = as.numeric(ilc1_expressed+ilc2_expressed+ilc3_expressed+lti_expressed > 0))  %>%
  separate("trait", c("cell_type", "gene"), remove=FALSE) %>% 
  mutate(cell_type = case_when(str_detect(trait, "ILC[1-3]_ILC[1-3]|ILC[1-3]_LTi") ~ trait,
                               TRUE ~ cell_type)) %>%
  mutate(trait_group = case_when(str_detect(trait, "ILC[1-3]_ILC[1-3]|ILC[1-3]_LTi") ~ "proportion",
                                 str_detect(trait, "topic") ~ "topic",
                                 str_detect(trait, "ILC[1-3]") ~ "expression",
                                 cell_type == "LTi" ~ "expression",
                                 TRUE ~ "cytokine")) %>%
  distinct(loci, cell_type, trait_group, count)




poly_vs_mono %>%
  filter(trait_group == "expression") %>%
  mutate(polygenic = case_when(count > 1 ~ "polygenic", TRUE ~ "mongenic")) %>%
  count(cell_type, polygenic) %>%
  ggplot(aes(cell_type, n, fill = cell_type, group = polygenic)) + 
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Monogenic and Polygenic eQTLs")


poly_vs_mono %>%
  filter(trait_group == "proportion") %>%
  distinct(loci, trait_group, count) %>%
  mutate(polygenic = case_when(count > 1 ~ "polygenic", TRUE ~ "mongenic")) %>%
  count(trait_group, polygenic) %>%
  ggplot(aes(trait_group, n, fill = trait_group, group = polygenic)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Monogenic and Polygenic Proportion")



poly_vs_mono %>%
  filter(trait_group == "topic") %>%
  distinct(loci, trait_group, count) %>%
  mutate(polygenic = case_when(count > 1 ~ "polygenic", TRUE ~ "mongenic")) %>%
  count(trait_group, polygenic) %>%
  ggplot(aes(trait_group, n, fill = trait_group, group = polygenic)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Monogenic and Polygenic Topics")

