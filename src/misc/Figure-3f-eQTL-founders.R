#' Figure 3E
#' 
#' 

# project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
# if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"
project_path <- my_path <-"./domice/"
project_path <- my_path <- "/home/kirk/Documents/projects/domice/"
figure_path <- paste0(project_path, "results/figures/")
data_path <- paste0(project_path, "data/")
source(paste0(figure_path, "helpers.R"))

load(paste0(data_path, "GM_snps.Rdata"))
## founder key
founderdf <- data.frame(founder_letter = LETTERS[1:8], 
                        founder_color = qtl2::CCcolors, 
                        founder = names(qtl2::CCcolors),
                        best_clean = c("A_J", "C57BL_6NJ", "129S1_SvImJ", "NOD_ShiLtJ",
                                       "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ"),
                        stringsAsFactors = FALSE)

topic_geno_df <- fread(paste0(project_path, "results/topic-qtl-genotypes.csv.gz"))

outcor <- sapply(colnames(topic_geno_df)[3:30], function(col) {
  # topic3
  # UNC4192456
  # dir(paste0(my_path, "data/"), pattern = "founder-coef.csv")
  var_col <- gsub("_.*", "", col)
  marker <- gsub(".*_", "", col)
  chr <- GM_snps$chr[GM_snps$marker == marker]
  
  tryCatch({
    coef <- read.csv(paste0(my_path, "results/founder-coefs/", var_col, "-chr", chr, "-founder-coef.csv"))
    coefdf <- coef[coef$X == marker, LETTERS[1:8]] %>% 
      t() %>% 
      data.frame() %>% 
      rownames_to_column("founder_letter") %>%
      left_join(founderdf) %>%
      dplyr::rename(founder_predicted = starts_with("X"))
    
    cat("Reading founders ... \n")
    f1obs <- fread(paste0(data_path, "founders1-louvain_labels-labels-transfer/obs.csv")) %>%
      mutate(best_clean2 = ifelse(best_clean=="AMB", SNG.1ST, best_clean),
             founder_topic = scale(.data[[var_col]]))
    
    f2obs <- fread(paste0(data_path, "founders2-louvain_labels-labels-transfer/obs.csv")) %>%
      mutate(best_clean2 = ifelse(best_clean=="AMB", SNG.1ST, best_clean),
             founder_topic = scale(.data[[var_col]]))
    
    fobs <- rbind(f1obs, f2obs)
    
    outdf <- fobs %>%
      filter(best_clean != "nan") %>%
      group_by(best_clean) %>%
      summarise(founder_topic = mean(founder_topic)) %>%
      left_join({
        coefdf %>%
          select(best_clean, founder_predicted)
      })
    cor(outdf$founder_topic, outdf$founder_predicted)
  }, error = function(e) 0)
  
}, simplify = TRUE, USE.NAMES = TRUE)


f1obs <- fread(paste0(data_path, "founders1-louvain_labels-labels-transfer/obs.csv")) %>%
  mutate(best_clean2 = ifelse(best_clean=="AMB", SNG.1ST, best_clean)) %>% 
  mutate(across(.cols = starts_with("topic"), .fns = scale,.names = "{.col}_scale"))

f2obs <- fread(paste0(data_path, "founders2-louvain_labels-labels-transfer/obs.csv")) %>%
  mutate(best_clean2 = ifelse(best_clean=="AMB", SNG.1ST, best_clean)) %>% 
  mutate(across(.cols = starts_with("topic"), .fns = scale,.names = "{.col}_scale"))


fobs <- rbind(f1obs, f2obs)


# fobs %>% 
#   filter(best_clean2 != "nan") %>%
#   group_by(best_clean2) %>%
#   summarise(topic3_scale = mean(topic3_scale)) %>%
#   ggplot(., aes(best_clean2, topic3_scale)) + geom_point()

fobs %>% 
  filter(best_clean2 != "nan") %>%
  ggplot(., aes(best_clean2, topic3_scale)) + geom_point()


topic_coefs %>%
  filter(topic == "topic3") %>%
  # mutate(predicted = scale(value)) %>%
  ggplot(., aes(best_clean, predicted, color = X)) + geom_point()






library(data.table)
load(paste0(data_path, "GM_snps.Rdata"))
founderdf <- data.frame(founder_letter = LETTERS[1:8], 
                        founder_color = qtl2::CCcolors, 
                        founder = names(qtl2::CCcolors),
                        best_clean = c("A_J", "C57BL_6NJ", "129S1_SvImJ", "NOD_ShiLtJ",
                                       "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ"),
                        stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
celltype <- as.character(args[1])
cat("Using celltype:", celltype, "...\n")


## RUN COEF QTL ###################

## eQTL #################
# eqtl_df <- fread(paste0(my_path, "results/eqtl-significant-loci.csv"), data.table = FALSE) %>%
#   # dplyr::filter(cell_type == celltype) %>% 
#   dplyr::group_by(cell_type, gene, marker_chr) %>% 
#   dplyr::arrange(desc(value_adj)) %>% 
#   mutate(min_rank = min_rank(value_adj), n = n()) %>%
#   dplyr::filter(min_rank == n)

vars <- read_csv(paste0(data_path, "vars.csv"))
eqtl_genes <-unique(c(vars$index[vars$ilc1_egenes_cv == 1],
                      vars$index[vars$ilc2_egenes_cv == 1],
                      vars$index[vars$ilc3_egenes_cv == 1],
                      vars$index[vars$lti_egenes_cv == 1]))

ccre <- fread(paste0(data_path, "GM_SNPS_Consequence_cCRE.csv"))
ilc1_ccre <- ccre$marker[ccre$ilc1_eqtl_loci_cv==1]
ilc2_ccre <- ccre$marker[ccre$ilc2_eqtl_loci_cv==1]
ilc3_ccre <- ccre$marker[ccre$ilc3_eqtl_loci_cv==1]
lti_ccre <- ccre$marker[ccre$lti_eqtl_loci_cv==1]
loci_markers <- unique(c(ilc1_ccre, ilc2_ccre, ilc3_ccre, lti_ccre))

  
outdf <- lapply(c("ILC1", "ILC2", "ILC3", "LTi"), function(celltype) {
  
  in_file <- paste0("zcat ", my_path, "results/eqtl-", celltype, "-founder-coefs.csv.gz")
  coef <- fread(in_file, data.table = FALSE)
  # coef <- read.csv(in_file) %>%
  #   dplyr::filter(marker %in% eqtl_df$marker[eqtl_df$cell_type == celltype])
  
  if(celltype == "ILC1"){
    coef <- coef %>% filter(marker %in% ilc1_ccre)
  }
  if(celltype == "ILC2"){
    coef <- coef %>% filter(marker %in% ilc2_ccre)
  }
  if(celltype == "ILC3"){
    coef <- coef %>% filter(marker %in% ilc3_ccre)
  }
  if(celltype == "LTi"){
    coef <- coef %>% filter(marker %in% lti_ccre)
  }
  coef
  
  
}) %>% bind_rows() %>% 
  dplyr::rename(lods = value) %>%
  tidyr::pivot_longer(cols = A:H) %>%
  dplyr::mutate(pred_exp = intercept + value) %>%
  dplyr::rename(founder_letter = name) %>%
  dplyr::left_join(founderdf)


f1_cnames <- readLines(paste0(data_path, "founders1-louvain_labels-labels-transfer/X.csv"), n=1)
f1_cnames <- unlist(strsplit(x = f1_cnames, split = ","))
f1_cnames <- intersect(eqtl_genes, f1_cnames)
f2_cnames <- readLines(paste0(data_path, "founders2-louvain_labels-labels-transfer/X.csv"), n=1)
f2_cnames <- unlist(strsplit(x = f2_cnames, split = ","))
f2_cnames <- intersect(eqtl_genes, f2_cnames)
genelist <- intersect(f1_cnames, f2_cnames)
genelist <- genelist[-c(1,2)]

f1X <- fread(paste0(data_path, "founders1-louvain_labels-labels-transfer/X.csv"),
             select = c("index", "best_clean", "SNG.1ST","called_cell_types_new", genelist)) %>%
  dplyr::filter(called_cell_types_new %in% c("ILC1", "ILC2", "ILC3", "ILC3(LTi-like)"))

f2X <- fread(paste0(data_path, "founders2-louvain_labels-labels-transfer/X.csv"),
             select = c("index", "best_clean", "SNG.1ST","called_cell_types_new", genelist)) %>%
  dplyr::filter(called_cell_types_new %in% c("ILC1", "ILC2", "ILC3", "ILC3(LTi-like)"))

fX <- rbind(f1X, f2X) %>%
  dplyr::mutate(best_clean2 = ifelse(best_clean=="AMB", SNG.1ST, best_clean),
                cell_type = str_replace(called_cell_types_new, "ILC3\\(LTi-like\\)", "LTi"))

fXsum <- fX %>%
  tidyr::pivot_longer(cols = Mrpl30:Uhrf2, names_to = "gene") %>%
  dplyr::group_by(best_clean, cell_type, gene) %>%
  dplyr::summarise(expression = mean(value, na.rm = TRUE))


merge_genes <- outdf %>%
  left_join(fXsum, by = c("best_clean", "gene", "cell_type")) %>%
  dplyr::mutate(expression = replace_na(expression, 0))



foundercor <- merge_genes %>%
  split(.$gene) %>%
  map_dbl(~cor(.x$pred_exp, .x$expression)) %>%
  as.data.frame() %>%
  dplyr::rename(correlation = .data[["."]]) %>%
  rownames_to_column("gene") %>% 
  left_join(outdf)

p1 <- foundercor %>%
  # filter(value_adj > 10) %>%
  # dplyr::select(cell_type, gene, marker, correlation, cis_effect) %>%
  left_join(GM_snps, by = "marker") %>%
  dplyr::group_by(gene, marker, pos, cis_effect)  %>%
  summarise(correlation = mean(correlation)) %>%
  unique %>%
  ggplot(., aes(x = correlation, fill = factor(cis_effect))) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "eQTL Prediction Correlation with Founder Strain", fill = "Cis Effect")

ggsave(filename = "results/figures/figure-3-eqtl-founderstrain.png",
       plot = p1,
       dpi = 330)

ggsave(filename = "results/figures/figure-3-eqtl-founderstrain.pdf",
       plot = p1,
       dpi = 330)


foundercor <- merge_genes %>%
  # split(.$gene) %>%
  split(list(.$cell_type, .$gene)) %>% 
  map_dbl(~cor(.x$pred_exp, .x$expression)) %>%
  as.data.frame() %>%
  dplyr::rename(correlation = .data[["."]]) %>%
  rownames_to_column("celltype_gene") %>% 
  separate(celltype_gene, c("cell_type", "gene")) %>%
  left_join(outdf)


p1 <- foundercor %>%
  # filter(value_adj > 10) %>%
  # dplyr::select(cell_type, gene, marker, correlation, cis_effect) %>%
  left_join(GM_snps, by = "marker") %>%
  dplyr::group_by(cell_type, gene, marker, pos, cis_effect)  %>%
  summarise(correlation = mean(correlation)) %>%
  unique %>%
  ggplot(., aes(x = correlation, fill = factor(cis_effect))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~cell_type) +
  theme_minimal() +
  labs(title = "Correlation with Founder Strain", fill = "Cis Effect")

p2 <- foundercor %>%
  filter(!is.na(cis_effect)) %>%
  dplyr::select(cell_type, gene, marker, correlation, value_adj, cis_effect) %>%
  unique %>%
  ggplot(., aes(correlation, value_adj, color = factor(cis_effect))) +
  geom_point() +
  facet_wrap(~cis_effect) +
  theme_minimal() +
  labs(title = "Correlation with Founder Strain", fill = "Cis Effect")


## top 20 correlated genes with founder 
# "Tmx1"   "Fgl2"   "Ywhah"  "Psma2"  "Rel"    "Mydgf"  "Ywhag"  "Grb2"   "Sptssa" "Mrps24" "Ppil2"  "Gpr171"
# [13] "Tecpr1" "Anp32e" "Hopx"   "Tyrobp" "Ostf1"  "Tmbim4" "Lgals3" "Msn" 


# Tmx1
p1 <- fX %>% 
  filter(cell_type=="LTi", best_clean != "nan", best_clean !="AMB") %>%
  # group_by(best_clean) %>%
  # summarise(Tmx1 =mean(Tmx1)) %>%
  ggplot(.,aes(best_clean, Tmx1)) + 
  geom_boxplot(fill = "lightgray") +
  geom_point() + 
  labs(x = "Founder")

p2 <- merge_genes %>% 
  filter(cell_type=="LTi", gene == "Tmx1") %>% 
  mutate(predicted= scale(pred_exp)) %>% 
  ggplot(.,aes(best_clean, predicted)) + 
  geom_point() + 
  labs(x = "Founder",
       y = "Predicted")

tmx1_plot <- p1 / p2

ggsave(filename = paste0("results/figures/pdf/figure-3f-founder-eQTL-validation-Tmx1.pdf"),
       plot = tmx1_plot,
       dpi = 330,
       width = 7,
       height = 5)


# Ywhah
p1 <- fX %>% 
  filter(cell_type=="ILC1", best_clean != "nan", best_clean !="AMB") %>%
  # group_by(best_clean) %>%
  # summarise(Tmx1 =mean(Tmx1)) %>%
  ggplot(.,aes(best_clean, Ywhah)) + 
  geom_boxplot(fill = "lightgray") +
  geom_point() + 
  labs(x = "Founder")

p2 <- merge_genes %>% 
  filter(cell_type=="ILC1", gene == "Ywhah") %>% 
  mutate(predicted= scale(pred_exp)) %>% 
  ggplot(.,aes(best_clean, predicted)) + 
  geom_point() + 
  labs(x = "Founder",
       y = "Predicted")

ywhah_plot <- p1 / p2

ggsave(filename = paste0("results/figures/pdf/figure-3f-founder-eQTL-validation-Ywhah.pdf"),
       plot = ywhah_plot,
       dpi = 330,
       width = 7,
       height = 5)


# Ywhag
p1 <- fX %>% 
  filter(cell_type=="LTi", best_clean != "nan", best_clean !="AMB") %>%
  # group_by(best_clean) %>%
  # summarise(Tmx1 =mean(Tmx1)) %>%
  ggplot(.,aes(best_clean, Ywhag)) + 
  geom_boxplot(fill = "lightgray") +
  geom_point() + 
  labs(x = "Founder")

p2 <- merge_genes %>% 
  filter(cell_type=="LTi", gene == "Ywhag") %>% 
  mutate(predicted= scale(pred_exp)) %>% 
  ggplot(.,aes(best_clean, predicted)) + 
  geom_point() + 
  labs(x = "Founder",
       y = "Predicted")

ywhag_plot <- p1 / p2

ggsave(filename = paste0("results/figures/pdf/figure-3f-founder-eQTL-validation-Ywhag.pdf"),
       plot = ywhag_plot,
       dpi = 330,
       width = 7,
       height = 5)


# Hopx
p1 <- fX %>% 
  filter(cell_type=="ILC1", best_clean != "nan", best_clean !="AMB") %>%
  group_by(best_clean) %>%
  summarise(Hopx =mean(Hopx)) %>%
  ggplot(.,aes(best_clean, Hopx)) + 
  geom_boxplot(fill = "lightgray") +
  geom_point() + 
  labs(x = "Founder")

p2 <- merge_genes %>% 
  filter(cell_type=="ILC1", gene == "Hopx") %>% 
  mutate(predicted= scale(pred_exp)) %>% 
  ggplot(.,aes(best_clean, predicted)) + 
  geom_point() + 
  labs(x = "Founder",
       y = "Predicted")

hopx_plot <- p1 / p2

ggsave(filename = paste0("results/figures/pdf/figure-3f-founder-eQTL-validation-Hopx.pdf"),
       plot = hopx_plot,
       dpi = 330,
       width = 7,
       height = 5)



## other relevant genes
# Aim2
p1 <- fX %>% 
  filter(cell_type=="ILC1", best_clean != "nan", best_clean !="AMB") %>% 
  group_by(best_clean) %>%
  summarise(Aim2 =mean(Aim2)) %>%
  ggplot(.,aes(best_clean, Aim2)) + 
  geom_boxplot(fill = "lightgray") +
  geom_point() + 
  labs(x = "Founder")

p2 <- merge_genes %>% 
  filter(cell_type=="ILC1", gene == "Aim2") %>% 
  mutate(predicted= scale(pred_exp)) %>% 
  ggplot(.,aes(best_clean, predicted)) + 
    geom_point() + 
    labs(x = "Founder",
         y = "Predicted")

aim2_plot <- p1 / p2

ggsave(filename = paste0("results/figures/pdf/figure-3f-founder-eQTL-validation-Aim2.pdf"),
       plot = aim2_plot,
       dpi = 330,
       width = 7,
       height = 5)



# Bcl2
p1 <- fX %>% 
  filter(cell_type=="ILC1", best_clean != "nan", best_clean !="AMB") %>% 
  ggplot(.,aes(best_clean, Bcl2l11)) + 
  geom_point() +
  geom_violin() +
  geom_jitter()

p2 <- merge_genes %>% 
  filter(cell_type=="ILC1", gene == "Bcl2l11") %>% 
  mutate(predicted= scale(pred_exp)) %>% 
  ggplot(.,aes(best_clean, predicted, color = marker)) + 
  geom_point() + theme(legend.position = "none")


# Cnn2
p1 <- fX %>% 
  filter(cell_type=="ILC1", best_clean != "nan", best_clean !="AMB") %>% 
  group_by(best_clean) %>%
  summarise(Cnn2 =mean(Cnn2)) %>%
  ggplot(.,aes(best_clean, Cnn2)) + 
  geom_boxplot(fill = "lightgray") +
  geom_point() + 
  labs(x = "Founder")

p2 <- merge_genes %>% 
  filter(cell_type=="ILC1", gene == "Cnn2", marker == "UNC4488618") %>% 
  mutate(predicted= scale(pred_exp)) %>% 
  ggplot(.,aes(best_clean, predicted)) + 
  geom_point() + 
  theme(legend.position = "none") + 
  labs(x = "Founder",
       y = "Predicted")

cnn2_plot <- p1 / p2

ggsave(filename = paste0("results/figures/pdf/figure-3f-founder-eQTL-validation-Cnn2.pdf"),
       plot = cnn2_plot,
       dpi = 330,
       width = 7,
       height = 5)




ccn2_plot <- fX %>% 
  filter(cell_type=="ILC1", best_clean != "nan", best_clean != "AMB") %>% 
  ggplot(.,aes(best_clean, Cnn2)) + 
  geom_boxplot(color = "gray") +
  labs(title = "Founder Gene Expression: Ccn2",
       x = "Founder Strain",
       y = "Expression")

ggsave(filename = paste0(data_path, "founder-expression-ccn2.pdf"),
       plot = ccn2_plot,
       dpi = 330)

areg_plot <- fX %>% 
  filter(cell_type=="ILC1", best_clean != "nan", best_clean != "AMB") %>% 
  ggplot(.,aes(best_clean, Areg)) + 
  geom_boxplot(color = "gray") +
  labs(title = "Founder Gene Expression: Ccn2",
       x = "Founder Strain",
       y = "Expression")

ggsave(filename = paste0(data_path, "founder-expression-areg.pdf"),
       plot = areg_plot,
       dpi = 330)


aim2_plot <- fX %>% 
  filter(cell_type=="ILC1", best_clean != "nan", best_clean != "AMB") %>% 
  ggplot(.,aes(best_clean, Aim2)) + 
  geom_boxplot(color = "gray") +
  labs(title = "Founder Gene Expression: Aim2",
       x = "Founder Strain",
       y = "Expression")

ggsave(filename = paste0(data_path, "founder-expression-aim2.pdf"),
       plot = aim2_plot,
       dpi = 330)




## Fgl2
plotout <-merge_genes %>% 
  filter(gene == "Fgl2") %>% 
  group_by(founder) %>% 
  summarise(expression = mean(expression),
            predicted_expression = mean(pred_exp))


ggplot(plotout, aes(founder, expression)) + 
  geom_point() + 
  geom_point(aes(founder, predicted_expression),color = "red") +
  labs(title = "Fgl2 Plots")



ggplot(plotout, aes(expression, predicted_expression)) + 
    geom_point() + 
    geom_label(aes(label = founder, color = founder_color)) +
    # geom_abline(slope=1, color = "red") +
    xlim(0,1) + ylim(0.5,1.5)




## Topics ###########################

plot_df <- fread(paste0(my_path, "data/manuscript-plot-data.csv"), data.table = FALSE)
ccre <- fread(paste0(my_path, "data/GM_SNPS_Consequence_cCRE.csv"))
ccre[, recomb_rate := round(cM / (pos / 1e6), 4)]

topic_lods <- fread(paste0(my_path, "results/topic-qtl-lods.csv.gz"), data.table = FALSE) %>%
  pivot_longer(topic0:topic19, names_to = "topic", values_to = "lods") %>%
  dplyr::select(-V1)

read_vars <- str_extract(string = dir(paste0(my_path, "results/founder-coefs"), "topic"), "topic[0-9]{1,2}-chr[0-9]{1,2}")

topic_coefs <- sapply(read_vars, function(f) {
  
  read.csv(dir(paste0(my_path, "results/founder-coefs"), f, full.names = TRUE))
  
}, USE.NAMES = TRUE, simplify = FALSE) %>%
  bind_rows(.,.id="topic_chr") %>%
  separate(topic_chr, c("topic", "chr"), "-") %>%
  tidyr::pivot_longer(cols = A:H) %>%
  dplyr::mutate(pred_exp = intercept + value) %>%
  dplyr::rename(founder_letter = name) %>%
  dplyr::left_join(founderdf)


f1obs <- fread(paste0(data_path, "founders1-louvain_labels-labels-transfer/obs.csv")) %>%
  dplyr::filter(called_cell_types_new %in% c("ILC1", "ILC2", "ILC3", "ILC3(LTi-like)"))

f2obs <- fread(paste0(data_path, "founders2-louvain_labels-labels-transfer/obs.csv")) %>%
  dplyr::filter(called_cell_types_new %in% c("ILC1", "ILC2", "ILC3", "ILC3(LTi-like)"))

fobs <- rbind(f1obs, f2obs) %>%
  dplyr::mutate(best_clean = ifelse(best_clean=="AMB", SNG.1ST, best_clean),
                cell_type = str_replace(called_cell_types_new, "ILC3\\(LTi-like\\)", "LTi")) %>%
  tidyr::pivot_longer(cols = topic0:topic19, names_to = "topic") %>%
  dplyr::group_by(best_clean, cell_type, topic) %>%
  dplyr::summarise(expression = mean(value, na.rm = TRUE))


merge_topics <- topic_coefs %>%
  dplyr::rename(marker = X) %>%
  left_join(fobs) %>%
  left_join(topic_lods, by = c("topic", "marker")) ## check if correct


founder_topic_cor <- merge_topics %>% 
  split(list(.$cell_type, .$topic)) %>%
  map_dbl(~cor(.x$pred_exp, .x$expression)) %>%
  as.data.frame() %>%
  dplyr::rename(correlation = .data[["."]]) %>%
  rownames_to_column("celltype_topic") %>% 
  separate(celltype_topic, c("cell_type", "topic")) %>%
  left_join(merge_topics)

founder_topic_cor %>%
  dplyr::select(cell_type, topic, marker, correlation) %>%
  unique %>%
  ggplot(., aes(x = correlation)) + 
  geom_density(alpha = 0.5) +
  facet_wrap(~cell_type)

founder_topic_cor %>%
  dplyr::select(cell_type, topic, marker, correlation, lods) %>%
  unique %>%
  ggplot(., aes(correlation, lods)) + 
  geom_point()




## No Celltype grouping ##################
founder_topic_cor <- merge_topics %>% 
  split(.$topic) %>%
  map_dbl(~cor(.x$pred_exp, .x$expression)) %>%
  as.data.frame() %>%
  dplyr::rename(correlation = .data[["."]]) %>%
  rownames_to_column("topic") %>% 
  # separate(celltype_topic, c("cell_type", "topic")) %>%
  left_join(merge_topics)



corplot_out <- founder_topic_cor %>%
  dplyr::select(cell_type, topic, marker, correlation) %>%
  unique %>%
  ggplot(., aes(x = correlation)) + 
  geom_density(adjust = 5, fill = "lightgray") +
  labs(title = "Correlation with Founder Strain",
       x = "Correlation",
       y = "Density")

ggsave(filename = "results/figures/figure-3g-topic-founder-correlation.pdf",
       plot = corplot_out,
       dpi = 330,
       height = 5,
       width = 7)


founder_topic_cor %>%
  dplyr::select(cell_type, topic, marker, correlation, lods) %>%
  unique %>%
  ggplot(., aes(correlation, lods)) + 
  geom_point()

