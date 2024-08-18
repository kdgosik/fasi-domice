#' Figure 5C
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



## eQTL #################
# eqtl_df <- fread(paste0(my_path, "results/eqtl-significant-loci.csv"), data.table = FALSE) %>%
#   # dplyr::filter(cell_type == celltype) %>% 
#   dplyr::group_by(cell_type, gene, marker_chr) %>% 
#   dplyr::arrange(desc(value_adj)) %>% 
#   mutate(min_rank = min_rank(value_adj), n = n()) %>%
#   dplyr::filter(min_rank == n)


library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)

vars <- fread(paste0(data_dir, "allchannels/vars.csv"), data.table=FALSE)
eqtl_genes <-unique(c(vars$index[vars$ilc1_egenes_cv == 1],
                      vars$index[vars$ilc2_egenes_cv == 1],
                      vars$index[vars$ilc3_egenes_cv == 1],
                      vars$index[vars$lti_egenes_cv == 1]))

ccre <- fread(paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"), data.table=FALSE)
ilc1_ccre <- ccre$marker[ccre$ilc1_eqtl_loci_cv==1]
ilc2_ccre <- ccre$marker[ccre$ilc2_eqtl_loci_cv==1]
ilc3_ccre <- ccre$marker[ccre$ilc3_eqtl_loci_cv==1]
lti_ccre <- ccre$marker[ccre$lti_eqtl_loci_cv==1]
loci_markers <- unique(c(ilc1_ccre, ilc2_ccre, ilc3_ccre, lti_ccre))

  
outdf <- lapply(c("ILC1", "ILC2", "ILC3", "LTi"), function(celltype) {
  
  in_file <- paste0("zcat ", domice_dir, "results/founders/eqtl-", celltype, "-founder-coefs.csv.gz")
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


f1_cnames <- readLines(paste0(data_dir, "founders1-louvain_labels-labels-transfer/X.csv"), n=1)
f1_cnames <- unlist(strsplit(x = f1_cnames, split = ","))
f1_cnames <- intersect(eqtl_genes, f1_cnames)
f2_cnames <- readLines(paste0(data_dir, "founders2-louvain_labels-labels-transfer/X.csv"), n=1)
f2_cnames <- unlist(strsplit(x = f2_cnames, split = ","))
f2_cnames <- intersect(eqtl_genes, f2_cnames)
genelist <- intersect(f1_cnames, f2_cnames)
genelist <- genelist[-c(1,2)]

f1X <- fread(paste0(data_dir, "founders1-louvain_labels-labels-transfer/X.csv"),
             select = c("index", "best_clean", "SNG.1ST","called_cell_types_new", genelist)) %>%
  dplyr::filter(called_cell_types_new %in% c("ILC1", "ILC2", "ILC3", "ILC3(LTi-like)"))

f2X <- fread(paste0(data_dir, "founders2-louvain_labels-labels-transfer/X.csv"),
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

ggsave(filename = "results/figures/figure-5C-eqtl-founderstrain.png",
       plot = p1,
       dpi = 330)

ggsave(filename = "results/figures/figure-5C-eqtl-founderstrain.pdf",
       plot = p1,
       dpi = 330)

