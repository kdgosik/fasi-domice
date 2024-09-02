#' @title Figure S4A Cytokine Heatmap
#' @author Kirk Gosik
#' @description
#'
#'


## Cytokine heatmap

cnames <- c("batch", "cage", "well", "IFNg", "IL5", "TNFa", "IL2", "IL6", "IL4", "IL10",
            "IL9", "IL17A", "IL17F", "IL22", "IL13", "coding", "mouse_id", "remove")
pheno_v1 <- lapply(1:4, function(i) {
  read_excel(path = paste0(pheno_path, "cytokine detection (version 1).xlsx"),
             sheet = i) %>% 
    na.omit
}) %>%
  do.call(rbind, .)
colnames(pheno_v1) <- cnames

# sample_map_df <- NULL
# ## batch3 - Broad_Inst_Xu_MURGIGV01_20180510 or Broad_Inst_Xu_MURGIGV01_20180607
batch3_samples <- pheno_v1$mouse_id[pheno_v1$batch == "Batch 3"]

# ## batch4 - Broad_Inst_Xu_MURGIGV01_20180706
batch4_samples <- pheno_v1$mouse_id[pheno_v1$batch == "Batch 4"]

pheno_bw <- read.csv(paste0(pheno_path, "pheno_bodyweight_processed.csv")) %>%
  dplyr::filter({ ID %in% batch3_samples & 
      genotyping_date %in% c("Broad_Inst_Xu_MURGIGV01_20180510", "Broad_Inst_Xu_MURGIGV01_20180607") } |
        {ID %in% batch4_samples & 
            genotyping_date %in% c("Broad_Inst_Xu_MURGIGV01_20180706")}) %>%
  dplyr::rename(mouse_id = ID)


sample_map_df <- readr::read_csv(paste0(geno_path, "reference-data/DODB_GM_inventory_20180816.csv")) 
sample_map_df_batch3 <- sample_map_df %>% 
  filter(Project == 'XBI_Set3', Original.Mouse.ID %in% batch3_samples)
sample_map_df_batch4 <- sample_map_df  %>% 
  filter(Original.Mouse.ID %in% batch4_samples)

sample_df <- rbind(sample_map_df_batch3, sample_map_df_batch4) %>%
  dplyr::select(mouse_id = Original.Mouse.ID, aprobs_sample = Unique.Sample.ID) %>%
  dplyr::left_join( pheno_bw , by = "mouse_id")

pheno <- pheno_v1 %>%
  dplyr::left_join(sample_df) %>%
  dplyr::filter(remove == 0, !is.na(aprobs_sample)) %>%
  dplyr::mutate(IL13 = ifelse(is.na(as.numeric(IL13)), 41, as.numeric(IL13))) %>%
  as.data.frame

# Load latest version of heatmap.3 function
devtools::source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

# pheno <- read.csv(paste0(pheno_path, "pheno_batch4_processed.csv"))
# pheno$Batches <- factor(pheno$Batches)
# cnames <- colnames(pheno)[c(3:14,16,20, 22:33)]
# options(na.action='na.pass')
# dat_na_removed <- model.matrix(object = as.formula(paste("~ 0 +", paste(cnames, collapse = "+"))), 
#                                  data = pheno,
#                                  contrasts.arg = list(Batches = contrasts(pheno$Batches, contrasts = F)))
# 
#   
# ## removing exluded mice
# pheno <- pheno %>% filter(excluded == 0)

# plotdata <- scale(pheno[, c("IgA_boxcox", "IgE_boxcox", "IgG1_boxcox", "IgG2_boxcox", "IgM_boxcox", "Mmcp1_boxcox")])
# rownames(plotdata) <- pheno$ID

plotdata <- scale(pheno[, 4:15])
plotdata <- pheno[, 4:15]
rownames(plotdata) <- pheno$mouse_id


# myheatmap(plotdata, "Boxcox Data")

bw_colors <- sort(unique(pheno$BW_1st))
names(bw_colors) <- colorRampPalette(brewer.pal(8, "Oranges"))(15)

mouse_color <- c("black","grey", "white")[pheno$Color]
mouse_batch <- c("Red", "Blue", "Green", "Purple")[pheno$Batches]
mouse_bw_colors <- names(bw_colors)[pheno$BW_1st - 12]

rowlabs <- t(cbind(mouse_color, mouse_batch, mouse_bw_colors))
rownames(rowlabs) <- c("Color", "Batch", "Weight")

rowlabs <- t(cbind(mouse_bw_colors))
rownames(rowlabs) <- c("Weight")

main_title = "Cytokine"
par(cex.main=1)
heatmap.3(x = plotdata, 
          na.rm = TRUE, 
          scale="none", 
          dendrogram = "both", 
          margins = c(6,12),
          Rowv = TRUE, 
          Colv = TRUE, 
          RowSideColors = rowlabs, 
          symbreaks = FALSE, 
          key = TRUE, 
          symkey = FALSE,
          density.info = "none", 
          trace = "none", 
          main = main_title, 
          col = rev(colorRampPalette(brewer.pal(8, "RdYlGn"))(25)), 
          labCol = cnames[4:15],
          labRow = FALSE,
          RowSideColorsSize = 2, 
          KeyValueName = "Cytokine")
legend("topright",
       legend = c("Batch1", "Batch2", "Batch3", "Batch4", "", "black mouse","white mouse", "grey mouse", "", "Low Weight", "Medium Weight", "High Weight"),
       fill = c("Red", "Blue", "Green", "Purple", "white", "black", "white", "grey", "white", names(bw_colors)[c(2,8,17)]), 
       border = FALSE, 
       bty = "n", 
       y.intersp = 0.7, 
       cex = 0.7)

