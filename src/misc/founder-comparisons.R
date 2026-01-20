
library(qtl2)
library(data.table)
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)

domice_dir <- "/workspace/fasi-domice/"
data_dir <- paste0(domice_dir, "data/")
results_dir <- paste0(domice_dir, "results/")
figure_dir <- paste0(domice_dir, "figures/")


load(paste0(data_dir, "genotype/alleleprobs.Rdata"))

chrom_names <- names(aprobs)


get_columns_above_threshold <- function(df, threshold = 0.4) {
  
  # For each row, return column names with value == 1 and mean > threshold
  apply(df, 1, function(row) {
    cols <- colnames(df)[row > threshold]
    if (length(cols) == 1) {
      paste0(rep(cols,2),  collapse="") }
    else paste0(cols,  collapse="")
    
  })
}


# helper function to process one marker
process_marker <- function(out_df, marker) {
  get_columns_above_threshold(out_df[,,marker]) %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::separate(col = ".", into = c("col1", "col2"), sep = 1) %>%
    tidyr::pivot_longer(cols = c(col1, col2), 
                        names_to = "position", 
                        values_to = "founder")
  
}

l <- map(chrom_names, function(chrom) {
  out_df <- aprobs[[chrom]]
  markers <- dimnames(out_df)[[3]]
  
  df <- map_dfr(markers, ~ process_marker(out_df, .x), .id = "marker")
  
  allele_df <- df %>%
    group_by(rowname, founder) %>%
    count()
  
  write.csv(allele_df, paste0("founder_allele_counts_chr",chrom,".csv"))
})
# names(l) <- chrom_names


# allele_df <- bind_rows(l) %>%
#   group_by(Var1, Var2) %>%
#   count()
# 
# write.csv(allele_df, "founder_allele_counts_chr11_12.csv")


df <- map_df(grep("founder_allele", dir("/workspace/fasi-domice/src/misc"), value = TRUE), read_csv)

