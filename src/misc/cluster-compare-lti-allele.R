library(data.table)
# library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(magrittr)


cat("Reading Utilities ... \n")
# fread of zip files function
fread_zip <- function(zipfile, skip = 9){
  
  # Create a name for the dir where we'll unzip
  zipdir <- tempfile()
  
  # Create the dir using that name
  dir.create(zipdir)
  
  # Unzip the file into the dir
  unzip(zipfile, exdir = zipdir)
  
  # Get a list of files in the dir
  files <- list.files(zipdir, full.names = TRUE)
  
  dat <- fread(files, skip = skip)
  gc(reset = TRUE)
  return(dat)
}

project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"

## assigning project paths
geno_path <- paste0(project_path, "genotype/")
pheno_path <- paste0(project_path, "phenotype/")
my_path <- paste0(project_path, "kirk/")
data_agg_path <- paste0(my_path, "data/aggregate_expression/")


## reading meta data
cat("Reading Meta Data file ... \n")
meta <- fread(paste0(my_path, "data/MetaDataAllChannels.csv"))


## read SNP consequence dataset
cat("Reading SNP Consequence file ... \n")
conseq <- fread(paste0(my_path, "data/reference/GM_SNPS_Consequence.csv"))

## read alleleprobs
load(paste0(geno_path, 'alleleprobs.Rdata'))
nms <- rownames(aprobs[[1]])

## subbset to DO mice for single cell data
samples <- grep("20180706", nms, value = T)
for( chr in names(aprobs) ) {
  aprobs[[chr]] <- aprobs[[chr]][samples, , ]
}


nms <- rownames(aprobs[[1]])
nms <- str_extract(nms, "C[0-9]{1,2}_[0-9]")

## rename sample ids to match phenotype data
for( chr in names(aprobs) ) {
  rownames(aprobs[[chr]]) <- nms
}


## Loading Marker Map
cat("Loading SNPs...", "\n")
## load(paste0(my_path, "data/GM_snps.Rdata"))
marker_map <- readRDS(paste0(geno_path, "Regev_map_20171221.rds"))
xpos <- qtl2:::map_to_xpos(marker_map, gap = 25)
xchrbound <- qtl2:::map_to_boundaries(marker_map, gap = 25)


## Make Marker Reference Dataset
MarkerRef <- data.frame(marker = names(xpos), xpos = xpos, stringsAsFactors = FALSE) %>%
  left_join(
    {conseq %>% select(marker, chr)}
  )


## get SNP report paths
report_paths <- list.files(path = geno_path, 
                           pattern = "FinalReport.zip", 
                           full.names = TRUE, 
                           recursive = TRUE)


## only need genotypes from G30 mice
cat("Reading SNP Data ... \n")
final <- fread_zip(report_paths[5])


## getting mice ids to use later
mice_ids <- unique(final$`Sample ID`)


## intial QC Filter
final <- final[`GC Score` > 0.5, ]
## assigning genotype with Reference Allele A
final[, genotype := as.numeric(`Allele1 - AB` == "A") + as.numeric(`Allele2 - AB` == "A")]
## assigning NA genotype
final[`Allele1 - AB` == "-", genotype := NA]


## making genotype dataframe
cat("Making genotype dataframe ... \n")
geno <- final[, .SD,.SDcols = c("SNP Name", "Sample ID", "genotype")] %>%
  tidyr::spread( `SNP Name`, genotype ) %>%
  rename( SampleID = `Sample ID` ) %>%
  as.data.frame()

rownames(geno) <- geno[,1]



## subsetting to LTi subset
## 5 would be stressed/activated
cat("Creating LTi dataset ... \n")
lti <- meta[louvain_labels %in% c(1,5), ][, .(nonactivated = sum(louvain_labels == 1),
                                              activated = sum(louvain_labels == 5), 
                                              total = .N), by = BestCall] %>%
  data.frame %>%
  # filter(BestCall != "AMB") %>%
  filter(BestCall %in% mice_ids) %>%
  rename(SampleID = BestCall)



## gettting a minor allele frequency of at least 20
maf <- colSums(geno[,-1], na.rm = TRUE)/(2*NROW(geno))
table(maf > 0.2)
geno_subset <- geno[, c(TRUE, maf > 0.2)]



## Running LTi QTL
lti_df <- lti %>% left_join(geno_subset)
var_cols <- apply(lti_df[,-c(1:4)], 2, function(x) {var(x, na.rm = TRUE) != 0})
lti_df <- lti_df[, c(rep(TRUE, 4), var_cols)]


## Creating PCA to add to analysis
to_pca <- lti_df[,5:NCOL(lti_df)]
to_pca[is.na(to_pca)] <- 0
lti_pca <- prcomp(to_pca, center = TRUE, scale = TRUE)
pc_df <- data.frame(lti_pca$x[,1:10])
pc_df$SampleID <- lti_df$SampleID


cat("Running LTi QTL Analysis ... \n")
fit_list_lti <- lapply(unique(conseq$chr), function(chrom) {
  
  markers_used <- intersect(colnames(lti_df)[-c(1:4)], dimnames(aprobs[[chrom]])[[3]])
  
  sapply(markers_used, function(g) {
    
    tmp_df <- lti_df[, c("SampleID", "nonactivated", "activated", "total")] %>%
      left_join(
        {rownames_to_column(data.frame(aprobs[[chrom]][,,g])) %>% 
            rename(SampleID = rowname)}
      ) %>%
      left_join(pc_df)
    
    summary(
      glm(cbind(nonactivated, activated)~0+A+B+C+D+E+F+G+H+PC1+PC2+PC3+PC4+PC5, 
          data = tmp_df, 
          family = binomial("logit"))
    )
    
  }, USE.NAMES = TRUE, simplify = FALSE)
  
})


rnames <- unlist(lapply(seq_along(fit_list_lti), function(i) names(fit_list_lti[[i]])))

lti_betas <- sapply(seq_along(fit_list_lti), function(i) {
  sapply(fit_list_lti[[i]], function(s) coef(s)[,"Estimate"][1:8], 
         USE.NAMES = TRUE, simplify = FALSE)
}, USE.NAMES = TRUE, simplify = FALSE) %>%
  do.call(rbind,.) %>%
  do.call(rbind,.) %>%
  as.data.frame() %>%
  mutate(marker = rnames)


lti_p_vals <- sapply(seq_along(fit_list_lti), function(i) {
  sapply(fit_list_lti[[i]], function(s) coef(s)[,"Pr(>|z|)"][1:8],
         USE.NAMES = TRUE, simplify = FALSE)
}, USE.NAMES = TRUE, simplify = FALSE) %>%
  do.call(rbind, .) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(marker = rnames)

lti_lods <- -log(lti_p_vals[,1:8]) %>%
  mutate(marker = rnames)


lti_output <- merge(lti_betas, lti_p_vals,  by = "marker", suffixes = c("_betas", "_pvals"))
lti_output <- merge(lti_output, lti_lods, by = "marker")

lti_output <- lti_output %>% left_join(conseq) %>% left_join(MarkerRef)

cat("Wrting LTi QTL Output to csv ... \n")
write.csv(lti_output, "LTi_stressed_vs_non_alleles_qtl.csv")


# cat("Plotting LTi LODs...", "\n")
# ggplot(lti_output, aes(x = xpos, y = lti_lods)) +
#   labs(x = "Chromosome of SNP", y = "LOD") +
#   geom_point(shape = 46) +
#   geom_line() +
#   scale_x_continuous(breaks = apply(xchrbound, 2, median), label = c(as.character(1:19), "X")) +
#   theme_minimal()
# ggsave(paste0(my_path, "LTi_stressed_vs_non_qtl.png"))

# library(vcfR)
vcf_file <- paste0(data_dir, "genotype/G30_scRNAseq_DOmice_SNPGrp1_QC0.50.vcf")
# vcf <- read.vcfR(vcf_file, verbose = FALSE )
vcf <- data.frame(t(fread(vcf_file, data.table = FALSE)))
colnames(vcf) <- as.character(make.unique(vcf[3,]))
vcf <- vcf[10:283,] %>%
  mutate(across(.col = everything(), function(x) str_count(str_remove(x, ":.*"),"1")))


geno <- readRDS('/workspace/fasi-domice/data/genotype/Regev_genoprobs_20171221.rds')


lti_stressed <- fread(paste0(results_dir, "proportions/LTi_stressed_vs_non_qtl.csv.gz"),
                      data.table = FALSE)

obs <- fread("/workspace/fasi-domice/data/allchannels/obs.csv", data.table=FALSE)

# Rbpj - rs229232435
gdf <- data.frame(aprobs[['5']][,,"UNCHS014395"])
gdf$BestCall <- str_extract(rownames(gdf), "C[0-9]{1,2}_[0-9]{1}")
gdf <- gdf %>% group_by(BestCall) %>% slice_head(n=1) %>% ungroup()


obs <- obs %>% dplyr::left_join(gdf)


obs %>% 
  filter(called_cell_types == "Lti ILC3") %>% 
  ggplot(aes(factor(louvain_labels), A)) + 
    geom_violin()
