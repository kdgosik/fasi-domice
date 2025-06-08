library(qtl2)
library(magrittr)
library(data.table)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)

# get arguments
args <- commandArgs(trailingOnly = TRUE)
# args <- 1

pheno_file <- args[1]
geno_file <- args[2]
marker_file <- args[3]
map_file <- args[4]
cores <- args[5]


## Project Path
project_path <- "/ahg/regevdata/projects/FASI_DOmice/"
if( length(dir(project_path)) == 0 ) project_path <- "/Volumes/ahg_regevdata/projects/FASI_DOmice/"

## assigning project paths
geno_path <- paste0(project_path, "genotype/")
pheno_path <- paste0(project_path, "phenotype/")
my_path <- paste0(project_path, "kirk/")



# ## Previous Processing Steps
# 
# cat("Reading phenotype data", paste0(my_path, "data/MeanExpressionCalls_Cluster", args[1], ".csv..."), "\n")
# mean_expression <- fread(paste0(my_path, "data/MeanExpressionCalls_Cluster", args[1], ".csv"), data.table = FALSE)
# rownames(mean_expression) <- mean_expression$V1
# mean_expression <- mean_expression[,-1]
# 
# variance_expression <- fread(paste0(my_path, "data/VarianceExpressionCalls_Cluster", args[1], ".csv"), data.table = FALSE)
# rownames(variance_expression) <- variance_expression$V1
# variance_expression <- variance_expression[,-1]
# 
# genes <- sqrt(variance_expression) / mean_expression
# genes[is.na(genes)] <- 0
# 
#   
# ## Using Just Mean Expression
# # genes <- fread(paste0(my_path, "data/MeanExpressionCalls_Cluster", args[1], ".csv"), data.table = FALSE)
# # rownames(genes) <- genes$V1
# # genes <- genes[,-1]
# 
# 
# ## keep columns with at least 70% expression
# keep_cols <- which(colMeans(genes > 0) > 0.7)




## Load Genotypes (Allele Probs)
cat("Reading allele probs data ...", "\n")
# load(paste0(geno_path, "alleleprobs.Rdata"))
load(geno_file)
nms <- rownames(aprobs[[1]])
#nms <- nms[grep("Broad_Inst_Xu_MURGIGV01_20180706", nms)]
nms <- str_extract(nms, "C[0-9]{1,2}_[0-9]")

## rename sample ids to match phenotype data
for( chr in names(aprobs) ) {
  rownames(aprobs[[chr]]) <- nms
}

## Load Phenotypes (all females)
# phenotype$ratio <- phenotype$IgE/phenotype$IgG1 
pheno <- genes[, keep_cols]
#rownames(pheno) <- phenotype$ID[phenotype$NOTES == "PASS"]
pheno <- read_csv(pheno_file)

## Sync sample names
samples <- intersect(rownames(pheno), rownames(aprobs[[1]]))
samples
length(samples)

pheno <- pheno[samples, ,drop = FALSE]
for( chr in names(aprobs) ) {
  aprobs[[chr]] <- aprobs[[chr]][samples, , ]
}


# Creating kinship matrix
K <- calc_kinship(probs = aprobs, type = "loco", cores = cores)


## Run QTL scan
cat("Running QTL analysis...", "\n")
qtl <- scan1(genoprobs = aprobs, 
             pheno = pheno, 
             kinship = K, 
             cores = cores)

lods <- as.data.frame(qtl)
lods$Marker <- rownames(lods)
lods <- lods[,c("Marker", colnames(lods)[-length(lods)])]
cat("Finished QTL Analysis ... \n")


## Plot Data Section ##################


## Loading Marker Map
cat("Loading SNPs...", "\n")
# load(paste0(my_path, "data/GM_snps.Rdata"))
load(marker_file)
# marker_map <- readRDS(paste0(geno_path, "Regev_map_20171221.rds"))
marker_map <- readRDS(map_file)

xpos <- qtl2:::map_to_xpos(marker_map, gap = 25)
xchrbound <- qtl2:::map_to_boundaries(marker_map, gap = 25)

cat("Loading Genes...", "\n")
ensembl <- fread(paste0(my_path, "data/ensembl.Mus_musculus.GRCm38.93.csv"))
ensembl <- ensembl[chr %in% c(as.character(1:19), "X"), .(start = min(start), chr = min(chr)), by = gene]

## Creating Gene Map
gene_map <- split(ensembl, f = ensembl$chr)
for( nm in names(gene_map) ) {
  vec <- gene_map[[nm]]$start
  names(vec) <- gene_map[[nm]]$gene
  gene_map[[nm]] <- vec
}
gene_map <- gene_map[c(as.character(1:19), "X")]
ypos <- qtl2:::map_to_xpos(gene_map, gap = 25)
ychrbound <- qtl2:::map_to_boundaries(gene_map, gap = 25)

## removing likely non-significant lods
cat("Removing low LOD Scores...", "\n")
lods[lods < 5] <- 0
keep_rows <- which(rowSums(lods[,-1]) != 0)
keep_cols <- c(1, 1 + which(colSums(lods[,-1]) != 0))

lods <- lods[keep_rows, keep_cols]
dim(lods)

cat("Transforming to Long Format ... \n")
lods_long <- melt(lods)
lods_long <- lods_long[lods_long$value > 0, ]
colnames(lods_long) <- c("Marker", "Gene", "value")


cat("Creating Position Reference Data ... \n")
MarkerRef <- data.frame(Marker = names(xpos), xpos = xpos, stringsAsFactors = FALSE) %>%
  left_join(
    {GM_snps %>% select(Marker = marker, MarkerChr = chr)}
  )
GeneRef <- data.frame(Gene = names(ypos), ypos = ypos, stringsAsFactors = FALSE) %>%
  left_join(
    {ensembl %>% select(Gene = gene, GeneChr = chr)}
  )

lods_long <- lods_long %>%
  left_join(MarkerRef) %>%
  left_join(GeneRef) %>%
  mutate(MarkerChr = as.numeric(ifelse(MarkerChr == "X", "20", MarkerChr)), 
         GeneChr = as.numeric(ifelse(GeneChr == "X", "20", GeneChr))) %>%
  filter(!is.na(MarkerChr) & !is.na(GeneChr))

dim(lods_long)

cat("\n", "Writing Plot Data to csv...", "\n")
write.csv(lods_long, paste0(my_path, "data/Cluster", args[1], "_CoefVar_PlotData.csv"), row.names = FALSE)


cat("Plotting LODs...", "\n")
ggplot(lods_long, aes(x = xpos, y = ypos)) +
  # geom_raster() +
  labs(x = "Chromosome of SNP", y = "Chromosome of Gene", title = paste0("GWAS Cluster ", args[1])) +
  geom_point(shape = 46) +
  scale_x_continuous(breaks = apply(xchrbound, 2, median), label = c(as.character(1:19), "X")) + 
  scale_y_continuous(breaks = apply(ychrbound, 2, median), label = c(as.character(1:19), "X")) + 
  theme_minimal() 
ggsave(paste0(my_path, "results/Cluster", args[1], "_SNP_Gene_CoefVar.png"))
