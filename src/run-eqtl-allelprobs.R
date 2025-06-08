library(qtl2)
library(magrittr)
library(data.table)
library(stringr)

# get arguments
args <- commandArgs(trailingOnly = TRUE)
# args <- 1
cat("Reading phenotype data", paste0("data/MeanExpressionCalls_Cluster", args[1], ".csv..."), "\n")
mean_expression <- fread(paste0("data/MeanExpressionCalls_Cluster", args[1], ".csv"), data.table = FALSE)
rownames(mean_expression) <- mean_expression$V1
mean_expression <- mean_expression[,-1]

variance_expression <- fread(paste0("data/VarianceExpressionCalls_Cluster", args[1], ".csv"), data.table = FALSE)
rownames(variance_expression) <- variance_expression$V1
variance_expression <- variance_expression[,-1]

genes <- 
  
## Using Just Mean Expression
# genes <- fread(paste0("data/MeanExpressionCalls_Cluster", args[1], ".csv"), data.table = FALSE)
# rownames(genes) <- genes$V1
# genes <- genes[,-1]


## keep columns with at least 70% expression
keep_cols <- which(colMeans(genes > 0) > 0.7)

## Load Genotypes (Allele Probs)
cat("Reading allele probs data ...", "\n")
load("/ahg/regevdata/projects/FASI_DOmice/genotype/alleleprobs.Rdata")
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


## Sync sample names
samples <- intersect(rownames(pheno), rownames(aprobs[[1]]))
samples
length(samples)

pheno <- pheno[samples, ,drop = FALSE]
for( chr in names(aprobs) ) {
  aprobs[[chr]] <- aprobs[[chr]][samples, , ]
}

# stopifnot(nrow(pheno) == nrow(aprobs[[1]]))
# stopifnot(rownames(pheno) == rownames(aprobs[[1]]))

# load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
load("data/GM_snps.Rdata")
head(GM_snps)
dim(GM_snps)

marker_names <- unlist(sapply(names(aprobs), function(n) dimnames(aprobs[[n]])[[3]]))
markers <- GM_snps[marker_names, ]

# stopifnot(rownames(markers) == marker_names)

map <- readRDS("/ahg/regevdata/projects/FASI_DOmice/genotype/Regev_map_20171221.rds")


# Creating kinship matrix
K <- calc_kinship(probs = aprobs, type = "loco", cores = 8)


## Run QTL scan
cat("Running qtl analysis...", "\n")
qtl <- scan1(genoprobs = aprobs, 
             pheno = pheno, 
             kinship = K, 
             cores = 8)


## save LOD scores as a csv file
out <- as.data.frame(qtl)
cat("Saving LOD Data...", "\n")
write.csv(out, paste0("results/Cluster", args[1], "_LOD_Output_alleleprobs.csv"))
cat("done!", "\n")
