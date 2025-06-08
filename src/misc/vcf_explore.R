# install.packages('vcfR')
library(vcfR)

vcf <- read.vcfR( paste0(data_dir, "genotype/G30_scRNAseq_DOmice_SNPGrp1_QC0.50.vcf"), verbose = FALSE )