# install.packages("googledrive")
library(googledrive)

project_dir <- "/Users/kirkgosik/Google Drive/My Drive/projects/"
domice_dir <- paste0(project_dir, "domice/paper of QTL study/Revised materials of ILC-QTL paper cell Science format/")

if( dir.exists("/workspace/fasi-domice/") ) {
  
  domice_dir <- "/workspace/fasi-domice/"
  
}

if( dir.exists("/Users/kirkgosik/Documents/projects/fasi-domice/") ) {
  
  domice_dir <- "/Users/kirkgosik/Documents/projects/fasi-domice/"
  
}

data_dir <- paste0(domice_dir, "data/")
results_dir <- paste0(data_dir, "results/")

dir.create(paste0(data_dir, "references"), showWarnings = FALSE)
dir.create(paste0(data_dir, "allchannels"), showWarnings = FALSE)
dir.create(paste0(data_dir, "eqtl"), showWarnings = FALSE)
dir.create(paste0(data_dir, "expression"), showWarnings = FALSE)
dir.create(paste0(data_dir, "cytokine"), showWarnings = FALSE)
dir.create(paste0(data_dir, "proportions"), showWarnings = FALSE)
dir.create(paste0(data_dir, "topics"), showWarnings = FALSE)
dir.create(paste0(data_dir, "genotype"), showWarnings = FALSE)


## load main files
drive_download(as_id("1UZ3KnuQ1OpsDqCefVH9CGhYyeeGiBfZI"), 
               path=paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"))

drive_download(as_id("1md2aNYx_XLE4r2mY7oncVfd5ZrqEnlCt"), 
               path=paste0(data_dir, "references/Mus_musculus.GRCm38.102.gtf"))

drive_download(as_id("1h6IShEYfcKSy3boXd-eSyLlf9Lj0XEmI"), 
               path=paste0(data_dir, "allchannels/vars.csv"))
# 






