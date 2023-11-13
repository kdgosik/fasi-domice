
## set initial paths
project_dir <- "/Users/kirkgosik/Google\\ Drive/My\\ Drive/projects/"
domice_dir <- paste0(project_dir, "domice/paper\\ of\\ QTL\\ study/Revised\\ materials\\ of\\ ILC-QTL\\ paper\\ cell\\ Science\\ format/")
data_dir <- paste0(domice_dir, "data/")
results_dir <- paste0(data_dir, "results/")

if( dir.exists("/workspace/fasi-domice/") ) {
  
  domice_dir <- "/workspace/fasi-domice/"
  results_dir <- paste0(domice_dir, "results/")
  
}

if( dir.exists("/Users/kirkgosik/Documents/projects/fasi-domice/") ) {
  
  domice_dir <- "/Users/kirkgosik/Documents/projects/fasi-domice/"
  results_dir <- paste0(domice_dir, "results/")
  
}

## set data dir
data_dir <- paste0(domice_dir, "data/")

## create subdirectories if they don't exist
dir.create(paste0(data_dir, "references"), showWarnings = FALSE)
dir.create(paste0(data_dir, "allchannels"), showWarnings = FALSE)
dir.create(paste0(data_dir, "eqtl"), showWarnings = FALSE)
dir.create(paste0(data_dir, "expression"), showWarnings = FALSE)
dir.create(paste0(data_dir, "cytokines"), showWarnings = FALSE)
dir.create(paste0(data_dir, "proportions"), showWarnings = FALSE)
dir.create(paste0(data_dir, "topics"), showWarnings = FALSE)
dir.create(paste0(data_dir, "genotype"), showWarnings = FALSE)
dir.create(paste0(data_dir, "founders1-louvain_labels-labels-transfer"), showWarnings = FALSE)
dir.create(paste0(data_dir, "founders2-louvain_labels-labels-transfer"), showWarnings = FALSE)


dir.create(paste0(results_dir), showWarnings = FALSE)
dir.create(paste0(results_dir, "cytokines"), showWarnings = FALSE)
dir.create(paste0(results_dir, "eqtl"), showWarnings = FALSE)
dir.create(paste0(results_dir, "proportions"), showWarnings = FALSE)
dir.create(paste0(results_dir, "topics"), showWarnings = FALSE)
dir.create(paste0(results_dir, "figures"), showWarnings = FALSE)
dir.create(paste0(results_dir, "founders"), showWarnings = FALSE)




