#############################. Gene Set 1. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes

library(data.table)
tabb = data.frame(fread("/data/deyk/kushal/Miao_DO/trait_by_loci_10kb_window.csv"))
idx = c(grep("_ILC1", tabb$trait), grep("_ILC2", tabb$trait), grep("_ILC3", tabb$trait), grep("_LTi", tabb$trait), grep("_activated", tabb$trait))
idx = unique(idx)
tabb = tabb[-idx, ]

celltypes = c("ILC1", "ILC2", "ILC3", "LTi")

ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  genes = unique(tabb2$eQTL_loci_gene_name)
  genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
  genes_human = genes_human[!is.na(genes_human)]
  write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL/GENESET1_HARBOR_eQTL_", celltypes[numc], ".txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  cat("We are at cell type:", numc, "\n")
}


#############################. Gene Set 2. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes and are also expressed in these cell types

library(data.table)
tabb = data.frame(fread("/data/deyk/kushal/Miao_DO/trait_by_loci_10kb_window.csv"))
idx = c(grep("_ILC1", tabb$trait), grep("_ILC2", tabb$trait), grep("_ILC3", tabb$trait), grep("_LTi", tabb$trait), grep("_activated", tabb$trait))
idx = unique(idx)
tabb = tabb[-idx, ]

celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
celltypes_idx = c(9, 10, 11, 12)

ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] == 1), ]
  genes = unique(tabb3$eQTL_loci_gene_name)
  genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
  genes_human = genes_human[!is.na(genes_human)]
  write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL/GENESET1_HARBOR_eQTL_expressed_", celltypes[numc], ".txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  cat("We are at cell type:", numc, "\n")
}

#############################. Gene Set 3. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes and are also specifically expressed in these cell types

library(data.table)
tabb = data.frame(fread("/data/deyk/kushal/Miao_DO/trait_by_loci_10kb_window.csv"))
idx = c(grep("_ILC1", tabb$trait), grep("_ILC2", tabb$trait), grep("_ILC3", tabb$trait), grep("_LTi", tabb$trait), grep("_activated", tabb$trait))
idx = unique(idx)
tabb = tabb[-idx, ]

celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
celltypes_idx = c(9, 10, 11, 12)

ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] == 1), ]
  xx = tabb3[, setdiff(celltypes_idx, celltypes_idx[numc])]
  xx[is.na(xx)] = 0
  idx2 = which(apply(xx, 1, max) == 0)
  tabb4 = tabb3[idx2, ]
  genes = unique(tabb4$eQTL_loci_gene_name)
  genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
  genes_human = genes_human[!is.na(genes_human)]
  write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL/GENESET2_HARBOR_eQTL_spec_expressed_", celltypes[numc], ".txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  cat("We are at cell type:", numc, "\n")
}


#############################. Gene Set 4. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes but not expressed in these cell types

library(data.table)
tabb = data.frame(fread("/data/deyk/kushal/Miao_DO/trait_by_loci_10kb_window.csv"))
idx = c(grep("_ILC1", tabb$trait), grep("_ILC2", tabb$trait), grep("_ILC3", tabb$trait), grep("_LTi", tabb$trait), grep("_activated", tabb$trait))
idx = unique(idx)
tabb = tabb[-idx, ]

celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
celltypes_idx = c(9, 10, 11, 12)

ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] != 1), ]
  genes = unique(tabb3$eQTL_loci_gene_name)
  genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
  genes_human = genes_human[!is.na(genes_human)]
  write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL/GENESET3_HARBOR_eQTL_notexpressed_", celltypes[numc], ".txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  cat("We are at cell type:", numc, "\n")
}


#############################. Gene Set 5. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes but not zero expressed in each cell type


library(data.table)
tabb = data.frame(fread("/data/deyk/kushal/Miao_DO/trait_by_loci_10kb_window.csv"))
idx = c(grep("_ILC1", tabb$trait), grep("_ILC2", tabb$trait), grep("_ILC3", tabb$trait), grep("_LTi", tabb$trait), grep("_activated", tabb$trait))
idx = unique(idx)
tabb = tabb[-idx, ]

celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
celltypes_idx = c(9, 10, 11, 12)

ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] == 0), ]
  genes = unique(tabb3$eQTL_loci_gene_name)
  genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
  genes_human = genes_human[!is.na(genes_human)]
  write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL/GENESET3_HARBOR_eQTL_zeroexpressed_", celltypes[numc], ".txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  cat("We are at cell type:", numc, "\n")
}

#############################. Gene Set 6. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes but not expressed in any cell type


library(data.table)
tabb = data.frame(fread("/data/deyk/kushal/Miao_DO/trait_by_loci_10kb_window.csv"))
idx = c(grep("_ILC1", tabb$trait), grep("_ILC2", tabb$trait), grep("_ILC3", tabb$trait), grep("_LTi", tabb$trait), grep("_activated", tabb$trait))
idx = unique(idx)
tabb = tabb[-idx, ]

celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
celltypes_idx = c(9, 10, 11, 12)

ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  xx = tabb2[, celltypes_idx]
  idx = which(apply(xx, 1, max) == 0)
  tabb3 = tabb2[idx, ]
  genes = unique(tabb3$eQTL_loci_gene_name)
  genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
  genes_human = genes_human[!is.na(genes_human)]
  write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL/GENESET4_HARBOR_eQTL_notexpressed_inanyILC_",
                                         celltypes[numc],  ".txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  cat("We are at cell type:", numc, "\n")
}


tabb = data.frame(fread("/data/deyk/kushal/Miao_DO/trait_by_loci_10kb_window.csv"))
celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
celltypes_idx = c(9, 10, 11, 12)
ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

tabb2 = tabb

genes_human_merged = c()
for(numc in 1:length(celltypes)){
  #xx = tabb2[grep(paste0(celltypes[numc], "_"), tabb$trait),celltypes_idx]
  xx = tabb2[,celltypes_idx]

  idx = which(apply(xx, 1, max) == 0)
  tabb3 = tabb2[idx, ]
  genes = unique(tabb3$eQTL_loci_gene_name)
  genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
  genes_human = genes_human[!is.na(genes_human)]
  genes_human_merged = c(genes_human_merged, genes_human)
}


ll = list.files("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL")
for(numl in 1:length(ll)){
  xx = data.frame(fread(paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL", "/", ll[numl]), header=F))[,1]
  outdff = cbind.data.frame(xx, 1)
  write.table(outdff, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL", "/",
                                    ll[numl]),
              row.names = F, col.names = F, sep = "\t", quote=F)
}


