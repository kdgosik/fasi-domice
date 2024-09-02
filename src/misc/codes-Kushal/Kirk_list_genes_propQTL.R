#############################. Gene Set 1. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes

tabb = data.frame(fread("/data/deyk/kushal/Miao_DO/trait_by_loci_10kb_window.csv"))

ll = list()
modes = c("ILC1_ILC2", "ILC1_ILC3", "ILC2_ILC3", "ILC1_LTi", "ILC2_LTi", "ILC3_LTi", "ILC3_activated", "LTi_activated")
for(numl in 1:length(modes)){
  idx = grep(modes[numl], tabb$trait)
  ll[[numl]] = idx
}
names(ll) = modes




ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

for(numl in 1:length(ll)){

  tabb2 = tabb[ll[[numl]], ]
  genes = unique(tabb2$eQTL_loci_gene_name)
  genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
  genes_human = genes_human[!is.na(genes_human)]
  if(length(genes_human) > 30){
    write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_propQTL/GENESET1_HARBOR_QTL_", modes[numl], ".txt"),
                row.names = F, col.names = F, sep = "\t", quote=F)
  }



  celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
  celltypes_idx = c(9, 10, 11, 12)

  for(numc in 1:length(celltypes)){
    tabb2 = tabb[ll[[numl]], ]
    tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] == 1), ]
    genes = unique(tabb3$eQTL_loci_gene_name)
    genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
    genes_human = genes_human[!is.na(genes_human)]
    if(length(genes_human) > 30){
      write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_propQTL/GENESET2_HARBOR_QTL_", modes[numl],
                                             "_expressed_", celltypes[numc], ".txt"),
                  row.names = F, col.names = F, sep = "\t", quote=F)
      cat("We are at cell type:", numc, "\n")
    }
  }

  celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
  celltypes_idx = c(9, 10, 11, 12)

  ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

  for(numc in 1:length(celltypes)){
    tabb2 = tabb[ll[[numl]], ]
    tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] == 1), ]
    xx = tabb3[, setdiff(celltypes_idx, celltypes_idx[numc])]
    xx[is.na(xx)] = 0
    idx2 = which(apply(xx, 1, max) == 0)
    tabb4 = tabb3[idx2, ]
    genes = unique(tabb4$eQTL_loci_gene_name)
    genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
    genes_human = genes_human[!is.na(genes_human)]
    if(length(genes_human) > 30){
      write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_propQTL/GENESET3_HARBOR_QTL_", modes[numl],
                                             "_spec_expressed_", celltypes[numc], ".txt"),
                  row.names = F, col.names = F, sep = "\t", quote=F)
      cat("We are at cell type:", numc, "\n")
    }
  }

  celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
  celltypes_idx = c(9, 10, 11, 12)

  ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

  for(numc in 1:length(celltypes)){
    tabb2 = tabb[ll[[numl]], ]
    tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] != 1), ]
    genes = unique(tabb3$eQTL_loci_gene_name)
    genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
    genes_human = genes_human[!is.na(genes_human)]
    if(length(genes_human) > 30){
      write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_propQTL/GENESET4_HARBOR_QTL_", modes[numl],
                                             "_notexpressed_", celltypes[numc], ".txt"),
                  row.names = F, col.names = F, sep = "\t", quote=F)
      cat("We are at cell type:", numc, "\n")
    }
  }

  celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
  celltypes_idx = c(9, 10, 11, 12)

  ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))

  tabb2 = tabb[ll[[numl]], ]
  xx = tabb2[, celltypes_idx]
  idx = which(apply(xx, 1, max) == 0)
  tabb3 = tabb2[idx, ]
  genes = unique(tabb3$eQTL_loci_gene_name)
  genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
  genes_human = genes_human[!is.na(genes_human)]

  if(length(genes_human) > 30){
    write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_propQTL/GENESET5_HARBOR_QTL_", modes[numl],
                                           "_notexpressed_inanyILC", ".txt"),
                row.names = F, col.names = F, sep = "\t", quote=F)
  }

  cat("We are at mode:", numl, "\n")

}


ll = list.files("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_propQTL/")
for(numl in 1:length(ll)){
  xx = data.frame(fread(paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_propQTL", "/", ll[numl]), header=F))[,1]
  outdff = cbind.data.frame(xx, 1)
  write.table(outdff, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_propQTL", "/",
                                    ll[numl]),
              row.names = F, col.names = F, sep = "\t", quote=F)
}
