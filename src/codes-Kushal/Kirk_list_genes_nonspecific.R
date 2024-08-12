

library(data.table)

ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))


tabb = data.frame(fread("/data/deyk/kushal/Miao_DO/trait_by_loci_10kb_window.csv"))

## prop QTL
idx1 = c(grep("_ILC1", tabb$trait), grep("_ILC2", tabb$trait), grep("_ILC3", tabb$trait), grep("_LTi", tabb$trait))
idx1 = unique(idx1)

## all QTL
idx2 = 1:nrow(tabb)

## topic QTL
idx3 = grep("topic", tabb$trait)

## expression QTL
idx_temp1 = c(grep("_ILC1", tabb$trait), grep("_ILC2", tabb$trait), grep("_ILC3", tabb$trait), grep("_LTi", tabb$trait), grep("_activated", tabb$trait))
idx_temp2 = grep("_", tabb$trait)
idx4 = setdiff(idx_temp2, idx_temp1)


ll = list()
ll[[1]] = idx4
ll[[2]] = idx1
ll[[3]] = idx3
ll[[4]] = idx2

names(ll) = c("expression", "proportion", "topic", "all")

modes = names(ll)

for(numl in 1:length(ll)){

    tabb2 = tabb[ll[[numl]], ]
    genes = unique(tabb2$eQTL_loci_gene_name)
    genes_human = ortholog_tabb$HGNC.symbol[match(genes, ortholog_tabb$MGI.symbol)]
    genes_human = genes_human[!is.na(genes_human)]
    if(length(genes_human) > 30){
    write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_Nonspecific/GENESET1_HARBOR_QTL_", modes[numl], ".txt"),
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
    write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_Nonspecific/GENESET2_HARBOR_QTL_", modes[numl],
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
    write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_Nonspecific/GENESET3_HARBOR_QTL_", modes[numl],
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
    write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_Nonspecific/GENESET4_HARBOR_QTL_", modes[numl],
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
    write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_Nonspecific/GENESET5_HARBOR_QTL_", modes[numl],
                                           "_notexpressed_inanyILC", ".txt"),
                row.names = F, col.names = F, sep = "\t", quote=F)
  }

  cat("We are at mode:", numl, "\n")

}


ll = list.files("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_Nonspecific/")
for(numl in 1:length(ll)){
  xx = data.frame(fread(paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_Nonspecific", "/", ll[numl]), header=F))[,1]
  outdff = cbind.data.frame(xx, 1)
  write.table(outdff, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_Nonspecific", "/",
                                    ll[numl]),
              row.names = F, col.names = F, sep = "\t", quote=F)
}
