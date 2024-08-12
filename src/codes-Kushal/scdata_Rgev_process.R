

##########################.   mouseDO_subsample.  ###############################################

library(data.table)
library(Seurat)

xx = get(load("/data/deyk/kushal/Miao_DO/mouseDO_subsample.Rda"))

celltypes1 = c("ILC1", "ILC2", "ILC3", "ILC3(LTi-like)")

celltypes = setdiff(levels(xx), celltypes1)

for(numc in 1:length(celltypes)){
  de.markers <- FindMarkers(xx, ident.1 = celltypes[numc], ident.2 = NULL, only.pos = TRUE)
  de.markers2 <- FindMarkers(xx, ident.1 = celltypes[numc], ident.2 = NULL)

  ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))
  genes_human1 = ortholog_tabb$HGNC.symbol[match(rownames(de.markers), ortholog_tabb$MGI.symbol)]
  genes_human2 = ortholog_tabb$HGNC.symbol[match(rownames(de.markers2), ortholog_tabb$MGI.symbol)]

  genes_human2 = genes_human2[!is.na(genes_human2)]
  genes_human1 = genes_human1[!is.na(genes_human1)]

  write.table(genes_human1, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores/scdata_mouseDO_subsample_", celltypes[numc],
                                         ".express.pos.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  write.table(genes_human2, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores/scdata_mouseDO_subsample_", celltypes[numc],
                                          ".express.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
}


##########################.   mouseENS.ss2.glia  ###############################################

library(data.table)
library(Seurat)

xx = get(load("/data/deyk/kushal/Miao_DO/mouseENS.ss2.glia.Rda"))

celltypes = levels(xx)

for(numc in 1:length(celltypes)){
  de.markers <- FindMarkers(xx, ident.1 = celltypes[numc], ident.2 = NULL, only.pos = TRUE)
  de.markers2 <- FindMarkers(xx, ident.1 = celltypes[numc], ident.2 = NULL)

  ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))
  genes_human1 = ortholog_tabb$HGNC.symbol[match(rownames(de.markers), ortholog_tabb$MGI.symbol)]
  genes_human2 = ortholog_tabb$HGNC.symbol[match(rownames(de.markers2), ortholog_tabb$MGI.symbol)]

  genes_human2 = genes_human2[!is.na(genes_human2)]
  genes_human1 = genes_human1[!is.na(genes_human1)]

  write.table(genes_human1, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores/scdata_mouseENS.ss2.glia_", celltypes[numc],
                                          ".express.pos.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  write.table(genes_human2, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores/scdata_mouseENS.ss2.glia_", celltypes[numc],
                                          ".express.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
}



##########################.   mouseENS.ss2.neur.Rda  ###############################################

library(data.table)
library(Seurat)

xx = get(load("/data/deyk/kushal/Miao_DO/mouseENS.ss2.neur.Rda"))

celltypes = levels(xx)

for(numc in 1:length(celltypes)){
  de.markers <- FindMarkers(xx, ident.1 = celltypes[numc], ident.2 = NULL, only.pos = TRUE)
  de.markers2 <- FindMarkers(xx, ident.1 = celltypes[numc], ident.2 = NULL)

  ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))
  genes_human1 = ortholog_tabb$HGNC.symbol[match(rownames(de.markers), ortholog_tabb$MGI.symbol)]
  genes_human2 = ortholog_tabb$HGNC.symbol[match(rownames(de.markers2), ortholog_tabb$MGI.symbol)]

  genes_human2 = genes_human2[!is.na(genes_human2)]
  genes_human1 = genes_human1[!is.na(genes_human1)]

  write.table(genes_human1, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores/scdata_mouseENS.ss2.neur_", celltypes[numc],
                                          ".express.pos.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  write.table(genes_human2, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores/scdata_mouseENS.ss2.neur_", celltypes[numc],
                                          ".express.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
}



##########################.   mouseintestineimmune.LP.Rda  ###############################################

library(data.table)
library(Seurat)

xx = get(load("/data/deyk/kushal/Miao_DO/mouseintestineimmune.PP.Rda"))

celltypes = levels(xx)

for(numc in 1:length(celltypes)){
  de.markers <- FindMarkers(xx, ident.1 = celltypes[numc], ident.2 = NULL, only.pos = TRUE)
  de.markers2 <- FindMarkers(xx, ident.1 = celltypes[numc], ident.2 = NULL)

  ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))
  genes_human1 = ortholog_tabb$HGNC.symbol[match(rownames(de.markers), ortholog_tabb$MGI.symbol)]
  genes_human2 = ortholog_tabb$HGNC.symbol[match(rownames(de.markers2), ortholog_tabb$MGI.symbol)]

  genes_human2 = genes_human2[!is.na(genes_human2)]
  genes_human1 = genes_human1[!is.na(genes_human1)]

  write.table(genes_human1, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores/scdata_mouseintestineimmune.LP_", celltypes[numc],
                                          ".express.pos.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  write.table(genes_human2, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores/scdata_mouseintestineimmune.LP_", celltypes[numc],
                                          ".express.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
}

