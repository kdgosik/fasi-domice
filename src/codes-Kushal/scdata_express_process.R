

##########################.   mouseDO_subsample.  ###############################################

library(data.table)
library(Seurat)

xx = get(load("/data/deyk/kushal/Miao_DO/mouseDO_subsample.Rda"))

celltypes1 = c("ILC1", "ILC2", "ILC3", "ILC3(LTi-like)")

celltypes = setdiff(levels(xx), celltypes1)

idx_all = c()

for(numc in 1:length(celltypes1)){
  arr = xx@assays$RNA
  idx1 = which(xx@meta.data$called_cell_types_new == celltypes1[numc])
  arr2 = arr[, idx1]
  barr2 = arr2
  barr2[barr2 > 0] = 1

  genes_mouse = rownames(arr2)[which(rowMeans(barr2) > 0.2)]

  ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))
  genes_human = ortholog_tabb$HGNC.symbol[match(genes_mouse, ortholog_tabb$MGI.symbol)]
  genes_human = genes_human[!is.na(genes_human)]

  write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores2/scdata_mouseDO_subsample_", celltypes1[numc],
                                          ".isexpress.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  idx_all = c(idx_all, idx1)

}


arr = xx@assays$RNA
idx2 = setdiff(1:ncol(arr), idx_all)
arr2 = arr[, idx2]
barr2 = arr2
barr2[barr2 > 0] = 1
genes_mouse = rownames(arr2)[which(rowMeans(barr2) > 0.2)]

ortholog_tabb = data.frame(fread("/data/deyk/extras/Orthologs_Yoshida_Mouse_Human.txt"))
genes_human = ortholog_tabb$HGNC.symbol[match(genes_mouse, ortholog_tabb$MGI.symbol)]
genes_human = genes_human[!is.na(genes_human)]

write.table(genes_human, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores2/scdata_mouseDO_subsample_", "non_ILC",
                                       ".isexpress.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)




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

  write.table(genes_human1, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores2/scdata_mouseDO_subsample_", celltypes[numc],
                                          ".DE.pos.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  write.table(genes_human2, file = paste0("/data/deyk/kushal/Miao_DO/Gene_Scores2/scdata_mouseDO_subsample_", celltypes[numc],
                                          ".DE.txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
}


ll = list.files("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Orr")
for(numl in 1:length(ll)){
  xx = data.frame(fread(paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Orr", "/", ll[numl]), header=F))[,1]
  outdff = cbind.data.frame(xx, 1)
  write.table(outdff, file = paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Orr", "/",
                                    ll[numl]),
              row.names = F, col.names = F, sep = "\t", quote=F)
}

Determining that the study has a maximization of benefits and a minimization of risks.
Respect for Persons, Beneficence, Justice

Can qualify as an activity “preparatory to research,” at least for the initial contact, but data should not leave the covered entity.
Supplement those of the Common Rule and FDA.
Identifiable health information that is created or held by covered entities and their business associates.
must be more detailed for disclosures that involve fewer than 50 subject records.
To all human subjects research that uses PHI without an authorization from the data subject.





