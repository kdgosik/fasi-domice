

ll = list.files("/data/deyk/kushal/Miao_DO/data/MAGMA/SETS_MIAO_eQTL", pattern = ".gsa.out")
pvals_mat = c()
for(numl in 1:length(ll)){
  tabb = read.table(paste0("/data/deyk/kushal/Miao_DO/data/MAGMA/SETS_MIAO_eQTL", "/", ll[numl]), header=T)
  pvals_mat = cbind(pvals_mat, tabb$P)
  set_names = tabb$FULL_NAME
}
rownames(pvals_mat) = set_names
colnames(pvals_mat) = ll

write.table(pvals_mat, file = "/data/deyk/kushal/Miao_DO/data/MAGMA/Miao_eQTL_celltype_MAGMA_gsa.out",
            col.names = T, row.names = T, sep = "\t", quote=F)

cutoff = max(as.numeric(pvals_mat)[which(p.adjust(as.numeric(pvals_mat), method = "BH") < 0.05)])
idx = which(pvals_mat < 0.005, arr.ind=T)
xx= c()
for(mm in 1:nrow(idx)){
  xx = rbind(xx, c(rownames(pvals_mat)[idx[mm, 1]], colnames(pvals_mat)[idx[mm, 2]], pvals_mat[idx[mm, 1], idx[mm, 2]]))
}


dff = cbind.data.frame(colnames(pvals_mat), c("Ectopic Pregnancy", "Miscarriage/Stillbirth", "Num. miscsarriages", "Num.pregnancy.terminations",
  "Pregnancy hypertension", "Ectopic Pregnancy", "Gestational hypertension", "Hemorrhage early pregnancy",
  "Excessive vomiting. pregnancy", "Prolonged pregnancy", "Supervision.normal.pregnancy", "Supervision.highrisk.pregnancy"))

write.table(dff, file = "/data/deyk/kushal/Placenta/ext_data/trait_id_mapping.txt", row.names = F, col.names = F)


yy = dff[match(xx[,2], dff[,1]), 2]

xx = cbind.data.frame(xx, yy)



ll = list.files("/data/deyk/kushal/Miao_DO/data/MAGMA/SETS_MIAO_NONSPECIFIC", pattern = ".gsa.out")
pvals_mat = c()
for(numl in 1:length(ll)){
  tabb = read.table(paste0("/data/deyk/kushal/Miao_DO/data/MAGMA/SETS_MIAO_NONSPECIFIC", "/", ll[numl]), header=T)
  pvals_mat = cbind(pvals_mat, tabb$P)
  set_names = tabb$FULL_NAME
}
rownames(pvals_mat) = set_names
colnames(pvals_mat) = ll

write.table(pvals_mat, file = "/data/deyk/kushal/Miao_DO/data/MAGMA/Miao_Nonspecific_celltype_MAGMA_gsa.out",
            col.names = T, row.names = T, sep = "\t", quote=F)

cutoff = max(as.numeric(pvals_mat)[which(p.adjust(as.numeric(pvals_mat), method = "BH") < 0.05)])
idx = which(pvals_mat < 0.005, arr.ind=T)
xx= c()
for(mm in 1:nrow(idx)){
  xx = rbind(xx, c(rownames(pvals_mat)[idx[mm, 1]], colnames(pvals_mat)[idx[mm, 2]], pvals_mat[idx[mm, 1], idx[mm, 2]]))
}



ll = list.files("/data/deyk/kushal/Miao_DO/data/MAGMA/SETS_MIAO_propQTL", pattern = ".gsa.out")
pvals_mat = c()
for(numl in 1:length(ll)){
  tabb = read.table(paste0("/data/deyk/kushal/Miao_DO/data/MAGMA/SETS_MIAO_propQTL", "/", ll[numl]), header=T)
  pvals_mat = cbind(pvals_mat, tabb$P)
  set_names = tabb$FULL_NAME
}
rownames(pvals_mat) = set_names
colnames(pvals_mat) = ll

write.table(pvals_mat, file = "/data/deyk/kushal/Miao_DO/data/MAGMA/Miao_propQTL_celltype_MAGMA_gsa.out",
            col.names = T, row.names = T, sep = "\t", quote=F)

cutoff = max(as.numeric(pvals_mat)[which(p.adjust(as.numeric(pvals_mat), method = "BH") < 0.05)])
idx = which(pvals_mat < 0.005, arr.ind=T)
xx= c()
for(mm in 1:nrow(idx)){
  xx = rbind(xx, c(rownames(pvals_mat)[idx[mm, 1]], colnames(pvals_mat)[idx[mm, 2]], pvals_mat[idx[mm, 1], idx[mm, 2]]))
}
