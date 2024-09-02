
tt = list.files("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO/")
celltypes = as.character(sapply(tt, function(x) return(strsplit(x, ".txt")[[1]][1])))

ncbi_tabb = read.table("/data/deyk/kushal/extras/NCBI37.3.ensembl.gene.loc", header=T)

ll = list()
for(numl in 1:length(celltypes)){
  tabb = read.table(paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO/", "/", tt[numl]), header=F)
  genes1 = tabb[which(tabb[,2] > 0.8), 1]
  ncbi_genes1_human= ncbi_tabb$NCBI[match(genes1, ncbi_tabb$HGNC)]
  ncbi_genes1_human = ncbi_genes1_human[!is.na(ncbi_genes1_human)]
  ll[[numl]] = ncbi_genes1_human
  cat("We are at celltype:", numl, "\n")
}
names(ll) = celltypes

ll2 = lapply(ll, function(x) paste0(x, collapse = " "))
outdf = cbind(names(ll2), ll2)

write.table(outdf, file = "/data/deyk/kushal/Miao_DO/data/MAGMA/miao_celltypes.set",
            row.names = F, col.names = F, sep = "\t", quote=F)


tt = list.files("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL/")
celltypes = as.character(sapply(tt, function(x) return(strsplit(x, ".txt")[[1]][1])))

ncbi_tabb = read.table("/data/deyk/kushal/extras/NCBI37.3.ensembl.gene.loc", header=T)

ll = list()
for(numl in 1:length(celltypes)){
  tabb = read.table(paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_eQTL/", "/", tt[numl]), header=F)
  genes1 = tabb[which(tabb[,2] > 0.8), 1]
  ncbi_genes1_human= ncbi_tabb$NCBI[match(genes1, ncbi_tabb$HGNC)]
  ncbi_genes1_human = ncbi_genes1_human[!is.na(ncbi_genes1_human)]
  ll[[numl]] = ncbi_genes1_human
  cat("We are at celltype:", numl, "\n")
}
names(ll) = celltypes

ll2 = lapply(ll, function(x) paste0(x, collapse = " "))
outdf = cbind(names(ll2), ll2)

write.table(outdf, file = "/data/deyk/kushal/Miao_DO/data/MAGMA/miao_eQTL_celltypes.set",
            row.names = F, col.names = F, sep = "\t", quote=F)



tt = list.files("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_propQTL/")
celltypes = as.character(sapply(tt, function(x) return(strsplit(x, ".txt")[[1]][1])))

ncbi_tabb = read.table("/data/deyk/kushal/extras/NCBI37.3.ensembl.gene.loc", header=T)

ll = list()
for(numl in 1:length(celltypes)){
  tabb = read.table(paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_propQTL/", "/", tt[numl]), header=F)
  genes1 = tabb[which(tabb[,2] > 0.8), 1]
  ncbi_genes1_human= ncbi_tabb$NCBI[match(genes1, ncbi_tabb$HGNC)]
  ncbi_genes1_human = ncbi_genes1_human[!is.na(ncbi_genes1_human)]
  ll[[numl]] = ncbi_genes1_human
  cat("We are at celltype:", numl, "\n")
}
names(ll) = celltypes

ll2 = lapply(ll, function(x) paste0(x, collapse = " "))
outdf = cbind(names(ll2), ll2)

write.table(outdf, file = "/data/deyk/kushal/Miao_DO/data/MAGMA/miao_propQTL_celltypes.set",
            row.names = F, col.names = F, sep = "\t", quote=F)


tt = list.files("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_Nonspecific/")
celltypes = as.character(sapply(tt, function(x) return(strsplit(x, ".txt")[[1]][1])))

ncbi_tabb = read.table("/data/deyk/kushal/extras/NCBI37.3.ensembl.gene.loc", header=T)

ll = list()
for(numl in 1:length(celltypes)){
  tabb = read.table(paste0("/data/deyk/kushal/Miao_DO/data/Gene_Scores/Miao_DO_Nonspecific/", "/", tt[numl]), header=F)
  genes1 = tabb[which(tabb[,2] > 0.8), 1]
  ncbi_genes1_human= ncbi_tabb$NCBI[match(genes1, ncbi_tabb$HGNC)]
  ncbi_genes1_human = ncbi_genes1_human[!is.na(ncbi_genes1_human)]
  ll[[numl]] = ncbi_genes1_human
  cat("We are at celltype:", numl, "\n")
}
names(ll) = celltypes

ll2 = lapply(ll, function(x) paste0(x, collapse = " "))
outdf = cbind(names(ll2), ll2)

write.table(outdf, file = "/data/deyk/kushal/Miao_DO/data/MAGMA/miao_Nonspecific_celltypes.set",
            row.names = F, col.names = F, sep = "\t", quote=F)
