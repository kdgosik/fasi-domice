tissue = "Miao_DO_eQTL"
traits = as.character(read.delim("/data/deyk/kushal/LDSC/TASKFILES/sumstats.txt", header=F)[,1])
traits2 = unique(as.character(sapply(traits, function(x) return(strsplit(x, ".sumstats")[[1]][1]))))
annot_cell = paste0("/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/", tissue)
results_cell = paste0("/data/deyk/kushal/Miao_DO/data/LDSC_RESULTS_hg38/", tissue, "/baseline_Epi_hg38")
annot_modules = list.files(results_cell)

tempp = read.delim(paste0(results_cell, "/", annot_modules[1], "/",
                          traits2[1], "_ldsc_postprocess.txt"), header=T)

tt= array(0, c(length(annot_modules), length(traits2), nrow(tempp), ncol(tempp)))
for(aa in 1:length(annot_modules)){
  for(kk in 1:length(traits2)){
    temp = read.delim(paste0(results_cell, "/", annot_modules[aa], "/",
                             traits2[kk], "_ldsc_postprocess.txt"), header=T)
    outt = as.matrix(temp[, 1:8, drop=F])
    tt[aa, kk, , ] = outt
  }
}
dimnames(tt)[[1]] = annot_modules
dimnames(tt)[[2]] = traits2
dimnames(tt)[[3]] = rownames(temp)
dimnames(tt)[[4]] = c("taustar", "se-taustar", "p-taustar", "E", "se.E", "p.E", "tau.z", "ptau.z")

save(tt, file = paste0("/data/deyk/kushal/Miao_DO/output/sclinker_2023_", tissue, "_hg38_baselineEpi.rda"))



tissue = "Miao_DO_propQTL"
traits = as.character(read.delim("/data/deyk/kushal/LDSC/TASKFILES/sumstats.txt", header=F)[,1])
traits2 = unique(as.character(sapply(traits, function(x) return(strsplit(x, ".sumstats")[[1]][1]))))
annot_cell = paste0("/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/", tissue)
results_cell = paste0("/data/deyk/kushal/Miao_DO/data/LDSC_RESULTS_hg38/", tissue, "/baseline_Epi_hg38")
annot_modules = list.files(results_cell)

tempp = read.delim(paste0(results_cell, "/", annot_modules[1], "/",
                          traits2[1], "_ldsc_postprocess.txt"), header=T)

tt= array(0, c(length(annot_modules), length(traits2), nrow(tempp), ncol(tempp)))
for(aa in 1:length(annot_modules)){
  for(kk in 1:length(traits2)){
    temp = read.delim(paste0(results_cell, "/", annot_modules[aa], "/",
                             traits2[kk], "_ldsc_postprocess.txt"), header=T)
    outt = as.matrix(temp[, 1:8, drop=F])
    tt[aa, kk, , ] = outt
  }
}
dimnames(tt)[[1]] = annot_modules
dimnames(tt)[[2]] = traits2
dimnames(tt)[[3]] = rownames(temp)
dimnames(tt)[[4]] = c("taustar", "se-taustar", "p-taustar", "E", "se.E", "p.E", "tau.z", "ptau.z")

save(tt, file = paste0("/data/deyk/kushal/Miao_DO/output/sclinker_2023_", tissue, "_hg38_baselineEpi.rda"))


tissue = "Miao_DO_Nonspecific"
traits = as.character(read.delim("/data/deyk/kushal/LDSC/TASKFILES/sumstats.txt", header=F)[,1])
traits2 = unique(as.character(sapply(traits, function(x) return(strsplit(x, ".sumstats")[[1]][1]))))
annot_cell = paste0("/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/", tissue)
results_cell = paste0("/data/deyk/kushal/Miao_DO/data/LDSC_RESULTS_hg38/", tissue, "/baseline_Epi_hg38")
annot_modules = list.files(results_cell)

tempp = read.delim(paste0(results_cell, "/", annot_modules[1], "/",
                          traits2[1], "_ldsc_postprocess.txt"), header=T)

tt= array(0, c(length(annot_modules), length(traits2), nrow(tempp), ncol(tempp)))
for(aa in 1:length(annot_modules)){
  for(kk in 1:length(traits2)){
    temp = read.delim(paste0(results_cell, "/", annot_modules[aa], "/",
                             traits2[kk], "_ldsc_postprocess.txt"), header=T)
    outt = as.matrix(temp[, 1:8, drop=F])
    tt[aa, kk, , ] = outt
  }
}
dimnames(tt)[[1]] = annot_modules
dimnames(tt)[[2]] = traits2
dimnames(tt)[[3]] = rownames(temp)
dimnames(tt)[[4]] = c("taustar", "se-taustar", "p-taustar", "E", "se.E", "p.E", "tau.z", "ptau.z")

save(tt, file = paste0("/data/deyk/kushal/Miao_DO/output/sclinker_2023_", tissue, "_hg38_baselineEpi.rda"))


tissue = "Orr"
traits = as.character(read.delim("/data/deyk/kushal/LDSC/TASKFILES/sumstats.txt", header=F)[,1])
traits2 = unique(as.character(sapply(traits, function(x) return(strsplit(x, ".sumstats")[[1]][1]))))
annot_cell = paste0("/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/", tissue)
results_cell = paste0("/data/deyk/kushal/Miao_DO/data/LDSC_RESULTS_hg38/", tissue, "/baseline_Epi_hg38")
annot_modules = list.files(results_cell)

tempp = read.delim(paste0(results_cell, "/", annot_modules[1], "/",
                          traits2[1], "_ldsc_postprocess.txt"), header=T)

tt= array(0, c(length(annot_modules), length(traits2), nrow(tempp), ncol(tempp)))
for(aa in 1:length(annot_modules)){
  for(kk in 1:length(traits2)){
    temp = read.delim(paste0(results_cell, "/", annot_modules[aa], "/",
                             traits2[kk], "_ldsc_postprocess.txt"), header=T)
    outt = as.matrix(temp[, 1:8, drop=F])
    tt[aa, kk, , ] = outt
  }
}
dimnames(tt)[[1]] = annot_modules
dimnames(tt)[[2]] = traits2
dimnames(tt)[[3]] = rownames(temp)
dimnames(tt)[[4]] = c("taustar", "se-taustar", "p-taustar", "E", "se.E", "p.E", "tau.z", "ptau.z")

save(tt, file = paste0("/data/deyk/kushal/Miao_DO/output/sclinker_2023_", tissue, "_hg38_baselineEpi.rda"))

