
annot_dir = "/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Miao_DO_eQTL"
out_dir = "/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Joint_Model_eQTL"
ll= list.files(annot_dir)

for(numchr in 1:22){
  merged_annot = c()
  for(numl in 1:length(ll)){
    dff = data.frame(fread(paste0(annot_dir, "/", ll[numl], "/", "ABC_Road_GI_BLD", "/", "ABC_Road_GI_BLD", ".", numchr, ".annot.gz")))
    merged_annot = cbind(merged_annot, dff[,5])
  }
  outdff = cbind.data.frame(dff[,1:4], merged_annot)
  colnames(outdff) = c("CHR", "BP", "SNP", "CM", ll)

  if(!dir.exists(paste0(out_dir, "/", "FULL_COMBO"))){
    dir.create(paste0(out_dir, "/", "FULL_COMBO"))
  }

  write.table(outdff, file = gzfile(paste0(out_dir, "/",
                                           "FULL_COMBO", "/",
                                           "FULL_COMBO", ".",
                                           numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("we are at chr:", numchr, "\n")
}




ll= list.files(annot_dir)
ll = ll[grep("eQTL_expressed", ll)]

for(numchr in 1:22){
  merged_annot = c()
  for(numl in 1:length(ll)){
    dff = data.frame(fread(paste0(annot_dir, "/", ll[numl], "/", "ABC_Road_GI_BLD", "/", "ABC_Road_GI_BLD", ".", numchr, ".annot.gz")))
    merged_annot = cbind(merged_annot, dff[,5])
  }
  outdff = cbind.data.frame(dff[,1:4], merged_annot)
  colnames(outdff) = c("CHR", "BP", "SNP", "CM", ll)

  if(!dir.exists(paste0(out_dir, "/", "eQTL_expressed_COMBO"))){
    dir.create(paste0(out_dir, "/", "eQTL_expressed_COMBO"))
  }

  write.table(outdff, file = gzfile(paste0(out_dir, "/",
                                           "eQTL_expressed_COMBO", "/",
                                           "eQTL_expressed_COMBO", ".",
                                           numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("we are at chr:", numchr, "\n")
}


ll= list.files(annot_dir)
out_dir = "/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Joint_Model_eQTL"

ll = ll[grep("spec_expressed", ll)]

for(numchr in 1:22){
  merged_annot = c()
  for(numl in 1:length(ll)){
    dff = data.frame(fread(paste0(annot_dir, "/", ll[numl], "/", "ABC_Road_GI_BLD", "/", "ABC_Road_GI_BLD", ".", numchr, ".annot.gz")))
    merged_annot = cbind(merged_annot, dff[,5])
  }
  outdff = cbind.data.frame(dff[,1:4], merged_annot)
  colnames(outdff) = c("CHR", "BP", "SNP", "CM", ll)

  if(!dir.exists(paste0(out_dir, "/", "eQTL_spec_expressed_COMBO"))){
    dir.create(paste0(out_dir, "/", "eQTL_spec_expressed_COMBO"))
  }

  write.table(outdff, file = gzfile(paste0(out_dir, "/",
                                           "eQTL_spec_expressed_COMBO", "/",
                                           "eQTL_spec_expressed_COMBO", ".",
                                           numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("we are at chr:", numchr, "\n")
}


ll= list.files(annot_dir)
out_dir = "/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Joint_Model_eQTL"

ll = ll[grep("notexpressed_inanyILC", ll)]

for(numchr in 1:22){
  merged_annot = c()
  for(numl in 1:length(ll)){
    dff = data.frame(fread(paste0(annot_dir, "/", ll[numl], "/", "ABC_Road_GI_BLD", "/", "ABC_Road_GI_BLD", ".", numchr, ".annot.gz")))
    merged_annot = cbind(merged_annot, dff[,5])
  }
  outdff = cbind.data.frame(dff[,1:4], merged_annot)
  colnames(outdff) = c("CHR", "BP", "SNP", "CM", ll)

  if(!dir.exists(paste0(out_dir, "/", "eQTL_notexpressed_anyILC_COMBO"))){
    dir.create(paste0(out_dir, "/", "eQTL_notexpressed_anyILC_COMBO"))
  }

  write.table(outdff, file = gzfile(paste0(out_dir, "/",
                                           "eQTL_notexpressed_anyILC_COMBO", "/",
                                           "eQTL_notexpressed_anyILC_COMBO", ".",
                                           numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("we are at chr:", numchr, "\n")
}


ll= list.files(annot_dir)
out_dir = "/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Joint_Model_eQTL"

ll = ll[grep("zeroexpressed", ll)]

for(numchr in 1:22){
  merged_annot = c()
  for(numl in 1:length(ll)){
    dff = data.frame(fread(paste0(annot_dir, "/", ll[numl], "/", "ABC_Road_GI_BLD", "/", "ABC_Road_GI_BLD", ".", numchr, ".annot.gz")))
    merged_annot = cbind(merged_annot, dff[,5])
  }
  outdff = cbind.data.frame(dff[,1:4], merged_annot)
  colnames(outdff) = c("CHR", "BP", "SNP", "CM", ll)

  if(!dir.exists(paste0(out_dir, "/", "eQTL_zeroexpressed_COMBO"))){
    dir.create(paste0(out_dir, "/", "eQTL_zeroexpressed_COMBO"))
  }

  write.table(outdff, file = gzfile(paste0(out_dir, "/",
                                           "eQTL_zeroexpressed_COMBO", "/",
                                           "eQTL_zeroexpressed_COMBO", ".",
                                           numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("we are at chr:", numchr, "\n")
}







ll= list.files(annot_dir)
out_dir = "/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Joint_Model_eQTL"

ll = ll[grep("LTi", ll)][-4]

for(numchr in 1:22){
  merged_annot = c()
  for(numl in 1:length(ll)){
    dff = data.frame(fread(paste0(annot_dir, "/", ll[numl], "/", "ABC_Road_GI_BLD", "/", "ABC_Road_GI_BLD", ".", numchr, ".annot.gz")))
    merged_annot = cbind(merged_annot, dff[,5])
  }
  outdff = cbind.data.frame(dff[,1:4], merged_annot)
  colnames(outdff) = c("CHR", "BP", "SNP", "CM", ll)

  if(!dir.exists(paste0(out_dir, "/", "LTi_COMBO"))){
    dir.create(paste0(out_dir, "/", "LTi_COMBO"))
  }

  write.table(outdff, file = gzfile(paste0(out_dir, "/",
                                           "LTi_COMBO", "/",
                                           "LTi_COMBO", ".",
                                           numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("we are at chr:", numchr, "\n")
}
