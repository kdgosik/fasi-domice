
library(data.table)
tabb = data.frame(fread("/data/deyk/kushal/Miao_DO/trait_by_loci_10kb_window.csv"))
idx = c(grep("_ILC1", tabb$trait), grep("_ILC2", tabb$trait), grep("_ILC3", tabb$trait), grep("_LTi", tabb$trait), grep("_activated", tabb$trait))
idx = unique(idx)
tabb = tabb[-idx, ]

library(rtracklayer)
lecif_score_tabb = import("/data/deyk/kushal/Miao_DO/data/LECIF/mm10.LECIFv1.1.bw")

gr1 = GRanges(seqnames = Rle(paste0("chr", tabb$loci_chr)),
              ranges = IRanges(start=tabb$loci_start, end = tabb$loci_end))
gr2 = lecif_score_tabb

cc=findOverlaps(gr1, gr2)
full_score = gr2$score[unique(subjectHits(cc))]

back_mean = mean(full_score)
back_sd = sd(full_score)

celltypes = c("ILC1", "ILC2", "ILC3", "LTi")
celltypes_idx = c(9, 10, 11, 12)



#############################. Gene Set 1. ##########################################
## comprises of genes that harbor eQTLs in each of 4 ILC subtypes

score_vec = c()
score_vec_sd = c()
for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  gr1 = GRanges(seqnames = Rle(paste0("chr", tabb2$loci_chr)),
                ranges = IRanges(start=tabb2$loci_start, end = tabb2$loci_end))
  gr2 = lecif_score_tabb
  cc=findOverlaps(gr1, gr2)
  score1 = gr2$score[unique(subjectHits(cc))]
  score_vec = c(score_vec, mean(score1))
  score_vec_sd = c(score_vec_sd, sd(score1)/sqrt(length(score1)))
  cat("we are at cell type:", numc, "\n")
}

score_vec1 = score_vec
score_vec_sd1 = score_vec_sd

#############################. Gene Set 2. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes and are also expressed in these cell types

score_vec = c()
score_vec_sd = c()
for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] == 1), ]
  gr1 = GRanges(seqnames = Rle(paste0("chr", tabb3$loci_chr)),
                ranges = IRanges(start=tabb3$loci_start, end = tabb3$loci_end))
  gr2 = lecif_score_tabb
  cc=findOverlaps(gr1, gr2)
  score1 = gr2$score[unique(subjectHits(cc))]
  score_vec = c(score_vec, mean(score1))
  score_vec_sd = c(score_vec_sd, sd(score1)/sqrt(length(score1)))
  cat("we are at cell type:", numc, "\n")
}

score_vec2 = score_vec
score_vec_sd2 = score_vec_sd

#############################. Gene Set 3. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes and are also specifically expressed in these cell types

score_vec = c()
score_vec_sd = c()
for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] == 1), ]
  xx = tabb3[, setdiff(celltypes_idx, celltypes_idx[numc])]
  xx[is.na(xx)] = 0
  idx2 = which(apply(xx, 1, max) == 0)
  tabb4 = tabb3[idx2, ]
  gr1 = GRanges(seqnames = Rle(paste0("chr", tabb4$loci_chr)),
                ranges = IRanges(start=tabb4$loci_start, end = tabb4$loci_end))
  gr2 = lecif_score_tabb
  cc=findOverlaps(gr1, gr2)
  score1 = gr2$score[unique(subjectHits(cc))]
  score_vec = c(score_vec, mean(score1))
  score_vec_sd = c(score_vec_sd, sd(score1)/sqrt(length(score1)))
  cat("we are at cell type:", numc, "\n")
}

score_vec3 = score_vec
score_vec_sd3 = score_vec_sd

#############################. Gene Set 4. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes but not expressed in these cell types

score_vec = c()
score_vec_sd = c()
for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] != 1), ]
  gr1 = GRanges(seqnames = Rle(paste0("chr", tabb3$loci_chr)),
                ranges = IRanges(start=tabb3$loci_start, end = tabb3$loci_end))
  gr2 = lecif_score_tabb
  cc=findOverlaps(gr1, gr2)
  score1 = gr2$score[unique(subjectHits(cc))]
  score_vec = c(score_vec, mean(score1))
  score_vec_sd = c(score_vec_sd, sd(score1)/sqrt(length(score1)))
  cat("we are at cell type:", numc, "\n")
}

score_vec4 = score_vec
score_vec_sd4 = score_vec_sd

#############################. Gene Set 5. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes but not zero expressed in each cell type


score_vec = c()
score_vec_sd = c()
for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  tabb3 = tabb2[which(tabb2[, celltypes_idx[numc]] == 0), ]
  gr1 = GRanges(seqnames = Rle(paste0("chr", tabb3$loci_chr)),
                ranges = IRanges(start=tabb3$loci_start, end = tabb3$loci_end))
  gr2 = lecif_score_tabb
  cc=findOverlaps(gr1, gr2)
  score1 = gr2$score[unique(subjectHits(cc))]
  score_vec = c(score_vec, mean(score1))
  score_vec_sd = c(score_vec_sd, sd(score1)/sqrt(length(score1)))
  cat("we are at cell type:", numc, "\n")
}

score_vec5 = score_vec
score_vec_sd5 = score_vec_sd

#############################. Gene Set 6. ##########################################

## comprises of genes that harbor eQTLs in each of 4 ILC subtypes but not expressed in any cell type

score_vec = c()
score_vec_sd = c()
for(numc in 1:length(celltypes)){
  tabb2 = tabb[grep(paste0(celltypes[numc], "_"), tabb$trait), ]
  xx = tabb2[, celltypes_idx]
  idx = which(apply(xx, 1, max) == 0)
  tabb3 = tabb2[idx, ]
  gr1 = GRanges(seqnames = Rle(paste0("chr", tabb3$loci_chr)),
                ranges = IRanges(start=tabb3$loci_start, end = tabb3$loci_end))
  gr2 = lecif_score_tabb
  cc=findOverlaps(gr1, gr2)
  score1 = gr2$score[unique(subjectHits(cc))]
  score_vec = c(score_vec, mean(score1))
  score_vec_sd = c(score_vec_sd, sd(score1)/sqrt(length(score1)))
  cat("we are at cell type:", numc, "\n")
}

score_vec6 = score_vec
score_vec_sd6 = score_vec_sd

ff = rbind(score_vec1, score_vec2, score_vec3, score_vec4, score_vec5, score_vec6)
ff_sd = rbind(score_vec_sd1, score_vec_sd2, score_vec_sd3, score_vec_sd4, score_vec_sd5, score_vec_sd6)
ff = ff/back_mean
rownames(ff) = rownames(ff_sd) = c("eQTL", "eQTL.expressed", "eQTL.specexpresed", "eQTL.noexpressed", "eQTL.zeroexpressed", "eQTL_noexpressed.anyILC")
colnames(ff) = colnames(ff_sd) = celltypes

ll = list("LECIF.Enrich" = ff, "LECIF.Enrich.sd" = ff_sd)
save(ll, file = "/data/deyk/kushal/Miao_DO/data/LECIF_eQTL_enrichment_all.rda")
