
library(data.table)
library(GenomicRanges)

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

library(rtracklayer)
lecif_score_tabb = import("/data/deyk/kushal/Miao_DO/data/LECIF/mm10.LECIFv1.1.bw")

gr1 = GRanges(seqnames = Rle(paste0("chr", tabb$loci_chr)),
              ranges = IRanges(start=tabb$loci_start, end = tabb$loci_end))
gr2 = lecif_score_tabb

cc=findOverlaps(gr1, gr2)
full_score = gr2$score[unique(subjectHits(cc))]

back_mean = mean(full_score)
back_sd = sd(full_score)


modes = names(ll)

score_vec = c()
score_vec_sd = c()

for(numl in 1:length(ll)){
  tabb2 = tabb[ll[[numl]], ]
  gr1 = GRanges(seqnames = Rle(paste0("chr", tabb2$loci_chr)),
                ranges = IRanges(start=tabb2$loci_start, end = tabb2$loci_end))
  gr2 = lecif_score_tabb
  cc=findOverlaps(gr1, gr2)
  score1 = gr2$score[unique(subjectHits(cc))]
  score_vec = c(score_vec, mean(score1))
  score_vec_sd = c(score_vec_sd, sd(score1)/sqrt(22))
  cat("We are at mode:", numl, "\n")
}

ll = list("LECIF.Enrich" = score_vec/back_mean, "LECIF.Enrich.sd" = score_vec_sd)



save(ll, file = "/data/deyk/kushal/Miao_DO/data/LECIF_nonspecific_enrichment_all.rda")



