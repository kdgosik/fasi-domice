library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)

ccre <- read_csv(paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"))


ccrenona <- ccre %>% dplyr::filter(!is.na(chr))

ccregr <- GRanges(
  seqnames = ccrenona$chr,
  ranges = IRanges(ccrenona$pos, end = ccrenona$pos + 1, names = ccrenona$marker),
  strand = "+"
)
mm10 <- import("../fasi-domice/data/references/Mus_musculus.GRCm38.102.gtf")

df <- fread("/workspace/fasi-domice/results/eqtl/qtl-plot-lods-NCR1\\+\\ ILC3-cv.csv.gz",
            data.table = FALSE) %>%
  left_join(dplyr::select(ccre, marker, marker_chr = chr, pos, ensembl_gene)) %>%
  filter(!is.na(pos)) %>%
  mutate(strand = "+")

ilc3eqtlgr <- GRanges(
  seqnames = df$marker_chr,
  ranges = IRanges(df$pos, end = df$pos + 1, names = df$marker),
  strand = "+"
)

ilc3eqtlgrl <- makeGRangesListFromDataFrame(df,
                                     split.field = "gene",
                                     names.field = "marker",
                                     seqnames="marker_chr",
                                     start.field="pos",
                                     end.field="pos",
                                     strand = "strand")



## possible number of QTLs for GigaMUGA
reduce(resize(ccregr,1000)) # 125931
reduce(resize(ccregr,10000)) # 69645
reduce(resize(ccregr,100000)) # 961

## count number of qtls for the first gene (check what to use)
reduce(resize(sort(ilc3eqtlgrl[[1]]),1000))
reduce(resize(sort(ilc3eqtlgrl[[1]]),10000))

## find overlap
GenomicRanges::intersect(reduce(resize(ccregr,10000)), reduce(resize(sort(ilc3eqtlgrl[[1]]),10000)) )