BiocManager::install(c("Gviz", "GenVisR", "TxDb.Mmusculus.UCSC.mm10.knownGene", "rtracklayer", 
                       "GenomicFeatures", "GenomicRanges", "BSgenome.Mmusculus.UCSC.mm10",
                       "GenomicAlignments", "ComplexHeatmap"))

remotes::install_github("jokergoo/ComplexHeatmap")

install.packages(c("R.utils", "cowplot", "ggpubr", "circlize", "qtl2", "CMplot"))

