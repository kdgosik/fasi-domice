install.packages(c("BiocManager","R.utils", "cowplot", "ggpubr", "circlize", "qtl2", "CMplot"))

BiocManager::install(c("XVector","SparseArray","Gviz", "GenVisR", "TxDb.Mmusculus.UCSC.mm10.knownGene", "rtracklayer", 
                       "GenomicFeatures", "GenomicRanges", "BSgenome.Mmusculus.UCSC.mm10",
                       "GenomicAlignments", "ComplexHeatmap"))

remotes::install_github("jokergoo/ComplexHeatmap")


