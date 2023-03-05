# install.packages("googledrive")
library(googledrive)

## download main files to work with
drive_download(as_id("1UZ3KnuQ1OpsDqCefVH9CGhYyeeGiBfZI"), 
               path=paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"))

drive_download(as_id("1md2aNYx_XLE4r2mY7oncVfd5ZrqEnlCt"), 
               path=paste0(data_dir, "references/Mus_musculus.GRCm38.102.gtf"))

drive_download(as_id("1h6IShEYfcKSy3boXd-eSyLlf9Lj0XEmI"), 
               path=paste0(data_dir, "allchannels/vars.csv"))
# 
drive_download(as_id("1B7DT53luVwd2wJXxJqnKPVD7go8gn1RZ"), 
               path=paste0(data_dir, "manuscript-plot-data.csv.gz"))


#
drive_download(as_id("1muSHQ_5UlHz6KexUR4ZL-Zaup6tjUrNh"), 
               path=paste0(results_dir, "qtl-plot-lods-NCR1\\+\\ ILC3-cv.csv.gz"))
