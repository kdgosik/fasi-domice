# install.packages("googledrive")
library(googledrive)

## download main files to work with
drive_download(as_id("1UZ3KnuQ1OpsDqCefVH9CGhYyeeGiBfZI"), 
               path=paste0(data_dir, "references/GM_SNPS_Consequence_cCRE.csv"),
               overwrite = TRUE)

drive_download(as_id("1md2aNYx_XLE4r2mY7oncVfd5ZrqEnlCt"), 
               path=paste0(data_dir, "references/Mus_musculus.GRCm38.102.gtf"),
               overwrite = TRUE)

drive_download(as_id("1h6IShEYfcKSy3boXd-eSyLlf9Lj0XEmI"), 
               path=paste0(data_dir, "allchannels/vars.csv"))
# 
drive_download(as_id("1B7DT53luVwd2wJXxJqnKPVD7go8gn1RZ"), 
               path=paste0(data_dir, "manuscript-plot-data.csv.gz"))


#
# name                                 id                               
# 1 qtl-plot-lods-ILC2-cv.csv.gz         1jwea1SSyrrl3bj2CzfUSwWny87pFY3oH
# 2 qtl-plot-lods-NCR1+ ILC3-mean.csv.gz 1DXAhEmem-f-eSsp-3xoXsbiojIJ08ex7
# 3 qtl-plot-lods-ILC1-cv.csv.gz         1VmtM7DUdWzPhaJWBBL-TC2FXzvYELgsG
# 4 qtl-plot-lods-NCR1+ ILC3-cv.csv.gz   1YNaVunrFZLqhVl-FdWRog9GQbUqKOzJI
# 5 qtl-plot-lods-ILC2-mean.csv.gz       1vX_xveg-jkRRxW4w7fNR8uKBXPulI7QY
# 6 qtl-plot-lods-Lti ILC3-cv.csv.gz     1F4o0IWLpOXPO9LSxTt6H5oIjKat-mkHk
# 7 qtl-plot-lods-ILC1-mean.csv.gz       16Fa0ou_k9z51pcDlVdHNVicBqtJGd-Dt
# 8 qtl-plot-lods-ILC2-cv.csv.gz         1PV4oZtdZuKt1gGblMm_sGf6L1Gw-8Wko
# 9 qtl-plot-lods-Lti ILC3-mean.csv.gz   1fvMnhkl9qO-rk4ql2VlzSCcI_9OGdzDt
# 10 qtl-plot-lods-NCR1+ ILC3-mean.csv.gz 1dUalfUyGJtE2QBtU3ICd-wLppV9THTGU
# 11 qtl-plot-lods-Lti ILC3-mean.csv.gz   1Yp8UdUwgLY7CozxAm5nrT3_dVPJzbVp7
# 12 qtl-plot-lods-NCR1+ ILC3-cv.csv.gz   1muSHQ_5UlHz6KexUR4ZL-Zaup6tjUrNh

# https://drive.google.com/file/d/1DUZo14WdcNJSbfyA2CC1LqOZ8uTaUeTl/view?usp=share_link
drive_download(as_id("1DUZo14WdcNJSbfyA2CC1LqOZ8uTaUeTl"), 
               path=paste0(results_dir, "eqt/qtl-plot-lods-ILC1-cv.csv.gz"),
               overwrite = TRUE)

drive_download(as_id("1muSHQ_5UlHz6KexUR4ZL-Zaup6tjUrNh"), 
               path=paste0(results_dir, "eqt/qtl-plot-lods-ILC2-cv.csv.gz"),
               overwrite = TRUE)

drive_download(as_id("1muSHQ_5UlHz6KexUR4ZL-Zaup6tjUrNh"), 
               path=paste0(results_dir, "eqt/qtl-plot-lods-NCR1\\+\\ ILC3-cv.csv.gz"),
               overwrite = TRUE)

drive_download(as_id("1muSHQ_5UlHz6KexUR4ZL-Zaup6tjUrNh"), 
               path=paste0(results_dir, "eqt/qtl-plot-lods-Lti ILC3-cv.csv.gz"),
               overwrite = TRUE)




