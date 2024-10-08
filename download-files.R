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


# https://drive.google.com/file/d/18hjkwig-CYFzVmkmqqll8ujDNQEvqoCY/view?usp=drive_link
drive_download(as_id("18hjkwig-CYFzVmkmqqll8ujDNQEvqoCY"), 
               path=paste0(data_dir, "allchannels/obs.csv"))






# 
drive_download(as_id("1B7DT53luVwd2wJXxJqnKPVD7go8gn1RZ"), 
               path=paste0(data_dir, "manuscript-plot-data.csv.gz"))

## eQTL #####
# drivedf <- drive_find()
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
               path=paste0(results_dir, "eqtl/qtl-plot-lods-ILC1-cv.csv.gz"),
               overwrite = TRUE)

# qtl-lods-ILC1-cv.csv.gz 1LHAYLCaTVQYeyO7iMsFOvGxIEHdxBxZ4
drive_download(as_id("1LHAYLCaTVQYeyO7iMsFOvGxIEHdxBxZ4"), 
               path=paste0(data_dir, "eqtl/qtl-lods-ILC1-cv.csv.gz"),
               overwrite = TRUE)

# qtl-loci-by-genes-ILC1.csv
drive_download(as_id("1KpsZg4TzAyufVMxXitjAi4OmYhjiptDn"), 
               path=paste0(results_dir, "eqtl/qtl-loci-by-genes-ILC1.csv"),
               overwrite = TRUE)


# qtl-plot-lods-ILC2-cv.csv.gz
drive_download(as_id("1jwea1SSyrrl3bj2CzfUSwWny87pFY3oH"), 
               path=paste0(results_dir, "eqtl/qtl-plot-lods-ILC2-cv.csv.gz"),
               overwrite = TRUE)

# qtl-lods-ILC2-cv.csv.gz 1Rkv7EXmHKp-h-0u-_xfD7KtA1IvXmqz3
drive_download(as_id("1Rkv7EXmHKp-h-0u-_xfD7KtA1IvXmqz3"), 
               path=paste0(data_dir, "eqtl/qtl-lods-ILC2-cv.csv.gz"),
               overwrite = TRUE)

# qtl-loci-by-genes-ILC2.csv
drive_download(as_id("1mq4qtjUeomfDqSZZ8MgrCYKeevkfcIMx"), 
               path=paste0(results_dir, "eqtl/qtl-loci-by-genes-ILC2.csv"),
               overwrite = TRUE)

drive_download(as_id("1muSHQ_5UlHz6KexUR4ZL-Zaup6tjUrNh"), 
               path=paste0(results_dir, "eqtl/qtl-plot-lods-NCR1\\+\\ ILC3-cv.csv.gz"),
               overwrite = TRUE)

# qtl-lods-NCR1+ ILC3-cv.csv.gz        1H4W0zSkcqoA_xuX8v8ByD17_pWBTowfU
drive_download(as_id("1H4W0zSkcqoA_xuX8v8ByD17_pWBTowfU"), 
               path=paste0(data_dir, "eqtl/qtl-lods-NCR1+ ILC3-cv.csv.gz"),
               overwrite = TRUE)


# qtl-loci-by-genes-ILC3.csv
drive_download(as_id("15Yx9U3HiXWCYdk12aPUevvd7jkxm04ZC"), 
               path=paste0(results_dir, "eqtl/qtl-loci-by-genes-ILC3.csv"),
               overwrite = TRUE)


drive_download(as_id("1F4o0IWLpOXPO9LSxTt6H5oIjKat-mkHk"), 
               path=paste0(results_dir, "eqtl/qtl-plot-lods-Lti ILC3-cv.csv.gz"),
               overwrite = TRUE)


# qtl-lods-Lti ILC3-cv.csv.gz        1Kq-U1RAdbMGTvtjezwigUr7xj2DETJzD
drive_download(as_id("1Kq-U1RAdbMGTvtjezwigUr7xj2DETJzD"), 
               path=paste0(data_dir, "eqtl/qtl-lods-Lti ILC3-cv.csv.gz"),
               overwrite = TRUE)

# qtl-loci-by-genes-LTi.csv
drive_download(as_id("1jBRgqaCcMvWUtyG7lBteGDQBURhaJLCN"), 
               path=paste0(results_dir, "eqtl/qtl-loci-by-genes-LTi.csv"),
               overwrite = TRUE)




## Proportion #########
# name                          id                               
# <chr>                         <drv_id>                         
# 1 gwas-ilc1-ilc3-results.csv.gz 18SPTbXiqx_j-O3pDP84qu-rQVsyYuZjS
# 2 gwas-ilc3-lti-results.csv.gz  1mBLJtausS8iz4gEQF1K-dJ5ToPT5TuLX
# 3 gwas-ilc1-ilc2-results.csv.gz 1zYmd5PYT2jXoRqgJdwRXuj1wE1LuTuMN
# 4 gwas-ilc2-lti-results.csv.gz  1kNnX0hw_Te_inGfCE5IUpu324fPhvk60
# 5 gwas-ilc1-lti-results.csv.gz  1zu5EF_4i2B1LNYaI9tE8Wvrpu29csgzM
# 6 gwas-ilc2-ilc3-results.csv.gz 16xyStzpJ1GAoHn1kvcomV_yTpLoQLoWe



drive_download(as_id("1zYmd5PYT2jXoRqgJdwRXuj1wE1LuTuMN"), 
               path=paste0(results_dir, "proportions/gwas-ilc1-ilc2-results.csv.gz"),
               overwrite = TRUE)


drive_download(as_id("18SPTbXiqx_j-O3pDP84qu-rQVsyYuZjS"), 
               path=paste0(results_dir, "proportions/gwas-ilc1-ilc3-results.csv.gz"),
               overwrite = TRUE)


drive_download(as_id("1zu5EF_4i2B1LNYaI9tE8Wvrpu29csgzM"), 
               path=paste0(results_dir, "proportions/gwas-ilc1-lti-results.csv.gz"),
               overwrite = TRUE)


drive_download(as_id("16xyStzpJ1GAoHn1kvcomV_yTpLoQLoWe"), 
               path=paste0(results_dir, "proportions/gwas-ilc2-ilc3-results.csv.gz"),
               overwrite = TRUE)


drive_download(as_id("1kNnX0hw_Te_inGfCE5IUpu324fPhvk60"), 
               path=paste0(results_dir, "proportions/gwas-ilc2-lti-results.csv.gz"),
               overwrite = TRUE)


drive_download(as_id("1mBLJtausS8iz4gEQF1K-dJ5ToPT5TuLX"), 
               path=paste0(results_dir, "proportions/gwas-ilc3-lti-results.csv.gz"),
               overwrite = TRUE)

# name                            id                               
# <chr>                           <drv_id>                         
# 1 ILC3_stressed_vs_non_qtl.csv.gz 1onABMrSmcUKLEn3-1TMvJuqqV3_AU_Jf
# 2 LTi_stressed_vs_non_qtl.csv.gz  16CyTv0EwsNIXuAhsoXW4LRQtZzvWRQ3S

drive_download(as_id("1onABMrSmcUKLEn3-1TMvJuqqV3_AU_Jf"), 
               path=paste0(results_dir, "proportions/ILC3_stressed_vs_non_qtl.csv.gz"),
               overwrite = TRUE)


drive_download(as_id("16CyTv0EwsNIXuAhsoXW4LRQtZzvWRQ3S"), 
               path=paste0(results_dir, "proportions/LTi_stressed_vs_non_qtl.csv.gz"),
               overwrite = TRUE)


## Topic #########
# name                              id                               
# <chr>                             <drv_id>                         
#   1 qtl-topic-lods.csv.gz         1u5ZuaI6YDRFxoWeYF252jTTL7fgRhwYn

drive_download(as_id("1u5ZuaI6YDRFxoWeYF252jTTL7fgRhwYn"), 
               path=paste0(results_dir, "topics/qtl-topic-lods.csv.gz"),
               overwrite = TRUE)



## Cytokines ###########
# name                             id                               
# <chr>                            <drv_id>                         
#   1 qtl-cytokines-steady-lods.csv.gz 1T2GDRh9cTyG_DX3C8mVj-Bzfxm9aNRcf

drive_download(as_id("1T2GDRh9cTyG_DX3C8mVj-Bzfxm9aNRcf"), 
               path=paste0(results_dir, "cytokines/qtl-cytokines-steady-lods.csv.gz"),
               overwrite = TRUE)


## genotypes #########

# G30_scRNAseq_DOmice_SNPGrp1_QC0.50.vcf  1R4IAU1zCiSSWwOI11ivg2CaNHp1Vox85
# Batch4_SNPGrp5_QC0.60.vcf 1cgj0lGRSvjF1kZNkPJDa7Aoi0E51e7or
# Batch4_SNPGrp4_QC0.60.vcf 1UInfnKAIq9ngiO-dvyM0nO6yTWUVew9V
# Batch4_SNPGrp3_QC0.60.vcf 1w9ZmcyxusO9BUZ-pZUeRqQLNYTmy3tXA
# Batch4_SNPGrp2_QC0.60.vcf  1QFU4NTIvrqSzrZrunrJynHnBWuWnQYQd
# Batch4_SNPGrp1_QC0.60.vcf  1MDjvyZFfZzh17pTs5p2WfLT8WMb9kv6l

drive_download(as_id("1R4IAU1zCiSSWwOI11ivg2CaNHp1Vox85"), 
               path=paste0(data_dir, "genotype/G30_scRNAseq_DOmice_SNPGrp1_QC0.50.vcf"),
               overwrite = TRUE)


# GM_geno.rdata 1mRIV9KYsPvGTR5hZAg5zyDzrnA4N3PYL
drive_download(as_id("1mRIV9KYsPvGTR5hZAg5zyDzrnA4N3PYL"), 
               path=paste0(data_dir, "genotype/GM_geno.rdata"),
               overwrite = TRUE)

# alleleprobs.Rdata  1VYk5rr9j04TKTfx8aG0NR6ZBshYAWj-D
drive_download(as_id("1VYk5rr9j04TKTfx8aG0NR6ZBshYAWj-D"), 
               path=paste0(data_dir, "genotype/alleleprobs.Rdata"),
               overwrite = TRUE)


# Regev_genoprobs_20171221.rds -  1NPnWooZC_p5HOaX9LM8x_QWv6ROSsSrB
drive_download(as_id("1NPnWooZC_p5HOaX9LM8x_QWv6ROSsSrB"), 
               path=paste0(data_dir, "genotype/Regev_genoprobs_20171221.rds"),
               overwrite = TRUE)




## phenotype files ######

# 3_immune_phenotype_data_four_batches_v2.xlsx
drive_download(as_id("1tcH3xQKI2oDXyFJa39hYZqCUwNq8Gofk"), 
               path=paste0(data_dir, "phenotype/3_immune_phenotype_data_four_batches_v2.xlsx"),
               overwrite = TRUE)


# phenotype_data_batch1_2_v1.csv
drive_download(as_id("1v6naph8y1FtQwPF1gNQ7_O-Tb1qqFR6Q"), 
               path=paste0(data_dir, "phenotype/phenotype_data_batch1_2_v1.csv"),
               overwrite = TRUE)


# phenotype_data_four_batches_v1.xlsx
drive_download(as_id("1I6pHL-ImsPZswQtsbFdYihHN3jswG---"), 
               path=paste0(data_dir, "phenotype/phenotype_data_four_batches_v1.xlsx"),
               overwrite = TRUE)




## founders data ######
# paste0(data_dir, "founders1-louvain_labels-labels-transfer")
# X.csv
## https://drive.google.com/file/d/1LALaZ8IX7qIQJtJHKsZxhdbk-iOMrEsj/view?usp=drive_link
drive_download(as_id("1LALaZ8IX7qIQJtJHKsZxhdbk-iOMrEsj"), 
               path=paste0(data_dir, "founders1-louvain_labels-labels-transfer/X.csv"),
               overwrite = TRUE)


# paste0(data_dir, "founders2-louvain_labels-labels-transfer")
# X.csv:
## https://drive.google.com/file/d/14Lv8RkSuAD3du9o9Xcb4YxNFScyLZEBu/view?usp=drive_link
drive_download(as_id("14Lv8RkSuAD3du9o9Xcb4YxNFScyLZEBu"), 
               path=paste0(data_dir, "founders2-louvain_labels-labels-transfer/X.csv"),
               overwrite = TRUE)





# eqtl-ILC1-founder-coefs.csv.gz    1BpgljQaHHzMTuVIoCRPq5ePw3LHIJLmX
# eqtl-LTi-founder-coefs.csv.gz     1OFYUHwgniE-AnuHABP8v7Glif6iO33e_
# eqtl-ILC3-founder-coefs.csv.gz    1DE37preC92Hzb3lZYRDE8XX1zRYezSZv
# eqtl-ILC2-founder-coefs.csv.gz    19OuIh0BcffFuPhOI0YLe6fDVUlS-yq0b
## ILC1
drive_download(as_id("1BpgljQaHHzMTuVIoCRPq5ePw3LHIJLmX"), 
               path=paste0(results_dir, "founders/eqtl-ILC1-founder-coefs.csv.gz"),
               overwrite = TRUE)

## ILC2
drive_download(as_id("19OuIh0BcffFuPhOI0YLe6fDVUlS-yq0b"), 
               path=paste0(results_dir, "founders/eqtl-ILC2-founder-coefs.csv.gz"),
               overwrite = TRUE)

## ILC3
drive_download(as_id("1DE37preC92Hzb3lZYRDE8XX1zRYezSZv"), 
               path=paste0(results_dir, "founders/eqtl-ILC3-founder-coefs.csv.gz"),
               overwrite = TRUE)

## LTi
drive_download(as_id("1OFYUHwgniE-AnuHABP8v7Glif6iO33e_"), 
               path=paste0(results_dir, "founders/eqtl-LTi-founder-coefs.csv.gz"),
               overwrite = TRUE)



# trait_by_loci_no_window.csv
drive_download(as_id("1nN7nkDeBM_SIX51jYS82ioFWAZCQy6SL"), 
               path=paste0(results_dir, "trait_by_loci_no_window.csv"),
               overwrite = TRUE)

# trait_by_loci_10kb_window.csv
drive_download(as_id("1bn4LbBI613GgDhNvEIFMZqR69D_0Opuw"), 
               path=paste0(results_dir, "trait_by_loci_10kb_window.csv"),
               overwrite = TRUE)


# trait_by_loci_50kb_window.csv
drive_download(as_id("1Kol5vTqNqQBSlUJKFXshZGLXlThxXL8B"), 
               path=paste0(results_dir, "trait_by_loci_50kb_window.csv"),
               overwrite = TRUE)


# trait_by_loci_500kb_window.csv
drive_download(as_id("176ldRbzb3vOON_eUjE3CkVG5VBof1eNS"), 
               path=paste0(results_dir, "trait_by_loci_500kb_window.csv"),
               overwrite = TRUE)



