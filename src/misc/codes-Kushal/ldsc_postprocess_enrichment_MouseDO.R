library(data.table)
library(rmeta)
annot_cell = "/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Miao_DO_eQTL"
results_cell = "/data/deyk/kushal/Miao_DO/data/LDSC_RESULTS_hg38/Miao_DO_eQTL/baseline_Epi_hg38/"
annot_names = list.files(results_cell)
annot_idx = 1

source("/data/deyk/kushal/Placenta/codes/rmeta/R/meta.R")

all_traits = read.table("/n/groups/price/kushal/LDSC/ldsc/TASKFILES/sumstats_encode2.txt")[,1]
all_traits = as.character(sapply(all_traits, function(x) return(strsplit(x, ".sumstats")[[1]][1])))

autoimmune_traits = c("Liu2015-CD", "Liu2015-IBD", "Liu2015-UC",
                      "PASS_Celiac", "PASS_Lupus",  "PASS_Primary_biliary_cirrhosis", "PASS_Rheumatoid_Arthritis",
                      "PASS_Type_1_Diabetes", "UKB.AID_Combined.SAIGE", "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED",
                      "UKB.Asthma.SAIGE", "PASS_Multiple_sclerosis", "PASS_Alzheimers_Jansen2019")

blood_bio_traits = c("UKB.Lym.BOLT", "UKB.Eosino.BOLT", "UKB.RBC.BOLT", "UKB.Plt.BOLT",  "UKB.Mono.BOLT", "UKB.Baso.BOLT", "UKB.Neutro.BOLT")

blood_traits = read.table("/data/deyk/kushal/LDSC/TASKFILES/sumstats.txt")[,1]
blood_traits = as.character(sapply(blood_traits, function(x) return(strsplit(x, ".sumstats")[[1]][1])))


run_single_enrichment_analysis = function(annot_cell,
                                 results_cell,
                                 annotation,
                                 traits,
                                 index_in_results=1){
  enrich_table = matrix(0, length(index_in_results), 3)
  cell_path = paste0(annot_cell, "/", annotation, "/", "ABC_Road_GI_BLD")
  res = paste0(results_cell, "/", annotation,  "/", "ABC_Road_GI_BLD", "/", traits[1], ".sumstats.gz.results")
  tab2 = read.table(res,header=T)
  cat("Number of annotations together with baseline : ", nrow(tab2), "\n")
  annot_names = as.character(tab2$Category[index_in_results])
  Mref = 5961159
  for(id in 1:length(index_in_results)){
    meta_enr        = NULL;
    meta_enrstat    = NULL;
    for(trait_id in 1:length(traits)){
      result.file=paste0(results_cell, "/", annotation, "/", "ABC_Road_GI_BLD", "/", traits[trait_id], ".sumstats.gz.results")
      res=read.table(result.file,header=T)
      logfile = paste(results_cell, "/", annotation, "/", "ABC_Road_GI_BLD", "/", traits[trait_id],".sumstats.gz.log", sep="")
      log = read.table(logfile,h=F,fill=T)
      h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
      myenrstat    = (h2g/Mref)*((res[index_in_results[id],3]/res[index_in_results[id],2])-(1-res[index_in_results[id],3])/(1-res[index_in_results[id],2])) #step1
      myenrstat_z  = qnorm(res[index_in_results[id],7]/2) #step2
      myenrstat_sd = myenrstat/myenrstat_z #step3
      meta_enrstat = rbind(meta_enrstat   , c(myenrstat, myenrstat_sd));
      meta_enr     = rbind(meta_enr, c(res[index_in_results[id],5],res[index_in_results[id],6] ));
    }
    test_eni1=meta.summaries(meta_enr[,1], meta_enr[,2],method="random")
    test_eni2=meta.summaries(meta_enrstat[,1], meta_enrstat[,2],method="random")
    cat("Printing enrichment results for annotation:", annot_names[id], "\n")
    cat(test_eni1$summary, " ",test_eni1$se.summary, " ", 2*pnorm(-abs(test_eni2$summary/test_eni2$se.summary)), "\n");
    enrich_table[id, ] = c(test_eni1$summary, test_eni1$se.summary, 2*pnorm(-abs(test_eni2$summary/test_eni2$se.summary)))
  }
  rownames(enrich_table) = annot_names
  return(enrich_table)
}


out3 = c()
for(m in 1:length(annot_names)){
  out3 = rbind(out3, run_single_enrichment_analysis(annot_cell, results_cell, annotation = annot_names[m],
                                             traits = blood_traits, index_in_results = 1))
}

rownames(out3) = annot_names

out4 = c()
for(m in 1:length(annot_names)){
  out4 = rbind(out4, run_single_enrichment_analysis(annot_cell, results_cell, annotation = annot_names[m],
                                             traits = autoimmune_traits))
}

rownames(out4) = annot_names

out5 = c()
for(m in 1:length(annot_names)){
  out5 = rbind(out5, run_single_enrichment_analysis(annot_cell, results_cell, annotation = annot_names[m],
                                                    traits = blood_bio_traits))
}

rownames(out5) = annot_names


ll <- list()
ll[["Blood"]] = out3
ll[["BloodBio"]] = out5
ll[["Autoimmune"]] = out4

save(ll, file = "/data/deyk/kushal/Miao_DO/data/marginal_meta_enrichment_Miao_DO_eQTL.rda")

