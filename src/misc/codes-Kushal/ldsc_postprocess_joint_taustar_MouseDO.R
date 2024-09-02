

###########################  many annots tau star  ##################################

library(data.table)
library(rmeta)
annot_cell = "/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Joint_Model"
results_cell = "/data/deyk/kushal/Miao_DO/data/LDSC_RESULTS_hg38/Joint_Model/"
base_path = "/data/deyk/kushal/ENCODE_Flagship2023/ANNOTATIONS/Baselines/baseline_Epi_hg38"
annot_names = "LTi_COMBO"

autoimmune_traits = c("Liu2015-CD", "Liu2015-IBD", "Liu2015-UC",
                      "PASS_Celiac", "PASS_Lupus",  "PASS_Primary_biliary_cirrhosis", "PASS_Rheumatoid_Arthritis",
                      "PASS_Type_1_Diabetes", "UKB.AID_Combined.SAIGE", "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED",
                      "UKB.Asthma.SAIGE", "PASS_Multiple_sclerosis", "PASS_Alzheimers_Jansen2019")

blood_bio_traits = c("UKB.Lym.BOLT", "UKB.Eosino.BOLT", "UKB.RBC.BOLT", "UKB.Plt.BOLT",  "UKB.Mono.BOLT", "UKB.Baso.BOLT", "UKB.Neutro.BOLT")

blood_traits = read.table("/data/deyk/kushal/LDSC/TASKFILES/sumstats.txt")[,1]
blood_traits = as.character(sapply(blood_traits, function(x) return(strsplit(x, ".sumstats")[[1]][1])))



get_sd_annot = function(cell_path, annot_index = 1, base_path, flag=0){
  if(flag == 0){
    sd_annot = rep(0, length(annot_index))
    for(i in 1:length(annot_index)){
      if(file.exists(paste0(cell_path, "/", "sd_annot_", annot_index[i], ".rda"))){
        sd_annot[i] = as.numeric(get(load(paste0(cell_path, "/", "sd_annot_", annot_index[i], ".rda"))))
      }else{
        flag = 1
        break
      }
    }
  }

  if(flag == 1){
    num = rep(0, length(annot_index))
    den = rep(0, length(annot_index))
    ll <- list.files(cell_path, pattern = ".annot.gz")
    ordering = c(10:19, 1, 20:22, 2, 3:9)
    for(m in 1:length(ll)){
      dat <- data.frame(fread(paste0("zcat ", cell_path, "/", ll[m])))
      base <- data.frame(fread(paste0("zcat ", base_path, "/", "baselineLD.", ordering[m], ".annot.gz")))
      pooled_dat <- cbind(dat[,-(1:4)], base[,-(1:4)])
      num = num  + (nrow(pooled_dat)-1) * apply(pooled_dat[,annot_index], 2, var)
      den = den + (nrow(pooled_dat)-1)
      rm(pooled_dat)
    }
    sd_annot = sqrt(num/den)
    for(i in 1:length(annot_index)){
      temp = sd_annot[i]
      save(temp, file = paste0(cell_path, "/", "sd_annot_", annot_index[i], ".rda"))
    }
  }
  return(sd_annot)
}

run_many_tau_analysis = function(annot_cell,
                                 results_cell,
                                 base_path,
                                 annotation,
                                 traits,
                                 index_in_results=NULL,  ### else specify 1:5 for first 5 annotations
                                 base_index = NULL,  ### number of annotations in the baseline
                                 flag = 1){
  base <- data.frame(fread(paste0("zcat ", base_path, "/", "baselineLD.", 22, ".annot.gz")))
  if(is.null(base_index)){
    base_index = ncol(base) - 4
  }
  cell_path = paste0(annot_cell, "/", annotation)
  res = paste0(results_cell, "/", annotation, "/", traits[1], ".sumstats.gz.results")
  tab2 = read.table(res,header=T)
  cat("Number of annotations together with baseline : ", nrow(tab2) , "\n")
  if(is.null(index_in_results)){index_in_results = 1:(nrow(tab2) - base_index)}
  tau_star_table = matrix(0, length(index_in_results), 3)
  annot_names = as.character(tab2$Category[index_in_results])
  sd_annot = get_sd_annot(cell_path, annot_index=index_in_results, base_path, flag = flag)
  for(id in 1:length(index_in_results)){
    sd_annot1=sd_annot[id]
    Mref = 5961159
    df = c()
    for(trait_id in 1:length(traits)){
      result.file=paste0(results_cell, "/", annotation, "/", traits[trait_id], ".sumstats.gz.part_delete")
      new_table=read.table(result.file,header=F)
      sc=c()
      logfile = paste(results_cell, "/", annotation, "/", traits[trait_id],".sumstats.gz.log", sep="")
      log = read.table(logfile,h=F,fill=T)
      h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
      coef1=sd_annot1*Mref/h2g
      for(i in 1:dim(new_table)[1]){
        tau1=as.numeric(new_table[i,index_in_results[id]])
        taus1=tau1*coef1
        #taus1 = tau1
        sc=c(sc,taus1)
        #cat("Block ", i, "\n")
      }
      mean_sc=mean(sc)
      se_sc=sqrt(199**2/200*var(sc))
      df = rbind(df, c(mean_sc,se_sc))
    }
    test_tauj=meta.summaries(df[,1],df[,2],method="random")
    tau=test_tauj$summary
    tau_se=test_tauj$se.summary
    z=tau/tau_se
    cat("Printing results for annotation:", annot_names[id], "\n")
    cat(tau, " ", tau_se, " ", 2*pnorm(-abs(z)), "\n")
    tau_star_table[id, ] = c(tau, tau_se, 2*pnorm(-abs(z)))
  }
  rownames(tau_star_table) = annot_names
  return(tau_star_table)
}


##########################################     brain      ###############################################


out3 = run_many_tau_analysis(annot_cell, results_cell, base_path,
                             annotation = annot_names,
                             traits = blood_traits,
                             flag = 0, base_index = NULL)
out4 = run_many_tau_analysis(annot_cell, results_cell, base_path,
                             annotation = annot_names,
                             traits = autoimmune_traits,
                             flag = 0, base_index = NULL)
out5 = run_many_tau_analysis(annot_cell, results_cell, base_path,
                             annotation = annot_names,
                             traits = blood_bio_traits,
                             flag = 0, base_index = NULL)


ll <- list()
ll[["Blood"]] = out3
ll[["BloodBio"]] = out5
ll[["Autoimmune"]] = out4

save(ll, file = paste0("/data/deyk/kushal/Miao_DO/data/joint_meta_taustar_", annot_names, ".rda"))








ll$All[order(ll$All[,3], decreasing = FALSE), ]



out1[which(out1[,3] < 0.1),]

save(ll, file = paste0("/n/groups/price/kushal/LDSC-Average/output/", annot_names, ".rda"))

