library(data.table)
library(rmeta)
annot_cell = "/data/deyk/kushal/Miao_DO/data/ANNOTATIONS_hg38/Miao_DO"
results_cell = "/data/deyk/kushal/Miao_DO/data/LDSC_RESULTS_hg38/Miao_DO/baseline_Epi_hg38/"
annot_names = list.files(results_cell)
annot_idx = 1

source("/data/deyk/kushal/Placenta/codes/rmeta/R/meta.R")

autoimmune_traits = c("Liu2015-CD", "Liu2015-IBD", "Liu2015-UC",
                      "PASS_Celiac", "PASS_Lupus",  "PASS_Primary_biliary_cirrhosis", "PASS_Rheumatoid_Arthritis",
                      "PASS_Type_1_Diabetes", "UKB.AID_Combined.SAIGE", "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED",
                      "UKB.Asthma.SAIGE", "PASS_Multiple_sclerosis", "PASS_Alzheimers_Jansen2019")

blood_bio_traits = c("UKB.Lym.BOLT", "UKB.Eosino.BOLT", "UKB.RBC.BOLT", "UKB.Plt.BOLT",  "UKB.Mono.BOLT", "UKB.Baso.BOLT", "UKB.Neutro.BOLT")

blood_traits = read.table("/data/deyk/kushal/LDSC/TASKFILES/sumstats.txt")[,1]
blood_traits = as.character(sapply(blood_traits, function(x) return(strsplit(x, ".sumstats")[[1]][1])))



get_sd_annot = function(cell_path, annot_index = 1, flag=0){
  if(flag == 0){
    if(file.exists(paste0(cell_path, "/", "sd_annot_", annot_index, ".rda"))){
        sd_annot = get(load(paste0(cell_path, "/", "sd_annot_", annot_index, ".rda")))
        return(sd_annot)
    }else{
    	flag = 1
    }}

    if(flag == 1){
    	num = 0
        den = 0
        ll <- list.files(cell_path, pattern = ".annot.gz")
        for(m in 1:length(ll)){
            dat <- data.frame(fread(paste0("zcat ", cell_path, "/", ll[m])))
            num = num  + (nrow(dat)-1) * var(dat[,4+annot_index])
            den = den + (nrow(dat)-1)
            rm(dat)
       }
    }

  estd_sd_annot = sqrt(num/den)
  save(estd_sd_annot, file = paste0(cell_path, "/", "sd_annot_", annot_index, ".rda"))
  return(estd_sd_annot)
}

run_single_tau_analysis = function(annot_cell,
                                   results_cell,
                                   annotations,
                                   traits,
                                   index_in_results=1,
                                   base_index = NULL,
                                   flag = 1){
    if(is.null(base_index)){base_index = index_in_results}
    tau_star_table = matrix(0, length(annotations), 3)
    for(annot_id in 1:length(annotations)){
        cell_path = paste0(annot_cell, "/", annotations[annot_id], "/", "ABC_Road_GI_BLD")
        sd_annot1=get_sd_annot(cell_path, annot_index=index_in_results, flag = flag)
        Mref = 5961159
        df = c()
        for(trait_id in 1:length(traits)){
            result.file=paste0(results_cell, "/", annotations[annot_id], "/", "ABC_Road_GI_BLD", "/", traits[trait_id], ".sumstats.gz.part_delete")
            new_table=read.table(result.file,header=F)
            sc=c()
            logfile = paste(results_cell, "/", annotations[annot_id], "/", "ABC_Road_GI_BLD", "/", traits[trait_id],".sumstats.gz.log", sep="")
            log = read.table(logfile,h=F,fill=T)
            h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
            coef1=sd_annot1*Mref/h2g
            for(i in 1:dim(new_table)[1]){
                  tau1=as.numeric(new_table[i,base_index])
                  taus1=tau1*coef1
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
        cat("Printing results for annotation:", annotations[annot_id], "\n")
        cat(tau, " ", tau_se, " ", 2*pnorm(-abs(z)), "\n")
        tau_star_table[annot_id, ] = c(tau, tau_se, 2*pnorm(-abs(z)))
    }
    rownames(tau_star_table) = annotations
    return(tau_star_table)
}






out3 = run_single_tau_analysis(annot_cell, results_cell, annotations = annot_names, traits = blood_traits,
                                index_in_results = annot_idx, flag = 0)
out4 = run_single_tau_analysis(annot_cell, results_cell, annotations = annot_names, traits = autoimmune_traits,
                               index_in_results = annot_idx, flag = 0)
out5 = run_single_tau_analysis(annot_cell, results_cell, annotations = annot_names, traits = blood_bio_traits,
                               index_in_results = annot_idx, flag = 0)



ll <- list()
ll[["Blood"]] = out3
ll[["BloodBio"]] = out5
ll[["Autoimmune"]] = out4

save(ll, file = "/data/deyk/kushal/Miao_DO/data/marginal_meta_taustar_Miao_DO.rda")


out3 = run_single_tau_analysis(annot_cell, results_cell, annotations = annot_names[9],
                               traits = temp_trait,
                               index_in_results = annot_idx, flag = 1)

