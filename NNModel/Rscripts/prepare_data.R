library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggprism)
library(stringr)
library(dplyr)
theme_set(theme_pubclean())
library(EnvStats) # for adding sample size
# library(ggpomological) #theme
library(survival)
library(ggfortify)
library(ranger)
library("anytime")
library(chron)
library(ggsurvfit)
library(gtsummary)
library(coxphf)
library(survminer)

setwd("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts")
source("survivalFunctionR.R")

get_data = function() {

  data <- read.csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
  data$status <- data$survivalStatus
  data$BCFS <-  data$survivalTime
  
  data$ichorPrediction[data$ichorTF>=0.03]<-"Positive"
  data$ichorPrediction[data$ichorTF<0.03]<-"Negative" 
  data$ctDNADetected[data$ctDNADetected==0]<-"Negative"
  data$ctDNADetected[data$ctDNADetected==1]<-"Positive" 
  data$ctDNA[data$ctDNADetected=="Fail"]<-"Fail" 
  
  # predictions <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ctDNADetection_corrected2.csv", sep="\t")
  
  predictions <- read.csv("~/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ichorTF.csv", sep="\t")
  
  data_pred <- merge(data,predictions, by="PatientId" )
  
  data_pred$Predictedlabel = data_pred$Predictedlabel/100
  
  data_pred$ALFAssay[data_pred$Predictedlabel>=0.03]<-"Positive"
  data_pred$ALFAssay[data_pred$Predictedlabel<0.03]<-"Negative" 
  
  
  # pearl 
  fragle_pred <-  read.csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")
  
  fragle_pred$PatientId <- str_replace(fragle_pred$Sample_ID, ".markDup.GCtagged", "")
  
  data_pearl<-data_pred %>% filter(study=="Pearl")
  data_pearl$VAF_pers <- data_pearl$VAF/100
  
  data_pearl <- data_pearl %>% filter(!is.na(ichorTF))
  data_pearl <- data_pearl %>% filter(!is.na(BCFS))
  data_pearl <- data_pearl %>% filter(!is.na(ALFAssay))
  
  data_pearl <- merge(fragle_pred, data_pearl, by="PatientId")
  
  data_pearl$ctDNA_Burden_surv[data_pearl$ctDNA_Burden>=0.01]<-"HighctDNA" 
  data_pearl$ctDNA_Burden_surv[data_pearl$ctDNA_Burden<0.01]<-"LowctDNA"
  data_pearl$timepoint[data_pearl$timepoint=="C1"] <- "D15" 
  data_pearl$timepoint[data_pearl$timepoint=="Surgery"] <- "PD" 
  
  data_pearl$patient <- substr(data_pearl$PatientId, start = 1, stop = 12)
  
  data_pearl_bl <- data_pearl %>% filter(timepoint=="BL")
  data_pearl_d15 <- data_pearl %>% filter(timepoint=="D15")
  data_pearl_pd <- data_pearl %>% filter(timepoint=="PD")
  
  cor(data_pearl$ichorTF, data_pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
  cor(data_pearl$ctDNA_Burden, data_pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
  cor(data_pearl$VAF_pers, data_pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
  
  
  ## neorhea
  data_neorhea <- data_pred %>% filter(study=="Neorhea")
  fragle_pred_neo <-  read.csv("/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")
  
  fragle_pred_neo$PatientId <- str_replace(fragle_pred_neo$Sample_ID, ".markDup.GCtagged", "")
  
  data_neorhea <- merge(fragle_pred_neo, data_neorhea, by="PatientId")
  
  
  data_neorhea <- data_neorhea %>% filter(!is.na(ctDNADetected))
  data_neorhea <- data_neorhea %>% filter(!is.na(BCFS))
  data_neorhea <- data_neorhea %>% filter(!is.na(ALFAssay))
  
  
  data_neorhea$ctDNA_Burden_surv[data_neorhea$ctDNA_Burden>=0.01]<-"HighctDNA" 
  data_neorhea$ctDNA_Burden_surv[data_neorhea$ctDNA_Burden<0.01]<-"LowctDNA"
  data_neorhea$VAF_perc <- data_neorhea$VAF*100
  data_neorhea$patient <- substr(data_neorhea$PatientId, start = 1, stop = 12)
  
  data_neorhea_bl <- data_neorhea %>% filter(timepoint=="BL")
  data_neorhea_c1d28 <- data_neorhea %>% filter(timepoint=="C1")
  data_neorhea_surgery <- data_neorhea %>% filter(timepoint=="Surgery")
  
  
  # snr
  pred_snr <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_train_SNR_healthy_ichorTF.csv", sep="\t")
  data_snr <- data %>% filter(study=="Synergy")
  data_pred_snr <- merge(data_snr,pred_snr, by="PatientId")
  
  data_pred_snr$Predictedlabel <- round(data_pred_snr$Predictedlabel/100,6)
  data_pred_snr$Predictedlabel[data_pred_snr$Predictedlabel<0] <- 0
  data_pred_snr$ALFAssay[data_pred_snr$Predictedlabel>=0.03]<-"Positive"
  data_pred_snr$ALFAssay[data_pred_snr$Predictedlabel<0.03]<-"Negative" 
  data_pred_snr$patient <- substr(data_pred_snr$PatientId, start = 1, stop = 10)
  data_pred_snr$timepoint[data_pred_snr$timepoint=="BL"] <- "W1" 
  data_pred_snr$timepoint[data_pred_snr$timepoint=="C1"] <- "W3" 
  data_pred_snr$timepoint[data_pred_snr$timepoint=="Surgery"] <- "W13"
  
  fragle_pred_snr <-  read.csv("/Users/alexandra/PhD/SynergyStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")
  
  fragle_pred_snr$PatientId <- str_replace(fragle_pred_snr$Sample_ID, ".markDup.GCtagged", "")
  
  data_snr <- merge(fragle_pred_snr, data_pred_snr, by="PatientId")
  
  
  data_snr$ctDNA_Burden_surv[data_snr$ctDNA_Burden>=0.01]<-"HighctDNA" 
  data_snr$ctDNA_Burden_surv[data_snr$ctDNA_Burden<0.01]<-"LowctDNA"
  
  
  data_snr_bl <- data_snr %>% filter(timepoint=="BL")
  data_snr_c1d28 <- data_snr %>% filter(timepoint=="C1")
  data_snr_surgery <- data_snr %>% filter(timepoint=="Surgery")
  
  # neorhea = list()
  # neorhea[[1]] <- data_neorhea
  # neorhea[[2]] <- data_neorhea_bl
  # neorhea[[3]] <- data_neorhea_c1d28
  # neorhea[[4]] <- data_neorhea_surgery
  # names(neorhea) <- names(c("data_neorhea", "data_neorhea_bl", "data_neorhea_c1d28", "data_neorhea_surgery"))
  neorhea <- list (all_neorhea=data_neorhea, bl_neorhea =data_neorhea_bl, c1d28_neorhea=data_neorhea_c1d28,surhery_neorhea= data_neorhea_surgery)
  
  # pearl = list()
  # pearl[[1]] <- data_pearl
  # pearl[[2]] <- data_pearl_bl
  # pearl[[3]] <- data_pearl_d15
  # pearl[[4]] <- data_pearl_pd
  # names(pearl) <- names(c("data_pearl", "data_pearl_bl", "data_pearl_d15", "data_pearl_pd"))
  
  pearl <- list (all_pearl=data_pearl, bl_pearl =data_pearl_bl, d15_pearl=data_pearl_d15,pd_pearl= data_pearl_pd)
  
  # synergy = list()
  # synergy[[1]] <- data_snr
  # synergy[[2]] <- data_snr_bl
  # synergy[[3]] <- data_snr_c1d28
  # synergy[[4]] <- data_snr_surgery
  # names(synergy) <- names(c("data_snr", "data_snr_bl", "data_snr_c1d28", "data_snr_surgery"))
  synergy <- list (all_synergy=data_snr, w1_snr =data_snr_bl, w3_snr=data_snr_c1d28, w13_snr= data_snr_surgery)
  
  # data = list()
  # data[[1]] <- pearl
  # data[[2]] <- neorhea
  # data[[3]] <- synergy
  # 
  # names(data) <- names(c("pearl", "neorhea", "synergy"))
  
  data <- list (prl=pearl, nrh =neorhea, snr=synergy)
  
  return(data)
}





get_data_2 = function(try) {
  
  data <- read.csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
  data$status <- data$survivalStatus
  data$BCFS <-  data$survivalTime
  
  data$ichorPrediction[data$ichorTF>=0.03]<-"Positive"
  data$ichorPrediction[data$ichorTF<0.03]<-"Negative" 
  data$ctDNADetected[data$ctDNADetected==0]<-"Negative"
  data$ctDNADetected[data$ctDNADetected==1]<-"Positive" 
  data$ctDNA[data$ctDNADetected=="Fail"]<-"Fail" 
  
  # predictions <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ctDNADetection_corrected2.csv", sep="\t")
  
  predictions_validation <- read.csv(paste0("~/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF",try,".csv"), sep="\t")
  predictions_test <- read.csv(paste0("~/PhD/PyCharmProjects/ALFAssayNN/results/results_train_ichorTF",try,".csv"), sep="\t")

  # predictions_validation <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssay/NNModel/results/results_validation_NRH_PRL_healthy_survLabel.csv", sep="\t")
  # predictions_test <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssay/NNModel/results/results_train_SNR_healthy_survLabel.csv", sep="\t")
  # 
  # 
  
  prediction <- rbind(predictions_validation,predictions_test)
  
  data_pred <- merge(data,prediction, by="PatientId" )
  data_pred$study <- if_else(str_detect(data_pred$PatientId,"NIPT-PIJB"),"Healthy_synergy", data_pred$study)
  
  data_pred$Predictedlabel = data_pred$Predictedlabel/100
  data_pred$Predictedlabel[data_pred$Predictedlabel<0]<-0

  data_pred$ALFAssay[data_pred$Predictedlabel>=0.03]<-"Positive"
  data_pred$ALFAssay[data_pred$Predictedlabel<0.03]<-"Negative"
  
  
  # data_pred$Predictedlabel = data_pred$Predictedlabel
  # data_pred$Predictedlabel[data_pred$Predictedlabel<0]<-0
  # 
  # data_pred$ALFAssay[data_pred$Predictedlabel>=0.5]<-1
  # data_pred$ALFAssay[data_pred$Predictedlabel<0.5]<-0 
  
  
  # pearl 
  fragle_pred <-  read.csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")
  
  fragle_pred$PatientId <- str_replace(fragle_pred$Sample_ID, ".markDup.GCtagged", "")
  
  data_pearl<- data_pred %>% filter(study=="Pearl")
  data_pearl$VAF_perc <- data_pearl$VAF/100
  data_pearl$PFS <- data_pearl$BCFS
  
  # data_pearl <- data_pearl %>% filter(!is.na(ichorTF))
  # data_pearl <- data_pearl %>% filter(!is.na(BCFS))
  # data_pearl <- data_pearl %>% filter(!is.na(ALFAssay))
  
  data_pearl <- merge(fragle_pred, data_pearl, by="PatientId")
  
  data_pearl$FragleDetection[data_pearl$ctDNA_Burden>=0.01]<-"HighctDNA" 
  data_pearl$FragleDetection[data_pearl$ctDNA_Burden<0.01]<-"LowctDNA"
  data_pearl$timepoint[data_pearl$timepoint=="C1"] <- "D15" 
  data_pearl$timepoint[data_pearl$timepoint=="Surgery"] <- "PD" 
  
  data_pearl$patient <- substr(data_pearl$PatientId, start = 1, stop = 12)
  
  data_pearl_bl <- data_pearl %>% filter(timepoint=="BL")
  data_pearl_d15 <- data_pearl %>% filter(timepoint=="D15")
  data_pearl_pd <- data_pearl %>% filter(timepoint=="PD")
  
  
  ## neorhea
  data_neorhea <- data_pred %>% filter(study=="Neorhea")
  fragle_pred_neo <-  read.csv("/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")
  
  fragle_pred_neo$PatientId <- str_replace(fragle_pred_neo$Sample_ID, ".markDup.GCtagged", "")
  
  data_neorhea <- merge(fragle_pred_neo, data_neorhea, by="PatientId")
  data_neorhea$PFS <- data_neorhea$BCFS
  
  # data_neorhea <- data_neorhea %>% filter(!is.na(ctDNADetected))
  # data_neorhea <- data_neorhea %>% filter(!is.na(BCFS))
  # data_neorhea <- data_neorhea %>% filter(!is.na(ALFAssay))
  # 
  
  data_neorhea$FragleDetection[data_neorhea$ctDNA_Burden>=0.01]<-"HighctDNA" 
  data_neorhea$FragleDetection[data_neorhea$ctDNA_Burden<0.01]<-"LowctDNA"
  data_neorhea$VAF_perc <- data_neorhea$VAF*100
  data_neorhea$patient <- substr(data_neorhea$PatientId, start = 1, stop = 12)
  
  data_neorhea_bl <- data_neorhea %>% filter(timepoint=="BL")
  data_neorhea_c1d28 <- data_neorhea %>% filter(timepoint=="C1")
  data_neorhea_surgery <- data_neorhea %>% filter(timepoint=="Surgery")
  
  
  # snr
  # pred_snr <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_train_SNR_healthy_ichorTF.csv", sep="\t")
  data_pred_snr <- data_pred %>% filter(study=="Synergy")
  # data_pred_snr <- merge(data_snr,pred_snr, by="PatientId")
  
  # data_pred_snr$Predictedlabel <- round(data_pred_snr$Predictedlabel/100,6)
  # data_pred_snr$Predictedlabel[data_pred_snr$Predictedlabel<0] <- 0
  # data_pred_snr$ALFAssay[data_pred_snr$Predictedlabel>=0.03]<-"Positive"
  # data_pred_snr$ALFAssay[data_pred_snr$Predictedlabel<0.03]<-"Negative" 
  data_pred_snr$patient <- substr(data_pred_snr$PatientId, start = 1, stop = 10)
  data_pred_snr$timepoint[data_pred_snr$timepoint=="BL"] <- "W1" 
  data_pred_snr$timepoint[data_pred_snr$timepoint=="C1"] <- "W3" 
  data_pred_snr$timepoint[data_pred_snr$timepoint=="Surgery"] <- "W13"
  data_pred_snr$PFS <- data_pred_snr$BCFS
  fragle_pred_snr <-  read.csv("/Users/alexandra/PhD/SynergyStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")
  
  fragle_pred_snr$PatientId <- str_replace(fragle_pred_snr$Sample_ID, ".markDup.GCtagged", "")
  
  data_snr <- merge(fragle_pred_snr, data_pred_snr, by="PatientId")
  
  
  data_snr$FragleDetection[data_snr$ctDNA_Burden>=0.01]<-"HighctDNA" 
  data_snr$FragleDetection[data_snr$ctDNA_Burden<0.01]<-"LowctDNA"
  data_snr$VAF_perc<- 0
  
  
  data_snr_bl <- data_snr %>% filter(timepoint=="W1")
  data_snr_c1d28 <- data_snr %>% filter(timepoint=="W3")
  data_snr_surgery <- data_snr %>% filter(timepoint=="W13")
  
  # neorhea = list()
  # neorhea[[1]] <- data_neorhea
  # neorhea[[2]] <- data_neorhea_bl
  # neorhea[[3]] <- data_neorhea_c1d28
  # neorhea[[4]] <- data_neorhea_surgery
  # names(neorhea) <- names(c("data_neorhea", "data_neorhea_bl", "data_neorhea_c1d28", "data_neorhea_surgery"))
  neorhea <- list (all_neorhea=data_neorhea, bl_neorhea =data_neorhea_bl, c1d28_neorhea=data_neorhea_c1d28,surhery_neorhea= data_neorhea_surgery)
  
  # pearl = list()
  # pearl[[1]] <- data_pearl
  # pearl[[2]] <- data_pearl_bl
  # pearl[[3]] <- data_pearl_d15
  # pearl[[4]] <- data_pearl_pd
  # names(pearl) <- names(c("data_pearl", "data_pearl_bl", "data_pearl_d15", "data_pearl_pd"))
  
  pearl <- list (all_pearl=data_pearl, bl_pearl =data_pearl_bl, d15_pearl=data_pearl_d15,pd_pearl= data_pearl_pd)
  
  # synergy = list()
  # synergy[[1]] <- data_snr
  # synergy[[2]] <- data_snr_bl
  # synergy[[3]] <- data_snr_c1d28
  # synergy[[4]] <- data_snr_surgery
  # names(synergy) <- names(c("data_snr", "data_snr_bl", "data_snr_c1d28", "data_snr_surgery"))
  synergy <- list (all_synergy=data_snr, w1_snr =data_snr_bl, w3_snr=data_snr_c1d28, w13_snr= data_snr_surgery)
  
 #controls
  data_healthy <- data_pred %>% filter(study %in% c("Healthy_synergy","Healthy"))
  data_healthy$patient <- substr(data_healthy$PatientId, start = 1, stop = 16)
  
  ## add Fragle pred and ichorCNA pred
  
  fragle_healthy <- read.csv("~/PhD/healthy_sWGS/fragle_healthy.csv", sep=",")
  
  fragle_healthy$PatientId <- sub("\\.markDup\\.GCtagged$", "", fragle_healthy$Sample_ID)
  ichor_healthy = read.csv("/Users/alexandra/PhD/healthy_sWGS/ichor_CNA_results.csv", sep="\t")
  
  data_healthy <- merge(fragle_healthy,data_healthy, by="PatientId") 
  data_healthy <- merge(ichor_healthy,data_healthy, by="PatientId") 
  
  data_healthy$ichorTF.y <- NULL
  names(data_healthy)[names(data_healthy) == "ichorTF.x"] <- "ichorTF"
  
  data_healthy$FragleDetection[data_healthy$ctDNA_Burden>=0.01]<-"HighctDNA" 
  data_healthy$FragleDetection[data_healthy$ctDNA_Burden<0.01]<-"LowctDNA"
  data_healthy$VAF_perc <-0
  
  data_healthy$PFS <- data_healthy$BCFS
  
  data <- list (prl=pearl, nrh =neorhea, snr=synergy, ctrl=data_healthy)
  
  return(data)
}

get_data_for_table = function(try) {
  cols <- c(
    "PatientId", "Sample_ID", "ctDNA_Burden", "Label",
    "ichorTF", "timepoint", "VAF", "ctDNADetected",
    "survivalStatus", "survivalTime", "VAFg0p001", "study",
    "MedianStatus", "batch", "Process", "status",
    "BCFS", "ichorPrediction", "ctDNA", "Predictedlabel",
    "Truelabel", "ALFAssay", "FragleDetection", "VAF_perc",
    "patient"
  )
  healthy <- select(data$ctrl,cols)
  pearl <- select(data$prl$all_pearl,cols)
  neorhea <- select(data$nrh$all_neorhea,cols)
  synergy <- select(data$snr$all_synergy,cols)
  
  data <- get_data_2(try)
  full_meta <- rbind(healthy,pearl,synergy,neorhea)
  full_meta <- full_meta %>%
    dplyr::rename(
      OncoFollowDNA = ctDNADetected,
      ichorCNA      = ichorPrediction,
      Fragle        = FragleDetection,
      FragleTF = ctDNA_Burden,
      OncoDNAVAF = VAF_perc
    )
  
  full_meta$ichorTF <- as.numeric(full_meta$ichorTF)
  
  write.table(full_meta, "~/PhD/FragmentationPatterns/Data/MetaData/MetaDataAndPredictions.csv", quote = FALSE, row.names = FALSE,sep="\t")
}
