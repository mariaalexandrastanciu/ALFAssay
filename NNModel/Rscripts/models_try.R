library(dplyr)
library(ggplot2)
library(magrittr)
library(data.table)
library(stringr)
rm(list = ls())
library(gmodels)
# set directory w
# set directory where data are reposited
ichorDir <- "~/PhD/healthy_sWGS/ichorCNA/results"
# get TF, MAD and ploidy in a dataframe with samples at rows and parameters at columns
paramFiles <- list.files(ichorDir, pattern = "param", recursive = T, full.names = T)

ichorParams <- sapply(paramFiles, function(file){
  lol <- readr::read_delim(file)
  return(list(TF = lol[4,2] %>% as.numeric,
              MAD = lol[16,2] %>% as.numeric,
              ploidy = lol[5,2] %>% as.numeric))
})
# ichorParams
ichorParams %<>% t
rownames(ichorParams) <- paramFiles %>% basename %>% stringr::str_split("\\.") %>% sapply(`[`, 1)
ichorParams %<>% data.frame

tf <- do.call(rbind.data.frame, ichorParams$TF)
colnames(tf) <- "ichorCNATF"
# 
newdataset <- data.frame(PatientId = rownames(ichorParams),ichorTF=tf$ichorCNATF)
write.table(newdataset, "/Users/alexandra/PhD/healthy_sWGS/ichor_CNA_results.csv",
            sep = '\t',col.names = TRUE, quote = F, row.names = FALSE)

data <- read.csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
ichor_healthy = read.csv("/Users/alexandra/PhD/healthy_sWGS/ichor_CNA_results.csv", sep="\t")

healthy <- data %>% dplyr::filter(Label==0)
data_1 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ichorTF_t2.csv", sep="\t")
data_1$Predictedlabel <- data_1$Predictedlabel/100
data_1 <- merge(data_1, ichor_healthy, on = "PatientId")

data_1 %>% filter(PatientId %in% healthy$PatientId &Predictedlabel >=0.01)


data_2 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ichorTF.csv", sep="\t")
data_2$Predictedlabel <- data_2$Predictedlabel/100
data_2 <- merge(data_2, ichor_healthy, on = "PatientId")
data_2 %>% filter(PatientId %in% healthy$PatientId &Predictedlabel >=0.1)


data_3 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ichorTF_t3.csv", sep="\t")
data_3$Predictedlabel <- data_3$Predictedlabel/100
data_3 <- merge(data_3, ichor_healthy, on = "PatientId")
data_3 %>% filter(PatientId %in% healthy$PatientId &Predictedlabel >=0.01)


data_4 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ichorTF_t4.csv", sep="\t")
data_4$Predictedlabel <- data_4$Predictedlabel/1000
data_4 <- merge(data_4, ichor_healthy, on = "PatientId")
data_4 %>% filter(PatientId %in% healthy$PatientId &Predictedlabel >=0.1)


data_5 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ichorTF_t5.csv", sep="\t")
data_5$Predictedlabel <- data_5$Predictedlabel/1000
data_5 <- merge(data_5, ichor_healthy, on = "PatientId")
data_5 %>% filter(PatientId %in% healthy$PatientId &Predictedlabel >=0.1)

data_6 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ichorTF_t6.csv", sep="\t")
data_6$Predictedlabel <- data_6$Predictedlabel/100
data_6 <- merge(data_6, ichor_healthy, on = "PatientId")
data_6 %>% filter(PatientId %in% healthy$PatientId &Predictedlabel >=0.1)


data_7 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t7.csv", sep="\t")
data_7$Predictedlabel <- data_7$Predictedlabel/100
data_7 <- merge(data_7, ichor_healthy, on = "PatientId")
data_7 %>% filter(PatientId %in% healthy$PatientId &Predictedlabel >=0.1)


data_8 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t8.csv", sep="\t")
data_8$Predictedlabel <- data_8$Predictedlabel/100
data_8$Truelabel <- data_8$Truelabel/100
data_8$ichorDetection[data_8$Truelabel>=0.03]<-"Positive"
data_8$ichorDetection[data_8$Truelabel<0.03]<-"Negative"
data_8$ALFAssay[data_8$Predictedlabel>=0.03]<-"Positive"
data_8$ALFAssay[data_8$Predictedlabel<0.03]<-"Negative"

data_8_healthy <- merge(data_8, ichor_healthy, on = "PatientId")
data_8 %>% filter(PatientId %in% data_8_healthy$PatientId &Predictedlabel >=0.1)


data_9 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t9.csv", sep="\t")
data_9$Predictedlabel <- data_9$Predictedlabel/100
data_9$Truelabel <- data_9$Truelabel/100
data_9$ichorDetection[data_9$Truelabel>=0.03]<-"Positive"
data_9$ichorDetection[data_9$Truelabel<0.03]<-"Negative"
data_9$ALFAssay[data_9$Predictedlabel>=0.03]<-"Positive"
data_9$ALFAssay[data_9$Predictedlabel<0.03]<-"Negative"
with(data_9, CrossTable(ichorDetection, ALFAssay))

data_9_healthy <- merge(data_9, ichor_healthy, on = "PatientId")
data_9_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >=0.03)




data_10 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t10.csv", sep="\t")
data_10$Predictedlabel <- data_10$Predictedlabel/100
data_10$Truelabel <- data_10$Truelabel/100
data_10$ichorDetection[data_10$Truelabel>=0.03]<-"Positive"
data_10$ichorDetection[data_10$Truelabel<0.03]<-"Negative"
data_10$ALFAssay[data_10$Predictedlabel>=0.03]<-"Positive"
data_10$ALFAssay[data_10$Predictedlabel<0.03]<-"Negative"
with(data_10, CrossTable(ichorDetection, ALFAssay))

data_10_healthy <- merge(data_10, ichor_healthy, on = "PatientId")
data_10_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >=0.01)




data_11 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t11.csv", sep="\t")
data_11$Predictedlabel <- data_11$Predictedlabel/100
data_11$Truelabel <- data_11$Truelabel/100
data_11$ichorDetection[data_11$Truelabel>=0.03]<-"Positive"
data_11$ichorDetection[data_11$Truelabel<0.03]<-"Negative"
data_11$ALFAssay[data_11$Predictedlabel>=0.03]<-"Positive"
data_11$ALFAssay[data_11$Predictedlabel<0.03]<-"Negative"
with(data_11, CrossTable(ichorDetection, ALFAssay))

data_11_healthy <- merge(data_11, ichor_healthy, on = "PatientId")
data_11_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)

data_9 %>% filter(PatientId %in% c(
'Genome-IJB-HP-81_S186', 'Genome-IJB-HP-77_S182', 'Genome-IJB-HP-70_S175', 'Genome-IJB-HP-71_S176', 'Genome-IJB-HP-56_S161', 'Genome-IJB-HP-57_S162', 'Genome-IJB-HP-58_S163',
'Genome-IJB-HP-53_S158', 'Genome-IJB-HP-51_S156'))



ratio_short_healthy <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssay/test/testdata/frag_ratio_short_healthy.csv", sep="\t")

list_healthy <- read.csv("/Users/alexandra/PhD/healthy_sWGS/mosdepth/list_healthy.txt", header=FALSE)
coverage <- read.csv("/Users/alexandra/PhD/healthy_sWGS/mosdepth/total_depth_healthy.txt", header=FALSE)

list_healthy$PatientId <- str_replace(str_replace(list_healthy$V1,"/globalscratch/ulb/bctr/astanciu/healthy_sWGS/mosdepth/",""),".mosdepth.summary.txt","")
list_healthy$coverage <- coverage$V1

frag_ratio_coverage <- merge(ratio_short_healthy,list_healthy, by="PatientId")

cor(frag_ratio_coverage$ratio, frag_ratio_coverage$coverage)



 ggscatter(frag_ratio_coverage, x = "ratio", y = "coverage",
                                   add = "reg.line",  # Add regressin line
                                   add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                                   conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "Short fragment ratio", y = "Coverage")

 
 
 
 data_12 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t12.csv", sep="\t")
 data_12$Predictedlabel <- data_12$Predictedlabel/100
 data_12$Truelabel <- data_12$Truelabel/100
 data_12$ichorDetection[data_12$Truelabel>=0.03]<-"Positive"
 data_12$ichorDetection[data_12$Truelabel<0.03]<-"Negative"
 data_12$ALFAssay[data_12$Predictedlabel>=0.03]<-"Positive"
 data_12$ALFAssay[data_12$Predictedlabel<0.03]<-"Negative"
 with(data_12, CrossTable(ichorDetection, ALFAssay))
 
 data_12_healthy <- merge(data_12, ichor_healthy, on = "PatientId")
 data_12_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 
 
 
 
 data_13 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t13.csv", sep="\t")
 data_13$Predictedlabel <- data_13$Predictedlabel/100
 data_13$Truelabel <- data_13$Truelabel/100
 data_13$ichorDetection[data_13$Truelabel>=0.03]<-"Positive"
 data_13$ichorDetection[data_13$Truelabel<0.03]<-"Negative"
 data_13$ALFAssay[data_13$Predictedlabel>=0.03]<-"Positive"
 data_13$ALFAssay[data_13$Predictedlabel<0.03]<-"Negative"
 with(data_13, CrossTable(ichorDetection, ALFAssay))
 
 data_13_healthy <- merge(data_13, ichor_healthy, on = "PatientId")
 data_13_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 
 data_14 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t14.csv", sep="\t")
 data_14$Predictedlabel <- data_14$Predictedlabel/100
 data_14$Truelabel <- data_14$Truelabel/100
 data_14$ichorDetection[data_14$Truelabel>=0.03]<-"Positive"
 data_14$ichorDetection[data_14$Truelabel<0.03]<-"Negative"
 data_14$ALFAssay[data_14$Predictedlabel>=0.03]<-"Positive"
 data_14$ALFAssay[data_14$Predictedlabel<0.03]<-"Negative"
 with(data_14, CrossTable(ichorDetection, ALFAssay))
 
 data_14_healthy <- merge(data_14, ichor_healthy, on = "PatientId")
 data_14_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)

 
 
 
 data_15 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t14.csv", sep="\t")
 data_15$Predictedlabel <- data_15$Predictedlabel/100
 data_15$Truelabel <- data_15$Truelabel/100
 data_15$ichorDetection[data_15$Truelabel>=0.03]<-"Positive"
 data_15$ichorDetection[data_15$Truelabel<0.03]<-"Negative"
 data_15$ALFAssay[data_15$Predictedlabel>=0.03]<-"Positive"
 data_15$ALFAssay[data_15$Predictedlabel<0.03]<-"Negative"
 with(data_15, CrossTable(ichorDetection, ALFAssay))
 
 data_15_healthy <- merge(data_15, ichor_healthy, on = "PatientId")
 data_15_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 
 data_16 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t16.csv", sep="\t")
 data_16$Predictedlabel <- data_16$Predictedlabel/100
 data_16$Truelabel <- data_16$Truelabel/100
 data_16$ichorDetection[data_16$Truelabel>=0.03]<-"Positive"
 data_16$ichorDetection[data_16$Truelabel<0.03]<-"Negative"
 data_16$ALFAssay[data_16$Predictedlabel>=0.03]<-"Positive"
 data_16$ALFAssay[data_16$Predictedlabel<0.03]<-"Negative"
 with(data_16, CrossTable(ichorDetection, ALFAssay))
 
 data_16_healthy <- merge(data_16, ichor_healthy, on = "PatientId")
 data_16_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 
 data_17 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t17.csv", sep="\t")
 data_17$Predictedlabel <- data_17$Predictedlabel/100
 data_17$Truelabel <- data_17$Truelabel/100
 data_17$ichorDetection[data_17$Truelabel>=0.03]<-"Positive"
 data_17$ichorDetection[data_17$Truelabel<0.03]<-"Negative"
 data_17$ALFAssay[data_17$Predictedlabel>=0.03]<-"Positive"
 data_17$ALFAssay[data_17$Predictedlabel<0.03]<-"Negative"
 with(data_17, CrossTable(ichorDetection, ALFAssay))
 
 data_17_healthy <- merge(data_17, ichor_healthy, on = "PatientId")
 data_17_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 
 
 data_18 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t18.csv", sep="\t")
 data_18$Predictedlabel <- data_18$Predictedlabel/100
 data_18$Truelabel <- data_18$Truelabel/100
 data_18$ichorDetection[data_18$Truelabel>=0.03]<-"Positive"
 data_18$ichorDetection[data_18$Truelabel<0.03]<-"Negative"
 data_18$ALFAssay[data_18$Predictedlabel>=0.03]<-"Positive"
 data_18$ALFAssay[data_18$Predictedlabel<0.03]<-"Negative"
 with(data_18, CrossTable(ichorDetection, ALFAssay))
 
 data_18_healthy <- merge(data_18, ichor_healthy, on = "PatientId")
 data_18_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 data_19 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t19.csv", sep="\t")
 data_19$Predictedlabel <- data_19$Predictedlabel/100
 data_19$Truelabel <- data_19$Truelabel/100
 data_19$ichorDetection[data_19$Truelabel>=0.03]<-"Positive"
 data_19$ichorDetection[data_19$Truelabel<0.03]<-"Negative"
 data_19$ALFAssay[data_19$Predictedlabel>=0.03]<-"Positive"
 data_19$ALFAssay[data_19$Predictedlabel<0.03]<-"Negative"
 with(data_19, CrossTable(ichorDetection, ALFAssay))
 
 data_19_healthy <- merge(data_19, ichor_healthy, on = "PatientId")
 data_19_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 data_20 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t20.csv", sep="\t")
 data_20$Predictedlabel <- data_20$Predictedlabel/100
 data_20$Truelabel <- data_20$Truelabel/100
 data_20$ichorDetection[data_20$Truelabel>=0.03]<-"Positive"
 data_20$ichorDetection[data_20$Truelabel<0.03]<-"Negative"
 data_20$ALFAssay[data_20$Predictedlabel>=0.03]<-"Positive"
 data_20$ALFAssay[data_20$Predictedlabel<0.03]<-"Negative"
 with(data_20, CrossTable(ichorDetection, ALFAssay))
 
 data_20_healthy <- merge(data_20, ichor_healthy, on = "PatientId")
 data_20_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 
 data_21 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t21.csv", sep="\t")
 data_21$Predictedlabel <- data_21$Predictedlabel/100
 data_21$Truelabel <- data_21$Truelabel/100
 data_21$ichorDetection[data_21$Truelabel>=0.03]<-"Positive"
 data_21$ichorDetection[data_21$Truelabel<0.03]<-"Negative"
 data_21$ALFAssay[data_21$Predictedlabel>=0.03]<-"Positive"
 data_21$ALFAssay[data_21$Predictedlabel<0.03]<-"Negative"
 with(data_21, CrossTable(ichorDetection, ALFAssay))
 
 data_21_healthy <- merge(data_21, ichor_healthy, on = "PatientId")
 data_21_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 
 
 data_22 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t22.csv", sep="\t")
 data_22$Predictedlabel <- data_22$Predictedlabel/100
 data_22$Truelabel <- data_22$Truelabel/100
 data_22$ichorDetection[data_22$Truelabel>=0.03]<-"Positive"
 data_22$ichorDetection[data_22$Truelabel<0.03]<-"Negative"
 data_22$ALFAssay[data_22$Predictedlabel>=0.03]<-"Positive"
 data_22$ALFAssay[data_22$Predictedlabel<0.03]<-"Negative"
 with(data_22, CrossTable(ichorDetection, ALFAssay))
 
 data_22_healthy <- merge(data_22, ichor_healthy, on = "PatientId")
 data_22_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.03)
 
 data_22_training <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_train_ichorTF_t22.csv", sep="\t")
 
 
 data_23 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t23.csv", sep="\t")
 data_23$Predictedlabel <- data_23$Predictedlabel/100
 data_23$Truelabel <- data_23$Truelabel/100
 data_23$ichorDetection[data_23$Truelabel>=0.03]<-"Positive"
 data_23$ichorDetection[data_23$Truelabel<0.03]<-"Negative"
 data_23$ALFAssay[data_23$Predictedlabel>=0.03]<-"Positive"
 data_23$ALFAssay[data_23$Predictedlabel<0.03]<-"Negative"
 with(data_23, CrossTable(ichorDetection, ALFAssay))
 
 data_23_healthy <- merge(data_23, ichor_healthy, on = "PatientId")
 data_23_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 data_23_training <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_train_ichorTF_t23.csv", sep="\t")
 
 
 data_24 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t24.csv", sep="\t")
 data_24$Predictedlabel <- data_24$Predictedlabel/100
 data_24$Truelabel <- data_24$Truelabel/100
 data_24$ichorDetection[data_24$Truelabel>=0.03]<-"Positive"
 data_24$ichorDetection[data_24$Truelabel<0.03]<-"Negative"
 data_24$ALFAssay[data_24$Predictedlabel>=0.03]<-"Positive"
 data_24$ALFAssay[data_24$Predictedlabel<0.03]<-"Negative"
 with(data_24, CrossTable(ichorDetection, ALFAssay))
 
 data_24_healthy <- merge(data_24, ichor_healthy, on = "PatientId")
 data_24_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 
 
 data_25 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t25.csv", sep="\t")
 data_25$Predictedlabel <- data_25$Predictedlabel/100
 data_25$Truelabel <- data_25$Truelabel/100
 data_25$ichorDetection[data_25$Truelabel>=0.03]<-"Positive"
 data_25$ichorDetection[data_25$Truelabel<0.03]<-"Negative"
 data_25$ALFAssay[data_25$Predictedlabel>=0.03]<-"Positive"
 data_25$ALFAssay[data_25$Predictedlabel<0.03]<-"Negative"
 with(data_25, CrossTable(ichorDetection, ALFAssay))
 
 data_25_healthy <- merge(data_25, ichor_healthy, on = "PatientId")
 data_25_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 data_26 <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_ichorTF_t26.csv", sep="\t")
 data_26$Predictedlabel <- data_26$Predictedlabel/100
 data_26$Truelabel <- data_26$Truelabel/100
 data_26$ichorDetection[data_26$Truelabel>=0.03]<-"Positive"
 data_26$ichorDetection[data_26$Truelabel<0.03]<-"Negative"
 data_26$ALFAssay[data_26$Predictedlabel>=0.03]<-"Positive"
 data_26$ALFAssay[data_26$Predictedlabel<0.03]<-"Negative"
 with(data_26, CrossTable(ichorDetection, ALFAssay))
 
 data_26_healthy <- merge(data_26, ichor_healthy, on = "PatientId")
 data_26_healthy %>% filter(PatientId %in% ichor_healthy$PatientId &Predictedlabel >0.01)
 
 
 data_24_training <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_train_ichorTF_t24.csv", sep="\t")
 
 
 
 data_regression <- read.csv("~/PhD/PyCharmProjects/ALFAssay/validations/regression_results_short_t26.csv", sep="\t")
 data_regression_v <- data_regression %>% filter(Process=="Validation")
 
 data_regression_v$ichorDetection[data_regression_v$ichorTF>=0.03]<-"Positive"
 data_regression_v$ichorDetection[data_regression_v$ichorTF<0.03]<-"Negative"
 data_regression_v$ALFAssay[data_regression_v$prediction>=0.03]<-"Positive"
 data_regression_v$ALFAssay[data_regression_v$prediction<0.03]<-"Negative"
 with(data_regression_v, CrossTable(ichorDetection, ALFAssay))
 
 
 data_regression_v_healthy <- merge(data_regression_v, ichor_healthy, on = "PatientId")
 data_regression_v_healthy %>% filter(PatientId %in% ichor_healthy$PatientId & prediction >0.01)
 
 
 
 # Build the confusion matrix
 conf_matrix_2 <- table(GroundTruth = data_regression_v$ichorDetection, Prediction = data_regression_v$ALFAssay)
 conf_matrix_1 <- table(GroundTruth = data_26$ichorDetection, Prediction = data_26$ALFAssay)

 
 
 
 # helper function to compute metrics from a 2×2 confusion matrix
 calc_metrics <- function(cm) {
   # assuming
   #           Prediction
   # GroundTruth   Neg   Pos
   #        Neg   TN    FP
   #        Pos   FN    TP
   TN <- cm[1,1]
   FP <- cm[1,2]
   FN <- cm[2,1]
   TP <- cm[2,2]
   
   sens <- TP / (TP + FN)                    # Sensitivity = TP/(TP+FN)
   spec <- TN / (TN + FP)                    # Specificity = TN/(TN+FP)
   ppv  <- TP / (TP + FP)                    # PPV = TP/(TP+FP)
   npv  <- TN / (TN + FN)                    # NPV = TN/(TN+FN)
   
   data.frame(
     Sensitivity = sens,
     Specificity = spec,
     PPV         = ppv,
     NPV         = npv
   )
 }
 
 # compute for each model
 metrics_alfassay <- calc_metrics(conf_matrix_1)
 metrics_regression <- calc_metrics(conf_matrix_2)
 
 # add a column for model names
 metrics_alfassay$Model    <- "ALFAssay"
 metrics_regression$Model  <- "Basic regression model"
 
 # combine into one table
 results <- rbind(
   metrics_regression[, c("Model","Sensitivity","Specificity","PPV","NPV")],
   metrics_alfassay[, c("Model","Sensitivity","Specificity","PPV","NPV")]
   
 )
 
 # optionally, row-name reset
 rownames(results) <- NULL
 
 # view
 print(results)
 write.csv(results, paste0("plots",try, "/ALFAssay_Performance_Metrics.csv",sep=""), row.names = FALSE)
 
 data$snr$all_synergy %>% filter(PatientId=="NIPT-P0001-ArmPhaseI-SNRBE0025-W1-Screening_S118")
 snr <- data$snr$all_synergy 
 nrow( data$snr$all_synergy )
length(unique(data$snr$all_synergy$patient))



