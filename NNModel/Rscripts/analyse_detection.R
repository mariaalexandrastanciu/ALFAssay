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
source("prepare_data.R")

try="_t10"

data <- get_data_2(try)


healthy_neorhea <- data$ctrl %>% filter(study=="Healthy")

# cols = c("PatientId", "Label" , "ichorTF","Predictedlabel")


pearl <- select(data$prl$all_pearl,"PatientId", "Label" , "ichorTF","Predictedlabel", "study", "VAF")
healthy_neorhea_subset <- select(healthy_neorhea,"PatientId", "Label" , "ichorTF","Predictedlabel", "study", "VAF")

unseen_cohort <- rbind(pearl, healthy_neorhea_subset)

unseen_cohort$Category <- "Controls"
unseen_cohort$Category[unseen_cohort$study=="Pearl"]<-"Pearl"
# unseen_cohort$Category[unseen_cohort$study=="Neorhea"]<-"Early HR+/HER2- BreastCancer"

stat.test_unseen <- unseen_cohort %>% 
  wilcox_test(Predictedlabel~Category, p.adjust.method = 'BH') %>%
  add_xy_position(x = "Category")  #%>% filter(!(group1==0 & group2 %in% c(">0.01 and <=0.1",">0.1 and <=0.2", ">0.2"))) %>% filter(!(group1==">0 and <=0.01" & group2 %in% c(">0.1 and <=0.2", ">0.2"))) %>% filter(!(group1==">0.01 and <=0.1" & group2 %in% c(">0.2")))


prediction_over_unseen_cohorts <- ggplot(unseen_cohort,
       aes(x = Category, y = Predictedlabel)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  # geom_hline(yintercept = 0.01, color = "red", lty = 2) +
  labs(x = "Cohort", y = "ALFAssay ctDNA fraction prediction") +
  ggtitle("ALFAssay ctDNA fraction prediction on VALIDATION") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_pvalue_manual(stat.test_unseen,  label = 'p')


healthy_synergy <- data$ctrl %>% filter(study=="Healthy_synergy")

# cols = c("PatientId", "Label" , "ichorTF","Predictedlabel")

neorhea <- select(data$nrh$all_neorhea,"PatientId", "Label" , "ichorTF","Predictedlabel", "study", "VAF")
synergy <- select(data$snr$all_synergy,"PatientId", "Label" , "ichorTF","Predictedlabel", "study", "VAF")
healthy_synergy_subset <- select(healthy_synergy,"PatientId", "Label" , "ichorTF","Predictedlabel", "study", "VAF")

test_cohorts<- rbind(neorhea,synergy,healthy_synergy_subset)

test_cohorts$Category <- "Controls"
test_cohorts$Category[test_cohorts$study=="Synergy"]<-"Synergy"
test_cohorts$Category[test_cohorts$study=="Neorhea"]<-"Neorhea"


stat.test <- test_cohorts %>% 
  wilcox_test(Predictedlabel~Category, p.adjust.method = 'BH') %>%
  add_xy_position(x = "Category")  #%>% filter(!(group1==0 & group2 %in% c(">0.01 and <=0.1",">0.1 and <=0.2", ">0.2"))) %>% filter(!(group1==">0 and <=0.01" & group2 %in% c(">0.1 and <=0.2", ">0.2"))) %>% filter(!(group1==">0.01 and <=0.1" & group2 %in% c(">0.2")))


prediction_over_test_cohorts <- ggplot(test_cohorts,
                                  aes(x = Category, y = Predictedlabel)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  # geom_hline(yintercept = 0.01, color = "red", lty = 2) +
  labs(x = "Cohort", y = "ALFAssay ctDNA fraction prediction") +
  ggtitle("ALFAssay ctDNA fraction prediction on TRAINING") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_pvalue_manual(stat.test,  label = 'p')



prediction_alfassay <- cowplot::plot_grid(prediction_over_test_cohorts,   
                                                 prediction_over_unseen_cohorts,
                                                 labels = c("A", "B"),
                                                 ncol = 2, nrow=1)

ggsave(paste0("plots", try, "/plots_progress/PredictedTFOverBC.pdf", sep=""),prediction_alfassay,width=16, height=8)

# title <- cowplot::ggdraw() + cowplot::draw_label("ALFAssay tumor fraction prediction on breast cancer cohorts", fontface='bold')
# prediction_alfassay  <- cowplot::plot_grid(title, prediction_alfassay, ncol=1, rel_heights=c(0.1, 1))

### vaf categories unseen data - pearl
data_pred_complete <- unseen_cohort[complete.cases(unseen_cohort[ , 6:7]),]

data_pred_complete <- unseen_cohort

data_pred_complete$VAF[data_pred_complete$study=="Pearl"] <- data_pred_complete$VAF[data_pred_complete$study=="Pearl"]/100

data_pred_complete$VAF_category <- 0
data_pred_complete$VAF_category[data_pred_complete$ichorTF>0 & data_pred_complete$ichorTF<=0.01] <- ">0 and <=0.01"
data_pred_complete$VAF_category[data_pred_complete$ichorTF>0.01 & data_pred_complete$ichorTF<=0.1] <- ">0.01 and <=0.1"
data_pred_complete$VAF_category[data_pred_complete$ichorTF>0.1 & data_pred_complete$ichorTF<=0.2] <- ">0.1 and <=0.2"
data_pred_complete$VAF_category[data_pred_complete$ichorTF>0.2 ] <- ">0.2"
data_pred_complete$VAF_category <- factor(data_pred_complete$VAF_category, levels=c("0", ">0 and <=0.01", ">0.01 and <=0.1", ">0.1 and <=0.2", ">0.2"))

stat.test <- data_pred_complete %>% 
  wilcox_test(Predictedlabel~VAF_category, p.adjust.method = 'BH') %>%
  add_xy_position(x = "VAF_category") %>% filter(!(group1==0 & group2 %in% c(">0.01 and <=0.1",">0.1 and <=0.2", ">0.2"))) %>% filter(!(group1==">0 and <=0.01" & group2 %in% c(">0.1 and <=0.2", ">0.2"))) %>% filter(!(group1==">0.01 and <=0.1" & group2 %in% c(">0.2")))

prediction_over_vaf_pearl_unseen <- ggplot(data_pred_complete,
                                  aes(x = VAF_category, y = Predictedlabel)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  # geom_hline(yintercept = 0.01, color = "red", lty = 2) +
  labs(x = "expected TF", y = "ALFAssay ctDNA fraction prediction") +
  ggtitle("ALFAssay ctDNA fraction prediction on expected TF for unseen cohorts(Healthy, Pearl)") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_pvalue_manual(stat.test,  label = 'p')


## vaf categories on test data

data_pred_complete_test <- test_cohorts[complete.cases(test_cohorts[ , 6:7]),]
data_pred_complete_test <- test_cohorts

data_pred_complete_test$VAF_category <- 0
data_pred_complete_test$VAF_category[data_pred_complete_test$ichorTF>0 & data_pred_complete_test$ichorTF<=0.01] <- ">0 and <=0.01"
data_pred_complete_test$VAF_category[data_pred_complete_test$ichorTF>0.01 & data_pred_complete_test$ichorTF<=0.1] <- ">0.01 and <=0.1"
data_pred_complete_test$VAF_category[data_pred_complete_test$ichorTF>0.1 & data_pred_complete_test$ichorTF<=0.2] <- ">0.1 and <=0.2"
data_pred_complete_test$VAF_category[data_pred_complete_test$ichorTF>0.2 ] <- ">0.2"
data_pred_complete_test$VAF_category <- factor(data_pred_complete_test$VAF_category, levels=c("0", ">0 and <=0.01", ">0.01 and <=0.1", ">0.1 and <=0.2", ">0.2"))

stat.test <- data_pred_complete_test %>% 
  wilcox_test(Predictedlabel~VAF_category, p.adjust.method = 'BH') %>%
  add_xy_position(x = "VAF_category") %>% filter(!(group1==0 & group2 %in% c(">0.01 and <=0.1",">0.1 and <=0.2", ">0.2"))) %>% filter(!(group1==">0 and <=0.01" & group2 %in% c(">0.1 and <=0.2", ">0.2"))) %>% filter(!(group1==">0.01 and <=0.1" & group2 %in% c(">0.2")))

prediction_over_vaf_test <- ggplot(data_pred_complete_test,
                                         aes(x = VAF_category, y = Predictedlabel)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  # geom_hline(yintercept = 0.01, color = "red", lty = 2) +
  labs(x = "expected TF", y = "ALFAssay ctDNA fraction prediction") +
  ggtitle("ALFAssay ctDNA fraction prediction on expected TF for test cohorts(Healthy,Synergy, Neorhea)") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_pvalue_manual(stat.test,  label = 'p')



prediction_over_vaf_levels <- cowplot::plot_grid(prediction_over_test_cohorts,   
                                                 prediction_over_unseen_cohorts,
                                                 prediction_over_vaf_test,
                                                 prediction_over_vaf_pearl_unseen,
                                          labels = c("A", "B", "C", "D"),
                                          ncol = 2, nrow=2)

title <- cowplot::ggdraw() + cowplot::draw_label("ALFAssay tumor fraction analysis on breast cancer cohorts and on expected VAF categories", fontface='bold')
prediction_over_vaf_levels <- cowplot::plot_grid(title, prediction_over_vaf_levels, ncol=1, rel_heights=c(0.1, 1))

ggsave(paste0("plots", try, "/PredictedTFOverBC.pdf", sep=""),prediction_over_vaf_levels,width=16, height=10)
