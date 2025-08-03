setwd("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts")
source("survivalFunctionR.R")

ctDNADetection <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/denoising/plots/denoising_try1/predictions_PRL_gc_corrected.csv", sep="\t" )

with(ctDNADetection, CrossTable(Prediction, Truth))

data <- read.csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
data$status <- data$survivalStatus
data$BCFS <-  data$survivalTime


data_pred <- merge(data,ctDNADetection, by="PatientId" )


# data_pred <- data_pred %>% filter(!is.na(ichorTF))
data_pred <- data_pred %>% filter(!is.na(BCFS))
# data_pred <- data_pred %>% filter(!is.na(ALFAssay))

data_pearl_bl <- data_pred %>% filter(timepoint=="BL" & study=="Pearl")
data_pearl_d15 <- data_pred %>% filter(timepoint=="C1" & study=="Pearl")
data_pearl_pd <- data_pred %>% filter(timepoint=="Surgery" & study=="Pearl")


ichor <- list()

ichor[[1]] <- survival_plot(data_pearl_bl, "BCFS", "Prediction","Baseline")
ichor[[2]] <- survival_plot(data_pearl_d15, "BCFS", "Prediction","D15")
ichor[[3]] <- survival_plot(data_pearl_pd, "BCFS", "Prediction","PD")

library(ggtext)
caption_text <- "KM on Pearl data using denoised ALFAssay on ctDNADetection results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3

rey="_t10"
ggsave(paste0("plots", try, "/PearlALFAssayDenoisedCtDNADetectionSurvival.pdf",sep=""), p3, width=12, height=9)



##### TF

ichorTF <- read.csv("~/PhD/PyCharmProjects/ALFAssayNN/denoising/plots/denoising_try1/predictions_PRL_gc_corrected_tf.csv",sep="\t")

ichorTF$Prediction <- ichorTF$Prediction/100
ichorTF$Truth <- ichorTF$Truth/100
ichorTF$ALFAssay[ichorTF$Prediction>=0.03]<-"Positive" 
ichorTF$ALFAssay[ichorTF$Prediction<0.03]<-"Negative" 
ichorTF$ichorCNA[ichorTF$Truth>=0.03]<-"Positive" 
ichorTF$ichorCNA[ichorTF$Truth<0.03]<-"Negative" 

data <- read.csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
data$status <- data$survivalStatus
data$BCFS <-  data$survivalTime


data_pred <- merge(data,ichorTF, by="PatientId" )
with(data_pred, CrossTable(ALFAssay, ichorCNA))

# data_pred <- data_pred %>% filter(!is.na(ichorTF))
data_pred <- data_pred %>% filter(!is.na(BCFS))
# data_pred <- data_pred %>% filter(!is.na(ALFAssay))

data_pearl_bl <- data_pred %>% filter(timepoint=="BL" & study=="Pearl")
data_pearl_d15 <- data_pred %>% filter(timepoint=="C1" & study=="Pearl")
data_pearl_pd <- data_pred %>% filter(timepoint=="Surgery" & study=="Pearl")


ichor <- list()

ichor[[1]] <- survival_plot(data_pearl_bl, "BCFS", "ALFAssay","Baseline")
ichor[[2]] <- survival_plot(data_pearl_d15, "BCFS", "ALFAssay","D15")
ichor[[3]] <- survival_plot(data_pearl_pd, "BCFS", "ALFAssay","PD")

library(ggtext)
caption_text <- "KM on Pearl data using ALFAssay no denoising on TF results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3

rey="_t10"
ggsave(paste0("plots", try, "/PearlALFAssayNoDenoisedTFSurvival.pdf",sep=""), p3, width=12, height=9)




meta <- read.csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep ="\t")

meta_healtthy <- meta %>% filter(study =="Healthy")

meta_healtthy$PatientId
gsub("HP-84", "", "Genome-IJB-HP-84_S189")
input_string<-"Genome-IJB-HP-84_S189"
output_string <- gsub("-{2,}", "-", gsub("HP-84", "", input_string))  # Removes double hyphens
output_string <- gsub("_+", "_", output_string)  # Cleans repeated underscores
output_string <- gsub("-_", "_", output_string) 

extracted <- regmatches(input_string, regexpr("HP-\\d+", input_string)
meta_healtthy$Patient <-  regmatches(meta_healtthy$PatientId, regexpr("HP-\\d+", meta_healtthy$PatientId))

meta_healtthy$matchedStudy <- ifelse(grepl("NIPT-PIJB", meta_healtthy$PatientId), "Synergy", "Neorhea")

meta_healtthy$common <- ifelse(, "Synergy", "Neorhea")

healthy_common <- meta_healtthy %>% select(PatientId,Patient, matchedStudy)

write.table(healthy_common, "/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/healthy_common.csv", row.names=FALSE)


Patients_with_Common_Column<- read.csv("~/Downloads/patients_with_common_cleaned.csv", sep =",")

Patients_with_Common_Column %>% filter(matchedStudy=="Neorhea" & common=="no")

