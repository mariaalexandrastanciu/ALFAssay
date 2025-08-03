
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
try="_t24"
data = get_data_2(try)
pearl <- data$prl$all_pearl

pearl$timepoint <- factor(pearl$timepoint ,
                    levels = c("BL","D15", "PD"),ordered = TRUE)
# plot only TF by timepoint with lines connecting patients
pearl_alfassay <- ggplot(pearl,
            aes(x = timepoint, y = Predictedlabel)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.03, color = "red", lty = 2) +
  labs(x = "Time point", y = "ALFAssay log10 TF") +
  ggtitle("ALFAssay results on Pearl") +
  theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5))

pearl_ichor <- ggplot(pearl,
                         aes(x = timepoint, y = ichorTF)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.03, color = "red", lty = 2) +
  labs(x = "Time point", y = "ichorCNA log10 TF") +
  ggtitle("ichorCNA results on Pearl") +
  theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(pearl, aes(x = factor(timepoint), y = ichorTF)) +
  # draw individual patient trajectories
  geom_line(aes(group = patient),
            color = "grey70",
            alpha = 0.5,
            size = 0.5) +
  # then overlay the points with a bit of jitter
  geom_jitter(width  = 0.2,
              height = 0,
              size   = 2,
              alpha  = 0.7,
              color  = "steelblue") +
  geom_hline(yintercept = 0.03,
             color      = "red",
             linetype   = "dashed") +
  labs(
    x     = "Time point",
    y     = "ichorCNA log10 TF",
    title = "ichorCNA Results on Pearl"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )




pearl_fragle <- ggplot(pearl,
                      aes(x = timepoint, y = ctDNA_Burden)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.01, color = "red", lty = 2) +
  labs(x = "Time point", y = "Fragle log10 TF") +
  ggtitle("Fragle results on Pearl") +
  theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5))

pearl_patients_complete <- pearl[complete.cases(pearl[ , 21:22]),] %>% group_by(patient) %>% summarise(count = n() ) %>% filter(count >2)
pearl_onco<- pearl %>% filter(patient %in% pearl_patients_complete$patient)

pearl_oncofollow <- ggplot(pearl,
                       aes(x = timepoint, y = VAF_pers)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.01, color = "red", lty = 2) +
  labs(x = "Time point", y = "OncoFollow DNA TF") +
  ggtitle("OncoFollow results on Pearl") +
  theme_classic(base_size=10) +
  # coord_trans(y='log10')  +
  theme(plot.title = element_text(hjust = 0.5))


pearl_timepoints <- cowplot::plot_grid(pearl_alfassay,   pearl_ichor,pearl_fragle,pearl_oncofollow,
                                          labels = c("A", "B", "C", "D"),
                                          ncol = 2, nrow=2)

ggsave(paste0("plots", try, "/PearlTimepoints.pdf", sep=""),pearl_timepoints,width=15, height=8)


## neorhea

neorhea <- data$nrh$all_neorhea

neorhea$timepoint <- factor(neorhea$timepoint ,
                          levels = c("BL","C1", "Surgery"),ordered = TRUE)
# plot only TF by timepoint with lines connecting patients
neorhea_alfassay <- ggplot(neorhea,
                         aes(x = timepoint, y = Predictedlabel)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.03, color = "red", lty = 2) +
  labs(x = "Time point", y = "ALFAssay log10 TF") +
  ggtitle("ALFAssay results on Neorhea") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5))

neorhea_ichor <- ggplot(neorhea,
                      aes(x = timepoint, y = ichorTF)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.03, color = "red", lty = 2) +
  labs(x = "Time point", y = "ichorCNA TF") +
  ggtitle("ichorCNA results on Neorhea") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5))

neorhea_fragle <- ggplot(neorhea,
                       aes(x = timepoint, y = ctDNA_Burden)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.01, color = "red", lty = 2) +
  labs(x = "Time point", y = "Fragle log10 TF") +
  ggtitle("Fragle results on Neorhea") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5))

neorhea_patients_complete <- neorhea[complete.cases(neorhea[ , 6:7]),] %>% group_by(patient) %>% summarise(count = n() ) %>% filter(count >2)
neorhea_radar<- neorhea %>% filter(patient %in% neorhea_patients_complete$patient)
neorhea_radar$VAF_perc<- neorhea_radar$VAF*100
neorhea_radar <- ggplot(neorhea,
                           aes(x = timepoint, y = VAF_perc)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.01, color = "red", lty = 2) +
  labs(x = "Time point", y = "Radar DNA % TF") +
  ggtitle("Radar results on Neorhea") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')  +
  theme(plot.title = element_text(hjust = 0.5))


neorhea_timepoints <- cowplot::plot_grid(neorhea_alfassay,   neorhea_ichor,neorhea_fragle,neorhea_radar,
                                       labels = c("A", "B", "C", "D"),
                                       ncol = 2, nrow=2)

ggsave(paste0("plots", try, "/NeorheaTimepoints.pdf",sep=""),neorhea_timepoints,width=15, height=8)



## synergy

synergy <- data$snr$all_synergy

synergy$timepoint <- factor(synergy$timepoint ,
                            levels = c("W1","W3", "W13"),ordered = TRUE)
# plot only TF by timepoint with lines connecting patients
synergy_alfassay <- ggplot(synergy,
                           aes(x = timepoint, y = Predictedlabel)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.03, color = "red", lty = 2) +
  labs(x = "Time point", y = "ALFAssay  TF") +
  ggtitle("ALFAssay results on Synergy") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5))

synergy_ichor <- ggplot(synergy,
                        aes(x = timepoint, y = ichorTF)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.03, color = "red", lty = 2) +
  labs(x = "Time point", y = "ichorCNA TF") +
  ggtitle("ichorCNA results on Synergy") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5))

synergy_fragle <- ggplot(synergy,
                         aes(x = timepoint, y = ctDNA_Burden)) +
  geom_boxplot() +
  # geom_line(aes(group = patient), alpha = 0.1) +
  geom_hline(yintercept = 0.01, color = "red", lty = 2) +
  labs(x = "Time point", y = "Fragle  TF") +
  ggtitle("Fragle results on Synergy") +
  # theme_classic(base_size=10) +
  # coord_trans(y='log10')+
  theme(plot.title = element_text(hjust = 0.5))



synergy_timepoints <- cowplot::plot_grid(synergy_alfassay,   synergy_ichor,synergy_fragle,
                                         labels = c("A", "B", "C"),
                                         ncol = 2, nrow=2)

ggsave(paste0("plots", try, "/SynergyTimepoints.pdf", sep=""),synergy_timepoints,width=15, height=8)

