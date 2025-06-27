



# Statistical analysis


# 1.0 Load libraries and data----
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(DHARMa)
setwd("output/")


data<- read_csv("for_stats.csv")
data$study<- as.factor(as.character(data$study))
data$material2<- as.factor(as.character(data$material2))

immune_proteins<- read.csv("working_data/immuneproteins26062025.txt", comment.char="#", header = FALSE) %>% 
  rename(protein_name = "V1")

data2<- merge(immune_proteins, data, by = "protein_name")

# 2.0 Explore data----
# 2.1 Plot distribution of PSM across both study and material sampled.
ggplot(data %>% subset(psm < 1000), aes(sample_ID, psm, fill = material2)) +
  geom_boxplot()+
  ylab("") +
  xlab("\nStudy")+ 
  labs(fill = "Material sampled")+
  ggtitle("Distributions of sample- and tissue-specific peptide sequence matches (PSM) where PSM = >= 2\nProteins with PSM >= 50 labeled") +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black", angle = 0, vjust = 0.5, hjust = 0.5))+
  geom_text_repel(data = data %>% subset(psm >= 50 & psm < 1000), aes(sample_ID, psm, label = protein_name), 
                  col = "black", size = 2.5,
                  max.overlaps = 100000,
                  box.padding = 0.25, 
                  point.padding = 0.1, 
                  min.segment.length = 0.1)+
  facet_grid(. ~ study, scales = "free_x", space = "free_x")



# 2.2 Look at data without collagen proteins and PSM less than 200
data_noCollagen<- data %>% 
  filter(!grepl("Collagen", protein_name))

ggplot(data_noCollagen %>% subset(psm < 200), aes(sample_ID, psm, fill = material2)) +
  geom_boxplot()+
  ylab("") +
  xlab("\nStudy")+ 
  labs(fill = "Material sampled")+
  ggtitle("Distributions of sample- and tissue-specific peptide sequence matches (PSM) where PSM = >= 2\nProteins with PSM >= 50 labeled. Collagen excluded") +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black", angle = 0, vjust = 0.5, hjust = 0.5))+
  geom_text_repel(data = data_noCollagen %>% subset(psm >= 50 & psm < 200), aes(sample_ID, psm, label = protein_name), 
                  col = "black", size = 2.5,
                  max.overlaps = 100000,
                  box.padding = 0.25, 
                  point.padding = 0.1, 
                  min.segment.length = 0.1)+
  facet_grid(. ~ study, scales = "free_x", space = "free_x")



# 2.3 Plot distribution of PSM across both study and material sampled for immune proteins----
ggplot(data2, aes(sample_ID, psm, fill = material2)) +
  geom_boxplot()+
  ylab("") +
  xlab("\nStudy")+ 
  labs(fill = "Material sampled")+
  ggtitle("Distributions of sample- and tissue-specific peptide sequence matches (PSM) for immune proteins where PSM = >= 2\nProteins with PSM >= 25 labeled") +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black", angle = 0, vjust = 0.5, hjust = 0.5))+
  geom_text_repel(data = data2 %>% subset(psm >= 25 & psm < 1000), aes(sample_ID, psm, label = protein_name), 
                  col = "black", size = 2.5,
                  max.overlaps = 100000,
                  box.padding = 0.25,
                  point.padding = 0.1, 
                  min.segment.length = 0.1)+
  facet_grid(. ~ study, scales = "free_x", space = "free_x")




# 3.0 Statistics
summary(data)
data$protein_name<- as.factor(as.character(data$protein_name))
data$sample_ID<- as.character(data$sample_ID)
data$material<- as.factor(as.character(data$material))


# 3.1 Add immune protein information to main dataframe----
data3<- data %>% 
  mutate(is_immuno = if_else(protein_name %in% data2$protein_name, "immuno", NA))
summary(data3)


# 3.2 Due to material and study being partially confounded, create a compound variable from both----
table(data3$material, data3$study)
data3$StudyTissue <- interaction(data3$study, data3$material)
data3$StudyTissue<- factor(data3$StudyTissue)
table(data3$StudyTissue)


# 4.0 All proteins comparison----
data4<- data3 %>% 
  mutate(bone_softTissue = if_else(StudyTissue %in% "Seili.soft_tissue", "Soft tissue",
                           if_else(StudyTissue %in% "Seili.coffin_bedding", "Coffin bedding", "Bone (all)"))) %>% 
  
  mutate(arch_forensic_st = if_else(bone_softTissue %in% "Soft tissue", "Archaeological soft tissue",
                                    ifelse(bone_softTissue %in% "Coffin bedding", NA,
                                    if_else(study %in% "Mickleburgh", "Forensic bone", "Archaeological bone"))))


model_1<- glmmTMB(psm ~ StudyTissue + (1 | sample_ID),
                   data = data4[data4$material %in% "bone", ],
                   dispformula = ~ study, 
                   family = truncated_nbinom2("log"))

model_2<- glmmTMB(psm ~ StudyTissue + (1 | sample_ID),
                   data = data4[data4$material %in% "bone", ],
                   # dispformula = ~ study, 
                   family = truncated_nbinom2("log"))

model_3<- glmmTMB(psm ~ 1 + (1 | sample_ID),
                   data = data4[data4$material %in% "bone", ],
                   # dispformula = ~ study, 
                   family = truncated_nbinom2("log"))


anova(model_1, model_2, model_3)
summary(model_1)
simres_model_1<- simulateResiduals(fittedModel = model_1, plot = TRUE)
testUniformity(simres_model_1)
testOutliers(simres_model_1)


# 5.0 Immune proteins per study and tissue type----
data4<- data4 %>% 
  filter(is_immuno %in% "immuno") %>% 
  filter(!bone_softTissue %in% "Coffin bedding")

table(data4$material, data4$study)
table(data4$bone_softTissue, data4$arch_forensic_st)


model_1 <- glmmTMB(psm ~ bone_softTissue + (1 | sample_ID),
                   data = data4,
                   dispformula = ~ study, 
                   family = truncated_nbinom2("log"))

model_2 <- glmmTMB(psm ~ bone_softTissue + (1 | sample_ID),
                   data = data4,
                   # dispformula = ~ study, 
                   family = truncated_nbinom2("log"))

model_3 <- glmmTMB(psm ~ 1 + (1 | sample_ID),
                   data = data4,
                   # dispformula = ~ study, 
                   family = truncated_nbinom2("log"))


anova(model_1, model_2, model_3)
summary(model_2)
simres_model_2<- simulateResiduals(fittedModel = model_2, plot = TRUE)
testUniformity(simres_model_2)
testOutliers(simres_model_2)

emm1<- emmeans(model_2, ~ bone_softTissue, type = "response", comparisons = TRUE)
pairs(emm1, reverse = TRUE)

emm1_df <- as.data.frame(emm1) 

immuneproteins<- ggplot(emm1_df, aes(x = bone_softTissue, y = response)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, lwd = 1) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        title = element_text(size=11))+
  ylab("")+
  xlab("")



model_4 <- glmmTMB(psm ~ arch_forensic_st + (1 | sample_ID),
                   data = data4,
                   dispformula = ~ study, 
                   family = truncated_nbinom2("log"))

model_5 <- glmmTMB(psm ~ arch_forensic_st + (1 | sample_ID),
                   data = data4,
                   # dispformula = ~ study, 
                   family = truncated_nbinom2("log"))

model_6 <- glmmTMB(psm ~ 1 + (1 | sample_ID),
                   data = data4,
                   # dispformula = ~ study, 
                   family = truncated_nbinom2("log"))


anova(model_4, model_5, model_6)
summary(model_5)
simres_model_5<- simulateResiduals(fittedModel = model_5, plot = TRUE)
testUniformity(simres_model_5)
testOutliers(simres_model_5)


library(sjPlot)
model_2_output<- tab_model(model_2,
                           transform = NULL,
                           file = NULL,
                           pred.labels = c("FE intercept: Bone", "Soft tissue", "Intercept:random effect"))


model_5_output<- tab_model(model_5,
                           transform = NULL,
                           file = NULL,
                           pred.labels = c("FE intercept: Archaeological bone", "Archaeological soft tissue", "Forensic bone", "Intercept:random effect"))




emm2<- emmeans(model_5, ~ arch_forensic_st, type = "response", comparisons = TRUE)
pairs(emm2, reverse = TRUE)

emm2_df <- as.data.frame(emm2) 

immuneproteins2<- ggplot(emm2_df, aes(x = arch_forensic_st, y = response)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, lwd = 1) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        title = element_text(size=11))+
  ylab("")+
  xlab("")

grid.arrange(immuneproteins,immuneproteins2, nrow = 1)
