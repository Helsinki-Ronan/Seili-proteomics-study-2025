



# 1.0 Load required libraries and data.
library(tidyverse)
library(viridis)
library(UpSetR)
setwd("output/")
all_studies_data <- read_csv("all_studies_and_samples.csv")


# 2.0 Reorder the 'study' variable as required.
data2 <- all_studies_data %>%
  mutate(study = factor(study, levels = c("Seili", "Mickleburgh", "Sawafuji", "Ntasi")))

# 2.1.1 Sort by decreasing "psm" within each study.
data3 <- data2 %>%
  arrange(study, desc(psm)) %>%
  mutate(sample_ID = factor(sample_ID, levels = unique(sample_ID)))

# 2.1.2 Rename some variables. 
data3<- data3 %>% 
  mutate(sample_ID = if_else(sample_ID %in% "EE5973", "294",
                     if_else(sample_ID %in% "EE5974", "318",
                     if_else(sample_ID %in% "EE5975", "319",
                     if_else(sample_ID %in% "EE5977", "296",
                     if_else(sample_ID %in% "EE5978", "301",
                     if_else(sample_ID %in% "EE5979", "302",
                     if_else(sample_ID %in% "EE5980", "293",
                      
                     if_else(sample_ID %in% "ZH1656", "1-ZH1656",
                     if_else(sample_ID %in% "ZH1657", "1-ZH1657",
                     if_else(sample_ID %in% "ZH1658", "1-ZH1658",
                     if_else(sample_ID %in% "ZH1659", "1-ZH1659",
                     if_else(sample_ID %in% "ZH1670", "1-ZH1670",
                     if_else(sample_ID %in% "ZH1672", "2-ZH1672",
                     if_else(sample_ID %in% "ZH1673", "2-ZH1673",
                     if_else(sample_ID %in% "ZH1676", "3-ZH1676", sample_ID))))))))))))))))


data3$material<- if_else(data3$sample_ID %in% "1-ZH1656", "coffin_bedding",
                 if_else(data3$sample_ID %in% "1-ZH1657", "coffin_bedding",
                 if_else(data3$sample_ID %in% "1-ZH1658", "bone",
                 if_else(data3$sample_ID %in% "1-ZH1659", "soft_tissue",
                 if_else(data3$sample_ID %in% "1-ZH1670", "soft_tissue",
                 if_else(data3$sample_ID %in% "2-ZH1672", "bone",
                 if_else(data3$sample_ID %in% "2-ZH1673", "soft_tissue",
                 if_else(data3$sample_ID %in% "3-ZH1676", "bone", 
                 if_else(data3$sample_ID %in% c("MS","BS"), "buccal_swab",
                 if_else(data3$sample_ID %in% "BC", "cloth", "bone"))))))))))

data3$material2<- if_else(data3$sample_ID %in% "1-ZH1656", "Coffin\nbedding",
                  if_else(data3$sample_ID %in% "1-ZH1657", "Coffin\nbedding",
                  if_else(data3$sample_ID %in% "1-ZH1658", "Bone\n ",
                  if_else(data3$sample_ID %in% "1-ZH1659", "Soft\ntissue",
                  if_else(data3$sample_ID %in% "1-ZH1670", "Soft\ntissue",
                  if_else(data3$sample_ID %in% "2-ZH1672", "Bone\n ",
                  if_else(data3$sample_ID %in% "2-ZH1673", "Soft\ntissue",
                  if_else(data3$sample_ID %in% "3-ZH1676", "Bone\n ", 
                  if_else(data3$sample_ID %in% c("MS","BS"), "Buccal swab",
                  if_else(data3$sample_ID %in% "BC", "Cloth", "Bone\n "))))))))))



# 3.0 High-resolution figures with individual protein names for supplementary materials----
data3 <- data3 %>%
  mutate(x_label = paste0(sample_ID, "\n\n", material2))
#write.csv(data3, "for_stats.csv", row.names = FALSE)

all_PSM5_highres_for_supplementary<- 
  ggplot(data3 %>% subset(psm >= 5), aes(x_label, protein_name, fill = psm)) +
  geom_tile() +
  scale_fill_viridis(discrete = FALSE) +
  ylab("") +
  xlab("\nSample ID") + 
  ggtitle("Heatmap of study-specific peptide sequence matches (PSM), where PSM >=5.") +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 3, colour = "black"),
   #  axis.text.y = element_blank(),
  #   axis.ticks.y = element_blank(),
     axis.text.x = element_text(size = 6, colour = "black", angle = 0, vjust = 0.5, hjust = 0.5)) +
  facet_grid(. ~ study, scales = "free_x", space = "free_x")


# 4.0 Make an overall list of all proteins recovered and study-specific lists----
# 4.1 These lists are of all proteins recovered in each study----
write.table(data3 %>% filter(study %in% "Sawafuji"), "all_proteins_PSM2_Sawafuji.txt")
write.table(data3 %>% filter(study %in% "Seili"), "all_proteins_PSM2_Seili.txt")
write.table(data3 %>% filter(study %in% "Mickleburgh"), "all_proteins_PSM2_Mickleburgh.txt")
write.table(data3 %>% filter(study %in% "Ntasi"), "all_proteins_PSM2_Ntasi.txt")


# 4.2 These lists are of the unique proteins recovered in each study, i.e. only found in that particular study----
# 4.2.1 Without Mickleburgh----
unique_proteins_Sawafuji_noM <- data3 %>%
  filter(!study %in% "Mickleburgh") %>% 
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Sawafuji") %>% 
  distinct(protein_name)

unique_proteins_Seili_noM <- data3 %>%
  filter(!study %in% "Mickleburgh") %>% 
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Seili") %>%
  distinct(protein_name)

unique_proteins_Mickleburgh_noM <- data3 %>%
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Mickleburgh") %>%
  distinct(protein_name)

unique_proteins_Ntasi_noM <- data3 %>%
  filter(!study %in% "Mickleburgh") %>% 
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Ntasi") %>%
  distinct(protein_name)

write.table(unique_proteins_Sawafuji_noM, "unique_proteins_PSM2_Sawafuji_noM.txt")
write.table(unique_proteins_Seili_noM, "unique_proteins_PSM2_Seili_noM.txt")
write.table(unique_proteins_Mickleburgh_noM, "unique_proteins_PSM2_Mickleburgh_noM.txt")
write.table(unique_proteins_Ntasi_noM, "unique_proteins_PSM2_Ntasi_noM.txt")

# 4.2.2 With Mickleburgh----
unique_proteins_Sawafuji_withM <- data3 %>%
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Sawafuji") %>% 
  distinct(protein_name)

unique_proteins_Seili_withM <- data3 %>%
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Seili") %>%
  distinct(protein_name)

unique_proteins_Mickleburgh_withM <- data3 %>%
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Mickleburgh") %>%
  distinct(protein_name)

unique_proteins_Ntasi_withM <- data3 %>%
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Ntasi") %>%
  distinct(protein_name)

write.table(unique_proteins_Sawafuji_withM, "unique_proteins_PSM2_Sawafuji_withM.txt")
write.table(unique_proteins_Seili_withM, "unique_proteins_PSM2_Seili_withM.txt")
write.table(unique_proteins_Mickleburgh_withM, "unique_proteins_PSM2_Mickleburgh_withM.txt")
write.table(unique_proteins_Ntasi_withM, "unique_proteins_PSM2_Ntasi_withM.txt")


# 4.2.3 With Mickleburgh and PSM >= 5----
unique_proteins_Sawafuji_withM_PSM5<- data3 %>%
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Sawafuji") %>% 
  filter(psm >= 5) %>% 
  distinct(protein_name) 


unique_proteins_Seili_withM_PSM5 <- data3 %>%
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Seili") %>%
  filter(psm >= 5) %>% 
  distinct(protein_name)


unique_proteins_Mickleburgh_withM_PSM5 <- data3 %>%
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Mickleburgh") %>%
  filter(psm >= 5) %>% 
  distinct(protein_name)

unique_proteins_Ntasi_withM_PSM5 <- data3 %>%
  group_by(protein_name) %>%
  filter(n_distinct(study) %in% 1 & study %in% "Ntasi") %>%
  filter(psm >= 5) %>%
  distinct(protein_name)

write.table(unique_proteins_Sawafuji_withM_PSM5, "unique_proteins_Sawafuji_withM_PSM5.txt")
write.table(unique_proteins_Seili_withM_PSM5, "unique_proteins_Seili_withM_PSM5.txt")
write.table(unique_proteins_Mickleburgh_withM_PSM5, "unique_proteins_Mickleburgh_withM_PSM5.txt")
write.table(unique_proteins_Ntasi_withM_PSM5, "unique_proteins_Ntasi_withM_PSM5.txt")




# 5.0 Upset plots----
# 5.1 List unique protein names per study----
protein_lists <- all_studies_data %>%
  group_by(study) %>%
  summarise(proteins = list(unique(protein_name))) %>%
  deframe()

# 5.2 Get overlap between studies----
combn(names(protein_lists), 2, function(pair) {
  a <- pair[1]
  b <- pair[2]
  overlap <- intersect(protein_lists[[a]], protein_lists[[b]])
  data.frame(study1 = a, study2 = b, n_overlap = length(overlap))
}, simplify = FALSE) %>%
  bind_rows()


# 5.3 Create dataframe for upsetR---
all_proteins_upSet_plot<- c(
  # Number of unique protein names in each study
  Seili = 279,
  Sawafuji = 82,
  Ntasi = 30,
  Mickleburgh = 118,
  # Number from each overlap intersection
  "Mickleburgh&Ntasi" = 23,
  "Mickleburgh&Sawafuji" = 46,
  "Mickleburgh&Seili" = 71,
  "Ntasi&Sawafuji" = 22,
  "Ntasi&Seili" = 27,
  "Sawafuji&Seili" = 56
)

# 5.4 Plot----
all<- upset(fromExpression(all_proteins_upSet_plot), 
      nintersects = 10, 
      nsets = 4, 
      order.by = "freq", 
      decreasing = TRUE, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 2, 
      point.size = 5, 
      line.size = 2, 
      main.bar.color = "#f8766dff",
      matrix.color = "#e66101ff", 
      sets.bar.color = "#f8766dff", 
      shade.color = "#7fcdbbff",
      shade.alpha = 1,
      matrix.dot.alpha = 1)




# 6.0 Upset plot for immune proteins----
# 6.1 Load data and merge with full data set----
immune_proteins<- read.delim("working_data/immuneproteins25062025.txt", comment.char="#", header = FALSE) %>% 
  rename(protein_name = "V1")


data_immune<- merge(immune_proteins, all_studies_data, by = "protein_name")

# 6.2 List unique immune protein names per study----
protein_lists_immune <- data_immune %>%
  group_by(study) %>%
  summarise(proteins = list(unique(protein_name))) %>%
  deframe()

# 6.3 Get overlap between studies----
combn(names(protein_lists_immune), 2, function(pair) {
  a <- pair[1]
  b <- pair[2]
  overlap <- intersect(protein_lists_immune[[a]], protein_lists_immune[[b]])
  data.frame(study1 = a, study2 = b, n_overlap = length(overlap))
}, simplify = FALSE) %>%
  bind_rows()

# 6.4 Create dataframe for upsetR---
immune_proteins_upSet_plot<- c(
  # Number of unique protein names in each study
  Seili = 48,
  Sawafuji = 17,
  Mickleburgh = 33,
  # Number from each overlap intersection
  "Mickleburgh&Sawafuji" = 9,
  "Mickleburgh&Seili" = 19,
  "Sawafuji&Seili" = 9
)

# 6.5 Plot----
immune<- upset(fromExpression(immune_proteins_upSet_plot), 
      nintersects = 10, 
      nsets = 4, 
      order.by = "freq", 
      decreasing = TRUE, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 2, 
      point.size = 5, 
      line.size = 2, 
      main.bar.color = "#f8766dff",
      matrix.color = "#e66101ff", 
      sets.bar.color = "#f8766dff", 
      shade.color = "#7fcdbbff",
      shade.alpha = 1,
      matrix.dot.alpha = 1)


# 6.6 Names of immune proteins per study and overlap----
write.table(mickleburgh_immune<- protein_lists_immune[["Mickleburgh"]], "mickleburgh_immune.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(seili_immune<- protein_lists_immune[["Seili"]], "seili_immune.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(sawafuji_immune<- protein_lists_immune[["Sawafuji"]], "sawafuji_immune.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

overlap_immune_Mickleburgh_Seili<- write.table(intersect(protein_lists_immune[["Mickleburgh"]], protein_lists_immune[["Seili"]]), 
                                               "overlap_immune_Mickleburgh_Seili.txt",
                                               quote = FALSE, row.names = FALSE, col.names = FALSE)

overlap_immune_Mickleburgh_Sawafuji<- write.table(intersect(protein_lists_immune[["Mickleburgh"]], protein_lists_immune[["Sawafuji"]]), 
                                               "overlap_immune_Mickleburgh_Sawafuji.txt",
                                               quote = FALSE, row.names = FALSE, col.names = FALSE)

overlap_immune_Seili_Sawafuji<- write.table(intersect(protein_lists_immune[["Seili"]], protein_lists_immune[["Sawafuji"]]), 
                                                  "overlap_immune_Seili_Sawafuji.txt",
                                                  quote = FALSE, row.names = FALSE, col.names = FALSE)
