# META DATA FRAME GENERATION FOR NEW SEGMENTATIONS #

# SET WD
wd.filepath <- "Z:/Raredon_Lab_Personal_Folders/Vi/Image Quantification/SimpleSeg"
setwd(wd.filepath)

# LOAD PACKAGES
library(ggplot2)
library(dplyr)

# LOAD DATA FRAMES
load("NEW Assay 1 Norm Data Frame.Robj")
load("NEW Assay 2 Norm Data Frame.Robj")


#### ADD CELL ID ####
df.meta.1 <- df.new.1 %>%
  mutate(
    Assay = 1
  )

df.meta.2 <- df.new.2 %>%
  mutate(
    Assay = 2
  )

m1cols <- colnames(df.meta.1)[!(colnames(df.meta.1) %in% colnames(df.meta.2))]
m2cols <- colnames(df.meta.2)[!(colnames(df.meta.2) %in% colnames(df.meta.1))]

# df.meta.1long <- df.meta.1
df.meta.1[m2cols] <- NA
# df.meta.2long <- df.meta.2
df.meta.2[m1cols] <- NA

meta <- rbind(df.meta.1, df.meta.2)

meta$`Cell ID` <- paste(meta$imageID, meta$object_id)
# meta <- meta %>%
#   relocate(`Cell ID`)

#### EXTRACT PASSAGE DATA ####
for (i in 1:nrow(meta)){
  if (startsWith(meta$`Cell ID`[i],"P1"))
  {meta$Passage[i] = "P1"}
  
  else if (startsWith(meta$`Cell ID`[i],"P3"))
  {meta$Passage[i] = "P3"}
  
  else if (startsWith(meta$`Cell ID`[i],"P5"))
  {meta$Passage[i] = "P5"}
  
}

#### EXTRACT LINEAGE DATA ####
for (i in 1:nrow(meta)){
  if (grepl("#12", meta$imageID[i]))
  {meta$Lineage[i] = "12"}
  
  else if (grepl("#4", meta$imageID[i]))
  {meta$Lineage[i] = "4"}
  
  else if (grepl("#5", meta$imageID[i]))
  {meta$Lineage[i] = "5"}
  
  else if (grepl("6", meta$imageID[i]))
  {meta$Lineage[i] = "6"}
}

#### WRITE REPLICATE DATA ####
for (i in 1:nrow(meta)){
  if (meta$Lineage[i] == 5)
  {meta$Replicate[i] = "1"}
  
  else if (meta$Lineage[i] == 6)
  {meta$Replicate[i] = "2"}
  
  else
  {meta$Replicate[i] = "3"}
}

save(meta, file = "NEW Norm Complete Data Frame.Robj")


#### ARRANGE COLUMNS ####
meta <- meta %>%
  relocate(
    `Assay`,
    `Cell ID`,
    `Passage`,
    `Lineage`,
    `Replicate`
  ) %>%
  mutate(
    Feature = NA
  ) %>%
  relocate(
    Feature,
    .after = objectNum
  )

save(meta, file = "NEW Norm Complete Data Frame.Robj")

#### REFRAME META DATA FRAME ####

# Extract each feature
metaDAPI <- meta %>%
  filter(!is.na(meta$DAPI.Raw)) %>%
  select(Assay:Feature, DAPI.Raw) %>%
  mutate(Feature = 'DAPI') %>%
  rename(`Intensity` = DAPI.Raw)

metaEPCAM <- meta %>%
  filter(!is.na(meta$EPCAM.Raw)) %>%
  select(Assay:Feature, EPCAM.Raw) %>%
  mutate(Feature = 'EPCAM') %>%
  rename(`Intensity` = EPCAM.Raw)

metaAger <- meta %>%
  filter(!is.na(meta$Ager.Raw)) %>%
  select(Assay:Feature, Ager.Raw) %>%
  mutate(Feature = 'Ager') %>%
  rename(`Intensity` = Ager.Raw)

metaKrt5 <- meta %>%
  filter(!is.na(meta$Krt5.Raw)) %>%
  select(Assay:Feature, Krt5.Raw) %>%
  mutate(Feature = 'Krt5') %>%
  rename(`Intensity` = Krt5.Raw)

metaABCA3 <- meta %>%
  filter(!is.na(meta$ABCA3.Raw)) %>%
  select(Assay:Feature, ABCA3.Raw) %>%
  mutate(Feature = 'ABCA3') %>%
  rename(`Intensity` = ABCA3.Raw)

metaSox2 <- meta %>%
  filter(!is.na(meta$Sox2.Raw)) %>%
  select(Assay:Feature, Sox2.Raw) %>%
  mutate(Feature = 'Sox2') %>%
  rename(`Intensity` = Sox2.Raw)

metaSox9 <- meta %>%
  filter(!is.na(meta$Sox9.Raw)) %>%
  select(Assay:Feature, Sox9.Raw) %>%
  mutate(Feature = 'Sox9') %>%
  rename(`Intensity` = Sox9.Raw)

bigmeta <- rbind(metaDAPI, metaEPCAM, metaAger, metaKrt5, metaABCA3, metaSox2, metaSox9)
save(bigmeta, file = "NEW Meta Data Frame by Feature 06.21.2025.Robj")


## ASSIGN POS/NEG ####
#Initialize thresholds
dapi.thresh = 0.10
epcam.thresh = 0.05
ager.thresh = 0.06
krt5.thresh = 0.0275
abca3.thresh = 0.05
sox2.thresh = 0.04
sox9.thresh = 0.035

# Update dataframes
metaDAPI <- metaDAPI %>%
  mutate(`Pos/Neg` = (if_else(Intensity >= dapi.thresh, true = "Pos", false = "Neg")))

metaEPCAM <- metaEPCAM %>%
  mutate(`Pos/Neg` = (if_else(Intensity >= epcam.thresh, true = "Pos", false = "Neg")))

metaAger <- metaAger %>%
  mutate(`Pos/Neg` = (if_else(Intensity >= ager.thresh, true = "Pos", false = "Neg")))

metaKrt5 <- metaKrt5 %>%
  mutate(`Pos/Neg` = (if_else(Intensity >= krt5.thresh, true = "Pos", false = "Neg")))

metaABCA3 <- metaABCA3 %>%
  mutate(`Pos/Neg` = (if_else(Intensity >= abca3.thresh, true = "Pos", false = "Neg")))

metaSox2 <- metaSox2 %>%
  mutate(`Pos/Neg` = (if_else(Intensity >= sox2.thresh, true = "Pos", false = "Neg")))

metaSox9 <- metaSox9 %>%
  mutate(`Pos/Neg` = (if_else(Intensity >= sox9.thresh, true = "Pos", false = "Neg")))

bigmeta <- rbind(metaDAPI, metaEPCAM, metaAger, metaKrt5, metaABCA3, metaSox2, metaSox9)
save(bigmeta, file = "NEW Meta Data Frame by Feature 06.21.2025.Robj")

# Summarize Pos/Neg Ratios
CtsPcts <- bigmeta %>%
  filter(Feature != 'DAPI') %>%
  group_by(Feature, Passage) %>%
  summarise(PosCt = sum(`Pos/Neg` == 'Pos'),
            NegCt = sum(`Pos/Neg` == 'Neg'),
            TotalCt = PosCt + NegCt,
            PosPct = PosCt/TotalCt,
            NegPct = NegCt/TotalCt,
            TotalPct = PosPct + NegPct)

save(CtsPcts, file = "NEW Table Summarizing Pos.Neg Counts.Percents.Robj")

pctplot <- ggplot(CtsPcts,
                  aes(x = Passage, y = PosPct, fill = Passage)) +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent'),
        panel.border = element_rect(colour = "black", fill='transparent')) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#fb5607", "#8338ec","#76c893")) +
  ylim(0,1) +
  labs(x = 'Passage',
       y = 'Percent Positive') +
  facet_wrap(~Feature, nrow = 1)

png(filename = "NEW Percent Positive Feature Expression Bar Plots in a Row BW.png", width = 8, height = 1.85, units = 'in', res = 350)
pctplot
dev.off()


