##### Code for Generating All Marker Lists
Idents(merged.obj) <- merged.obj$Condition
######St. Dclk1
#####
StartingOverEmbryo<- FindMarkers(merged.obj, ident.1 = "St. Dclk1", ident.2 = "Embryogenesis", only.pos = T, logfc.threshold = 0.25)
StartingOverEmbryo$ratio <- StartingOverEmbryo$pct.1/StartingOverEmbryo$pct.2
StartingOverEmbryo$power <-StartingOverEmbryo$ratio*StartingOverEmbryo$avg_log2FC
StartingOverEmbryo<- StartingOverEmbryo %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
StartingOverEmbryo <- StartingOverEmbryo[order(-StartingOverEmbryo$power), ]

StartingOverEmbryoGenes <- rownames(StartingOverEmbryo)

StartingOverEmbryoGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartingOverEmbryoGenes, value = TRUE)
StartingOverEmbryoGenes.cleaned <- StartingOverEmbryoGenes[! StartingOverEmbryoGenes %in% StartingOverEmbryoGenes.toRemove]

StartingOverEmbryoGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartingOverEmbryoGenes.cleaned, value = TRUE)

StartingOverEmbryoGenes.cleaned.NoCC <- StartingOverEmbryoGenes.cleaned[! StartingOverEmbryoGenes.cleaned %in% StartingOverEmbryoGenes.CellCycle]

###StartOverP0
StartingOverP0<- FindMarkers(merged.obj, ident.1 = "St. Dclk1", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
StartingOverP0$ratio <- StartingOverP0$pct.1/StartingOverP0$pct.2
StartingOverP0$power <-StartingOverP0$ratio*StartingOverP0$avg_log2FC
StartingOverP0<- StartingOverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
StartingOverP0 <- StartingOverP0[order(-StartingOverP0$power), ]

StartingOverP0Genes <- rownames(StartingOverP0)

StartingOverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartingOverP0Genes, value = TRUE)
StartingOverP0Genes.cleaned <- StartingOverP0Genes[! StartingOverP0Genes %in% StartingOverP0Genes.toRemove]

StartingOverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartingOverP0Genes.cleaned, value = TRUE)

StartingOverP0Genes.cleaned.NoCC <- StartingOverP0Genes.cleaned[! StartingOverP0Genes.cleaned %in% StartingOverP0Genes.CellCycle]

###StartOverP3
StartingOverP3<- FindMarkers(merged.obj, ident.1 = "St. Dclk1", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
StartingOverP3$ratio <- StartingOverP3$pct.1/StartingOverP3$pct.2
StartingOverP3$power <-StartingOverP3$ratio*StartingOverP3$avg_log2FC
StartingOverP3<- StartingOverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
StartingOverP3 <- StartingOverP3[order(-StartingOverP3$power), ]

StartingOverP3Genes <- rownames(StartingOverP3)

StartingOverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartingOverP3Genes, value = TRUE)
StartingOverP3Genes.cleaned <- StartingOverP3Genes[! StartingOverP3Genes %in% StartingOverP3Genes.toRemove]

StartingOverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartingOverP3Genes.cleaned, value = TRUE)

StartingOverP3Genes.cleaned.NoCC <- StartingOverP3Genes.cleaned[! StartingOverP3Genes.cleaned %in% StartingOverP3Genes.CellCycle]

###StartOver24Hr
StartingOver24Hr<- FindMarkers(merged.obj, ident.1 = "St. Dclk1", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
StartingOver24Hr$ratio <- StartingOver24Hr$pct.1/StartingOver24Hr$pct.2
StartingOver24Hr$power <-StartingOver24Hr$ratio*StartingOver24Hr$avg_log2FC
StartingOver24Hr<- StartingOver24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
StartingOver24Hr <- StartingOver24Hr[order(-StartingOver24Hr$power), ]

StartingOver24HrGenes <- rownames(StartingOver24Hr)

StartingOver24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartingOver24HrGenes, value = TRUE)
StartingOver24HrGenes.cleaned <- StartingOver24HrGenes[! StartingOver24HrGenes %in% StartingOver24HrGenes.toRemove]

StartingOver24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartingOver24HrGenes.cleaned, value = TRUE)

StartingOver24HrGenes.cleaned.NoCC <- StartingOver24HrGenes.cleaned[! StartingOver24HrGenes.cleaned %in% StartingOver24HrGenes.CellCycle]


StartingOver5D<- FindMarkers(merged.obj, ident.1 = "St. Dclk1", ident.2 = "Engineered-FiveDay", only.pos = T, logfc.threshold = 0.25)
StartingOver5D$ratio <- StartingOver5D$pct.1/StartingOver5D$pct.2
StartingOver5D$power <-StartingOver5D$ratio*StartingOver5D$avg_log2FC
StartingOver5D<- StartingOver5D %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
StartingOver5D <- StartingOver5D[order(-StartingOver5D$power), ]

StartingOverFiveDayGenes <- rownames(StartingOver5D)

StartingOverFiveDayGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartingOverFiveDayGenes, value = TRUE)
StartingOverFiveDayGenes.cleaned <- StartingOverFiveDayGenes[! StartingOverFiveDayGenes %in% StartingOverFiveDayGenes.toRemove]

StartingOverFiveDayGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartingOverFiveDayGenes.cleaned, value = TRUE)

StartingOverFiveDayGenes.cleaned.NoCC <- StartingOverFiveDayGenes.cleaned[! StartingOverFiveDayGenes.cleaned %in% StartingOverFiveDayGenes.CellCycle]


StartInt1 <- intersect(StartingOverEmbryoGenes.cleaned.NoCC, 
                       StartingOverP0Genes.cleaned.NoCC)
StartInt2 <- intersect(StartInt1, StartingOverP3Genes.cleaned.NoCC)
StartInt3 <- intersect(StartInt2, StartingOver24HrGenes.cleaned.NoCC)
StartInt4 <- intersect(StartInt3, StartingOverFiveDayGenes.cleaned.NoCC)
#####


#####Engineered-FiveDay
#####
FiveDayOverEmbryo<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "Embryogenesis", only.pos = T, logfc.threshold = 0.25)
FiveDayOverEmbryo$ratio <- FiveDayOverEmbryo$pct.1/FiveDayOverEmbryo$pct.2
FiveDayOverEmbryo$power <-FiveDayOverEmbryo$ratio*FiveDayOverEmbryo$avg_log2FC
FiveDayOverEmbryo<- FiveDayOverEmbryo %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
FiveDayOverEmbryo <- FiveDayOverEmbryo[order(-FiveDayOverEmbryo$power), ]

FiveDayOverEmbryoGenes <- rownames(FiveDayOverEmbryo)

FiveDayOverEmbryoGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayOverEmbryoGenes, value = TRUE)
FiveDayOverEmbryoGenes.cleaned <- FiveDayOverEmbryoGenes[! FiveDayOverEmbryoGenes %in% FiveDayOverEmbryoGenes.toRemove]

FiveDayOverEmbryoGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayOverEmbryoGenes.cleaned, value = TRUE)

FiveDayOverEmbryoGenes.cleaned.NoCC <- FiveDayOverEmbryoGenes.cleaned[! FiveDayOverEmbryoGenes.cleaned %in% FiveDayOverEmbryoGenes.CellCycle]

###FiveDayOverP0
FiveDayOverP0<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
FiveDayOverP0$ratio <- FiveDayOverP0$pct.1/FiveDayOverP0$pct.2
FiveDayOverP0$power <-FiveDayOverP0$ratio*FiveDayOverP0$avg_log2FC
FiveDayOverP0<- FiveDayOverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
FiveDayOverP0 <- FiveDayOverP0[order(-FiveDayOverP0$power), ]

FiveDayOverP0Genes <- rownames(FiveDayOverP0)

FiveDayOverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayOverP0Genes, value = TRUE)
FiveDayOverP0Genes.cleaned <- FiveDayOverP0Genes[! FiveDayOverP0Genes %in% FiveDayOverP0Genes.toRemove]

FiveDayOverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayOverP0Genes.cleaned, value = TRUE)

FiveDayOverP0Genes.cleaned.NoCC <- FiveDayOverP0Genes.cleaned[! FiveDayOverP0Genes.cleaned %in% FiveDayOverP0Genes.CellCycle]

###FiveDayOverP3
FiveDayOverP3<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
FiveDayOverP3$ratio <- FiveDayOverP3$pct.1/FiveDayOverP3$pct.2
FiveDayOverP3$power <-FiveDayOverP3$ratio*FiveDayOverP3$avg_log2FC
FiveDayOverP3<- FiveDayOverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
FiveDayOverP3 <- FiveDayOverP3[order(-FiveDayOverP3$power), ]

FiveDayOverP3Genes <- rownames(FiveDayOverP3)

FiveDayOverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayOverP3Genes, value = TRUE)
FiveDayOverP3Genes.cleaned <- FiveDayOverP3Genes[! FiveDayOverP3Genes %in% FiveDayOverP3Genes.toRemove]

FiveDayOverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayOverP3Genes.cleaned, value = TRUE)

FiveDayOverP3Genes.cleaned.NoCC <- FiveDayOverP3Genes.cleaned[! FiveDayOverP3Genes.cleaned %in% FiveDayOverP3Genes.CellCycle]

###FiveDayOver24Hr
FiveDayOver24Hr<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
FiveDayOver24Hr$ratio <- FiveDayOver24Hr$pct.1/FiveDayOver24Hr$pct.2
FiveDayOver24Hr$power <-FiveDayOver24Hr$ratio*FiveDayOver24Hr$avg_log2FC
FiveDayOver24Hr<- FiveDayOver24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
FiveDayOver24Hr <- FiveDayOver24Hr[order(-FiveDayOver24Hr$power), ]

FiveDayOver24HrGenes <- rownames(FiveDayOver24Hr)

FiveDayOver24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayOver24HrGenes, value = TRUE)
FiveDayOver24HrGenes.cleaned <- FiveDayOver24HrGenes[! FiveDayOver24HrGenes %in% FiveDayOver24HrGenes.toRemove]

FiveDayOver24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayOver24HrGenes.cleaned, value = TRUE)

FiveDayOver24HrGenes.cleaned.NoCC <- FiveDayOver24HrGenes.cleaned[! FiveDayOver24HrGenes.cleaned %in% FiveDayOver24HrGenes.CellCycle]

##FiveDayOverStart
FiveDayOverStart<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "St. Dclk1", only.pos = T, logfc.threshold = 0.25)
FiveDayOverStart$ratio <- FiveDayOverStart$pct.1/FiveDayOverStart$pct.2
FiveDayOverStart$power <-FiveDayOverStart$ratio*FiveDayOverStart$avg_log2FC
FiveDayOverStart<- FiveDayOverStart %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
FiveDayOverStart <- FiveDayOverStart[order(-FiveDayOverStart$power), ]

FiveDayOverStartGenes <- rownames(FiveDayOverStart)

FiveDayOverStartGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayOverStartGenes, value = TRUE)
FiveDayOverStartGenes.cleaned <- FiveDayOverStartGenes[! FiveDayOverStartGenes %in% FiveDayOverStartGenes.toRemove]

FiveDayOverStartGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayOverStartGenes.cleaned, value = TRUE)

FiveDayOverStartGenes.cleaned.NoCC <- FiveDayOverStartGenes.cleaned[! FiveDayOverStartGenes.cleaned %in% FiveDayOverStartGenes.CellCycle]


FiveDayInt1 <- intersect(FiveDayOverEmbryoGenes.cleaned.NoCC, 
                         FiveDayOverP0Genes.cleaned.NoCC)
FiveDayInt2 <- intersect(FiveDayInt1, FiveDayOverP3Genes.cleaned.NoCC)
FiveDayInt3 <- intersect(FiveDayInt2, FiveDayOver24HrGenes.cleaned.NoCC)
FiveDayInt4 <- intersect(FiveDayInt3, FiveDayOverStartGenes.cleaned.NoCC)
#####

#####Embryogenesis
#####
EmbryoOverFiveDay<- FindMarkers(merged.obj, ident.1 = "Embryogenesis", ident.2 = "Engineered-FiveDay", only.pos = T, logfc.threshold = 0.25)
EmbryoOverFiveDay$ratio <- EmbryoOverFiveDay$pct.1/EmbryoOverFiveDay$pct.2
EmbryoOverFiveDay$power <-EmbryoOverFiveDay$ratio*EmbryoOverFiveDay$avg_log2FC
EmbryoOverFiveDay<- EmbryoOverFiveDay %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
EmbryoOverFiveDay <- EmbryoOverFiveDay[order(-EmbryoOverFiveDay$power), ]

EmbryoOverFiveDayGenes <- rownames(EmbryoOverFiveDay)

EmbryoOverFiveDayGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoOverFiveDayGenes, value = TRUE)
EmbryoOverFiveDayGenes.cleaned <- EmbryoOverFiveDayGenes[! EmbryoOverFiveDayGenes %in% EmbryoOverFiveDayGenes.toRemove]

EmbryoOverFiveDayGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoOverFiveDayGenes.cleaned, value = TRUE)

EmbryoOverFiveDayGenes.cleaned.NoCC <- EmbryoOverFiveDayGenes.cleaned[! EmbryoOverFiveDayGenes.cleaned %in% EmbryoOverFiveDayGenes.CellCycle]

###EmbryoOverP0
EmbryoOverP0<- FindMarkers(merged.obj, ident.1 = "Embryogenesis", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
EmbryoOverP0$ratio <- EmbryoOverP0$pct.1/EmbryoOverP0$pct.2
EmbryoOverP0$power <-EmbryoOverP0$ratio*EmbryoOverP0$avg_log2FC
EmbryoOverP0<- EmbryoOverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
EmbryoOverP0 <- EmbryoOverP0[order(-EmbryoOverP0$power), ]

EmbryoOverP0Genes <- rownames(EmbryoOverP0)

EmbryoOverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoOverP0Genes, value = TRUE)
EmbryoOverP0Genes.cleaned <- EmbryoOverP0Genes[! EmbryoOverP0Genes %in% EmbryoOverP0Genes.toRemove]

EmbryoOverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoOverP0Genes.cleaned, value = TRUE)

EmbryoOverP0Genes.cleaned.NoCC <- EmbryoOverP0Genes.cleaned[! EmbryoOverP0Genes.cleaned %in% EmbryoOverP0Genes.CellCycle]

###EmbryoOverP3
EmbryoOverP3<- FindMarkers(merged.obj, ident.1 = "Embryogenesis", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
EmbryoOverP3$ratio <- EmbryoOverP3$pct.1/EmbryoOverP3$pct.2
EmbryoOverP3$power <-EmbryoOverP3$ratio*EmbryoOverP3$avg_log2FC
EmbryoOverP3<- EmbryoOverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
EmbryoOverP3 <- EmbryoOverP3[order(-EmbryoOverP3$power), ]

EmbryoOverP3Genes <- rownames(EmbryoOverP3)

EmbryoOverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoOverP3Genes, value = TRUE)
EmbryoOverP3Genes.cleaned <- EmbryoOverP3Genes[! EmbryoOverP3Genes %in% EmbryoOverP3Genes.toRemove]

EmbryoOverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoOverP3Genes.cleaned, value = TRUE)

EmbryoOverP3Genes.cleaned.NoCC <- EmbryoOverP3Genes.cleaned[! EmbryoOverP3Genes.cleaned %in% EmbryoOverP3Genes.CellCycle]

###EmbryoOver24Hr
EmbryoOver24Hr<- FindMarkers(merged.obj, ident.1 = "Embryogenesis", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
EmbryoOver24Hr$ratio <- EmbryoOver24Hr$pct.1/EmbryoOver24Hr$pct.2
EmbryoOver24Hr$power <-EmbryoOver24Hr$ratio*EmbryoOver24Hr$avg_log2FC
EmbryoOver24Hr<- EmbryoOver24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
EmbryoOver24Hr <- EmbryoOver24Hr[order(-EmbryoOver24Hr$power), ]

EmbryoOver24HrGenes <- rownames(EmbryoOver24Hr)

EmbryoOver24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoOver24HrGenes, value = TRUE)
EmbryoOver24HrGenes.cleaned <- EmbryoOver24HrGenes[! EmbryoOver24HrGenes %in% EmbryoOver24HrGenes.toRemove]

EmbryoOver24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoOver24HrGenes.cleaned, value = TRUE)

EmbryoOver24HrGenes.cleaned.NoCC <- EmbryoOver24HrGenes.cleaned[! EmbryoOver24HrGenes.cleaned %in% EmbryoOver24HrGenes.CellCycle]

##EmbryoOverStart
EmbryoOverStart<- FindMarkers(merged.obj, ident.1 = "Embryogenesis", ident.2 = "St. Dclk1", only.pos = T, logfc.threshold = 0.25)
EmbryoOverStart$ratio <- EmbryoOverStart$pct.1/EmbryoOverStart$pct.2
EmbryoOverStart$power <-EmbryoOverStart$ratio*EmbryoOverStart$avg_log2FC
EmbryoOverStart<- EmbryoOverStart %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
EmbryoOverStart <- EmbryoOverStart[order(-EmbryoOverStart$power), ]

EmbryoOverStartGenes <- rownames(EmbryoOverStart)

EmbryoOverStartGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoOverStartGenes, value = TRUE)
EmbryoOverStartGenes.cleaned <- EmbryoOverStartGenes[! EmbryoOverStartGenes %in% EmbryoOverStartGenes.toRemove]

EmbryoOverStartGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoOverStartGenes.cleaned, value = TRUE)

EmbryoOverStartGenes.cleaned.NoCC <- EmbryoOverStartGenes.cleaned[! EmbryoOverStartGenes.cleaned %in% EmbryoOverStartGenes.CellCycle]


EmbryoInt1 <- intersect(EmbryoOverFiveDayGenes.cleaned.NoCC, 
                        EmbryoOverP0Genes.cleaned.NoCC)

EmbryoInt2 <- intersect(EmbryoInt1, EmbryoOverP3Genes.cleaned.NoCC)

EmbryoInt3 <- intersect(EmbryoInt2, EmbryoOver24HrGenes.cleaned.NoCC)
EmbryoInt4 <- intersect(EmbryoInt3, EmbryoOverStartGenes.cleaned.NoCC)
#####

#####24 Hr
#####
#24 Hr over Embryo
Eng24HrOverEmbryo<- FindMarkers(merged.obj, ident.1 = "Engineered-24Hr", ident.2 = "Embryogenesis", only.pos = T, logfc.threshold = 0.25)
Eng24HrOverEmbryo$ratio <- Eng24HrOverEmbryo$pct.1/Eng24HrOverEmbryo$pct.2
Eng24HrOverEmbryo$power <-Eng24HrOverEmbryo$ratio*Eng24HrOverEmbryo$avg_log2FC
Eng24HrOverEmbryo<- Eng24HrOverEmbryo %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
Eng24HrOverEmbryo <- Eng24HrOverEmbryo[order(-Eng24HrOverEmbryo$power), ]

Eng24HrOverEmbryoGenes <- rownames(Eng24HrOverEmbryo)

Eng24HrOverEmbryoGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), Eng24HrOverEmbryoGenes, value = TRUE)
Eng24HrOverEmbryoGenes.cleaned <- Eng24HrOverEmbryoGenes[! Eng24HrOverEmbryoGenes %in% Eng24HrOverEmbryoGenes.toRemove]

Eng24HrOverEmbryoGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), Eng24HrOverEmbryoGenes.cleaned, value = TRUE)

Eng24HrOverEmbryoGenes.cleaned.NoCC <- Eng24HrOverEmbryoGenes.cleaned[! Eng24HrOverEmbryoGenes.cleaned %in% Eng24HrOverEmbryoGenes.CellCycle]

#24 Hr Over Start
Eng24HrOverStart<- FindMarkers(merged.obj, ident.1 = "Engineered-24Hr", ident.2 = "St. Dclk1", only.pos = T, logfc.threshold = 0.25)
Eng24HrOverStart$ratio <- Eng24HrOverStart$pct.1/Eng24HrOverStart$pct.2
Eng24HrOverStart$power <-Eng24HrOverStart$ratio*Eng24HrOverStart$avg_log2FC
Eng24HrOverStart<- Eng24HrOverStart %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
Eng24HrOverStart <- Eng24HrOverStart[order(-Eng24HrOverStart$power), ]

Eng24HrOverStartGenes <- rownames(Eng24HrOverStart)

Eng24HrOverStartGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), Eng24HrOverStartGenes, value = TRUE)
Eng24HrOverStartGenes.cleaned <- Eng24HrOverStartGenes[! Eng24HrOverStartGenes %in% Eng24HrOverStartGenes.toRemove]

Eng24HrOverStartGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), Eng24HrOverStartGenes.cleaned, value = TRUE)

Eng24HrOverStartGenes.cleaned.NoCC <- Eng24HrOverStartGenes.cleaned[! Eng24HrOverStartGenes.cleaned %in% Eng24HrOverStartGenes.CellCycle]

#24 Hr Over P0
Eng24HrOverP0<- FindMarkers(merged.obj, ident.1 = "Engineered-24Hr", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
Eng24HrOverP0$ratio <- Eng24HrOverP0$pct.1/Eng24HrOverP0$pct.2
Eng24HrOverP0$power <-Eng24HrOverP0$ratio*Eng24HrOverP0$avg_log2FC
Eng24HrOverP0<- Eng24HrOverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
Eng24HrOverP0 <- Eng24HrOverP0[order(-Eng24HrOverP0$power), ]

Eng24HrOverP0Genes <- rownames(Eng24HrOverP0)

Eng24HrOverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), Eng24HrOverP0Genes, value = TRUE)
Eng24HrOverP0Genes.cleaned <- Eng24HrOverP0Genes[! Eng24HrOverP0Genes %in% Eng24HrOverP0Genes.toRemove]

Eng24HrOverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), Eng24HrOverP0Genes.cleaned, value = TRUE)

Eng24HrOverP0Genes.cleaned.NoCC <- Eng24HrOverP0Genes.cleaned[! Eng24HrOverP0Genes.cleaned %in% Eng24HrOverP0Genes.CellCycle]

#24 Hr Over P3
Eng24HrOverP3<- FindMarkers(merged.obj, ident.1 = "Engineered-24Hr", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
Eng24HrOverP3$ratio <- Eng24HrOverP3$pct.1/Eng24HrOverP3$pct.2
Eng24HrOverP3$power <-Eng24HrOverP3$ratio*Eng24HrOverP3$avg_log2FC
Eng24HrOverP3<- Eng24HrOverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
Eng24HrOverP3 <- Eng24HrOverP3[order(-Eng24HrOverP3$power), ]

Eng24HrOverP3Genes <- rownames(Eng24HrOverP3)

Eng24HrOverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), Eng24HrOverP3Genes, value = TRUE)
Eng24HrOverP3Genes.cleaned <- Eng24HrOverP3Genes[! Eng24HrOverP3Genes %in% Eng24HrOverP3Genes.toRemove]

Eng24HrOverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), Eng24HrOverP3Genes.cleaned, value = TRUE)

Eng24HrOverP3Genes.cleaned.NoCC <- Eng24HrOverP3Genes.cleaned[! Eng24HrOverP3Genes.cleaned %in% Eng24HrOverP3Genes.CellCycle]

#24 Hr Over FiveDay
Eng24HrOver5D<- FindMarkers(merged.obj, ident.1 = "Engineered-24Hr", ident.2 = "Engineered-FiveDay", only.pos = T, logfc.threshold = 0.25)
Eng24HrOver5D$ratio <- Eng24HrOver5D$pct.1/Eng24HrOver5D$pct.2
Eng24HrOver5D$power <-Eng24HrOver5D$ratio*Eng24HrOver5D$avg_log2FC
Eng24HrOver5D<- Eng24HrOver5D %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
Eng24HrOver5D <- Eng24HrOver5D[order(-Eng24HrOver5D$power), ]

Eng24HrOverFiveDayGenes <- rownames(Eng24HrOver5D)

Eng24HrOverFiveDayGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), Eng24HrOverFiveDayGenes, value = TRUE)
Eng24HrOverFiveDayGenes.cleaned <- Eng24HrOverFiveDayGenes[! Eng24HrOverFiveDayGenes %in% Eng24HrOverFiveDayGenes.toRemove]

Eng24HrOverFiveDayGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), Eng24HrOverFiveDayGenes.cleaned, value = TRUE)

Eng24HrOverFiveDayGenes.cleaned.NoCC <- Eng24HrOverFiveDayGenes.cleaned[! Eng24HrOverFiveDayGenes.cleaned %in% Eng24HrOverFiveDayGenes.CellCycle]


Eng24HrInt1 <- intersect(Eng24HrOverFiveDayGenes.cleaned.NoCC, 
                         Eng24HrOverP0Genes.cleaned.NoCC)

Eng24HrInt2 <- intersect(Eng24HrInt1, Eng24HrOverP3Genes.cleaned.NoCC)

Eng24HrInt3 <- intersect(Eng24HrInt2, Eng24HrOverEmbryoGenes.cleaned.NoCC)
Eng24HrInt4 <- intersect(Eng24HrInt3, Eng24HrOverStartGenes.cleaned.NoCC)

#####

#####P0-2D
#####
#P0OVerEmbryo
P0OverEmbryo<- FindMarkers(merged.obj, ident.1 = "P0-2D", ident.2 = "Embryogenesis", only.pos = T, logfc.threshold = 0.25)
P0OverEmbryo$ratio <- P0OverEmbryo$pct.1/P0OverEmbryo$pct.2
P0OverEmbryo$power <-P0OverEmbryo$ratio*P0OverEmbryo$avg_log2FC
P0OverEmbryo<- P0OverEmbryo %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
P0OverEmbryo <- P0OverEmbryo[order(-P0OverEmbryo$power), ]

P0OverEmbryoGenes <- rownames(P0OverEmbryo)

P0OverEmbryoGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), P0OverEmbryoGenes, value = TRUE)
P0OverEmbryoGenes.cleaned <- P0OverEmbryoGenes[! P0OverEmbryoGenes %in% P0OverEmbryoGenes.toRemove]

P0OverEmbryoGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), P0OverEmbryoGenes.cleaned, value = TRUE)

P0OverEmbryoGenes.cleaned.NoCC <- P0OverEmbryoGenes.cleaned[! P0OverEmbryoGenes.cleaned %in% P0OverEmbryoGenes.CellCycle]
#P0OverStart
P0OverStart<- FindMarkers(merged.obj, ident.1 = "P0-2D", ident.2 = "St. Dclk1", only.pos = T, logfc.threshold = 0.25)
P0OverStart$ratio <- P0OverStart$pct.1/P0OverStart$pct.2
P0OverStart$power <-P0OverStart$ratio*P0OverStart$avg_log2FC
P0OverStart<- P0OverStart %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
P0OverStart <- P0OverStart[order(-P0OverStart$power), ]

P0OverStartGenes <- rownames(P0OverStart)

P0OverStartGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), P0OverStartGenes, value = TRUE)
P0OverStartGenes.cleaned <- P0OverStartGenes[! P0OverStartGenes %in% P0OverStartGenes.toRemove]

P0OverStartGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), P0OverStartGenes.cleaned, value = TRUE)

P0OverStartGenes.cleaned.NoCC <- P0OverStartGenes.cleaned[! P0OverStartGenes.cleaned %in% P0OverStartGenes.CellCycle]

#P0OverP3
P0OverP3<- FindMarkers(merged.obj, ident.1 = "P0-2D", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
P0OverP3$ratio <- P0OverP3$pct.1/P0OverP3$pct.2
P0OverP3$power <-P0OverP3$ratio*P0OverP3$avg_log2FC
P0OverP3<- P0OverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
P0OverP3 <- P0OverP3[order(-P0OverP3$power), ]

P0OverP3Genes <- rownames(P0OverP3)

P0OverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), P0OverP3Genes, value = TRUE)
P0OverP3Genes.cleaned <- P0OverP3Genes[! P0OverP3Genes %in% P0OverP3Genes.toRemove]

P0OverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), P0OverP3Genes.cleaned, value = TRUE)

P0OverP3Genes.cleaned.NoCC <- P0OverP3Genes.cleaned[! P0OverP3Genes.cleaned %in% P0OverP3Genes.CellCycle]

#P0Over24Hr
P0Over24Hr<- FindMarkers(merged.obj, ident.1 = "P0-2D", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
P0Over24Hr$ratio <- P0Over24Hr$pct.1/P0Over24Hr$pct.2
P0Over24Hr$power <-P0Over24Hr$ratio*P0Over24Hr$avg_log2FC
P0Over24Hr<- P0Over24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
P0Over24Hr <- P0Over24Hr[order(-P0Over24Hr$power), ]

P0Over24HrGenes <- rownames(P0Over24Hr)

P0Over24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), P0Over24HrGenes, value = TRUE)
P0Over24HrGenes.cleaned <- P0Over24HrGenes[! P0Over24HrGenes %in% P0Over24HrGenes.toRemove]

P0Over24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), P0Over24HrGenes.cleaned, value = TRUE)

P0Over24HrGenes.cleaned.NoCC <- P0Over24HrGenes.cleaned[! P0Over24HrGenes.cleaned %in% P0Over24HrGenes.CellCycle]

#P0OverFiveDay
P0OverFiveDay<- FindMarkers(merged.obj, ident.1 = "P0-2D", ident.2 = "Engineered-FiveDay", only.pos = T, logfc.threshold = 0.25)
P0OverFiveDay$ratio <- P0OverFiveDay$pct.1/P0OverFiveDay$pct.2
P0OverFiveDay$power <-P0OverFiveDay$ratio*P0OverFiveDay$avg_log2FC
P0OverFiveDay<- P0OverFiveDay %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
P0OverFiveDay <- P0OverFiveDay[order(-P0OverFiveDay$power), ]

P0OverFiveDayGenes <- rownames(P0OverFiveDay)

P0OverFiveDayGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), P0OverFiveDayGenes, value = TRUE)
P0OverFiveDayGenes.cleaned <- P0OverFiveDayGenes[! P0OverFiveDayGenes %in% P0OverFiveDayGenes.toRemove]

P0OverFiveDayGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), P0OverFiveDayGenes.cleaned, value = TRUE)

P0OverFiveDayGenes.cleaned.NoCC <- P0OverFiveDayGenes.cleaned[! P0OverFiveDayGenes.cleaned %in% P0OverFiveDayGenes.CellCycle]


P0Int1 <- intersect(P0OverEmbryoGenes.cleaned.NoCC, 
                    P0OverFiveDayGenes.cleaned.NoCC)

P0Int2 <- intersect(P0Int1, P0OverP3Genes.cleaned.NoCC)
P0Int3 <- intersect(P0Int2, P0Over24HrGenes.cleaned.NoCC)
P0Int4 <- intersect(P0Int3, P0OverStartGenes.cleaned.NoCC)
#####

#####P3-2D
#####
#P3OverEmbryo
P3OverEmbryo<- FindMarkers(merged.obj, ident.1 = "P3-2D", ident.2 = "Embryogenesis", only.pos = T, logfc.threshold = 0.25)
P3OverEmbryo$ratio <- P3OverEmbryo$pct.1/P3OverEmbryo$pct.2
P3OverEmbryo$power <-P3OverEmbryo$ratio*P3OverEmbryo$avg_log2FC
P3OverEmbryo<- P3OverEmbryo %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
P3OverEmbryo <- P3OverEmbryo[order(-P3OverEmbryo$power), ]

P3OverEmbryoGenes <- rownames(P3OverEmbryo)

P3OverEmbryoGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), P3OverEmbryoGenes, value = TRUE)
P3OverEmbryoGenes.cleaned <- P3OverEmbryoGenes[! P3OverEmbryoGenes %in% P3OverEmbryoGenes.toRemove]

P3OverEmbryoGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), P3OverEmbryoGenes.cleaned, value = TRUE)

P3OverEmbryoGenes.cleaned.NoCC <- P3OverEmbryoGenes.cleaned[! P3OverEmbryoGenes.cleaned %in% P3OverEmbryoGenes.CellCycle]
#P3OverStart
P3OverStart<- FindMarkers(merged.obj, ident.1 = "P3-2D", ident.2 = "St. Dclk1", only.pos = T, logfc.threshold = 0.25)
P3OverStart$ratio <- P3OverStart$pct.1/P3OverStart$pct.2
P3OverStart$power <-P3OverStart$ratio*P3OverStart$avg_log2FC
P3OverStart<- P3OverStart %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
P3OverStart <- P3OverStart[order(-P3OverStart$power), ]

P3OverStartGenes <- rownames(P3OverStart)

P3OverStartGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), P3OverStartGenes, value = TRUE)
P3OverStartGenes.cleaned <- P3OverStartGenes[! P3OverStartGenes %in% P3OverStartGenes.toRemove]

P3OverStartGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), P3OverStartGenes.cleaned, value = TRUE)

P3OverStartGenes.cleaned.NoCC <- P3OverStartGenes.cleaned[! P3OverStartGenes.cleaned %in% P3OverStartGenes.CellCycle]

#P3OverP0
P3OverP0<- FindMarkers(merged.obj, ident.1 = "P3-2D", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
P3OverP0$ratio <- P3OverP0$pct.1/P3OverP0$pct.2
P3OverP0$power <-P3OverP0$ratio*P3OverP0$avg_log2FC
P3OverP0<- P3OverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
P3OverP0 <- P3OverP0[order(-P3OverP0$power), ]

P3OverP0Genes <- rownames(P3OverP0)

P3OverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), P3OverP0Genes, value = TRUE)
P3OverP0Genes.cleaned <- P3OverP0Genes[! P3OverP0Genes %in% P3OverP0Genes.toRemove]

P3OverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), P3OverP0Genes.cleaned, value = TRUE)

P3OverP0Genes.cleaned.NoCC <- P3OverP0Genes.cleaned[! P3OverP0Genes.cleaned %in% P3OverP0Genes.CellCycle]

#P3Over24Hr
P3Over24Hr<- FindMarkers(merged.obj, ident.1 = "P3-2D", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
P3Over24Hr$ratio <- P3Over24Hr$pct.1/P3Over24Hr$pct.2
P3Over24Hr$power <-P3Over24Hr$ratio*P3Over24Hr$avg_log2FC
P3Over24Hr<- P3Over24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
P3Over24Hr <- P3Over24Hr[order(-P3Over24Hr$power), ]

P3Over24HrGenes <- rownames(P3Over24Hr)

P3Over24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), P3Over24HrGenes, value = TRUE)
P3Over24HrGenes.cleaned <- P3Over24HrGenes[! P3Over24HrGenes %in% P3Over24HrGenes.toRemove]

P3Over24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), P3Over24HrGenes.cleaned, value = TRUE)

P3Over24HrGenes.cleaned.NoCC <- P3Over24HrGenes.cleaned[! P3Over24HrGenes.cleaned %in% P3Over24HrGenes.CellCycle]

#P3OverFiveDay
P3OverFiveDay<- FindMarkers(merged.obj, ident.1 = "P3-2D", ident.2 = "Engineered-FiveDay", only.pos = T, logfc.threshold = 0.25)
P3OverFiveDay$ratio <- P3OverFiveDay$pct.1/P3OverFiveDay$pct.2
P3OverFiveDay$power <-P3OverFiveDay$ratio*P3OverFiveDay$avg_log2FC
P3OverFiveDay<- P3OverFiveDay %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(2000, power)
P3OverFiveDay <- P3OverFiveDay[order(-P3OverFiveDay$power), ]

P3OverFiveDayGenes <- rownames(P3OverFiveDay)

P3OverFiveDayGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), P3OverFiveDayGenes, value = TRUE)
P3OverFiveDayGenes.cleaned <- P3OverFiveDayGenes[! P3OverFiveDayGenes %in% P3OverFiveDayGenes.toRemove]

P3OverFiveDayGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), P3OverFiveDayGenes.cleaned, value = TRUE)

P3OverFiveDayGenes.cleaned.NoCC <- P3OverFiveDayGenes.cleaned[! P3OverFiveDayGenes.cleaned %in% P3OverFiveDayGenes.CellCycle]


P3Int1 <- intersect(P3OverEmbryoGenes.cleaned.NoCC, 
                    P3OverFiveDayGenes.cleaned.NoCC)
P3Int2 <- intersect(P3Int1, P3OverP0Genes.cleaned.NoCC)
P3Int3 <- intersect(P3Int2, P0Over24HrGenes.cleaned.NoCC)
P3Int4 <- intersect(P3Int3, P0OverStartGenes.cleaned.NoCC)
#####


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.1.30/Heatmaps_(2000Feats)")
#####Heatmaps#####
#####Embryo
#####
png(file = 'EmbryoOnly.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = EmbryoInt4[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoInt4[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'EmbryoOnly.2.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = EmbryoInt4[51:100], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoInt4[51:100],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'EmbryoOnly.3.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = EmbryoInt4[101:150], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoInt4[101:150],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'EmbryoOnly.4.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = EmbryoInt4[151:200], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoInt4[151:200],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

#####

#####Start
#####
png(file = 'StartOnly.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = StartInt4[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartInt4[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'StartOnly.2.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = StartInt4[51:100], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartInt4[51:100],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'StartOnly.3.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = StartInt4[101:150], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartInt4[101:150],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'StartOnly.4.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = StartInt4[151:200], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartInt4[151:200],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()
#####

#####P0 
#####
png(file = 'P0Only.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = P0Int4[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=P0Int4[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'P0Only.2.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = P0Int4[51:81], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=P0Int4[51:81],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

# png(file = 'P0Only.3.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
# CustomHeatmapOnly3(object = downsampled,
#                    data.type = 'RNA', # sets the assay to pull data from
#                    primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
#                    secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
#                    #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
#                    #quarternary = 'orig.ident' , # not used in this specific function
#                    primary.cols = col.pal$Condition, # Needs to be a named list of colors
#                    secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
#                    #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
#                    #quarternary.cols = NULL, # not used in this specific function
#                    features = P0Int4[101:150], # defined above, customize
#                    labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
#                    selected.row.anotations=P0Int4[101:150],
#                    selected.label.size = 10,
#                    use.scale.data = T,
#                    range.frac = 0.35)
# dev.off()
# 
# png(file = 'P0Only.4.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
# CustomHeatmapOnly3(object = downsampled,
#                    data.type = 'RNA', # sets the assay to pull data from
#                    primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
#                    secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
#                    #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
#                    #quarternary = 'orig.ident' , # not used in this specific function
#                    primary.cols = col.pal$Condition, # Needs to be a named list of colors
#                    secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
#                    #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
#                    #quarternary.cols = NULL, # not used in this specific function
#                    features = P0Int4[151:200], # defined above, customize
#                    labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
#                    selected.row.anotations=P0Int4[151:200],
#                    selected.label.size = 10,
#                    use.scale.data = T,
#                    range.frac = 0.35)
# dev.off()
#####

#####P3 
#####
png(file = 'P3Only.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = P3Int4[1:6], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=P3Int4[1:6],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

# png(file = 'P3Only.2.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
# CustomHeatmapOnly3(object = downsampled,
#                    data.type = 'RNA', # sets the assay to pull data from
#                    primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
#                    secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
#                    #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
#                    #quarternary = 'orig.ident' , # not used in this specific function
#                    primary.cols = col.pal$Condition, # Needs to be a named list of colors
#                    secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
#                    #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
#                    #quarternary.cols = NULL, # not used in this specific function
#                    features = P3Int2[51:100], # defined above, customize
#                    labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
#                    selected.row.anotations=P3Int2[51:100],
#                    selected.label.size = 10,
#                    use.scale.data = T,
#                    range.frac = 0.35)
# dev.off()
# 
# png(file = 'P3Only.3.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
# CustomHeatmapOnly3(object = downsampled,
#                    data.type = 'RNA', # sets the assay to pull data from
#                    primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
#                    secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
#                    #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
#                    #quarternary = 'orig.ident' , # not used in this specific function
#                    primary.cols = col.pal$Condition, # Needs to be a named list of colors
#                    secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
#                    #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
#                    #quarternary.cols = NULL, # not used in this specific function
#                    features = P3Int2[101:126], # defined above, customize
#                    labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
#                    selected.row.anotations=P3Int2[101:126],
#                    selected.label.size = 10,
#                    use.scale.data = T,
#                    range.frac = 0.35)
# dev.off()
# 
# png(file = 'P3Only.4.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
# CustomHeatmapOnly3(object = downsampled,
#                    data.type = 'RNA', # sets the assay to pull data from
#                    primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
#                    secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
#                    #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
#                    #quarternary = 'orig.ident' , # not used in this specific function
#                    primary.cols = col.pal$Condition, # Needs to be a named list of colors
#                    secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
#                    #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
#                    #quarternary.cols = NULL, # not used in this specific function
#                    features = P3Int4[151:200], # defined above, customize
#                    labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
#                    selected.row.anotations=P3Int4[151:200],
#                    selected.label.size = 10,
#                    use.scale.data = T,
#                    range.frac = 0.35)
# dev.off()

#####

#####24 Hr
#####
png(file = '24HrOnly.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = Eng24HrInt4[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=Eng24HrInt4[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = '24HrOnly.2.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = Eng24HrInt4[51:100], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=Eng24HrInt4[51:100],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = '24HrOnly.3.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = Eng24HrInt4[101:150], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=Eng24HrInt4[101:150],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = '24HrOnly.4.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = Eng24HrInt4[151:179], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=Eng24HrInt4[151:179],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

#####FiveDay
#####
png(file = 'FiveDayOnly.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = FiveDayInt4[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FiveDayInt4[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'FiveDayOnly.2.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = FiveDayInt4[51:100], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FiveDayInt4[51:100],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'FiveDayOnly.3.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = FiveDayInt4[101:150], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FiveDayInt4[101:150],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'FiveDayOnly.4.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = FiveDayInt4[151:200], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FiveDayInt4[151:200],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()