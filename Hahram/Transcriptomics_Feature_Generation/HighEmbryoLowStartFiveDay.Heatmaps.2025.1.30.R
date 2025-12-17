# Features high in  embryogeneis, low in adulthood and five day, re-activate during either P0, P3, or BDL series [‘Regenerative programming that is reduced during ATI-differentiation’]


##EMBRYO AND P0, Over Start, P3, BDL, FiveDay
#####

merged.obj$EmbryoPlusP0 <- merged.obj$Condition

Idents(merged.obj) <- merged.obj$EmbryoPlusP0

merged.obj <- RenameIdents(merged.obj,
                           "Embryogenesis" = "Focus",
                           "P0-2D" = "Focus"
)

merged.obj$EmbryoPlusP0 <- Idents(merged.obj)
Idents(merged.obj) <- merged.obj$EmbryoPlusP0
#OverStart
EmbryoP0FocusOverStart<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "St. Dclk1", only.pos = T, logfc.threshold = 0.25)
EmbryoP0FocusOverStart$ratio <- EmbryoP0FocusOverStart$pct.1/EmbryoP0FocusOverStart$pct.2
EmbryoP0FocusOverStart$power <-EmbryoP0FocusOverStart$ratio*EmbryoP0FocusOverStart$avg_log2FC
EmbryoP0FocusOverStart<- EmbryoP0FocusOverStart %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoP0FocusOverStartGenes <- rownames(EmbryoP0FocusOverStart)
EmbryoP0FocusOverStartGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoP0FocusOverStartGenes, value = TRUE)
EmbryoP0FocusOverStartGenes.cleaned <- EmbryoP0FocusOverStartGenes[! EmbryoP0FocusOverStartGenes %in% EmbryoP0FocusOverStartGenes.toRemove]
EmbryoP0FocusOverStartGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoP0FocusOverStartGenes.cleaned, value = TRUE)
EmbryoP0FocusOverStartGenes.cleaned.NoCC <- EmbryoP0FocusOverStartGenes.cleaned[! EmbryoP0FocusOverStartGenes.cleaned %in% EmbryoP0FocusOverStartGenes.CellCycle]

#OverFiveDay
EmbryoP0FocusOverFiveDay<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "Engineered-FiveDay", only.pos = T, logfc.threshold = 0.25)
EmbryoP0FocusOverFiveDay$ratio <- EmbryoP0FocusOverFiveDay$pct.1/EmbryoP0FocusOverFiveDay$pct.2
EmbryoP0FocusOverFiveDay$power <-EmbryoP0FocusOverFiveDay$ratio*EmbryoP0FocusOverFiveDay$avg_log2FC
EmbryoP0FocusOverFiveDay<- EmbryoP0FocusOverFiveDay %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoP0FocusOverFiveDayGenes <- rownames(EmbryoP0FocusOverFiveDay)
EmbryoP0FocusOverFiveDayGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoP0FocusOverFiveDayGenes, value = TRUE)
EmbryoP0FocusOverFiveDayGenes.cleaned <- EmbryoP0FocusOverFiveDayGenes[! EmbryoP0FocusOverFiveDayGenes %in% EmbryoP0FocusOverFiveDayGenes.toRemove]
EmbryoP0FocusOverFiveDayGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoP0FocusOverFiveDayGenes.cleaned, value = TRUE)
EmbryoP0FocusOverFiveDayGenes.cleaned.NoCC <- EmbryoP0FocusOverFiveDayGenes.cleaned[! EmbryoP0FocusOverFiveDayGenes.cleaned %in% EmbryoP0FocusOverFiveDayGenes.CellCycle]

###OverP3
EmbryoP0FocusOverP3<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
EmbryoP0FocusOverP3$ratio <- EmbryoP0FocusOverP3$pct.1/EmbryoP0FocusOverP3$pct.2
EmbryoP0FocusOverP3$power <-EmbryoP0FocusOverP3$ratio*EmbryoP0FocusOverP3$avg_log2FC
EmbryoP0FocusOverP3<- EmbryoP0FocusOverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoP0FocusOverP3Genes <- rownames(EmbryoP0FocusOverP3)
EmbryoP0FocusOverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoP0FocusOverP3Genes, value = TRUE)
EmbryoP0FocusOverP3Genes.cleaned <- EmbryoP0FocusOverP3Genes[! EmbryoP0FocusOverP3Genes %in% EmbryoP0FocusOverP3Genes.toRemove]
EmbryoP0FocusOverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoP0FocusOverP3Genes.cleaned, value = TRUE)
EmbryoP0FocusOverP3Genes.cleaned.NoCC <- EmbryoP0FocusOverP3Genes.cleaned[! EmbryoP0FocusOverP3Genes.cleaned %in% EmbryoP0FocusOverP3Genes.CellCycle]

###Over24Hr
EmbryoP0FocusOver24Hr<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
EmbryoP0FocusOver24Hr$ratio <- EmbryoP0FocusOver24Hr$pct.1/EmbryoP0FocusOver24Hr$pct.2
EmbryoP0FocusOver24Hr$power <-EmbryoP0FocusOver24Hr$ratio*EmbryoP0FocusOver24Hr$avg_log2FC
EmbryoP0FocusOver24Hr<- EmbryoP0FocusOver24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoP0FocusOver24HrGenes <- rownames(EmbryoP0FocusOver24Hr)
EmbryoP0FocusOver24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoP0FocusOver24HrGenes, value = TRUE)
EmbryoP0FocusOver24HrGenes.cleaned <- EmbryoP0FocusOver24HrGenes[! EmbryoP0FocusOver24HrGenes %in% EmbryoP0FocusOver24HrGenes.toRemove]
EmbryoP0FocusOver24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoP0FocusOver24HrGenes.cleaned, value = TRUE)
EmbryoP0FocusOver24HrGenes.cleaned.NoCC <- EmbryoP0FocusOver24HrGenes.cleaned[! EmbryoP0FocusOver24HrGenes.cleaned %in% EmbryoP0FocusOver24HrGenes.CellCycle]

inst1<- intersect(EmbryoP0FocusOverStartGenes.cleaned.NoCC, EmbryoP0FocusOverFiveDayGenes.cleaned.NoCC)
inst2 <- intersect(inst1, EmbryoP0FocusOverP3Genes.cleaned.NoCC)
inst3 <- intersect(inst2, EmbryoP0FocusOver24HrGenes.cleaned.NoCC)
EmbryoAndP0.OverRest <- inst3
#####

##EMBRYO AND P3, over Start, P0, BDL, FiveDay
#####
merged.obj$EmbryoPlusP3 <- merged.obj$Condition

Idents(merged.obj) <- merged.obj$EmbryoPlusP3

merged.obj <- RenameIdents(merged.obj,
                           "Embryogenesis" = "Focus",
                           "P3-2D" = "Focus"
)

merged.obj$EmbryoPlusP3 <- Idents(merged.obj)


#OverStart
EmbryoP3FocusOverStart<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "St. Dclk1", only.pos = T, logfc.threshold = 0.25)
EmbryoP3FocusOverStart$ratio <- EmbryoP3FocusOverStart$pct.1/EmbryoP3FocusOverStart$pct.2
EmbryoP3FocusOverStart$power <-EmbryoP3FocusOverStart$ratio*EmbryoP3FocusOverStart$avg_log2FC
EmbryoP3FocusOverStart<- EmbryoP3FocusOverStart %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoP3FocusOverStartGenes <- rownames(EmbryoP3FocusOverStart)
EmbryoP3FocusOverStartGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoP3FocusOverStartGenes, value = TRUE)
EmbryoP3FocusOverStartGenes.cleaned <- EmbryoP3FocusOverStartGenes[! EmbryoP3FocusOverStartGenes %in% EmbryoP3FocusOverStartGenes.toRemove]
EmbryoP3FocusOverStartGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoP3FocusOverStartGenes.cleaned, value = TRUE)
EmbryoP3FocusOverStartGenes.cleaned.NoCC <- EmbryoP3FocusOverStartGenes.cleaned[! EmbryoP3FocusOverStartGenes.cleaned %in% EmbryoP3FocusOverStartGenes.CellCycle]

#OverFiveDay
EmbryoP3FocusOverFiveDay<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "Engineered-FiveDay", only.pos = T, logfc.threshold = 0.25)
EmbryoP3FocusOverFiveDay$ratio <- EmbryoP3FocusOverFiveDay$pct.1/EmbryoP3FocusOverFiveDay$pct.2
EmbryoP3FocusOverFiveDay$power <-EmbryoP3FocusOverFiveDay$ratio*EmbryoP3FocusOverFiveDay$avg_log2FC
EmbryoP3FocusOverFiveDay<- EmbryoP3FocusOverFiveDay %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoP3FocusOverFiveDayGenes <- rownames(EmbryoP3FocusOverFiveDay)
EmbryoP3FocusOverFiveDayGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoP3FocusOverFiveDayGenes, value = TRUE)
EmbryoP3FocusOverFiveDayGenes.cleaned <- EmbryoP3FocusOverFiveDayGenes[! EmbryoP3FocusOverFiveDayGenes %in% EmbryoP3FocusOverFiveDayGenes.toRemove]
EmbryoP3FocusOverFiveDayGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoP3FocusOverFiveDayGenes.cleaned, value = TRUE)
EmbryoP3FocusOverFiveDayGenes.cleaned.NoCC <- EmbryoP3FocusOverFiveDayGenes.cleaned[! EmbryoP3FocusOverFiveDayGenes.cleaned %in% EmbryoP3FocusOverFiveDayGenes.CellCycle]

###OverP0
EmbryoP3FocusOverP0<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
EmbryoP3FocusOverP0$ratio <- EmbryoP3FocusOverP0$pct.1/EmbryoP3FocusOverP0$pct.2
EmbryoP3FocusOverP0$power <-EmbryoP3FocusOverP0$ratio*EmbryoP3FocusOverP0$avg_log2FC
EmbryoP3FocusOverP0<- EmbryoP3FocusOverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoP3FocusOverP0Genes <- rownames(EmbryoP3FocusOverP0)
EmbryoP3FocusOverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoP3FocusOverP0Genes, value = TRUE)
EmbryoP3FocusOverP0Genes.cleaned <- EmbryoP3FocusOverP0Genes[! EmbryoP3FocusOverP0Genes %in% EmbryoP3FocusOverP0Genes.toRemove]
EmbryoP3FocusOverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoP3FocusOverP0Genes.cleaned, value = TRUE)
EmbryoP3FocusOverP0Genes.cleaned.NoCC <- EmbryoP3FocusOverP0Genes.cleaned[! EmbryoP3FocusOverP0Genes.cleaned %in% EmbryoP3FocusOverP0Genes.CellCycle]

###Over24Hr
EmbryoP3FocusOver24Hr<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
EmbryoP3FocusOver24Hr$ratio <- EmbryoP3FocusOver24Hr$pct.1/EmbryoP3FocusOver24Hr$pct.2
EmbryoP3FocusOver24Hr$power <-EmbryoP3FocusOver24Hr$ratio*EmbryoP3FocusOver24Hr$avg_log2FC
EmbryoP3FocusOver24Hr<- EmbryoP3FocusOver24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoP3FocusOver24HrGenes <- rownames(EmbryoP3FocusOver24Hr)
EmbryoP3FocusOver24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoP3FocusOver24HrGenes, value = TRUE)
EmbryoP3FocusOver24HrGenes.cleaned <- EmbryoP3FocusOver24HrGenes[! EmbryoP3FocusOver24HrGenes %in% EmbryoP3FocusOver24HrGenes.toRemove]
EmbryoP3FocusOver24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoP3FocusOver24HrGenes.cleaned, value = TRUE)
EmbryoP3FocusOver24HrGenes.cleaned.NoCC <- EmbryoP3FocusOver24HrGenes.cleaned[! EmbryoP3FocusOver24HrGenes.cleaned %in% EmbryoP3FocusOver24HrGenes.CellCycle]

inst1<- intersect(EmbryoP3FocusOverStartGenes.cleaned.NoCC, EmbryoP3FocusOverFiveDayGenes.cleaned.NoCC)
inst2 <- intersect(inst1, EmbryoP3FocusOverP0Genes.cleaned.NoCC)
inst3 <- intersect(inst2, EmbryoP3FocusOver24HrGenes.cleaned.NoCC)
EmbryoAndP3.OverRest <- inst3
#####

##EMBRYO AND BDL, over Start, P0, P3, FiveDaY
#####
merged.obj$EmbryoPlus24Hr <- merged.obj$Condition

Idents(merged.obj) <- merged.obj$EmbryoPlus24Hr

merged.obj <- RenameIdents(merged.obj,
                           "Embryogenesis" = "Focus",
                           "Engineered-24Hr" = "Focus"
)

merged.obj$EmbryoPlus24Hr <- Idents(merged.obj)


#OverStart
Embryo24HrFocusOverStart<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "St. Dclk1", only.pos = T, logfc.threshold = 0.25)
Embryo24HrFocusOverStart$ratio <- Embryo24HrFocusOverStart$pct.1/Embryo24HrFocusOverStart$pct.2
Embryo24HrFocusOverStart$power <-Embryo24HrFocusOverStart$ratio*Embryo24HrFocusOverStart$avg_log2FC
Embryo24HrFocusOverStart<- Embryo24HrFocusOverStart %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
Embryo24HrFocusOverStartGenes <- rownames(Embryo24HrFocusOverStart)
Embryo24HrFocusOverStartGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), Embryo24HrFocusOverStartGenes, value = TRUE)
Embryo24HrFocusOverStartGenes.cleaned <- Embryo24HrFocusOverStartGenes[! Embryo24HrFocusOverStartGenes %in% Embryo24HrFocusOverStartGenes.toRemove]
Embryo24HrFocusOverStartGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), Embryo24HrFocusOverStartGenes.cleaned, value = TRUE)
Embryo24HrFocusOverStartGenes.cleaned.NoCC <- Embryo24HrFocusOverStartGenes.cleaned[! Embryo24HrFocusOverStartGenes.cleaned %in% Embryo24HrFocusOverStartGenes.CellCycle]

#OverFiveDay
Embryo24HrFocusOverFiveDay<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "Engineered-FiveDay", only.pos = T, logfc.threshold = 0.25)
Embryo24HrFocusOverFiveDay$ratio <- Embryo24HrFocusOverFiveDay$pct.1/Embryo24HrFocusOverFiveDay$pct.2
Embryo24HrFocusOverFiveDay$power <-Embryo24HrFocusOverFiveDay$ratio*Embryo24HrFocusOverFiveDay$avg_log2FC
Embryo24HrFocusOverFiveDay<- Embryo24HrFocusOverFiveDay %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
Embryo24HrFocusOverFiveDayGenes <- rownames(Embryo24HrFocusOverFiveDay)
Embryo24HrFocusOverFiveDayGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), Embryo24HrFocusOverFiveDayGenes, value = TRUE)
Embryo24HrFocusOverFiveDayGenes.cleaned <- Embryo24HrFocusOverFiveDayGenes[! Embryo24HrFocusOverFiveDayGenes %in% Embryo24HrFocusOverFiveDayGenes.toRemove]
Embryo24HrFocusOverFiveDayGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), Embryo24HrFocusOverFiveDayGenes.cleaned, value = TRUE)
Embryo24HrFocusOverFiveDayGenes.cleaned.NoCC <- Embryo24HrFocusOverFiveDayGenes.cleaned[! Embryo24HrFocusOverFiveDayGenes.cleaned %in% Embryo24HrFocusOverFiveDayGenes.CellCycle]

###OverP0
Embryo24HrFocusOverP0<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
Embryo24HrFocusOverP0$ratio <- Embryo24HrFocusOverP0$pct.1/Embryo24HrFocusOverP0$pct.2
Embryo24HrFocusOverP0$power <-Embryo24HrFocusOverP0$ratio*Embryo24HrFocusOverP0$avg_log2FC
Embryo24HrFocusOverP0<- Embryo24HrFocusOverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
Embryo24HrFocusOverP0Genes <- rownames(Embryo24HrFocusOverP0)
Embryo24HrFocusOverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), Embryo24HrFocusOverP0Genes, value = TRUE)
Embryo24HrFocusOverP0Genes.cleaned <- Embryo24HrFocusOverP0Genes[! Embryo24HrFocusOverP0Genes %in% Embryo24HrFocusOverP0Genes.toRemove]
Embryo24HrFocusOverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), Embryo24HrFocusOverP0Genes.cleaned, value = TRUE)
Embryo24HrFocusOverP0Genes.cleaned.NoCC <- Embryo24HrFocusOverP0Genes.cleaned[! Embryo24HrFocusOverP0Genes.cleaned %in% Embryo24HrFocusOverP0Genes.CellCycle]

###OverP3
Embryo24HrFocusOverP3<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
Embryo24HrFocusOverP3$ratio <- Embryo24HrFocusOverP3$pct.1/Embryo24HrFocusOverP3$pct.2
Embryo24HrFocusOverP3$power <-Embryo24HrFocusOverP3$ratio*Embryo24HrFocusOverP3$avg_log2FC
Embryo24HrFocusOverP3<- Embryo24HrFocusOverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
Embryo24HrFocusOverP3Genes <- rownames(Embryo24HrFocusOverP3)
Embryo24HrFocusOverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), Embryo24HrFocusOverP3Genes, value = TRUE)
Embryo24HrFocusOverP3Genes.cleaned <- Embryo24HrFocusOverP3Genes[! Embryo24HrFocusOverP3Genes %in% Embryo24HrFocusOverP3Genes.toRemove]
Embryo24HrFocusOverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), Embryo24HrFocusOverP3Genes.cleaned, value = TRUE)
Embryo24HrFocusOverP3Genes.cleaned.NoCC <- Embryo24HrFocusOverP3Genes.cleaned[! Embryo24HrFocusOverP3Genes.cleaned %in% Embryo24HrFocusOverP3Genes.CellCycle]

inst1<- intersect(Embryo24HrFocusOverStartGenes.cleaned.NoCC, Embryo24HrFocusOverFiveDayGenes.cleaned.NoCC)
inst2 <- intersect(inst1, Embryo24HrFocusOverP0Genes.cleaned.NoCC)
inst3 <- intersect(inst2, Embryo24HrFocusOverP3Genes.cleaned.NoCC)
EmbryoAnd24Hr.OverRest <- inst3
#####

##EMBRYO AND FiveDay, over Start, P0, P3, BDL
#####
merged.obj$EmbryoPlusFiveDay <- merged.obj$Condition

Idents(merged.obj) <- merged.obj$EmbryoPlusFiveDay

merged.obj <- RenameIdents(merged.obj,
                           "Embryogenesis" = "Focus",
                           "Engineered-FiveDay" = "Focus"
)

merged.obj$EmbryoPlusFiveDay <- Idents(merged.obj)


#OverStart
EmbryoFiveDayFocusOverStart<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "St. Dclk1", only.pos = T, logfc.threshold = 0.25)
EmbryoFiveDayFocusOverStart$ratio <- EmbryoFiveDayFocusOverStart$pct.1/EmbryoFiveDayFocusOverStart$pct.2
EmbryoFiveDayFocusOverStart$power <-EmbryoFiveDayFocusOverStart$ratio*EmbryoFiveDayFocusOverStart$avg_log2FC
EmbryoFiveDayFocusOverStart<- EmbryoFiveDayFocusOverStart %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoFiveDayFocusOverStartGenes <- rownames(EmbryoFiveDayFocusOverStart)
EmbryoFiveDayFocusOverStartGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoFiveDayFocusOverStartGenes, value = TRUE)
EmbryoFiveDayFocusOverStartGenes.cleaned <- EmbryoFiveDayFocusOverStartGenes[! EmbryoFiveDayFocusOverStartGenes %in% EmbryoFiveDayFocusOverStartGenes.toRemove]
EmbryoFiveDayFocusOverStartGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoFiveDayFocusOverStartGenes.cleaned, value = TRUE)
EmbryoFiveDayFocusOverStartGenes.cleaned.NoCC <- EmbryoFiveDayFocusOverStartGenes.cleaned[! EmbryoFiveDayFocusOverStartGenes.cleaned %in% EmbryoFiveDayFocusOverStartGenes.CellCycle]

#Over24Hr
EmbryoFiveDayFocusOver24Hr<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
EmbryoFiveDayFocusOver24Hr$ratio <- EmbryoFiveDayFocusOver24Hr$pct.1/EmbryoFiveDayFocusOver24Hr$pct.2
EmbryoFiveDayFocusOver24Hr$power <-EmbryoFiveDayFocusOver24Hr$ratio*EmbryoFiveDayFocusOver24Hr$avg_log2FC
EmbryoFiveDayFocusOver24Hr<- EmbryoFiveDayFocusOver24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoFiveDayFocusOver24HrGenes <- rownames(EmbryoFiveDayFocusOver24Hr)
EmbryoFiveDayFocusOver24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoFiveDayFocusOver24HrGenes, value = TRUE)
EmbryoFiveDayFocusOver24HrGenes.cleaned <- EmbryoFiveDayFocusOver24HrGenes[! EmbryoFiveDayFocusOver24HrGenes %in% EmbryoFiveDayFocusOver24HrGenes.toRemove]
EmbryoFiveDayFocusOver24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoFiveDayFocusOver24HrGenes.cleaned, value = TRUE)
EmbryoFiveDayFocusOver24HrGenes.cleaned.NoCC <- EmbryoFiveDayFocusOver24HrGenes.cleaned[! EmbryoFiveDayFocusOver24HrGenes.cleaned %in% EmbryoFiveDayFocusOver24HrGenes.CellCycle]

###OverP0
EmbryoFiveDayFocusOverP0<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
EmbryoFiveDayFocusOverP0$ratio <- EmbryoFiveDayFocusOverP0$pct.1/EmbryoFiveDayFocusOverP0$pct.2
EmbryoFiveDayFocusOverP0$power <-EmbryoFiveDayFocusOverP0$ratio*EmbryoFiveDayFocusOverP0$avg_log2FC
EmbryoFiveDayFocusOverP0<- EmbryoFiveDayFocusOverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoFiveDayFocusOverP0Genes <- rownames(EmbryoFiveDayFocusOverP0)
EmbryoFiveDayFocusOverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoFiveDayFocusOverP0Genes, value = TRUE)
EmbryoFiveDayFocusOverP0Genes.cleaned <- EmbryoFiveDayFocusOverP0Genes[! EmbryoFiveDayFocusOverP0Genes %in% EmbryoFiveDayFocusOverP0Genes.toRemove]
EmbryoFiveDayFocusOverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoFiveDayFocusOverP0Genes.cleaned, value = TRUE)
EmbryoFiveDayFocusOverP0Genes.cleaned.NoCC <- EmbryoFiveDayFocusOverP0Genes.cleaned[! EmbryoFiveDayFocusOverP0Genes.cleaned %in% EmbryoFiveDayFocusOverP0Genes.CellCycle]

###OverP3
EmbryoFiveDayFocusOverP3<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
EmbryoFiveDayFocusOverP3$ratio <- EmbryoFiveDayFocusOverP3$pct.1/EmbryoFiveDayFocusOverP3$pct.2
EmbryoFiveDayFocusOverP3$power <-EmbryoFiveDayFocusOverP3$ratio*EmbryoFiveDayFocusOverP3$avg_log2FC
EmbryoFiveDayFocusOverP3<- EmbryoFiveDayFocusOverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
EmbryoFiveDayFocusOverP3Genes <- rownames(EmbryoFiveDayFocusOverP3)
EmbryoFiveDayFocusOverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), EmbryoFiveDayFocusOverP3Genes, value = TRUE)
EmbryoFiveDayFocusOverP3Genes.cleaned <- EmbryoFiveDayFocusOverP3Genes[! EmbryoFiveDayFocusOverP3Genes %in% EmbryoFiveDayFocusOverP3Genes.toRemove]
EmbryoFiveDayFocusOverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), EmbryoFiveDayFocusOverP3Genes.cleaned, value = TRUE)
EmbryoFiveDayFocusOverP3Genes.cleaned.NoCC <- EmbryoFiveDayFocusOverP3Genes.cleaned[! EmbryoFiveDayFocusOverP3Genes.cleaned %in% EmbryoFiveDayFocusOverP3Genes.CellCycle]

inst1<- intersect(EmbryoFiveDayFocusOverStartGenes.cleaned.NoCC, EmbryoFiveDayFocusOver24HrGenes.cleaned.NoCC)
inst2 <- intersect(inst1, EmbryoFiveDayFocusOverP0Genes.cleaned.NoCC)
inst3 <- intersect(inst2, EmbryoFiveDayFocusOverP3Genes.cleaned.NoCC)
EmbryoAndFiveDay.OverRest <- inst3
#####

#Heatmaps
#####
#Embryo and P0 Over All
png(file = 'EmbryoAndP0OverRest.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = EmbryoAndP0.OverRest[1:61], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoAndP0.OverRest[1:61],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

#Embryo and P3 Over All
png(file = 'EmbryoAndP3OverRest.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = EmbryoAndP3.OverRest[1:38], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoAndP3.OverRest[1:38],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

png(file = 'EmbryoAnd24HrOverRest.Heatmap.1.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = EmbryoAnd24Hr.OverRest[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoAnd24Hr.OverRest[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

png(file = 'EmbryoAnd24HrOverRest.Heatmap.2.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = EmbryoAnd24Hr.OverRest[51:100], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoAnd24Hr.OverRest[51:100],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

png(file = 'EmbryoAnd24HrOverRest.Heatmap.3.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = EmbryoAnd24Hr.OverRest[101:140], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoAnd24Hr.OverRest[101:140],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

#Embryo and FiveDay Over All

png(file = 'EmbryoAndFiveDayOverRest.Heatmap.1.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = EmbryoAndFiveDay.OverRest[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoAndFiveDay.OverRest[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

png(file = 'EmbryoAndFiveDayOverRest.Heatmap.2.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = EmbryoAndFiveDay.OverRest[51:100], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoAndFiveDay.OverRest[51:100],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()
png(file = 'EmbryoAndFiveDayOverRest.Heatmap.3.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = EmbryoAndFiveDay.OverRest[101:160], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoAndFiveDay.OverRest[101:160],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()
