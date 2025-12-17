####What are the other heatmaps that Sam wants to see?


###Low during embryogenesis, up during adult, down during P0 to 24Hr, and then up in 5D
### How do we do this?
setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BDL_DoOver")

Idents(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate) <- Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate$Condition
Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate <- RenameIdents(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate,
                           "Engineered-3D" = "Engineered-24Hr",
                           "Engineered-P0" = "P0-2D",
                           "Engineered-P3" = "P3-2D"
)

Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate$Condition <- Idents(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate)

save(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, file = "Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate.Robj")
merged.obj <- Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate

###Two approaches: either we try the "Othermetada approach, or the iterative intersection

merged.obj$StartAnd5D <- merged.obj$Condition
Idents(merged.obj) <- merged.obj$StartAnd5D

merged.obj <- RenameIdents( 
                          merged.obj,
                          "Embryogenesis" = "Other",
                          "P0-2D" = "Other",
                          "P3-2D" = "Other",
                          "Engineered-24Hr" = "Other"
  
)
merged.obj$StartAnd5D<- Idents(merged.obj)
Idents(merged.obj) <- merged.obj$StartAnd5D
merged.obj <- JoinLayers(merged.obj)

FiveDay<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "Other", only.pos = T, logfc.threshold = 0.25)
FiveDay$ratio <- FiveDay$pct.1/FiveDay$pct.2
FiveDay$power <-FiveDay$ratio*FiveDay$avg_log2FC
FiveDay<- FiveDay %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
FiveDayGenes <- rownames(FiveDay)

FiveDayGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayGenes, value = TRUE)
FiveDayGenes.cleaned <- FiveDayGenes[! FiveDayGenes %in% FiveDayGenes.toRemove]

FiveDayGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayGenes.cleaned, value = TRUE)

FiveDayGenes.cleaned.NoCC <- FiveDayGenes.cleaned[! FiveDayGenes.cleaned %in% FiveDayGenes.CellCycle]


Starting<- FindMarkers(merged.obj, ident.1 = "St. Dclk1", ident.2 = "Other", only.pos = T, logfc.threshold = 0.25)
Starting$ratio <- Starting$pct.1/Starting$pct.2
Starting$power <-Starting$ratio*Starting$avg_log2FC
Starting<- Starting %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
StartingGenes <- rownames(Starting)

StartingGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartingGenes, value = TRUE)
StartingGenes.cleaned <- StartingGenes[! StartingGenes %in% StartingGenes.toRemove]

StartingGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartingGenes.cleaned, value = TRUE)

StartingGenes.cleaned.NoCC <- StartingGenes.cleaned[! StartingGenes.cleaned %in% StartingGenes.CellCycle]


FiveDay.Start<- intersect(StartingGenes.cleaned.NoCC, FiveDayGenes.cleaned.NoCC)

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
                   features = FiveDay.Start[1:48], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FiveDay.Start[1:48],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.2)

###In theory the above should have given us everything that is in 5day and starting but not embryo, p0, or p3

###The other way to do this would be to create multiple combos and use iterative  intersect()

Idents(merged.obj) <- merged.obj$Condition
FiveDay<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "Other", only.pos = T, logfc.threshold = 0.25)
FiveDay$ratio <- FiveDay$pct.1/FiveDay$pct.2
FiveDay$power <-FiveDay$ratio*FiveDay$avg_log2FC
FiveDay<- FiveDay %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
FiveDayGenes <- rownames(FiveDay)

FiveDayGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayGenes, value = TRUE)
FiveDayGenes.cleaned <- FiveDayGenes[! FiveDayGenes %in% FiveDayGenes.toRemove]

FiveDayGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayGenes.cleaned, value = TRUE)

FiveDayGenes.cleaned.NoCC <- FiveDayGenes.cleaned[! FiveDayGenes.cleaned %in% FiveDayGenes.CellCycle]

StartingOverEmbryo<- FindMarkers(merged.obj, ident.1 = "St. Dclk1", ident.2 = "Embryogenesis", only.pos = T, logfc.threshold = 0.25)
StartingOverEmbryo$ratio <- StartingOverEmbryo$pct.1/StartingOverEmbryo$pct.2
StartingOverEmbryo$power <-StartingOverEmbryo$ratio*StartingOverEmbryo$avg_log2FC
StartingOverEmbryo<- StartingOverEmbryo %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
StartingOverEmbryoGenes <- rownames(StartingOverEmbryo)

StartingOverEmbryoGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartingOverEmbryoGenes, value = TRUE)
StartingOverEmbryoGenes.cleaned <- StartingOverEmbryoGenes[! StartingOverEmbryoGenes %in% StartingOverEmbryoGenes.toRemove]

StartingOverEmbryoGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartingOverEmbryoGenes.cleaned, value = TRUE)

StartingOverEmbryoGenes.cleaned.NoCC <- StartingOverEmbryoGenes.cleaned[! StartingOverEmbryoGenes.cleaned %in% StartingOverEmbryoGenes.CellCycle]

###StartOverP0
StartingOverP0<- FindMarkers(merged.obj, ident.1 = "St. Dclk1", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
StartingOverP0$ratio <- StartingOverP0$pct.1/StartingOverP0$pct.2
StartingOverP0$power <-StartingOverP0$ratio*StartingOverP0$avg_log2FC
StartingOverP0<- StartingOverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
StartingOverP0Genes <- rownames(StartingOverP0)

StartingOverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartingOverP0Genes, value = TRUE)
StartingOverP0Genes.cleaned <- StartingOverP0Genes[! StartingOverP0Genes %in% StartingOverP0Genes.toRemove]

StartingOverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartingOverP0Genes.cleaned, value = TRUE)

StartingOverP0Genes.cleaned.NoCC <- StartingOverP0Genes.cleaned[! StartingOverP0Genes.cleaned %in% StartingOverP0Genes.CellCycle]

###StartOverP3
StartingOverP3<- FindMarkers(merged.obj, ident.1 = "St. Dclk1", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
StartingOverP3$ratio <- StartingOverP3$pct.1/StartingOverP3$pct.2
StartingOverP3$power <-StartingOverP3$ratio*StartingOverP3$avg_log2FC
StartingOverP3<- StartingOverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
StartingOverP3Genes <- rownames(StartingOverP3)

StartingOverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartingOverP3Genes, value = TRUE)
StartingOverP3Genes.cleaned <- StartingOverP3Genes[! StartingOverP3Genes %in% StartingOverP3Genes.toRemove]

StartingOverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartingOverP3Genes.cleaned, value = TRUE)

StartingOverP3Genes.cleaned.NoCC <- StartingOverP3Genes.cleaned[! StartingOverP3Genes.cleaned %in% StartingOverP3Genes.CellCycle]

###StartOver24Hr
StartingOver24Hr<- FindMarkers(merged.obj, ident.1 = "St. Dclk1", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
StartingOver24Hr$ratio <- StartingOver24Hr$pct.1/StartingOver24Hr$pct.2
StartingOver24Hr$power <-StartingOver24Hr$ratio*StartingOver24Hr$avg_log2FC
StartingOver24Hr<- StartingOver24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
StartingOver24HrGenes <- rownames(StartingOver24Hr)

StartingOver24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartingOver24HrGenes, value = TRUE)
StartingOver24HrGenes.cleaned <- StartingOver24HrGenes[! StartingOver24HrGenes %in% StartingOver24HrGenes.toRemove]

StartingOver24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartingOver24HrGenes.cleaned, value = TRUE)

StartingOver24HrGenes.cleaned.NoCC <- StartingOver24HrGenes.cleaned[! StartingOver24HrGenes.cleaned %in% StartingOver24HrGenes.CellCycle]


intersect1 <- intersect(StartingOverEmbryoGenes.cleaned.NoCC, 
                        StartingOverP0Genes.cleaned.NoCC)

intersect2 <- intersect(intersect1, StartingOverP3Genes.cleaned.NoCC)

intersect3 <- intersect(intersect2, StartingOver24HrGenes.cleaned.NoCC)

FiveDayOverEmbryo<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "Embryogenesis", only.pos = T, logfc.threshold = 0.25)
FiveDayOverEmbryo$ratio <- FiveDayOverEmbryo$pct.1/FiveDayOverEmbryo$pct.2
FiveDayOverEmbryo$power <-FiveDayOverEmbryo$ratio*FiveDayOverEmbryo$avg_log2FC
FiveDayOverEmbryo<- FiveDayOverEmbryo %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
FiveDayOverEmbryoGenes <- rownames(FiveDayOverEmbryo)

FiveDayOverEmbryoGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayOverEmbryoGenes, value = TRUE)
FiveDayOverEmbryoGenes.cleaned <- FiveDayOverEmbryoGenes[! FiveDayOverEmbryoGenes %in% FiveDayOverEmbryoGenes.toRemove]

FiveDayOverEmbryoGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayOverEmbryoGenes.cleaned, value = TRUE)

FiveDayOverEmbryoGenes.cleaned.NoCC <- FiveDayOverEmbryoGenes.cleaned[! FiveDayOverEmbryoGenes.cleaned %in% FiveDayOverEmbryoGenes.CellCycle]

###FiveDayOverP0
FiveDayOverP0<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
FiveDayOverP0$ratio <- FiveDayOverP0$pct.1/FiveDayOverP0$pct.2
FiveDayOverP0$power <-FiveDayOverP0$ratio*FiveDayOverP0$avg_log2FC
FiveDayOverP0<- FiveDayOverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
FiveDayOverP0Genes <- rownames(FiveDayOverP0)

FiveDayOverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayOverP0Genes, value = TRUE)
FiveDayOverP0Genes.cleaned <- FiveDayOverP0Genes[! FiveDayOverP0Genes %in% FiveDayOverP0Genes.toRemove]

FiveDayOverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayOverP0Genes.cleaned, value = TRUE)

FiveDayOverP0Genes.cleaned.NoCC <- FiveDayOverP0Genes.cleaned[! FiveDayOverP0Genes.cleaned %in% FiveDayOverP0Genes.CellCycle]

###FiveDayOverP3
FiveDayOverP3<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
FiveDayOverP3$ratio <- FiveDayOverP3$pct.1/FiveDayOverP3$pct.2
FiveDayOverP3$power <-FiveDayOverP3$ratio*FiveDayOverP3$avg_log2FC
FiveDayOverP3<- FiveDayOverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
FiveDayOverP3Genes <- rownames(FiveDayOverP3)

FiveDayOverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayOverP3Genes, value = TRUE)
FiveDayOverP3Genes.cleaned <- FiveDayOverP3Genes[! FiveDayOverP3Genes %in% FiveDayOverP3Genes.toRemove]

FiveDayOverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayOverP3Genes.cleaned, value = TRUE)

FiveDayOverP3Genes.cleaned.NoCC <- FiveDayOverP3Genes.cleaned[! FiveDayOverP3Genes.cleaned %in% FiveDayOverP3Genes.CellCycle]

###FiveDayOver24Hr
FiveDayOver24Hr<- FindMarkers(merged.obj, ident.1 = "Engineered-FiveDay", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
FiveDayOver24Hr$ratio <- FiveDayOver24Hr$pct.1/FiveDayOver24Hr$pct.2
FiveDayOver24Hr$power <-FiveDayOver24Hr$ratio*FiveDayOver24Hr$avg_log2FC
FiveDayOver24Hr<- FiveDayOver24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
FiveDayOver24HrGenes <- rownames(FiveDayOver24Hr)

FiveDayOver24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), FiveDayOver24HrGenes, value = TRUE)
FiveDayOver24HrGenes.cleaned <- FiveDayOver24HrGenes[! FiveDayOver24HrGenes %in% FiveDayOver24HrGenes.toRemove]

FiveDayOver24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), FiveDayOver24HrGenes.cleaned, value = TRUE)

FiveDayOver24HrGenes.cleaned.NoCC <- FiveDayOver24HrGenes.cleaned[! FiveDayOver24HrGenes.cleaned %in% FiveDayOver24HrGenes.CellCycle]


intersectA <- intersect(FiveDayOverEmbryoGenes.cleaned.NoCC, 
                        FiveDayOverP0Genes.cleaned.NoCC)

intersectB <- intersect(intersectA, FiveDayOverP3Genes.cleaned.NoCC)

intersectC <- intersect(intersectB, FiveDayOver24HrGenes.cleaned.NoCC)

FiveDayStart <- intersect(intersectC, intersect3)
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
                   features = intersect1[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=intersect1[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)


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
                   features = testaaGenes.cleaned.NoCC[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=testaaGenes.cleaned.NoCC[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.45)
