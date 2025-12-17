#####From what we have seen in our attempt to generate very cluster specific features, there are 
#####Certain conditions where feature overlap cannot be avoided such as the later 



# Features that are high during embryogeneis, low during adulthood, re-activate during either P0, P3, BDL series, or 5-day series. [‘Regenerative programming’]
# Features that are high during embryogeneis, low during adulthood, re-activate during either P0, P3, or BDL series, and then de-activate during 5-day series. [‘Regenerative programming that is reduced during ATI-differentiation’]
# Features that are low during embryogenesis, high during adulthood, de-activate during Po, P3, or BDL series, and then re-activate during 5-day series. [‘Adult differentiation programs induced during ATI-differentiation’]
# ATI-associated genes / gene module [ask Satoshi for this]

###Lets see if we can bring about features that are high during start and five day


###combined attempt

###Five Day and Starting

merged.obj$StartPlusFiveDay <- merged.obj$Condition

Idents(merged.obj) <- merged.obj$StartPlusFiveDay

merged.obj <- RenameIdents(merged.obj,
                           "St. Dclk1" = "Focus",
                           "Engineered-FiveDay" = "Focus"
                           )

merged.obj$StartPlusFiveDay <- Idents(merged.obj)

Idents(merged.obj) <- merged.obj$StartPlusFiveDay
StartFiveDayFocusOverEmbryo<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "Embryogenesis", only.pos = T, logfc.threshold = 0.25)
StartFiveDayFocusOverEmbryo$ratio <- StartFiveDayFocusOverEmbryo$pct.1/StartFiveDayFocusOverEmbryo$pct.2
StartFiveDayFocusOverEmbryo$power <-StartFiveDayFocusOverEmbryo$ratio*StartFiveDayFocusOverEmbryo$avg_log2FC
StartFiveDayFocusOverEmbryo<- StartFiveDayFocusOverEmbryo %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
StartFiveDayFocusOverEmbryoGenes <- rownames(StartFiveDayFocusOverEmbryo)

StartFiveDayFocusOverEmbryoGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartFiveDayFocusOverEmbryoGenes, value = TRUE)
StartFiveDayFocusOverEmbryoGenes.cleaned <- StartFiveDayFocusOverEmbryoGenes[! StartFiveDayFocusOverEmbryoGenes %in% StartFiveDayFocusOverEmbryoGenes.toRemove]

StartFiveDayFocusOverEmbryoGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartFiveDayFocusOverEmbryoGenes.cleaned, value = TRUE)

StartFiveDayFocusOverEmbryoGenes.cleaned.NoCC <- StartFiveDayFocusOverEmbryoGenes.cleaned[! StartFiveDayFocusOverEmbryoGenes.cleaned %in% StartFiveDayFocusOverEmbryoGenes.CellCycle]

###StartFiveDayFocusOverP0
StartFiveDayFocusOverP0<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "P0-2D", only.pos = T, logfc.threshold = 0.25)
StartFiveDayFocusOverP0$ratio <- StartFiveDayFocusOverP0$pct.1/StartFiveDayFocusOverP0$pct.2
StartFiveDayFocusOverP0$power <-StartFiveDayFocusOverP0$ratio*StartFiveDayFocusOverP0$avg_log2FC
StartFiveDayFocusOverP0<- StartFiveDayFocusOverP0 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
StartFiveDayFocusOverP0Genes <- rownames(StartFiveDayFocusOverP0)

StartFiveDayFocusOverP0Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartFiveDayFocusOverP0Genes, value = TRUE)
StartFiveDayFocusOverP0Genes.cleaned <- StartFiveDayFocusOverP0Genes[! StartFiveDayFocusOverP0Genes %in% StartFiveDayFocusOverP0Genes.toRemove]

StartFiveDayFocusOverP0Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartFiveDayFocusOverP0Genes.cleaned, value = TRUE)

StartFiveDayFocusOverP0Genes.cleaned.NoCC <- StartFiveDayFocusOverP0Genes.cleaned[! StartFiveDayFocusOverP0Genes.cleaned %in% StartFiveDayFocusOverP0Genes.CellCycle]

###StartFiveDayFocusOverP3
StartFiveDayFocusOverP3<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "P3-2D", only.pos = T, logfc.threshold = 0.25)
StartFiveDayFocusOverP3$ratio <- StartFiveDayFocusOverP3$pct.1/StartFiveDayFocusOverP3$pct.2
StartFiveDayFocusOverP3$power <-StartFiveDayFocusOverP3$ratio*StartFiveDayFocusOverP3$avg_log2FC
StartFiveDayFocusOverP3<- StartFiveDayFocusOverP3 %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
StartFiveDayFocusOverP3Genes <- rownames(StartFiveDayFocusOverP3)

StartFiveDayFocusOverP3Genes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartFiveDayFocusOverP3Genes, value = TRUE)
StartFiveDayFocusOverP3Genes.cleaned <- StartFiveDayFocusOverP3Genes[! StartFiveDayFocusOverP3Genes %in% StartFiveDayFocusOverP3Genes.toRemove]

StartFiveDayFocusOverP3Genes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartFiveDayFocusOverP3Genes.cleaned, value = TRUE)

StartFiveDayFocusOverP3Genes.cleaned.NoCC <- StartFiveDayFocusOverP3Genes.cleaned[! StartFiveDayFocusOverP3Genes.cleaned %in% StartFiveDayFocusOverP3Genes.CellCycle]

###StartFiveDayFocusOver24Hr
StartFiveDayFocusOver24Hr<- FindMarkers(merged.obj, ident.1 = "Focus", ident.2 = "Engineered-24Hr", only.pos = T, logfc.threshold = 0.25)
StartFiveDayFocusOver24Hr$ratio <- StartFiveDayFocusOver24Hr$pct.1/StartFiveDayFocusOver24Hr$pct.2
StartFiveDayFocusOver24Hr$power <-StartFiveDayFocusOver24Hr$ratio*StartFiveDayFocusOver24Hr$avg_log2FC
StartFiveDayFocusOver24Hr<- StartFiveDayFocusOver24Hr %>% filter (p_val_adj < 0.05) %>% filter (pct.1 >= 0.1) %>%  top_n(900, power)
StartFiveDayFocusOver24HrGenes <- rownames(StartFiveDayFocusOver24Hr)

StartFiveDayFocusOver24HrGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), StartFiveDayFocusOver24HrGenes, value = TRUE)
StartFiveDayFocusOver24HrGenes.cleaned <- StartFiveDayFocusOver24HrGenes[! StartFiveDayFocusOver24HrGenes %in% StartFiveDayFocusOver24HrGenes.toRemove]

StartFiveDayFocusOver24HrGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), collapse = "|"), StartFiveDayFocusOver24HrGenes.cleaned, value = TRUE)

StartFiveDayFocusOver24HrGenes.cleaned.NoCC <- StartFiveDayFocusOver24HrGenes.cleaned[! StartFiveDayFocusOver24HrGenes.cleaned %in% StartFiveDayFocusOver24HrGenes.CellCycle]


StartFiveDayInt1 <- intersect(StartFiveDayFocusOverEmbryoGenes.cleaned.NoCC, 
                         StartFiveDayFocusOverP0Genes.cleaned.NoCC)
StartFiveDayInt2 <- intersect(StartFiveDayInt1, StartFiveDayFocusOverP3Genes.cleaned.NoCC)
StartFiveDayInt3 <- intersect(StartFiveDayInt2, StartFiveDayFocusOver24HrGenes.cleaned.NoCC)

png(file = 'StartFiveDayFocusTest.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = StartFiveDayInt3[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartFiveDayInt3[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

###List of genes
StartFiveDayFocusOverEmbryoGenes.cleaned.NoCC
StartFiveDayFocusOverP3Genes.cleaned.NoCC
StartFiveDayFocusOver24HrGenes.cleaned.NoCC
StartFiveDayFocusOverP0Genes.cleaned.NoCC

#####Features that are upregulated in Starting and Fiveday, Downregulated In Embryo and P0 
StartFiveDayOverEmbryoP0<- intersect(StartFiveDayFocusOverEmbryoGenes.cleaned.NoCC, StartFiveDayFocusOverP0Genes.cleaned.NoCC)
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
                   features = StartFiveDayOverEmbryoP0[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartFiveDayOverEmbryoP0[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
StartFiveDayOverEmbryoP0.v2 <- setdiff(StartFiveDayOverEmbryoP0, EmbryoInt4)
StartFiveDayOverEmbryoP0.v3 <- setdiff(StartFiveDayOverEmbryoP0, P0Int4)

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
                   features = StartFiveDayOverEmbryoP0.v2[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartFiveDayOverEmbryoP0.v2[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)

StartFiveDayOverEmbryoP0.v2 <- StartFiveDayOverEmbryoP0[! StartFiveDayOverEmbryoP0 %in% EmbryoInt4]

#####Features that are upregulated in Starting and Fiveday, Downregulated In Embryo and P3 
StartFiveDayOverEmbryoP3 <- intersect(StartFiveDayFocusOverEmbryoGenes.cleaned.NoCC, StartFiveDayFocusOverP3Genes.cleaned.NoCC)
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
                   features = StartFiveDayOverEmbryoP3[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartFiveDayOverEmbryoP3[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)


StartFiveDayOverEmbryoP3.v2 <- setdiff(StartFiveDayOverEmbryoP3, EmbryoInt4)


#####Features that are upregulated in Starting and Fiveday, Downregulated In Embryo and 24Hr
StartFiveDayOverEmbryo24Hr <- intersect(StartFiveDayFocusOverEmbryoGenes.cleaned.NoCC, StartFiveDayFocusOver24HrGenes.cleaned.NoCC)
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
                   features = StartFiveDayOverEmbryo24Hr[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartFiveDayOverEmbryo24Hr[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.4)


# Features high in  embryogeneis, low in adulthood and five day, re-activate during either P0, P3, or BDL series, and then de-activate during 5-day series. [‘Regenerative programming that is reduced during ATI-differentiation’]

