
#Finding Basal Markers For ST

ST.downregVsBasal <- c("Goblet", "Tuft", "ATI", "Secretory", 
                  "Ciliated", "ATII", "Cell_Cycle")

for(i in ST.downregVsBasal) {
  BasalOverOther<- FindMarkers(merged.obj.St.Dclk1.Integ, ident.1 = "Basal", 
                                  ident.2 = i, only.pos = T, logfc.threshold = 0.25)
  BasalOverOther$ratio <- BasalOverOther$pct.1/BasalOverOther$pct.2
  BasalOverOther$power <-BasalOverOther$ratio*BasalOverOther$avg_log2FC
  
  BasalOverOther <- BasalOverOther %>% 
    filter (p_val_adj < 0.05) %>% 
    filter (pct.1 >= 0.1) %>%  
    top_n(2000, power)
  
  BasalOverOther <- BasalOverOther[order(-BasalOverOther$power), ]
  
  BasalOverOtherGenes <- rownames(BasalOverOther)
  
  BasalOverOtherGenes.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                          BasalOverOtherGenes, value = TRUE)
  BasalOverOtherGenes.cleaned <- BasalOverOtherGenes[! BasalOverOtherGenes %in% BasalOverOtherGenes.toRemove]
  
  BasalOverOtherGenes.CellCycle <- grep(paste(c("Cc", "Cen", "Cdc", "Cep", "Kif", "Kntc", "Ndc", "Tubb", "Top2a"), 
                                                 collapse = "|"), BasalOverOtherGenes.cleaned, value = TRUE)
  
  BasalOverOtherGenes.cleaned.NoCC <- BasalOverOtherGenes.cleaned[! BasalOverOtherGenes.cleaned %in% BasalOverOtherGenes.CellCycle]
  
  assign(paste0("BasalOver", i), BasalOverOther)
  assign(paste0("BasalOver", i, "cleaned.NoCC"), BasalOverOtherGenes.cleaned.NoCC)
  rm(BasalOverOtherGenes.cleaned.NoCC)
}
###########

FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Ecm1", reduction = "umap.test")
########Basal
Starting.Basal_Features <- intersect(intersect(intersect(intersect(intersect(intersect(BasalOverATIcleaned.NoCC, BasalOverATIIcleaned.NoCC), BasalOverCell_Cyclecleaned.NoCC), BasalOverCiliatedcleaned.NoCC), BasalOverGobletcleaned.NoCC), BasalOverSecretorycleaned.NoCC), BasalOverTuftcleaned.NoCC) 

P0.Basal_Features <- intersect(intersect(intersect(intersect(intersect(BasalOverAT1_Likecleaned.NoCC, BasalOverCiliatedcleaned.NoCC), BasalOverCycling_Epithelialcleaned.NoCC), BasalOverHillock_Basalcleaned.NoCC), BasalOverSecretorycleaned.NoCC), BasalOverCycling_Basalcleaned.NoCC)

P3.Basal_Features <- intersect(intersect(intersect(BasalOverAT1_Likecleaned.NoCC,BasalOverCycling_Epithelialcleaned.NoCC), BasalOverHillock_Basalcleaned.NoCC), BasalOverSecretorycleaned.NoCC)

Eng24Hr.Basal_Features <- intersect(intersect(BasalOverAT1_Likecleaned.NoCC, BasalOverHillock_Likecleaned.NoCC), BasalOverSecretorycleaned.NoCC)

#########Hillock
P0.Hillock_Features <- intersect(intersect(intersect(intersect(intersect(HillockOverAT1_Likecleaned.NoCC, HillockOverCiliatedcleaned.NoCC), HillockOverCycling_Epithelialcleaned.NoCC), HillockOverBasal_Likecleaned.NoCC), HillockOverSecretorycleaned.NoCC), HillockOverCycling_Basalcleaned.NoCC)

P3.Hillock_Features <- intersect(intersect(intersect(intersect(HillockOverAT1_Likecleaned.NoCC,HillockOverCycling_Epithelialcleaned.NoCC), HillockOverSecretorycleaned.NoCC), HillockOverCycling_Basalcleaned.NoCC), HillockOverBasal_Likecleaned.NoCC)

Eng24Hr.Hillock_Features <- intersect(intersect(intersect(HillockOverAT1_Likecleaned.NoCC, HillockOverBasal_Likecleaned.NoCC), HillockOverSecretorycleaned.NoCC), HillockOverCycling_Basalcleaned.NoCC)

FiveDay.Hillock_Features <- intersect(HillockOverCycling_Epithelialcleaned.NoCC, HillockOverTransitionalcleaned.NoCC)
