###ATI Intersect Testing
#######
a.intersect.test <- intersect(ATI_genes, ATII_ATI_genes)
#"Ajuba" "Ydjc"  "Hspb2" "Scx"  shared with each other
a.intersect.test <- intersect(ATI_genes, ATII_genes.revised)
#all unique
a.intersect.test <- intersect(ATI_genes, BASC_genes)
#all unique
a.intersect.test <- intersect(ATI_genes, Ciliated_genes)
#all unique
a.intersect.test <- intersect(ATI_genes, Secretory_genes)
#all unique
a.intersect.test <- intersect(ATI_genes, Tuft_genes)
#all unique
a.intersect.test <- intersect(ATI_genes, Hillock_genes)
#all unique
a.intersect.test <- intersect(ATI_genes, Basal_genes)

remove.fromATI <- c("Ajuba", "Ydjc",  "Hspb2", "Scx")
ATI_genes.revised <- setdiff(ATI_genes, remove.fromATI)
#all unique
######

###ATII Intersect Testing
######
a.intersect.test <- intersect(ATII_genes.revised, ATII_ATI_genes)
#all unique
a.intersect.test <- intersect(ATII_genes.revised, ATI_genes)
#all unique
a.intersect.test <- intersect(ATII_genes.revised, BASC_genes)
#all unique
a.intersect.test <- intersect(ATII_genes.revised, Ciliated_genes)
#all unique
a.intersect.test <- intersect(ATII_genes.revised, Secretory_genes)
#all unique
a.intersect.test <- intersect(ATII_genes.revised, Tuft_genes)
#all unique
a.intersect.test <- intersect(ATII_genes.revised, Hillock_genes)
#all unique
a.intersect.test <- intersect(ATII_genes.revised, Basal_genes)
#all unique
######

###BASC Intersect Testing
######
a.intersect.test <- intersect(BASC_genes, ATII_ATI_genes)
#all unique
a.intersect.test <- intersect(BASC_genes, ATI_genes)
#all unique
a.intersect.test <- intersect(BASC_genes, ATII_genes.revised)
#all unique
a.intersect.test <- intersect(BASC_genes, Ciliated_genes)
#all unique
a.intersect.test <- intersect(BASC_genes, Secretory_genes)
#all unique
a.intersect.test <- intersect(BASC_genes, Tuft_genes)
#all unique
a.intersect.test <- intersect(BASC_genes, Hillock_genes)
#all unique
a.intersect.test <- intersect(BASC_genes, Basal_genes)
#all unique
######

###Ciliated Intersect Testing
######
a.intersect.test <- intersect(Ciliated_genes, ATII_ATI_genes)
#all unique
a.intersect.test <- intersect(Ciliated_genes, ATI_genes)
#all unique
a.intersect.test <- intersect(Ciliated_genes, ATII_genes.revised)
#all unique
a.intersect.test <- intersect(Ciliated_genes, BASC_genes)
#all unique
a.intersect.test <- intersect(Ciliated_genes, Secretory_genes)
#all unique
a.intersect.test <- intersect(Ciliated_genes, Tuft_genes)
#all unique
a.intersect.test <- intersect(Ciliated_genes, Hillock_genes)
#all unique
a.intersect.test <- intersect(Ciliated_genes, Basal_genes)
#all unique
######

###Secretory Intersect Testing
######
a.intersect.test <- intersect(Secretory_genes, ATII_ATI_genes)
#all unique
a.intersect.test <- intersect(Secretory_genes, ATI_genes)
#all unique
a.intersect.test <- intersect(Secretory_genes, ATII_genes.revised)
#all unique
a.intersect.test <- intersect(Secretory_genes, BASC_genes)
#all unique
a.intersect.test <- intersect(Secretory_genes, Ciliated_genes)
#all unique
a.intersect.test <- intersect(Secretory_genes, Tuft_genes)
#all unique
a.intersect.test <- intersect(Secretory_genes, Hillock_genes)
#all unique
a.intersect.test <- intersect(Secretory_genes, Basal_genes)
#all unique
######

###Tuft Intersect Testing
######
a.intersect.test <- intersect(Tuft_genes, ATII_ATI_genes)
#all unique
a.intersect.test <- intersect(Tuft_genes, ATI_genes)
#all unique
a.intersect.test <- intersect(Tuft_genes, ATII_genes.revised)
#all unique
a.intersect.test <- intersect(Tuft_genes, BASC_genes)
#all unique
a.intersect.test <- intersect(Tuft_genes, Ciliated_genes)
#all unique
a.intersect.test <- intersect(Tuft_genes, Secretory_genes)
#all unique
a.intersect.test <- intersect(Tuft_genes, Hillock_genes)
#all unique
a.intersect.test <- intersect(Tuft_genes, Basal_genes)
#all unique
######

###Intersect Testing Basal
######
a.intersect.test <- intersect(Basal_genes, ATI_genes)
#all unique
a.intersect.test <- intersect(Basal_genes, ATII_genes.revised)
#all unique
a.intersect.test <- intersect(Basal_genes, ATII_ATI_genes)
#all unique
a.intersect.test <- intersect(Basal_genes, BASC_genes)
#all unique
a.intersect.test <- intersect(Basal_genes, Ciliated_genes)
#all unique
a.intersect.test <- intersect(Basal_genes, Secretory_genes)
#all unique
a.intersect.test <- intersect(Basal_genes, Tuft_genes)
#all unique
a.intersect.test <- intersect(Basal_genes, Hillock_genes)
#all unique
######

###Intersect Testing Hillock
######
a.intersect.test <- intersect(Hillock_genes, ATI_genes)
#Tacstd2
a.intersect.test <- intersect(Hillock_genes, ATII_genes.revised)
#all unique
a.intersect.test <- intersect(Hillock_genes, ATII_ATI_genes)
#all unique
a.intersect.test <- intersect(Hillock_genes, BASC_genes)
#all unique
a.intersect.test <- intersect(Hillock_genes, Ciliated_genes)
#all unique
a.intersect.test <- intersect(Hillock_genes, Secretory_genes)
#all unique
a.intersect.test <- intersect(Hillock_genes, Tuft_genes)
#all unique
a.intersect.test <- intersect(Hillock_genes, Basal_genes)
#all unique
######

VlnPlot(merged.obj.P0.Integ, features = "Tacstd2")

VlnPlot(lung.combined, features = "Tacstd2")
VlnPlot(merged.obj.24Hr.Integ, features = "Tacstd2")

removing.gene <- c("Tacstd2")


ATI_genes <- setdiff(ATI_genes, removing.gene)
a.intersect.test <- intersect(Hillock_genes, ATI_genes)
#all unique



lung.combined <- AddModuleScore(lung.combined,
                                           features = list(ATI_genes),
                                           name="ATI_Module")

lung.combined <- AddModuleScore(lung.combined,
                                      features = list(ATII_genes),
                                      name="ATII_Module")


lung.combined <- AddModuleScore(lung.combined,
                                      features = list(Secretory_genes),
                                      name="Secretory_Module")

lung.combined <- AddModuleScore(lung.combined,
                                      features = list(Ciliated_genes),
                                      name="Ciliated_Module")

lung.combined <- AddModuleScore(lung.combined,
                                      features = list(Tuft_genes),
                                      name="Tuft_Module")


lung.combined <- AddModuleScore(lung.combined,
                                      features = list(BASC_genes),
                                      name="BASC_Module")


lung.combined <- AddModuleScore(lung.combined,
                                      features = list(ATII_ATI_genes),
                                      name="ATII_ATI_Module")



png(file = "ATII.ModuleScore.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
                  features = "ATII_Module1", label = TRUE, repel = TRUE, order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

png(file = "BASC.ModuleScore.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
                 features = "BASC_Module1", label = TRUE, repel = TRUE, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

png(file = "Secretory.ModuleScore.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
                 features = "Secretory_Module1", label = TRUE, repel = TRUE, order = F) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()


png(file = "Ciliated.ModuleScore.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
            features = "Ciliated_Module1", label = TRUE, repel = TRUE, order = F) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

png(file = "ATI.ModuleScore.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
                 features = "ATI_Module1", label = TRUE, repel = TRUE, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()


png(file = "ATII_ATI.ModuleScore.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
            features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()


png(file = "Tuft.ModuleScore.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
            features = "Tuft_Module1", label = TRUE, repel = TRUE, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()




####Overlap Areas:
##Our Secretory Genes bleed into our Ciliated
# "Scgb3a2" "Bpifa5"  "Bpifb1"  "Scgb1a1" "Tbx1"    "Ptpn13"  "Cpe"     "Gstm7"   "Cldn4"   "Pir"    
VlnPlot(lung.combined, features = c("Pir"), idents = c("Ciliated", "Secretory"))
FeaturePlot(lung.combined, features = "Scgb1a1", label = T, pt.size = 2, order = T)



###Secretory Over Ciliated
SecretoryOverCiliated <- FindMarkers(lung.combined, 
                          ident.1 = "Secretory", 
                          ident.2 = "Ciliated", 
                          only.pos = T, 
                          logfc.threshold = 0.1)
SecretoryOverCiliated$ratio <- SecretoryOverCiliated$pct.1/SecretoryOverCiliated$pct.2
SecretoryOverCiliated$power <-SecretoryOverCiliated$ratio*SecretoryOverCiliated$avg_log2FC

SecretoryOverCiliated <- SecretoryOverCiliated %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

SecretoryOverCiliated <- SecretoryOverCiliated %>%  top_n(200, power)
SecretoryOverCiliated <- rownames(SecretoryOverCiliated)


SecretoryOverCiliated.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                       SecretoryOverCiliated, value = TRUE)
SecretoryOverCiliated <- SecretoryOverCiliated[! SecretoryOverCiliated %in% SecretoryOverCiliated.toRemove]
Secretory_genes.2 <- intersect(Secretory_genes, SecretoryOverCiliated)
###Preliminary finding indicates there isnt much we can do about secretory or ciliated without removing almost all our features



###ATII Over ATII_ATI
ATIIOverATII_ATI <- FindMarkers(lung.combined, 
                                     ident.1 = "ATII", 
                                     ident.2 = "ATII-ATI", 
                                     only.pos = T, 
                                     logfc.threshold = 0.1)
ATIIOverATII_ATI$ratio <- ATIIOverATII_ATI$pct.1/ATIIOverATII_ATI$pct.2
ATIIOverATII_ATI$power <-ATIIOverATII_ATI$ratio*ATIIOverATII_ATI$avg_log2FC

ATIIOverATII_ATI <- ATIIOverATII_ATI %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

ATIIOverATII_ATI <- ATIIOverATII_ATI %>%  top_n(200, power)
ATIIOverATII_ATI <- rownames(ATIIOverATII_ATI)


ATIIOverATII_ATI.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                  ATIIOverATII_ATI, value = TRUE)
ATIIOverATII_ATI <- ATIIOverATII_ATI[! ATIIOverATII_ATI %in% ATIIOverATII_ATI.toRemove]
ATII_genes.2 <- intersect(ATII_genes, ATIIOverATII_ATI)


lung.combined <- AddModuleScore(lung.combined,
                                features = list(ATII_genes.2),
                                name="ATII_genes.2_Module")



png(file = "ATII.ModuleScore.Adjust.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
            features = "ATII_genes.2_Module1", label = TRUE, repel = TRUE, order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()


###ATII Over BASC
ATIIOverBASC <- FindMarkers(lung.combined, 
                                ident.1 = "ATII", 
                                ident.2 = "BASC", 
                                only.pos = T, 
                                logfc.threshold = 0.1)
ATIIOverBASC$ratio <- ATIIOverBASC$pct.1/ATIIOverBASC$pct.2
ATIIOverBASC$power <-ATIIOverBASC$ratio*ATIIOverBASC$avg_log2FC

ATIIOverBASC <- ATIIOverBASC %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

ATIIOverBASC <- ATIIOverBASC %>%  top_n(200, power)
ATIIOverBASC <- rownames(ATIIOverBASC)


ATIIOverBASC.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                              ATIIOverBASC, value = TRUE)
ATIIOverBASC <- ATIIOverBASC[! ATIIOverBASC %in% ATIIOverBASC.toRemove]
ATII_genes.3 <- intersect(ATII_genes, ATIIOverBASC)


lung.combined <- AddModuleScore(lung.combined,
                                features = list(ATII_genes.3),
                                name="ATII_genes.3_Module")

png(file = "ATII.ModuleScore.Adjust.2.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
            features = "ATII_genes.3_Module1", label = TRUE, repel = TRUE, order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

ATII_genes.revised <- intersect(ATII_genes.2, ATII_genes.3)

lung.combined <- AddModuleScore(lung.combined,
                                features = list(ATII_genes.revised),
                                name="ATII_genes.revised_Module")

png(file = "ATII.ModuleScore.Adjust.Revised.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
            features = "ATII_genes.revised_Module1", label = TRUE, repel = TRUE, order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

ATII_genes.revised
# "Defb3"   "Lamp3"   "S100g"   "Slc34a2" "Defb4"   "Sftpb"   "Napsa"   "Fabp5" 