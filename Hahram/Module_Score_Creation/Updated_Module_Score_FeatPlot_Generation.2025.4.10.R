##Running setdiff tests

setdiff_test <- setdiff(ATI_genes, ATII_genes)
a.setdiff_test <- intersect(ATI_genes, ATII_genes)
a.setdiff_test <- intersect(ATI_genes, Tuft_genes)
#ALL unique confirmed ATI

a.setdiff_test <- intersect(Ciliated_genes, ATI_genes)
#ALL unique confirmed Ciliated

a.setdiff_test <- intersect(Secretory_genes, ATI_genes)
#ALL unique confirmed 

a.setdiff_test <- intersect(ATII_ATI_genes, ATI_genes)
#Only ATII_ATI shared, but to be expected

a.setdiff_test <- intersect(Tuft_genes, BASC_genes)
##Tuft and BASC share 23 genes. May give them to Tuft

remove.from.BASC <- intersect(BASC_genes, Tuft_genes)

BASC_genes <- BASC_genes[! BASC_genes %in% remove.from.BASC]


###Now for the actual Module score
Basal_Genes<- Basal_Module_Score
Hillock_Genes <- Hillock_Module_Score

DimPlot(merged.obj.P0.Integ, reduction = "umap.test.2")
DimPlot(merged.obj.P3.Integ, reduction = "umap.test")
DimPlot(merged.obj.24Hr.Integ, reduction = "umap.test")
DimPlot(merged.obj.FiveDay.Integ, reduction = "umap.test")


##ATI_Module
#####
merged.obj.P0.Integ <- AddModuleScore(merged.obj.P0.Integ,
                                      features = list(ATI_genes),
                                      name="ATI_Module")
merged.obj.P3.Integ <- AddModuleScore(merged.obj.P3.Integ,
                                      features = list(ATI_genes),
                                      name="ATI_Module")

merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                        features = list(ATI_genes),
                                        name="ATI_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                           features = list(ATI_genes),
                                           name="ATI_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                           features = list(ATI_genes),
                                           name="ATI_Module")

png(file = "AT1.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.St.Dclk1.Integ,
                features = "ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2<-FeaturePlot(merged.obj.P0.Integ,
                features = "ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<-FeaturePlot(merged.obj.P3.Integ,
                features = "ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<-FeaturePlot(merged.obj.24Hr.Integ,
                features = "ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a5<-FeaturePlot(merged.obj.FiveDay.Integ,
                features = "ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



cowplot::plot_grid(a1, a2, a3, a4, a5, nrow = 2)
dev.off()
#####

##ATII_Module

merged.obj.P0.Integ <- AddModuleScore(merged.obj.P0.Integ,
                                      features = list(ATII_genes),
                                      name="ATII_Module")

merged.obj.P3.Integ <- AddModuleScore(merged.obj.P3.Integ,
                                      features = list(ATII_genes),
                                      name="ATII_Module")

merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                        features = list(ATII_genes),
                                        name="ATII_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                           features = list(ATII_genes),
                                           name="ATII_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(ATII_genes),
                                            name="ATII_Module")

png(file = "AT2.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")

a1<-FeaturePlot(merged.obj.St.Dclk1.Integ,
                features = "ATII_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2<-FeaturePlot(merged.obj.P0.Integ,
                features = "ATII_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<-FeaturePlot(merged.obj.P3.Integ,
                features = "ATII_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<-FeaturePlot(merged.obj.24Hr.Integ,
                features = "ATII_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a5<-FeaturePlot(merged.obj.FiveDay.Integ,
                features = "ATII_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


cowplot::plot_grid(a1, a2, a3, a4, a5, nrow = 2)
dev.off()
##Secretory_Module


merged.obj.P0.Integ <- AddModuleScore(merged.obj.P0.Integ,
                                      features = list(Secretory_genes),
                                      name="Secretory_Module")

merged.obj.P3.Integ <- AddModuleScore(merged.obj.P3.Integ,
                                      features = list(Secretory_genes),
                                      name="Secretory_Module")

merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                        features = list(Secretory_genes),
                                        name="Secretory_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                           features = list(Secretory_genes),
                                           name="Secretory_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(Secretory_genes),
                                            name="Secretory_Module")

png(file = "Secretory.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")

a1<-FeaturePlot(merged.obj.St.Dclk1.Integ,
                features = "Secretory_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


a2<-FeaturePlot(merged.obj.P0.Integ,
                features = "Secretory_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


a3<-FeaturePlot(merged.obj.P3.Integ,
                features = "Secretory_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<-FeaturePlot(merged.obj.24Hr.Integ,
                features = "Secretory_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a5<-FeaturePlot(merged.obj.FiveDay.Integ,
                features = "Secretory_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


cowplot::plot_grid(a1, a2, a3, a4, a5, nrow = 2)
dev.off()

##Ciliated
#######
merged.obj.P0.Integ <- AddModuleScore(merged.obj.P0.Integ,
                                      features = list(Ciliated_genes),
                                      name="Ciliated_Module")

merged.obj.P3.Integ <- AddModuleScore(merged.obj.P3.Integ,
                                      features = list(Ciliated_genes),
                                      name="Ciliated_Module")

merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                        features = list(Ciliated_genes),
                                        name="Ciliated_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                           features = list(Ciliated_genes),
                                           name="Ciliated_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                           features = list(Ciliated_genes),
                                           name="Ciliated_Module")

png(file = "Ciliated.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")

a1<-FeaturePlot(merged.obj.St.Dclk1.Integ,
                features = "Ciliated_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2<-FeaturePlot(merged.obj.P0.Integ,
                features = "Ciliated_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<-FeaturePlot(merged.obj.P3.Integ,
                features = "Ciliated_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<-FeaturePlot(merged.obj.24Hr.Integ,
                features = "Ciliated_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a5<-FeaturePlot(merged.obj.FiveDay.Integ,
                features = "Ciliated_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


cowplot::plot_grid(a1, a2, a3, a4, a5, nrow = 2)
dev.off()
################


###Tuft

merged.obj.P0.Integ <- AddModuleScore(merged.obj.P0.Integ,
                                      features = list(Tuft_genes),
                                      name="Tuft_Module")

merged.obj.P3.Integ <- AddModuleScore(merged.obj.P3.Integ,
                                      features = list(Tuft_genes),
                                      name="Tuft_Module")

merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                        features = list(Tuft_genes),
                                        name="Tuft_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                           features = list(Tuft_genes),
                                           name="Tuft_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                           features = list(Tuft_genes),
                                           name="Tuft_Module")

png(file = "Tuft.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")

a1<- FeaturePlot(merged.obj.St.Dclk1.Integ,
                 features = "Tuft_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
cowplot::plot_grid(a1, a2, a3, a4, a5, nrow = 2)

a2 <- FeaturePlot(merged.obj.P0.Integ,
                  features = "Tuft_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<- FeaturePlot(merged.obj.P3.Integ,
                 features = "Tuft_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<- FeaturePlot(merged.obj.24Hr.Integ,
                 features = "Tuft_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a5<- FeaturePlot(merged.obj.FiveDay.Integ,
                 features = "Tuft_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
cowplot::plot_grid(a1, a2, a3, a4, a5, nrow = 2)


dev.off()



###BASC

merged.obj.P0.Integ <- AddModuleScore(merged.obj.P0.Integ,
                                      features = list(BASC_genes),
                                      name="BASC_Module")

merged.obj.P3.Integ <- AddModuleScore(merged.obj.P3.Integ,
                                      features = list(BASC_genes),
                                      name="BASC_Module")

merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                        features = list(BASC_genes),
                                        name="BASC_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                           features = list(BASC_genes),
                                           name="BASC_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                           features = list(BASC_genes),
                                           name="BASC_Module")


png(file = "BASC.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")

a1<- FeaturePlot(merged.obj.St.Dclk1.Integ,
                 features = "BASC_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2 <- FeaturePlot(merged.obj.P0.Integ,
                  features = "BASC_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<- FeaturePlot(merged.obj.P3.Integ,
                 features = "BASC_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4 <- FeaturePlot(merged.obj.24Hr.Integ,
                 features = "BASC_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a5<- FeaturePlot(merged.obj.FiveDay.Integ,
                 features = "BASC_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))




cowplot::plot_grid(a1, a2, a3, a4, a5, nrow = 2)
dev.off()


###ATII_ATI

merged.obj.P0.Integ <- AddModuleScore(merged.obj.P0.Integ,
                                      features = list(ATII_ATI_genes),
                                      name="ATII_ATI_Module")

merged.obj.P3.Integ <- AddModuleScore(merged.obj.P3.Integ,
                                      features = list(ATII_ATI_genes),
                                      name="ATII_ATI_Module")

merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                        features = list(ATII_ATI_genes),
                                        name="ATII_ATI_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                           features = list(ATII_ATI_genes),
                                           name="ATII_ATI_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                           features = list(ATII_ATI_genes),
                                           name="ATII_ATI_Module")


png(file = "ATII_ATI.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")

a1<- FeaturePlot(merged.obj.St.Dclk1.Integ,
                 features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2 <- FeaturePlot(merged.obj.P0.Integ,
                  features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<- FeaturePlot(merged.obj.P3.Integ,
                 features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<- FeaturePlot(merged.obj.24Hr.Integ,
                 features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a5<- FeaturePlot(merged.obj.FiveDay.Integ,
                 features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



cowplot::plot_grid(a1, a2, a3, a4, a5, nrow = 2)
dev.off()


###Basal

merged.obj.P0.Integ <- AddModuleScore(merged.obj.P0.Integ,
                                      features = list(Basal_Genes),
                                      name="Basal_Module")

merged.obj.P3.Integ <- AddModuleScore(merged.obj.P3.Integ,
                                      features = list(Basal_Genes),
                                      name="Basal_Module")

merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                        features = list(Basal_Genes),
                                        name="Basal_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                           features = list(Basal_Genes),
                                           name="Basal_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(Basal_Genes),
                                            name="Basal_Module")


png(file = "Basal.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<- FeaturePlot(merged.obj.St.Dclk1.Integ,
                 features = "Basal_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2 <- FeaturePlot(merged.obj.P0.Integ,
                  features = "Basal_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<- FeaturePlot(merged.obj.P3.Integ,
                 features = "Basal_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<- FeaturePlot(merged.obj.24Hr.Integ,
                 features = "Basal_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a5<- FeaturePlot(merged.obj.FiveDay.Integ,
                 features = "Basal_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


cowplot::plot_grid(a1, a2, a3, a4, a5, nrow = 2)
dev.off()

png(file = "Eng24HrUMAP.png", width=16, height=12, units = 'in',res=300, type = "cairo")
DimPlot(merged.obj.24Hr.Integ, reduction = "umap.test") + ggtitle("Integrated 24 Hour Object")
dev.off()

png(file = "FiveDayUMAP.png", width=16, height=12, units = 'in',res=300, type = "cairo")
DimPlot(merged.obj.FiveDay.Integ, reduction = "umap.test") + ggtitle("Integrated FiveDay Object")
dev.off()

###Hillock

merged.obj.P0.Integ <- AddModuleScore(merged.obj.P0.Integ,
                                      features = list(Hillock_Genes),
                                      name="Hillock_Module")

merged.obj.P3.Integ <- AddModuleScore(merged.obj.P3.Integ,
                                      features = list(Hillock_Genes),
                                      name="Hillock_Module")

merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                        features = list(Hillock_Genes),
                                        name="Hillock_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                           features = list(Hillock_Genes),
                                           name="Hillock_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(Hillock_Genes),
                                            name="Hillock_Module")


png(file = "Hillock.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<- FeaturePlot(merged.obj.St.Dclk1.Integ,
                 features = "Hillock_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
a2 <- FeaturePlot(merged.obj.P0.Integ,
                  features = "Hillock_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<- FeaturePlot(merged.obj.P3.Integ,
                 features = "Hillock_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<- FeaturePlot(merged.obj.24Hr.Integ,
                 features = "Hillock_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a5<- FeaturePlot(merged.obj.FiveDay.Integ,
                 features = "Hillock_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


cowplot::plot_grid(a1, a2, a3, a4, a5, nrow = 2)
dev.off()



lung.combined <- AddModuleScore(lung.combined,
                                features = list(Basal_Genes),
                                name="Basal_Module")


lung.combined <- AddModuleScore(lung.combined,
                                features = list(Hillock_Genes),
                                name="Hillock_Module")



png(file = "Basal.ModuleScore.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
            features = "Basal_Module1", label = TRUE, repel = TRUE, order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

png(file = "Hillock.ModuleScore.Full_Atlas.png", width=16, height=12, units = 'in',res=300, type = "cairo")
FeaturePlot(lung.combined,
            features = "Hillock_Module1", label = TRUE, repel = TRUE, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()




fix.sc <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-2, 2))

######UNIVERSAL SCALES
png(file = "Basal.ModuleScore.Full_Atlas.UniversalScale.PreAdjustments.png", width=16, height=12, units = 'in',res=300, type = "cairo")
Basal_Module_Atlas <- FeaturePlot(lung.combined, features = "Basal_Module1", label = T, repel = T, combine = FALSE)
Basal_Module_Atlas.2 <- lapply(Basal_Module_Atlas, function (x) x + fix.sc)
CombinePlots(Basal_Module_Atlas.2)
dev.off()

png(file = "Hillock.ModuleScore.Full_Atlas.UniversalScale.PreAdjustments.png", width=16, height=12, units = 'in',res=300, type = "cairo")
Hillock_Module_Atlas <- FeaturePlot(lung.combined, features = "Hillock_Module1", label = T, repel = T, combine = FALSE)
Hillock_Module_Atlas.2 <- lapply(Hillock_Module_Atlas, function (x) x + fix.sc)
CombinePlots(Hillock_Module_Atlas.2)
dev.off()


png(file = "StartingUMAP.png", width=16, height=12, units = 'in',res=300, type = "cairo")
DimPlot(merged.obj.St.Dclk1.Integ, reduction = "umap.test", shuffle = T, label = T) + ggtitle("Integrated St. Dclk1 Object")
dev.off()


png(file = "Starting.Basal_Integrated.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Tp63", label = T, reduction = "umap.test") 
a2<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Krt14", label = T, reduction = "umap.test")
a3<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Krt5", label = T, reduction = "umap.test")
a4<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Itga6", label = T, reduction = "umap.test")
cowplot::plot_grid(a1, a2, a3, a4, nrow = 2)
dev.off()

png(file = "Starting.ATI_Integrated.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Aqp5", label = T, reduction = "umap.test")
a2<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Pdpn", label = T, reduction = "umap.test")
a3<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Ager", label = T, reduction = "umap.test")
cowplot::plot_grid(a1, a2, a3, nrow = 2)
dev.off()


png(file = "Starting.ATII_Integrated.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Sftpc", label = T, reduction = "umap.test")
a2<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Lamp3", label = T, reduction = "umap.test")
a3<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Napsa", label = T, reduction = "umap.test")
cowplot::plot_grid(a1, a2, a3, nrow = 2)
dev.off()


png(file = "Starting.Hillock_Integrated.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Krt13", label = T, reduction = "umap.test") 
a2<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Ecm1", label = T, reduction = "umap.test")
a3<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Cldn3", label = T, reduction = "umap.test")
cowplot::plot_grid(a1, a2, a3, nrow = 2)
dev.off()
png(file = "Starting.Ciliated_Integrated.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Cfap161", label = T, reduction = "umap.test") 
a2<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Cfap44", label = T, reduction = "umap.test")
a3<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Dnali1", label = T, reduction = "umap.test")
a4<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Dnai1", label = T, reduction = "umap.test")
cowplot::plot_grid(a1, a2, a3, a4, nrow = 2)
dev.off()


png(file = "Secretory.ATI_Integrated.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Bpifa5", label = T, reduction = "umap.test")
a2<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Scgb3a2", label = T, reduction = "umap.test")
a3<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Muc5b", label = T, reduction = "umap.test")
cowplot::plot_grid(a1, a2, a3, nrow = 2)
dev.off()

png(file = "Starting.Tuft_Integrated.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Sox9", label = T, reduction = "umap.test")
a2<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Dclk1", label = T, reduction = "umap.test")
a3<-FeaturePlot(merged.obj.St.Dclk1.Integ, features = "Alox5ap", label = T, reduction = "umap.test")
cowplot::plot_grid(a1, a2, a3, nrow = 2)
dev.off()