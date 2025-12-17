##Running setdiff tests to confirm that each group is orthogonal

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

##Finding all the genes that are shared between BASC and Tuft, and thus removing them from BASC for exclusivity towards Tuft
remove.from.BASC <- intersect(BASC_genes, Tuft_genes)
  
BASC_genes <- BASC_genes[! BASC_genes %in% remove.from.BASC]


###Now for the actual Module scores


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


png(file = "AT1.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.P0.Integ,
            features = "ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2<-FeaturePlot(merged.obj.P3.Integ,
            features = "ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<-FeaturePlot(merged.obj.24Hr.Integ,
            features = "ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<-FeaturePlot(merged.obj.FiveDay.Integ,
            features = "ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

cowplot::plot_grid(a1, a2, a3, a4, nrow = 2)
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

png(file = "AT2.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.P0.Integ,
            features = "ATII_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2<-FeaturePlot(merged.obj.P3.Integ,
            features = "ATII_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<-FeaturePlot(merged.obj.24Hr.Integ,
            features = "ATII_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<-FeaturePlot(merged.obj.FiveDay.Integ,
            features = "ATII_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

cowplot::plot_grid(a1, a2, a3, a4, nrow = 2)
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

png(file = "Secretory.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.P0.Integ,
            features = "Secretory_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


a2<-FeaturePlot(merged.obj.P3.Integ,
            features = "Secretory_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<-FeaturePlot(merged.obj.24Hr.Integ,
            features = "Secretory_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<-FeaturePlot(merged.obj.FiveDay.Integ,
            features = "Secretory_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

cowplot::plot_grid(a1, a2, a3, a4, nrow = 2)
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


png(file = "Ciliated.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1<-FeaturePlot(merged.obj.P0.Integ,
            features = "Ciliated_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2<-FeaturePlot(merged.obj.P3.Integ,
            features = "Ciliated_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<-FeaturePlot(merged.obj.24Hr.Integ,
            features = "Ciliated_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<-FeaturePlot(merged.obj.FiveDay.Integ,
            features = "Ciliated_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
cowplot::plot_grid(a1, a2, a3, a4, nrow = 2)
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

png(file = "Tuft.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1 <- FeaturePlot(merged.obj.P0.Integ,
            features = "Tuft_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2<- FeaturePlot(merged.obj.P3.Integ,
            features = "Tuft_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<- FeaturePlot(merged.obj.24Hr.Integ,
            features = "Tuft_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<- FeaturePlot(merged.obj.FiveDay.Integ,
            features = "Tuft_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
cowplot::plot_grid(a1, a2, a3, a4, nrow = 2)
dev.off()



###BASC

merged.obj.P0.Integ <- AddModuleScore(merged.obj.P0.Integ,
                                      features = list(BASC_genes.adjust),
                                      name="BASC_Module")

merged.obj.P3.Integ <- AddModuleScore(merged.obj.P3.Integ,
                                      features = list(BASC_genes.adjust),
                                      name="BASC_Module")

merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                        features = list(BASC_genes.adjust),
                                        name="BASC_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                           features = list(BASC_genes.adjust),
                                           name="BASC_Module")

png(file = "BASC.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1 <- FeaturePlot(merged.obj.P0.Integ,
                  features = "BASC_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2<- FeaturePlot(merged.obj.P3.Integ,
                 features = "BASC_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<- FeaturePlot(merged.obj.24Hr.Integ,
                 features = "BASC_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<- FeaturePlot(merged.obj.FiveDay.Integ,
                 features = "BASC_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
cowplot::plot_grid(a1, a2, a3, a4, nrow = 2)
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

png(file = "ATII_ATI.ModuleScore.png", width=16, height=12, units = 'in',res=300, type = "cairo")
a1 <- FeaturePlot(merged.obj.P0.Integ,
                  features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test.2", order =T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a2<- FeaturePlot(merged.obj.P3.Integ,
                 features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a3<- FeaturePlot(merged.obj.24Hr.Integ,
                 features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

a4<- FeaturePlot(merged.obj.FiveDay.Integ,
                 features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
cowplot::plot_grid(a1, a2, a3, a4, nrow = 2)
dev.off()




png(file = "Eng24HrUMAP.png", width=16, height=12, units = 'in',res=300, type = "cairo")
DimPlot(merged.obj.24Hr.Integ, reduction = "umap.test") + ggtitle("Integrated 24 Hour Object")
dev.off()

png(file = "FiveDayUMAP.png", width=16, height=12, units = 'in',res=300, type = "cairo")
DimPlot(merged.obj.FiveDay.Integ, reduction = "umap.test") + ggtitle("Integrated FiveDay Object")
dev.off()


