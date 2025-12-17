Transition_genes <- intersect()


Abb_Basal_Genes<-c(
                      "Mdm2",
                      "Cdh2",
                      "Camk2n1",
                      "Mmp7",
                      "Ephb2",
                      "Spink1",
                      "Ptgs2",
                      "Cdkn2a",
                      "Cdkn2b",
                      "Hmga2",
                      "Vim",
                      "Fn1",
                      "Col1a1",
                      "Tnc",
                      "Vcan",
                      "Pcp4",
                      "Cux2",
                      "Tgfbi",
                      "Ociad2",
                      "Kcnn4",
                      "Hnf4a",
                      "Tp63",
                      "Sox9",
                      "Sox4",
                      "Sox11"
                      )
FeaturePlot(merged.obj.FiveDay.Integ, features = "Mmp7", order = T)
merged.obj.24Hr.Integ <- AddModuleScore(merged.obj.24Hr.Integ,
                                           features = list(Abb_Basal_Genes),
                                           name="Abb_Basal_Module")

merged.obj.FiveDay.Integ <- AddModuleScore(merged.obj.FiveDay.Integ,
                                            features = list(Abb_Basal_Genes),
                                            name="Abb_Basal_Module")

FeaturePlot(merged.obj.24Hr.Integ,
            features = "Abb_Basal_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(merged.obj.FiveDay.Integ,
            features = "Abb_Basal_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



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
