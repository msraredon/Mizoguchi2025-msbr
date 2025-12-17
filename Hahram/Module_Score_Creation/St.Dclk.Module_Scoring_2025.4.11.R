merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(ATI_genes),
                                            name="ATI_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(ATII_genes),
                                            name="ATII_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(Ciliated_genes),
                                            name="Ciliated_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(Secretory_genes),
                                            name="Secretory_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(Tuft_genes),
                                            name="Tuft_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(BASC_genes),
                                            name="BASC_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(ATII_ATI_genes),
                                            name="ATII_ATI_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(Basal_Module_Score),
                                            name="Basal_Module")

merged.obj.St.Dclk1.Integ <- AddModuleScore(merged.obj.St.Dclk1.Integ,
                                            features = list(Hillock_Module_Score),
                                            name="Hillock_Module")


FeaturePlot(merged.obj.St.Dclk1.Integ,
            features = "ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(merged.obj.St.Dclk1.Integ,
            features = "ATII_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(merged.obj.St.Dclk1.Integ,
            features = "Ciliated_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(merged.obj.St.Dclk1.Integ,
            features = "Secretory_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(merged.obj.St.Dclk1.Integ,
            features = "Tuft_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(merged.obj.St.Dclk1.Integ,
            features = "BASC_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(merged.obj.St.Dclk1.Integ,
            features = "ATII_ATI_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(merged.obj.St.Dclk1.Integ,
            features = "Basal_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(merged.obj.St.Dclk1.Integ,
            features = "Hillock_Module1", label = TRUE, repel = TRUE, reduction = "umap.test", order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
