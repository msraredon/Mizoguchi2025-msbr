testerfiles <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5

AbbasalGenes.Curated <- c("Cdh2", "Cdkn2a", "Ptgs2", "Hmga2", "Tp63", "Lamb3", "Mmp7", "Ephb2", "Prss2", "Gdf15", "Pxdn")
a.intersect.test <- intersect(AbbasalGenes.Curated, ATI_genes.revised)
#Remove nothing from ab_basal
a.intersect.test <- intersect(AbbasalGenes.Curated, ATII_genes.revised)
#Remove nothing from ab_basal
a.intersect.test <- intersect(AbbasalGenes.Curated, Basal.Genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(AbbasalGenes.Curated, BASC_genes)
#Remove Map2 and Cldn6 from ab_basal
a.intersect.test <- intersect(AbbasalGenes.Curated, ATII_ATI_genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(AbbasalGenes.Curated, Secretory_genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(AbbasalGenes.Curated, Tuft_genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(AbbasalGenes.Curated, Hillock.Genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(AbbasalGenes.Curated, Ciliated_genes)
#Remove nothing from ab_basal


lung.combined <- AddModuleScore(lung.combined, features = list(AbbasalGenes.Curated), name="AbBasal_Module")

fix.sc.ab_basaloid <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.2, 0.35), oob = scales::squish)


png(file = "Ab_Basal.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
AbBasaloid_Module_Atlas <- FeaturePlot(lung.combined, features = "AbBasal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
AbBasaloid_Module_Atlas.2 <- lapply(AbBasaloid_Module_Atlas, function (x) x + fix.sc.ab_basaloid)
CombinePlots(AbBasaloid_Module_Atlas.2) + NoAxes() + ggtitle("Aberrant Basal Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

testerfiles <- AddModuleScore(testerfiles, features = list(AbbasalGenes.Curated), name="AbBasal_Module")

fix.sc.ab_basaloid <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.25, 0.75), oob = scales::squish)

png(file = "Ab_Basal.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
AbBasaloid_Module_Atlas <- FeaturePlot(testerfiles, features = "AbBasal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
AbBasaloid_Module_Atlas.2 <- lapply(AbBasaloid_Module_Atlas, function (x) x + fix.sc.ab_basaloid)
CombinePlots(AbBasaloid_Module_Atlas.2) + NoAxes() + ggtitle("Aberrant Basal Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()

ab_basal_genes.v2.human <- toupper(AbbasalGenes.Curated)
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(ab_basal_genes.v2.human), name="AbBasal_Module")

fix.sc.ab_basaloid <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.2, 2), oob = scales::squish)


png(file = "Ab_Basal.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
AbBasaloid_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "AbBasal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
AbBasaloid_Module_Atlas.2 <- lapply(AbBasaloid_Module_Atlas, function (x) x + fix.sc.ab_basaloid)
CombinePlots(AbBasaloid_Module_Atlas.2) + NoAxes() + ggtitle("Aberrant Basal Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

Hillock.Genes <- c("Itga3", "Itga2", "Itga6", "Itgb1", "Itgb4", "Dsg3", "Cwh43", "Ecm1", "Cav1", "Serpinb2", "Krt6a", "Krt13", "Igfbp3", "Muc4")
a.intersect.test <- intersect(Hillock.Genes, ATI_genes.revised)
#Remove Itga3, Cav1 from ab_basal
a.intersect.test <- intersect(Hillock.Genes, ATII_genes.revised)
#Remove nothing from ab_basal
a.intersect.test <- intersect(Hillock.Genes, Basal.Genes)
#Remove Dsg3 from ab_basal
a.intersect.test <- intersect(Hillock.Genes, BASC_genes)
#Remove Cwh43 from ab_basal
a.intersect.test <- intersect(Hillock.Genes, ATII_ATI_genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(Hillock.Genes, Secretory_genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(Hillock.Genes, Tuft_genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(AbbasalGenes.Curated, Hillock.Genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(Hillock.Genes, Ciliated_genes)
#Remove nothing from ab_basal
Hillock.Genes <- c("Itga2", "Itga6", "Itgb1", "Itgb4", "Ecm1", "Serpinb2", "Krt6a", "Krt13", "Igfbp3", "Muc4")
hillock.human <- toupper(Hillock.Genes)

lung.combined <- AddModuleScore(lung.combined, features = list(Hillock.Genes), name="Hillock_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(hillock.human), name="Hillock_Module")
testerfiles <- AddModuleScore(testerfiles, features = list(Hillock.Genes), name="Hillock_Module")

fix.sc.hillock <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.5, 2), oob = scales::squish)

png(file = "Hillock.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Hillock_Module_Atlas <- FeaturePlot(testerfiles, features = "Hillock_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Hillock_Module_Atlas.2 <- lapply(Hillock_Module_Atlas, function (x) x + fix.sc.hillock)
CombinePlots(Hillock_Module_Atlas.2) + NoAxes() + ggtitle("Hillock Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()

fix.sc.hillock <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.1, .75), oob = scales::squish)

png(file = "Hillock.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Hillock_Module_Atlas <- FeaturePlot(lung.combined, features = "Hillock_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Hillock_Module_Atlas.2 <- lapply(Hillock_Module_Atlas, function (x) x + fix.sc.hillock)
CombinePlots(Hillock_Module_Atlas.2) + NoAxes() + ggtitle("Hillock Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

fix.sc.hillock <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.5, 4), oob = scales::squish)

png(file = "Hillock.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Hillock_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "Hillock_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Hillock_Module_Atlas.2 <- lapply(Hillock_Module_Atlas, function (x) x + fix.sc.hillock)
CombinePlots(Hillock_Module_Atlas.2) + NoAxes() + ggtitle("Hillock Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()




VlnPlot(ipf.control.merge, features = "Hillock_Module1")


VlnPlot(ipf.control.merge, features = "IGFBP2", pt.size = 2.5) + NoLegend()
Basal.Genes <- c("Krt5", "Krt14", "Aqp3", "Krt17", "Col17a1", "Sprr1b", "Krt16", "Dsg3", "Dlk2", "Lama3", "Serpinb5")
Basal.Genes.human <- toupper(Basal.Genes)

a.intersect.test <- intersect(Basal.Genes, ATI_genes.revised)
#Remove Itga3, Cav1 from ab_basal
a.intersect.test <- intersect(Basal.Genes, ATII_genes.revised)
#Remove nothing from ab_basal
a.intersect.test <- intersect(Hillock.Genes, Basal.Genes)
#Remove Dsg3 from ab_basal
a.intersect.test <- intersect(Basal.Genes, BASC_genes)
#Remove Cwh43 from ab_basal
a.intersect.test <- intersect(Basal.Genes, ATII_ATI_genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(Basal.Genes, Secretory_genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(Basal.Genes, Tuft_genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(AbbasalGenes.Curated, Basal.Genes)
#Remove nothing from ab_basal
a.intersect.test <- intersect(Basal.Genes, Ciliated_genes)
#Remove nothing from ab_basal
lung.combined <- AddModuleScore(lung.combined, features = list(Basal.Genes), name="Basal_Module")

ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(Basal.Genes.human), name="Basal_Module")

testerfiles <- AddModuleScore(testerfiles, features = list(Basal.Genes), name="Basal_Module")


fix.sc.basal <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-1, 1.5), oob = scales::squish)

png(file = "Basal.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Basal_Module_Atlas <- FeaturePlot(testerfiles, features = "Basal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Basal_Module_Atlas.2 <- lapply(Basal_Module_Atlas, function (x) x + fix.sc.basal)
CombinePlots(Basal_Module_Atlas.2) + NoAxes() + ggtitle("Basal Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()


fix.sc.basal <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-1, 4.5), oob = scales::squish)

png(file = "Basal.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Basal_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "Basal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Basal_Module_Atlas.2 <- lapply(Basal_Module_Atlas, function (x) x + fix.sc.basal)
CombinePlots(Basal_Module_Atlas.2) + NoAxes() + ggtitle("Basal Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()


fix.sc.basal <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, .6), oob = scales::squish)

png(file = "Basal.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Basal_Module_Atlas <- FeaturePlot(lung.combined, features = "Basal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Basal_Module_Atlas.2 <- lapply(Basal_Module_Atlas, function (x) x + fix.sc.basal)
CombinePlots(Basal_Module_Atlas.2) + NoAxes() + ggtitle("Basal Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()



AbbasalGenes.Curated.Edited <- c("Cdh2", "Cdkn2a", "Ptgs2", "Hmga2", "Itgb6", "Gpr87", "Lamb3", "Mmp7", "Ephb2", "Prss2", "Gdf15", "Pxdn")
lung.combined <- AddModuleScore(lung.combined, features = list(AbbasalGenes.Curated.Edited), name="AbBasal_Module")
testerfiles <- AddModuleScore(testerfiles, features = list(AbbasalGenes.Curated.Edited), name="AbBasal_Module")
ab_basal_genes.v2.human <- toupper(AbbasalGenes.Curated.Edited)
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(ab_basal_genes.v2.human), name="AbBasal_Module")


png(file = "Krt17_Expression.png", width=18, height=12, units = 'in',res=300, type = "cairo")
VlnPlot(ipf.control.merge, features = "KRT17") + NoLegend()
dev.off()

VlnPlot(ipf.control.merge, features = c("KRT17", "TP63"), pt.size = 0) + NoLegend()

krt17.avg <-AverageExpression(ipf.control.merge, features = c("KRT17", "TP63"), return.seurat = T)
View(krt17.avg$RNA)
krt17.avg <- as.data.frame(krt17.avg)
