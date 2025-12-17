lung.combined <- AddModuleScore(lung.combined, features = list(Basal_genes.revised), name="Basal_Module")
lung.combined <- AddModuleScore(lung.combined, features = list(Hillock_genes.revised), name="Hillock_Module")
lung.combined <- AddModuleScore(lung.combined, features = list(ATI_genes.revised), name="ATI_Module")
lung.combined <- AddModuleScore(lung.combined, features = list(ATII_genes.revised), name="ATII_Module")
lung.combined <- AddModuleScore(lung.combined, features = list(Tuft_genes), name="Tuft_Module")
lung.combined <- AddModuleScore(lung.combined,  features = list(Ciliated_genes), name="Ciliated_Module")
lung.combined <- AddModuleScore(lung.combined,  features = list(Secretory_genes), name="Secretory_Module")
lung.combined <- AddModuleScore(lung.combined, features = list(BASC_genes), name="BASC_Module")
lung.combined <- AddModuleScore(lung.combined, features = list(ATII_ATI_genes), name="ATII_ATI_Module")
lung.combined <- AddModuleScore(lung.combined, features = list(Hillock_genes), name="Hillock_Module")
lung.combined <- AddModuleScore(lung.combined, features = list(AbBasal_genes.revised), name="AbBasal_Module")

fix.sc.basal <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, .6), oob = scales::squish)

png(file = "Basal.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Basal_Module_Atlas <- FeaturePlot(lung.combined, features = "Basal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Basal_Module_Atlas.2 <- lapply(Basal_Module_Atlas, function (x) x + fix.sc.basal)
CombinePlots(Basal_Module_Atlas.2) + NoAxes() + ggtitle("Basal Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

fix.sc.hillock <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.1, .75), oob = scales::squish)

png(file = "Hillock.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Hillock_Module_Atlas <- FeaturePlot(lung.combined, features = "Hillock_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Hillock_Module_Atlas.2 <- lapply(Hillock_Module_Atlas, function (x) x + fix.sc.hillock)
CombinePlots(Hillock_Module_Atlas.2) + NoAxes() + ggtitle("Hillock Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

fix.sc.at2 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-1, 4), oob = scales::squish)

png(file = "AT2.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
AT2_Module_Atlas <- FeaturePlot(lung.combined, features = "ATII_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
AT2_Module_Atlas.2 <- lapply(AT2_Module_Atlas, function (x) x + fix.sc.at2)
CombinePlots(AT2_Module_Atlas.2) + NoAxes() + ggtitle("ATII Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

fix.sc.at1 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.5, 1.5), oob = scales::squish)

png(file = "ATI.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ATI_Module_Atlas <- FeaturePlot(lung.combined, features = "ATI_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
ATI_Module_Atlas.2 <- lapply(ATI_Module_Atlas, function (x) x + fix.sc.at1)
CombinePlots(ATI_Module_Atlas.2) + NoAxes() + ggtitle("ATI Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()


fix.sc.secretory <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.5, 2.5), oob = scales::squish)

png(file = "Secretory.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Secretory_Module_Atlas <- FeaturePlot(lung.combined, features = "Secretory_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Secretory_Module_Atlas.2 <- lapply(Secretory_Module_Atlas, function (x) x + fix.sc.secretory)
CombinePlots(Secretory_Module_Atlas.2) + NoAxes() + ggtitle("Secretory Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()


fix.sc.ciliated <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, .75), oob = scales::squish)

png(file = "Ciliated.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Ciliated_Module_Atlas <- FeaturePlot(lung.combined, features = "Ciliated_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Ciliated_Module_Atlas.2 <- lapply(Ciliated_Module_Atlas, function (x) x + fix.sc.ciliated)
CombinePlots(Ciliated_Module_Atlas.2) + NoAxes() + ggtitle("Ciliated Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")+ theme(title = element_text(size = 25))
dev.off() 


fix.sc.tuft <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.5, 2.5), oob = scales::squish)

png(file = "Tuft.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Tuft_Module_Atlas <- FeaturePlot(lung.combined, features = "Tuft_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Tuft_Module_Atlas.2 <- lapply(Tuft_Module_Atlas, function (x) x + fix.sc.tuft)
CombinePlots(Tuft_Module_Atlas.2) + NoAxes() + ggtitle("Tuft Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()


fix.sc.basc <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.1, 0.5), oob = scales::squish)

png(file = "BASC.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
BASC_Module_Atlas <- FeaturePlot(lung.combined, features = "BASC_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
BASC_Module_Atlas.2 <- lapply(BASC_Module_Atlas, function (x) x + fix.sc.basc)
CombinePlots(BASC_Module_Atlas.2) + NoAxes() + ggtitle("BASC Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()


fix.sc.at2_at1 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.2, .75), oob = scales::squish)

png(file = "ATII_ATI_ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ATII_ATI_Module_Atlas <- FeaturePlot(lung.combined, features = "ATII_ATI_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
ATII_ATI_Module_Atlas.2 <- lapply(ATII_ATI_Module_Atlas, function (x) x + fix.sc.at2_at1)
CombinePlots(ATII_ATI_Module_Atlas.2) + NoAxes() + ggtitle("ATII-ATI Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

fix.sc.ab_basaloid <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.2, 0.35), oob = scales::squish)


png(file = "Ab_Basal.ModuleScore.LungCombined.png", width=18, height=12, units = 'in',res=300, type = "cairo")
AbBasaloid_Module_Atlas <- FeaturePlot(lung.combined, features = "AbBasal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
AbBasaloid_Module_Atlas.2 <- lapply(AbBasaloid_Module_Atlas, function (x) x + fix.sc.ab_basaloid)
CombinePlots(AbBasaloid_Module_Atlas.2) + NoAxes() + ggtitle("Aberrant Basal Module Score Feature Plot In Rat Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()




Basal_genes.human <- toupper(Basal_genes.revised)
Hillock_genes.human <- toupper(Hillock_genes.revised)
ATI_genes.revised.human <- toupper(ATI_genes.revised)
ATII_genes.revised.human <- toupper(ATII_genes.revised)
Tuft_genes.human <- toupper(Tuft_genes)
Ciliated_genes.human <- toupper(Ciliated_genes)
Secretory_genes.human <- toupper(Secretory_genes)
BASC_genes.human <- toupper(BASC_genes)
ATII_ATI_genes.human <- toupper(ATII_ATI_genes)
Hillock_genes.human <- toupper(Hillock_genes)
ab_basal_genes.v2.human <- toupper(AbBasal_genes.revised)

ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(Basal_genes.human), name="Basal_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(Hillock_genes.human), name="Hillock_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(ATI_genes.revised.human), name="ATI_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(ATII_genes.revised.human), name="ATII_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(Tuft_genes.human), name="Tuft_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge,  features = list(Ciliated_genes.human), name="Ciliated_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge,  features = list(Secretory_genes.human), name="Secretory_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(BASC_genes.human), name="BASC_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(ATII_ATI_genes.human), name="ATII_ATI_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(Hillock_genes.human), name="Hillock_Module")
ipf.control.merge <- AddModuleScore(ipf.control.merge, features = list(ab_basal_genes.v2.human), name="AbBasal_Module")


fix.sc.basal <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-1, 4.5), oob = scales::squish)

png(file = "Basal.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Basal_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "Basal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Basal_Module_Atlas.2 <- lapply(Basal_Module_Atlas, function (x) x + fix.sc.basal)
CombinePlots(Basal_Module_Atlas.2) + NoAxes() + ggtitle("Basal Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

fix.sc.hillock <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.5, 4), oob = scales::squish)

png(file = "Hillock.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Hillock_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "Hillock_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Hillock_Module_Atlas.2 <- lapply(Hillock_Module_Atlas, function (x) x + fix.sc.hillock)
CombinePlots(Hillock_Module_Atlas.2) + NoAxes() + ggtitle("Hillock Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

fix.sc.at2 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-2, 10), oob = scales::squish)

png(file = "AT2.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
AT2_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "ATII_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
AT2_Module_Atlas.2 <- lapply(AT2_Module_Atlas, function (x) x + fix.sc.at2)
CombinePlots(AT2_Module_Atlas.2) + NoAxes() + ggtitle("ATII Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

fix.sc.at1 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.5, 2), oob = scales::squish)

png(file = "ATI.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ATI_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "ATI_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
ATI_Module_Atlas.2 <- lapply(ATI_Module_Atlas, function (x) x + fix.sc.at1)
CombinePlots(ATI_Module_Atlas.2) + NoAxes() + ggtitle("ATI Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()


fix.sc.secretory <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-2, 8), oob = scales::squish)

png(file = "Secretory.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Secretory_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "Secretory_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Secretory_Module_Atlas.2 <- lapply(Secretory_Module_Atlas, function (x) x + fix.sc.secretory)
CombinePlots(Secretory_Module_Atlas.2) + NoAxes() + ggtitle("Secretory Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()


fix.sc.ciliated <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.5, 2), oob = scales::squish)

png(file = "Ciliated.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Ciliated_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "Ciliated_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Ciliated_Module_Atlas.2 <- lapply(Ciliated_Module_Atlas, function (x) x + fix.sc.ciliated)
CombinePlots(Ciliated_Module_Atlas.2) + NoAxes() + ggtitle("Ciliated Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")+ theme(title = element_text(size = 25))
dev.off() 


fix.sc.tuft <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.5, 1), oob = scales::squish)

png(file = "Tuft.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Tuft_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "Tuft_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Tuft_Module_Atlas.2 <- lapply(Tuft_Module_Atlas, function (x) x + fix.sc.tuft)
CombinePlots(Tuft_Module_Atlas.2) + NoAxes() + ggtitle("Tuft Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()


fix.sc.basc <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.1, 0.35), oob = scales::squish)

png(file = "BASC.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
BASC_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "BASC_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
BASC_Module_Atlas.2 <- lapply(BASC_Module_Atlas, function (x) x + fix.sc.basc)
CombinePlots(BASC_Module_Atlas.2) + NoAxes() + ggtitle("BASC Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()


fix.sc.at2_at1 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, 1.5), oob = scales::squish)

png(file = "ATII_ATI_ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ATII_ATI_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "ATII_ATI_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
ATII_ATI_Module_Atlas.2 <- lapply(ATII_ATI_Module_Atlas, function (x) x + fix.sc.at2_at1)
CombinePlots(ATII_ATI_Module_Atlas.2) + NoAxes() + ggtitle("ATII-ATI Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()

fix.sc.ab_basaloid <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.2, 2), oob = scales::squish)


png(file = "Ab_Basal.ModuleScore.IPFAtlas.png", width=18, height=12, units = 'in',res=300, type = "cairo")
AbBasaloid_Module_Atlas <- FeaturePlot(ipf.control.merge, features = "AbBasal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
AbBasaloid_Module_Atlas.2 <- lapply(AbBasaloid_Module_Atlas, function (x) x + fix.sc.ab_basaloid)
CombinePlots(AbBasaloid_Module_Atlas.2) + NoAxes() + ggtitle("Aberrant Basal Module Score Feature Plot In Human IPF Atlas") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold") + theme(title = element_text(size = 25))
dev.off()


