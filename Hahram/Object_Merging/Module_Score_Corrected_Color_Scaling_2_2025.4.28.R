load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.4.16/merged.obj.Embryo.Integ.updated.2025.4.16.Robj")

load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.3.27/merged.obj.P0.Integ.2025.3.27.updated.Robj")

load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.3.27/merged.obj.P3.Integ.2025.3.27.updated.Robj")

load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.3.28/merged.obj.24Hr.Integ.updated.2025.3.28.Robj")

DimPlot(merged.obj.Embryo.Integ, reduction = "umap.test")
DimPlot(merged.obj.P0.Integ, reduction = "umap.test")
DimPlot(merged.obj.P3.Integ, reduction = "umap.test")

png(file = 'Figure3.UMAP.Embryo.Split.png', width=9, height=5, units = 'in',res=300, type = "cairo")
DimPlot(merged.obj.Embryo.Integ, split.by = "Sample", cols = c("#ECA72C",  "#FF47DA",  "#40F99B",  "#462255",  "#8DE4FF",  "#FFA69E",  "#BAFF29",  "#8D0801"),reduction = "umap.test") + NoAxes() + ggtitle("Embryogenesis UMAP Split By Sample") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

png(file = 'Figure3.UMAP.Starting.Split.png', width=9, height=5, units = 'in',res=300, type = "cairo")
DimPlot(merged.obj.St.Dclk1.Integ, split.by = "Sample", cols = c("#ECA72C",  "#FF47DA",  "#40F99B",  "#462255",  "#8DE4FF",  "#FFA69E",  "#BAFF29",  "#8D0801"),reduction = "umap.test") + NoAxes() + ggtitle("St. Dclk1 UMAP Split By Sample") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

png(file = 'Figure3.UMAP.Starting.png', width=16, height=12, units = 'in',res=300, type = "cairo")
DimPlot(merged.obj.St.Dclk1.Integ, cols = c("#ECA72C",  "#FF47DA",  "#40F99B",  "#462255",  "#8DE4FF",  "#FFA69E",  "#BAFF29",  "#8D0801"),reduction = "umap.test") + NoAxes() + ggtitle("UMAP of Integrated St. Dclk1 Population") + theme(plot.title = element_text(hjust = 0.5))
dev.off()


png(file = 'Figure3.UMAP.P0.Split.png', width=9, height=5, units = 'in',res=300, type = "cairo")
DimPlot(merged.obj.P0.Integ, group.by= "AnnotationOrdered", split.by = "Sample", cols = c("#ECA72C",  "#FF47DA",  "#40F99B",  "#462255",  "#8DE4FF",  "#FFA69E",  "#BAFF29",  "#8D0801"),reduction = "umap.test") + NoAxes() + ggtitle("P0 UMAP Split By Sample") + theme(plot.title = element_text(hjust = 0.5))
dev.off()


png(file = 'Figure3.UMAP.P3.Split.png', width=9, height=5, units = 'in',res=300, type = "cairo")
DimPlot(merged.obj.P3.Integ, group.by= "AnnotationOrdered", split.by = "Sample", cols = c("#ECA72C",  "#FF47DA",  "#40F99B",  "#462255",  "#8DE4FF",  "#FFA69E",  "#BAFF29",  "#8D0801"), reduction = "umap.test") + NoAxes() + ggtitle("P3 UMAP Split By Sample") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

DimPlot(merged.obj.24Hr.Integ, group.by = "prelim.annotations.2025.3.28", reduction = "umap.test")
DimPlot(merged.obj.Embryo.Integ, group.by = "annotations.2025.4.16", reduction = "umap.test")
DimPlot(merged.obj.FiveDay.Integ, group.by = "annotations.updated.2025.3.28", reduction = "umap.test")
DimPlot(merged.obj.P0.Integ, group.by = "prelim.annotations.2025.3.26", reduction = "umap.test.2")
DimPlot(merged.obj.P3.Integ, group.by = "cell_population", reduction = "umap.test")
DimPlot(merged.obj.St.Dclk1.Integ, group.by = "new.annotations.2025.4.10", reduction = "umap.test")


DimPlot(merged.obj.24Hr.Integ, reduction = "umap.test")
DimPlot(merged.obj.Embryo.Integ, reduction = "umap.test")
DimPlot(merged.obj.FiveDay.Integ, reduction = "umap.test")
DimPlot(merged.obj.P0.Integ, reduction = "umap.test.2")
DimPlot(merged.obj.P3.Integ, reduction = "umap.test")
DimPlot(merged.obj.St.Dclk1.Integ, reduction = "umap.test")

DimPlot(tissue_engineering_epithelial_merged.dataset, group.by = "annotations.updated.2025.4.28")

tissue_engineering_epithelial_merged.dataset$annotations.updated.2025.4.28 <- c(tissue_engineering_epithelial_merged.dataset$prelim.annotations.2025.3.28, tissue_engineering_epithelial_merged.dataset$annotations.2025.4.16, tissue_engineering_epithelial_merged.dataset$annotations.updated.2025.3.28, tissue_engineering_epithelial_merged.dataset$prelim.annotations.2025.3.26, tissue_engineering_epithelial_merged.dataset$new.annotations.2025.4.10)


tissue_engineering_epithelial_merged.dataset$annotations.updated.2025.4.28 <- plyr::mapvalues(
  x = tissue_engineering_epithelial_merged.dataset$orig.ident, 
  from = c("sample1_lane1", "sample1_lane2", "sample2_lane1", "sample2_lane2"), 
  to = c("Sample1", "Sample1", "Sample2", "Sample2")
)


DimPlot(tissue_engineering_epithelial_merged.dataset, group.by = "annotations.updated.2025.3.28")
DimPlot(tissue_engineering_epithelial_merged.dataset, group.by = "prelim.annotations.2025.3.28")
DimPlot(tissue_engineering_epithelial_merged.dataset, group.by = "prelim.annotations.2025.3.26")
DimPlot(tissue_engineering_epithelial_merged.dataset, group.by = "prelim.annotations.2025.3.26")
DimPlot(tissue_engineering_epithelial_merged.dataset, group.by = "prelim.annotations.2025.3.26")
DimPlot(tissue_engineering_epithelial_merged.dataset, group.by = "annotations.2025.4.16")
DimPlot(tissue_engineering_epithelial_merged.dataset, group.by = "cell_population")


HK_tissue_engineering_epithelial_merged.dataset.2025.4.28 <- tissue_engineering_epithelial_merged.dataset
HK_tissue_engineering_epithelial_merged.dataset.2025.4.28$Cell_Populations.Updated.2025.4.28 <-
  HK_tissue_engineering_epithelial_merged.dataset.2025.4.28$cell_population
save(HK_tissue_engineering_epithelial_merged.dataset.2025.4.28, file = "HK_tissue_engineering_epithelial_merged.dataset.2025.4.28.Robj")


Idents(tissue_engineering_epithelial_merged.dataset) <- tissue_engineering_epithelial_merged.dataset$Condition
DimPlot(tissue_engineering_epithelial_merged.dataset)

fix.sc <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.5, 2), oob = scales::squish)

fix.sc.basal <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-1, 2), oob = scales::squish)

png(file = "Basal.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Basal_Module_Atlas <- FeaturePlot(tissue_engineering_epithelial_merged.dataset, features = "Basal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE)
Basal_Module_Atlas.2 <- lapply(Basal_Module_Atlas, function (x) x + fix.sc.basal)
CombinePlots(Basal_Module_Atlas.2)
dev.off()

fix.sc.hillock <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.5, 2), oob = scales::squish)

png(file = "Hillock.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Hillock_Module_Atlas <- FeaturePlot(tissue_engineering_epithelial_merged.dataset, features = "Hillock_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE)
Hillock_Module_Atlas.2 <- lapply(Hillock_Module_Atlas, function (x) x + fix.sc.hillock)
CombinePlots(Hillock_Module_Atlas.2)
dev.off()

fix.sc.at2 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.5, 2), oob = scales::squish)

png(file = "AT2.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
AT2_Module_Atlas <- FeaturePlot(tissue_engineering_epithelial_merged.dataset, features = "ATII_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE)
AT2_Module_Atlas.2 <- lapply(AT2_Module_Atlas, function (x) x + fix.sc.at2)
CombinePlots(AT2_Module_Atlas.2)
dev.off()

fix.sc.at1 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, .75), oob = scales::squish)

png(file = "ATI.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ATI_Module_Atlas <- FeaturePlot(tissue_engineering_epithelial_merged.dataset, features = "ATI_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE)
ATI_Module_Atlas.2 <- lapply(ATI_Module_Atlas, function (x) x + fix.sc.at1)
CombinePlots(ATI_Module_Atlas.2)
dev.off()


fix.sc.secretory <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.5, 2), oob = scales::squish)

png(file = "Secretory.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Secretory_Module_Atlas <- FeaturePlot(tissue_engineering_epithelial_merged.dataset, features = "Secretory_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE)
Secretory_Module_Atlas.2 <- lapply(Secretory_Module_Atlas, function (x) x + fix.sc.secretory)
CombinePlots(Secretory_Module_Atlas.2)
dev.off()


fix.sc.ciliated <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, 1), oob = scales::squish)

png(file = "Ciliated.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Ciliated_Module_Atlas <- FeaturePlot(tissue_engineering_epithelial_merged.dataset, features = "Ciliated_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE)
Ciliated_Module_Atlas.2 <- lapply(Ciliated_Module_Atlas, function (x) x + fix.sc.ciliated)
CombinePlots(Ciliated_Module_Atlas.2)
dev.off()


fix.sc.tuft <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, .75), oob = scales::squish)

png(file = "Tuft.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Tuft_Module_Atlas <- FeaturePlot(tissue_engineering_epithelial_merged.dataset, features = "Tuft_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE)
Tuft_Module_Atlas.2 <- lapply(Tuft_Module_Atlas, function (x) x + fix.sc.tuft)
CombinePlots(Tuft_Module_Atlas.2)
dev.off()


fix.sc.basc <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.25, 0.75), oob = scales::squish)

png(file = "BASC.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
BASC_Module_Atlas <- FeaturePlot(tissue_engineering_epithelial_merged.dataset, features = "BASC_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE)
BASC_Module_Atlas.2 <- lapply(BASC_Module_Atlas, function (x) x + fix.sc.basc)
CombinePlots(BASC_Module_Atlas.2)
dev.off()


fix.sc.at2_at1 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, 1), oob = scales::squish)

png(file = "ATII_ATI_ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ATII_ATI_Module_Atlas <- FeaturePlot(tissue_engineering_epithelial_merged.dataset, features = "ATII_ATI_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE)
ATII_ATI_Module_Atlas.2 <- lapply(ATII_ATI_Module_Atlas, function (x) x + fix.sc.at2_at1)
CombinePlots(ATII_ATI_Module_Atlas.2)
dev.off()


FeaturePlot(tissue_engineering_epithelial_merged.dataset, features = "Dppa2")
