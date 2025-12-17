####Generating Heatmaps of Our Module Scores To Examine Specificity

png(file = "Rat_Lung.png", width=16, height=12, units = 'in',res=300, type = "cairo")
DimPlot(lung.combined) + NoAxes() + ggtitle("Rat Lung Atlas Grouped By Cell Class")
dev.off()

png(file = "Rat_Lung.CellType.png", width=16, height=12, units = 'in',res=300, type = "cairo")
DimPlot(lung.combined, group.by = "CellType_Final") +NoAxes() + ggtitle("Rat Lung Atlas Grouped By Annotations")
dev.off()



library(Seurat)
library(ggplot2)

Idents(lung.combined) <- lung.combined$CellClass_Final
lung.combined.epithelial <- subset(lung.combined, idents = "Epithelium")

lung.combined.epithelial <-FindVariableFeatures(lung.combined.epithelial)
lung.combined.epithelial <-ScaleData(lung.combined.epithelial)
lung.combined.epithelial <-RunPCA(lung.combined.epithelial)

ElbowPlot(lung.combined.epithelial, ndims = 50)
PCHeatmap(lung.combined.epithelial, cells = 500, balanced = T, dims =1:9)
PCHeatmap(lung.combined.epithelial, cells = 500, balanced = T, dims =10:18)
PCHeatmap(lung.combined.epithelial, cells = 500, balanced = T, dims =19:27)
PCHeatmap(lung.combined.epithelial, cells = 500, balanced = T, dims =28:36)
lung.combined.epithelial <- RunUMAP(lung.combined.epithelial, dims = c(1:3, 5, 10, 13))
DimPlot(lung.combined.epithelial)

lung.combined.epithelial <- FindNeighbors(lung.combined.epithelial,  dims = c(1:3, 5, 10, 13))
lung.combined.epithelial <- FindClusters(lung.combined.epithelial, resolution = .1)
DimPlot(lung.combined.epithelial, group.by = "CellType_Final")
DimPlot(lung.combined.epithelial, group.by = "CellType.combined.Integrated", reduction = "umap.test")

DimPlot(lung.combined, group.by = "orig.ident", split.by = "CellType_Final")
DimPlot(lung.combined)
table(lung.combined.epithelial$orig.ident)
DimPlot(lung.combined.epithelial)

table(lung.combined.epithelial$CellType_Final)
Idents(lung.combined.epithelial) <- lung.combined.epithelial$CellType_Final
downsampled.lung.combined.epithelial <- subset(lung.combined.epithelial, downsample = 201)
downsampled.lung.combined.epithelial <- ScaleData(downsampled.lung.combined.epithelial,
                                                  features = rownames(downsampled.lung.combined.epithelial))

color.palette <- list()

color.palette$Sample <- c("#A3C3D9", "#FBAF00","#574B60", "green","purple","orange", 
                          "#FF47DA", "red4", "#8DE4FF", "#D81E5B", "#87A878", "#DBF9B8", "#090909", "#FFCDBC")
names(color.palette$Sample) <- c("2-12", "2-14", "2-16", "2-28", "2-7", "2-8", "3-1", 
                                 "8-10", "fRat", "mRat", "P0-f", "P0-m",  
                                 "rat1A", "rat811")

color.palette$Annotation <- c("red", "blue", "pink", "yellow", "purple", "orange", "green", "brown")
DimPlot(lung.combined.epithelial, cols = color.palette$Annotation)
unique(downsampled.lung.combined.epithelial$CellType_Final)

png(file = "Rat_Lung_Epithelial.png", width=16, height=12, units = 'in',res=300, type = "cairo")
DimPlot(lung.combined.epithelial, group.by = "CellType_Final", cols = color.palette$Annotation, label = T, repel = T) +NoAxes() + ggtitle("Epithelial Subset of Rat Cell Atlas")
dev.off()

names(color.palette$Annotation) <- c("ATI", "ATII", "ATII-ATI", "BASC", "Ciliated", "Secretory", "Tuft")


downsampled.lung.combined.epithelial$reordered_anno <- factor(downsampled.lung.combined.epithelial$CellType_Final, levels = c("ATI", "ATII", "ATII-ATI", "BASC", "Ciliated", "Secretory", "Tuft"))

png(file = "ModScore_Heatmap.Tuft.png", width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.lung.combined.epithelial,
                   data.type = 'RNA', 
                   primary = 'reordered_anno', 
                   secondary = 'orig.ident',
                   #tertiary = 'Orig_ID', 
                   primary.cols = color.palette$Annotation, 
                   secondary.cols = color.palette$Sample, 
                   #tertiary.cols = color.palette$Sample, 
                   features = Tuft_genes, 
                   labels = c('Annotation', 'Sample'), 
                   selected.row.anotations=Tuft_genes,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()


png(file = "ModScore_Heatmap.ATII.png", width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.lung.combined.epithelial,
                   data.type = 'RNA', 
                   primary = 'reordered_anno', 
                   secondary = 'orig.ident',
                   #tertiary = 'Orig_ID', 
                   primary.cols = color.palette$Annotation, 
                   secondary.cols = color.palette$Sample, 
                   #tertiary.cols = color.palette$Sample, 
                   features = ATII_genes[1:50], 
                   labels = c('Annotation', 'Sample'), 
                   selected.row.anotations=ATII_genes[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

png(file = "ModScore_Heatmap.ATII.2.png", width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.lung.combined.epithelial,
                   data.type = 'RNA', 
                   primary = 'reordered_anno', 
                   secondary = 'orig.ident',
                   #tertiary = 'Orig_ID', 
                   primary.cols = color.palette$Annotation, 
                   secondary.cols = color.palette$Sample, 
                   #tertiary.cols = color.palette$Sample, 
                   features = ATII_genes[51:77], 
                   labels = c('Annotation', 'Sample'), 
                   selected.row.anotations=ATII_genes[51:77],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

png(file = "ModScore_Heatmap.Secretory.png", width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.lung.combined.epithelial,
                   data.type = 'RNA', 
                   primary = 'reordered_anno', 
                   secondary = 'orig.ident',
                   #tertiary = 'Orig_ID', 
                   primary.cols = color.palette$Annotation, 
                   secondary.cols = color.palette$Sample, 
                   #tertiary.cols = color.palette$Sample, 
                   features = Secretory_genes, 
                   labels = c('Annotation', 'Sample'), 
                   selected.row.anotations=Secretory_genes,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.2)
dev.off()

png(file = "ModScore_Heatmap.Ciliated.png", width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.lung.combined.epithelial,
                   data.type = 'RNA', 
                   primary = 'reordered_anno', 
                   secondary = 'orig.ident',
                   #tertiary = 'Orig_ID', 
                   primary.cols = color.palette$Annotation, 
                   secondary.cols = color.palette$Sample, 
                   #tertiary.cols = color.palette$Sample, 
                   features = Ciliated_genes[1:55], 
                   labels = c('Annotation', 'Sample'), 
                   selected.row.anotations=Ciliated_genes[1:55],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

png(file = "ModScore_Heatmap.ATI.png", width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.lung.combined.epithelial,
                   data.type = 'RNA', 
                   primary = 'reordered_anno', 
                   secondary = 'orig.ident',
                   #tertiary = 'Orig_ID', 
                   primary.cols = color.palette$Annotation, 
                   secondary.cols = color.palette$Sample, 
                   #tertiary.cols = color.palette$Sample, 
                   features = ATI_genes[1:50], 
                   labels = c('Annotation', 'Sample'), 
                   selected.row.anotations=ATI_genes[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()

png(file = "ModScore_Heatmap.ATI.2.png", width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.lung.combined.epithelial,
                   data.type = 'RNA', 
                   primary = 'reordered_anno', 
                   secondary = 'orig.ident',
                   #tertiary = 'Orig_ID', 
                   primary.cols = color.palette$Annotation, 
                   secondary.cols = color.palette$Sample, 
                   #tertiary.cols = color.palette$Sample, 
                   features = ATI_genes[51:89], 
                   labels = c('Annotation', 'Sample'), 
                   selected.row.anotations=ATI_genes[51:89],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()
png(file = "ModScore_Heatmap.BASC.png", width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.lung.combined.epithelial,
                   data.type = 'RNA', 
                   primary = 'reordered_anno', 
                   secondary = 'orig.ident',
                   #tertiary = 'Orig_ID', 
                   primary.cols = color.palette$Annotation, 
                   secondary.cols = color.palette$Sample, 
                   #tertiary.cols = color.palette$Sample, 
                   features = BASC_genes.adjust[1:50], 
                   labels = c('Annotation', 'Sample'), 
                   selected.row.anotations=BASC_genes.adjust[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()
png(file = "ModScore_Heatmap.BASC.2.png", width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.lung.combined.epithelial,
                   data.type = 'RNA', 
                   primary = 'reordered_anno', 
                   secondary = 'orig.ident',
                   #tertiary = 'Orig_ID', 
                   primary.cols = color.palette$Annotation, 
                   secondary.cols = color.palette$Sample, 
                   #tertiary.cols = color.palette$Sample, 
                   features = BASC_genes.adjust[51:77], 
                   labels = c('Annotation', 'Sample'), 
                   selected.row.anotations=BASC_genes.adjust[51:77],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()
png(file = "ModScore_Heatmap.ATII_ATI.png", width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.lung.combined.epithelial,
                   data.type = 'RNA', 
                   primary = 'reordered_anno', 
                   secondary = 'orig.ident',
                   #tertiary = 'Orig_ID', 
                   primary.cols = color.palette$Annotation, 
                   secondary.cols = color.palette$Sample, 
                   #tertiary.cols = color.palette$Sample, 
                   features = ATII_ATI_genes, 
                   labels = c('Annotation', 'Sample'), 
                   selected.row.anotations=ATII_ATI_genes,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.3)
dev.off()