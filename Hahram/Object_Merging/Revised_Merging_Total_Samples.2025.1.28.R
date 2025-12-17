BDL.merged.2025.1.28 <-merge(BDL1, y = c(BDL2, BDL3), add.cell.ids = c("BDL1", "BDL2", "BDL3"), project = "BDL.Merged", merge.data = TRUE)

BDL.merged.2025.1.28 <- FindVariableFeatures(BDL.merged.2025.1.28, verbose = FALSE, nfeatures = 2000)
BDL.merged.2025.1.28 <- ScaleData(BDL.merged.2025.1.28, features = VariableFeatures(BDL.merged.2025.1.28))


BDL.merged.2025.1.28 <- RunPCA(BDL.merged.2025.1.28, 
                               features = VariableFeatures(BDL.merged.2025.1.28), 
                               npcs = 100)

ElbowPlot(BDL.merged.2025.1.28, ndims = 100)

BDL.merged.2025.1.28 <- FindNeighbors(BDL.merged.2025.1.28, dims = 1:45)
BDL.merged.2025.1.28 <- FindClusters(BDL.merged.2025.1.28, resolution = .5)
BDL.merged.2025.1.28 <- RunUMAP(BDL.merged.2025.1.28, dims = 1:45)

DimPlot(BDL.merged.2025.1.28)

DimPlot(BDL.merged.2025.1.28, group.by = "prelim.annotations")
DimPlot(BDL.merged.2025.1.28, group.by = "Condition")
DimPlot(BDL.merged.2025.1.28, group.by = "Sample")

FeaturePlot(BDL.merged.2025.1.28, features = "Top2a")
FeaturePlot(BDL.merged.2025.1.28, features = "Krt13")

BDL.merged.2025.1.28$annotations <- BDL.merged.2025.1.28$prelim.annotations
save(BDL.merged.2025.1.28, file = "BDL.merged.2025.1.28.Robj")


#####Now we can merge with everything
load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2024.11.22/BASC_St.P0.P3_combined.updated.annotations.2024.11.22.Robj")

load("~/project/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2024.12.3/alltimepoints_allclasses (1).Robj")

Idents(frl.all) <-frl.all$Condition


frl.all.normal <- subset(frl.all, idents = "Normal")

frl.all.normal <- FindVariableFeatures(frl.all.normal, verbose = FALSE, nfeatures = 2000)
frl.all.normal <- ScaleData(frl.all.normal, features = VariableFeatures(frl.all.normal))
frl.all.normal <- RunPCA(frl.all.normal,features = VariableFeatures(frl.all.normal), npcs = 100)
DimHeatmap(frl.all.normal, dims = 37:45, cells = 500, balanced = TRUE)
DimHeatmap(frl.all.normal, dims = 46:52, cells = 500, balanced = TRUE)
DimHeatmap(frl.all.normal, dims = 37:45, cells = 500, balanced = TRUE)

ElbowPlot(frl.all.normal, ndims = 100)
frl.all.normal <- FindNeighbors(frl.all.normal, dims = 1:45)
frl.all.normal <- FindClusters(frl.all.normal, resolution = .5)
frl.all.normal <- RunUMAP(frl.all.normal, dims = 1:45)

DimPlot(frl.all.normal, label = TRUE)
DimPlot(frl.all.normal, label = TRUE, group.by = "Sample", repel = T)
DimPlot(frl.all.normal, label = TRUE, group.by = "CellClass", repel = T)

unique(frl.all.normal$GrossCellType)
frl.all.normal$annotations <- frl.all.normal$GrossCellType
frl.all.normal$Condition <- "Embryogenesis"

Idents(frl.all.normal) <- frl.all.normal$Sample
frl.all.normal.noNRL21 <- subset(frl.all.normal, idents = "NRL21", invert = T)

frl.all.normal.noNRL21 <- FindVariableFeatures(frl.all.normal.noNRL21, verbose = FALSE, nfeatures = 2000)
frl.all.normal.noNRL21 <- ScaleData(frl.all.normal.noNRL21, features = VariableFeatures(frl.all.normal.noNRL21))
frl.all.normal.noNRL21 <- RunPCA(frl.all.normal.noNRL21, 
                                 features = VariableFeatures(frl.all.normal.noNRL21), 
                                 npcs = 100)


ElbowPlot(frl.all.normal.noNRL21, ndims = 100)
frl.all.normal.noNRL21 <- FindNeighbors(frl.all.normal.noNRL21, dims = 1:45)
frl.all.normal.noNRL21 <- FindClusters(frl.all.normal.noNRL21, resolution = .5)
frl.all.normal.noNRL21 <- RunUMAP(frl.all.normal.noNRL21, dims = 1:45)

DimPlot(frl.all.normal.noNRL21, group.by = "Condition")

frl.all.normal.noNRL21$Condition <- "Embryogenesis"
DimPlot(frl.all.normal.noNRL21, group.by = "CellClass")
frl.all.normal.noNRL21$annotations <- frl.all.normal.noNRL21$GrossCellType

####################Combining Embryogenesis, Start, P0, P3, and 3D
Satoshi.Rachel.Combined.2025.1.28 <-merge(BDL.merged.2025.1.28, y = c(store, frl.all.normal.noNRL21), add.cell.ids = c("BDL", "Epi_2D", "Embryo"), project = "Engineered.Merged", merge.data = TRUE)

Satoshi.Rachel.Combined.2025.1.28 <- FindVariableFeatures(Satoshi.Rachel.Combined.2025.1.28, verbose = FALSE, nfeatures = 2000)
Satoshi.Rachel.Combined.2025.1.28 <- ScaleData(Satoshi.Rachel.Combined.2025.1.28, features = VariableFeatures(Satoshi.Rachel.Combined.2025.1.28))


Satoshi.Rachel.Combined.2025.1.28 <- RunPCA(Satoshi.Rachel.Combined.2025.1.28, 
                               features = VariableFeatures(Satoshi.Rachel.Combined.2025.1.28), 
                               npcs = 100)

DimHeatmap(Satoshi.Rachel.Combined.2025.1.28, dims = 28:36, cells = 500, balanced = TRUE)
DimHeatmap(Satoshi.Rachel.Combined.2025.1.28, dims = 37:45, cells = 500, balanced = TRUE)
DimHeatmap(Satoshi.Rachel.Combined.2025.1.28, dims = 46:54, cells = 500, balanced = TRUE)
DimHeatmap(Satoshi.Rachel.Combined.2025.1.28, dims = 55:63, cells = 500, balanced = TRUE)
DimHeatmap(Satoshi.Rachel.Combined.2025.1.28, dims = 64:72, cells = 500, balanced = TRUE)

ElbowPlot(Satoshi.Rachel.Combined.2025.1.28, ndims = 100)

Satoshi.Rachel.Combined.2025.1.28 <- FindNeighbors(Satoshi.Rachel.Combined.2025.1.28, dims = 1:45)
Satoshi.Rachel.Combined.2025.1.28 <- FindClusters(Satoshi.Rachel.Combined.2025.1.28, resolution = .5)
Satoshi.Rachel.Combined.2025.1.28 <- RunUMAP(Satoshi.Rachel.Combined.2025.1.28, dims = 1:45)

DimPlot(Satoshi.Rachel.Combined.2025.1.28, group.by = "Sample")

#####Removing all non-epithelial
Idents(Satoshi.Rachel.Combined.2025.1.28) <- Satoshi.Rachel.Combined.2025.1.28$CellClass
Satoshi.Rachel.Combined.2025.1.28.EpiOnly <- subset(Satoshi.Rachel.Combined.2025.1.28, idents = c("Endothelial", "Endothelium", "Immune", "Mesenchymal", "Mesenchyme"), invert = T)

Satoshi.Rachel.Combined.2025.1.28.EpiOnly <- FindVariableFeatures(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, verbose = FALSE, nfeatures = 2000)
Satoshi.Rachel.Combined.2025.1.28.EpiOnly <- ScaleData(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, features = VariableFeatures(Satoshi.Rachel.Combined.2025.1.28.EpiOnly))
Satoshi.Rachel.Combined.2025.1.28.EpiOnly <- RunPCA(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, 
                                 features = VariableFeatures(Satoshi.Rachel.Combined.2025.1.28.EpiOnly), 
                                 npcs = 100)


ElbowPlot(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, ndims = 100)
Satoshi.Rachel.Combined.2025.1.28.EpiOnly <- FindNeighbors(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, dims = 1:45)
Satoshi.Rachel.Combined.2025.1.28.EpiOnly <- FindClusters(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, resolution = .5)
Satoshi.Rachel.Combined.2025.1.28.EpiOnly <- RunUMAP(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, dims = 1:45)

DimPlot(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, group.by = "Condition")
DimPlot(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, group.by = "CellClass")
DimPlot(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, group.by = "annotations")

Satoshi.Rachel.Combined.2025.1.28.EpiOnly$CellClassUpdate <- "Epithelial"
save(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, file = "Satoshi.Rachel.Combined.2025.1.28.EpiOnly.Robj")

LEP_DL.merged.2025.1.28$Condition <- "Engineered-FiveDay"
########Merging With Combined 5-DAY
Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate <-merge(Satoshi.Rachel.Combined.2025.1.28.EpiOnly, y = LEP_DL.merged.2025.1.28, add.cell.ids = c("Embryo.Start.P0.P3.3D", "Five-Day"), project = "Satoshi.Rachel.Merged.Updated", merge.data = TRUE)

Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate <- FindVariableFeatures(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, verbose = FALSE, nfeatures = 2000)
Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate <- ScaleData(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, features = VariableFeatures(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate))


Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate <- RunPCA(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, 
                               features = VariableFeatures(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate), 
                               npcs = 100)

ElbowPlot(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, ndims = 100)

Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate <- FindNeighbors(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, dims = 1:45)
Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate <- FindClusters(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, resolution = .5)
Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate <- RunUMAP(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, dims = 1:45)

DimPlot(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate)


DimPlot(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, group.by = "Condition", label = T)
DimPlot(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, group.by = "Sample", label = T)
DimPlot(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, group.by = "annotations", label = T)

save(Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate, file = "Satoshi.Rachel.Combined.2025.1.28.EpiOnly.5dayUpdate.Robj")
