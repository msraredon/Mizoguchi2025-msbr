####Merging LEP_DL1-3

load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/LEP_DL/LEP_DL1/LEP_DL1.updated.2025.1.16.Robj")
load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/LEP_DL/LEP_DL2/LEP_DL2.cleaned.2025.1.13.Robj")
load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/LEP_DL/LEP_DL3/LEP_DL3.IT5.2025.1.23.Robj")
DimPlot(LEP_DL1)
DimPlot(LEP_DL2.IT4)
DimPlot(LEP_DL3.IT5)

LEP_DL.merged.2025.1.24 <-merge(LEP_DL1, y = c(LEP_DL2.IT4, LEP_DL3.IT5), add.cell.ids = c("LEP_DL1", "LEP_DL2", "LEP_DL3"), project = "LEP.Merged", merge.data = TRUE)

LEP_DL.merged.2025.1.24 <- FindVariableFeatures(LEP_DL.merged.2025.1.24, verbose = FALSE, nfeatures = 2000)
LEP_DL.merged.2025.1.24 <- ScaleData(LEP_DL.merged.2025.1.24, features = VariableFeatures(LEP_DL.merged.2025.1.24))


LEP_DL.merged.2025.1.24 <- RunPCA(LEP_DL.merged.2025.1.24, npcs = 100)

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/LEP_DL/LEP_DL_Merging/PCA")
png(filename = 'Satoshi.PCA_1.2025.1.21.png',width = 15,height = 15,units = 'in',res=300, type = 'cairo')
DimHeatmap(LEP_DL.merged.2025.1.24, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'Satoshi.PCA_2.2025.1.24.png',width = 15,height = 15,units = 'in',res=300, type = 'cairo')
DimHeatmap(LEP_DL.merged.2025.1.24, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'Satoshi.PCA_3.2025.1.24.png',width = 15,height = 15,units = 'in',res=300, type = 'cairo')
DimHeatmap(LEP_DL.merged.2025.1.24, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'Satoshi.PCA_4.2025.1.24.png',width = 15,height = 15,units = 'in',res=300, type = 'cairo')
DimHeatmap(LEP_DL.merged.2025.1.24, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'Satoshi.PCA_5.2025.1.24.png',width = 15,height = 15,units = 'in',res=300, type = 'cairo')
DimHeatmap(LEP_DL.merged.2025.1.24, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'Satoshi.PCA_6.2025.1.24.png',width = 15,height = 15,units = 'in',res=300, type = 'cairo')
DimHeatmap(LEP_DL.merged.2025.1.24, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'Satoshi.ElbowPlot.2025.1.24.png',width = 15,height = 15,units = 'in',res=300, type = 'cairo')
ElbowPlot(LEP_DL.merged.2025.1.24, ndims = 100)
dev.off()



LEP_DL.merged.2025.1.24 <- FindNeighbors(LEP_DL.merged.2025.1.24, dims = 1:30)
LEP_DL.merged.2025.1.24 <- FindClusters(LEP_DL.merged.2025.1.24, resolution = .5)
LEP_DL.merged.2025.1.24 <- RunUMAP(LEP_DL.merged.2025.1.24, dims = 1:30)

DimPlot(LEP_DL.merged.2025.1.24, label = TRUE)
DimPlot(LEP_DL.merged.2025.1.24, group.by = "Sample", label = TRUE)
FeaturePlot(LEP_DL.merged.2025.1.24, features = "Elf3")
LEP_DL.merged.2025.1.24$updated.cluster.numbers <- Idents(LEP_DL.merged.2025.1.24)

LEP_DL.merged.2025.1.24[['RNA']] <- JoinLayers(LEP_DL.merged.2025.1.24[['RNA']])
