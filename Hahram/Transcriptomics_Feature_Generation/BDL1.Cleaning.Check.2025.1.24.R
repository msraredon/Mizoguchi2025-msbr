###BDL1 Revised Cleaning and Annotation
###In Hindsight, I am not satisfied with how my initial cleaning went on my BDL samples, and have learned
###There are a number of practices that did not follow. I will maintain those in mind when I redo BDL1, 2, and 3


###Generation of the Seurat Object 

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BDL1_Check_2025.1.24")

library(dplyr)
library(Seurat)
library(patchwork)

BDL1.2025.1.24.data <- Read10X(data.dir = "/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BDL1_Check_2025.1.24/Raw_Data/")

BDL1.2025.1.24 <- CreateSeuratObject(counts = BDL1.2025.1.24.data)

# An object of class Seurat 
# 32883 features across 2609068 samples within 1 assay 
# Active assay: RNA (32883 features, 0 variable features)
# 1 layer present: counts


##just to double-check what we're looking at 

BDL1.2025.1.24[["percent.mt"]] <- PercentageFeatureSet(BDL1.2025.1.24, pattern = "^Mt-")

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BDL1_Check_2025.1.24/Plots/")
png("BDL1.2025.1.24_Pre_Thresholding.VlnPlot.png", res=300, unit="in", height=8, width=11)
VlnPlot(BDL1.2025.1.24, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), raster = FALSE)
dev.off()

png("BDL1.2025.1.24_Pre_Thresholding.VlnPlot_Log.png", res=300, unit="in", height=8, width=11)
VlnPlot(BDL1.2025.1.24, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), log = TRUE, raster = FALSE)
dev.off()

plot1 <- FeatureScatter(BDL1.2025.1.24, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BDL1.2025.1.24, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

###Looks fairly trim. We will eliminate any nCount below 500.
### We will also get rid of any mt percent above 25
BDL1.2025.1.24<-subset(BDL1.2025.1.24, subset = nCount_RNA > 500 & percent.mt < 25)

# An object of class Seurat 
# 32883 features across 45435 samples within 1 assay 
# Active assay: RNA (32883 features, 0 variable features)
# 1 layer present: counts

png("BDL1.2025.1.24_VlnPlotPostSubset.png", res=300, unit="in", height=8, width=11)
VlnPlot(BDL1.2025.1.24, features = c("nCount_RNA", "nFeature_RNA"))
dev.off()

png("BDL1.2025.1.24_VlnPlot_LogPostSubset.png", res=300, unit="in", height=8, width=11)
VlnPlot(BDL1.2025.1.24, features = c("nCount_RNA", "nFeature_RNA"), log = TRUE)
dev.off()

## Now proceed to the rest of the cleaning--normalize, scale, and run PCA
BDL1.2025.1.24<- NormalizeData(BDL1.2025.1.24, verbose=FALSE)

BDL1.2025.1.24<- FindVariableFeatures(BDL1.2025.1.24, verbose=FALSE, nfeatures=2000)
BDL1.2025.1.24<- ScaleData(BDL1.2025.1.24, features = VariableFeatures(BDL1.2025.1.24))
BDL1.2025.1.24<- RunPCA(BDL1.2025.1.24, features = VariableFeatures(BDL1.2025.1.24), npcs =100)

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BDL1_Check_2025.1.24/PCA/Iteration_1")

png("BDL1.2025.1.24_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("BDL1.2025.1.24_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()

png(filename = 'BDL1.2025.1.24.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(BDL1.2025.1.24, ndims = 100)
dev.off()

## We will use PCs 

dims.use <- c(1:19, 21, 22, 24, 28:36)



BDL1.2025.1.24 <- FindNeighbors(BDL1.2025.1.24,
                                reduction="pca",
                                k.param = 20,
                                dims=dims.use
                                
)

BDL1.2025.1.24 <- FindClusters(BDL1.2025.1.24,
                               resolution= 2
                               
)

BDL1.2025.1.24<- RunUMAP(BDL1.2025.1.24, 
                         dims=dims.use,
                         min.dist = .2
)
DimPlot(BDL1.2025.1.24, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE)

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BDL1_Check_2025.1.24/Plots/UMAP/Iteration_1/")

png("BDL1.2025.1.24_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(BDL1.2025.1.24, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("BDL1.2025.1.24_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BDL1.2025.1.24, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BDL1.2025.1.24_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BDL1.2025.1.24, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BDL1.2025.1.24_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BDL1.2025.1.24, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BDL1_Check_2025.1.24/Plots/Iteration_1/")

png("BDL1.2025.1.24_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BDL1.2025.1.24, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("BDL1.2025.1.24_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BDL1.2025.1.24, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("BDL1.2025.1.24_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(BDL1.2025.1.24, features = "percent.mt", split.by = "seurat_clusters")
dev.off()

BDL1.2025.1.24.Markers <- FindAllMarkers(BDL1.2025.1.24, logfc.threshold = 0.25, only.pos = T)
BDL1.2025.1.24.Markers$ratio <- BDL1.2025.1.24.Markers$pct.1/BDL1.2025.1.24.Markers$pct.2
BDL1.2025.1.24.Markers$power <- BDL1.2025.1.24.Markers$ratio*BDL1.2025.1.24.Markers$avg_log2FC

######RESUME HERE########
BDL1.2025.1.24.IT2<- subset(BDL1.2025.1.24, idents = c(0, 1, 2, 3, 4, 6, 7), invert = T)

# An object of class Seurat 
# 32883 features across 19979 samples within 1 assay 
# Active assay: RNA (32883 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

BDL1.2025.1.24.IT2<- FindVariableFeatures(BDL1.2025.1.24.IT2, verbose=FALSE, nfeatures=2000)
BDL1.2025.1.24.IT2<- ScaleData(BDL1.2025.1.24.IT2, features = VariableFeatures(BDL1.2025.1.24.IT2))
BDL1.2025.1.24.IT2<- RunPCA(BDL1.2025.1.24.IT2, features = VariableFeatures(BDL1.2025.1.24.IT2), npcs = 100)

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BDL1_Check_2025.1.24/PCA/Iteration_2")

png("BDL1.2025.1.24_PCA1.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24.IT2, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("BDL1.2025.1.24_PCA2.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24.IT2, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA3.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24.IT2, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA4.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24.IT2, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA5.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24.IT2, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA6.png", res=300, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24.IT2, dims = 46:54, cells = 500, balanced = TRUE)
dev.off()


png(filename = 'BDL1.2025.1.24.ElbowPlot.png', res=300, unit="in", height=8, width=11)
ElbowPlot(BDL1.2025.1.24.IT2, ndims = 100)
dev.off()

## We will use dimensions 1:14, 16:18, 20, 21, 23:30, 32:34
dims.use <- c(1:14, 16:18, 20, 21, 23:30, 32:34)



BDL1.2025.1.24.IT2 <- FindNeighbors(BDL1.2025.1.24.IT2,
                                    reduction="pca",
                                    k.param = 20,
                                    dims=dims.use
                                    
)

BDL1.2025.1.24.IT2 <- FindClusters(BDL1.2025.1.24.IT2,
                                   resolution= 2
                                   
)
BDL1.2025.1.24.IT2<- RunUMAP(BDL1.2025.1.24.IT2, 
                             dims=dims.use,
                             min.dist = .2
)
DimPlot(BDL1.2025.1.24.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=FALSE)


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BDL1_Check_2025.1.24/Plots/UMAP/Iteration_2/")

png("BDL1.2025.1.24_UMAP_1.png", res=300, unit="in", height=8, width=11)
DimPlot(BDL1.2025.1.24.IT2, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE, repel = TRUE) + NoLegend()
dev.off()

png("BDL1.2025.1.24_nCount_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BDL1.2025.1.24.IT2, features = "nCount_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BDL1.2025.1.24_nFeature_RNA.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BDL1.2025.1.24.IT2, features = "nFeature_RNA", label = TRUE, repel = TRUE)
dev.off()

png("BDL1.2025.1.24_percent.mt.png", res = 300, unit = "in", height=8, width=11)
FeaturePlot(BDL1.2025.1.24.IT2, features = "percent.mt", label = TRUE, repel = TRUE)
dev.off()

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BDL1_Check_2025.1.24/Plots/Iteration_2/")

png("BDL1.2025.1.24_cluster_VlnPlot.nCount_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BDL1.2025.1.24.IT2, features = "nCount_RNA", split.by = "seurat_clusters")
dev.off()

png("BDL1.2025.1.24_cluster_VlnPlot.nFeature_RNA.png", res=300, unit="in", height=8, width=11)
VlnPlot(BDL1.2025.1.24.IT2, features = "nFeature_RNA", split.by = "seurat_clusters")
dev.off()

png("BDL1.2025.1.24_cluster_VlnPlot.percent.mt.png", res=300, unit="in", height=8, width=11)
VlnPlot(BDL1.2025.1.24.IT2, features = "percent.mt", split.by = "seurat_clusters")
dev.off()

BDL1.2025.1.24.IT2.Markers <- FindAllMarkers(BDL1.2025.1.24.IT2, logfc.threshold = 0.25, only.pos = T)
BDL1.2025.1.24.IT2.Markers$ratio <- BDL1.2025.1.24.IT2.Markers$pct.1/BDL1.2025.1.24.IT2.Markers$pct.2
BDL1.2025.1.24.IT2.Markers$power <- BDL1.2025.1.24.IT2.Markers$ratio*BDL1.2025.1.24.IT2.Markers$avg_log2FC


####Iteration 3

BDL1.2025.1.24<- subset(BDL1.2025.1.24, idents = c("12", "17", "19"), invert = T)

## The object is currently sitting at 32883 features across 7246 samples within 1 assay 

BDL1.2025.1.24<- NormalizeData(BDL1.2025.1.24, verbose=FALSE)
BDL1.2025.1.24<- FindVariableFeatures(BDL1.2025.1.24, verbose=FALSE, nfeatures=2000)
all.genes <- rownames(BDL1.2025.1.24)
BDL1.2025.1.24<- ScaleData(BDL1.2025.1.24, features = VariableFeatures(BDL1.2025.1.24))
BDL1.2025.1.24<- RunPCA(BDL1.2025.1.24, features = VariableFeatures(BDL1.2025.1.24))

setwd("/home/hk738/Satoshi_Data/BASC_Ventilated_Retest/PCA/PCA_3/")

png("BDL1.2025.1.24_PCA1.png", res=200, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("BDL1.2025.1.24_PCA2.png", res=200, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA3.png", res=200, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA4.png", res=200, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA5.png", res=200, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

dims.use <- c(1:7, 9:10, 12, 14:16, 18, 20:25)

BDL1.2025.1.24<- RunUMAP(BDL1.2025.1.24, 
                         dims=dims.use,
                         min.dist = .2
)

DimPlot(BDL1.2025.1.24, shuffle=TRUE, raster=FALSE, pt.size=1, label=FALSE)

BDL1.2025.1.24 <- FindNeighbors(BDL1.2025.1.24,
                                reduction="pca",
                                k.param = 20,
                                dims=dims.use
                                
)

BDL1.2025.1.24 <- FindClusters(BDL1.2025.1.24,
                               resolution= 2
                               
)

setwd("/home/hk738/Satoshi_Data/BASC_Ventilated_Retest/UMAP/UMAP_3/")

png("BDL1.2025.1.24_UMAP_3.png", res=200, unit="in", height=8, width=11)
DimPlot(BDL1.2025.1.24, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE) + NoLegend()
dev.off()

setwd("/home/hk738/Satoshi_Data/BASC_Ventilated_Retest/Plots/Iteration_3/")

png("BDL1.2025.1.24_cluster_VlnPlot.png", res=200, unit="in", height=8, width=11)
vplot1 <- VlnPlot(BDL1.2025.1.24, features = "nCount_RNA", split.by = "seurat_clusters")
vplot2  <- VlnPlot(BDL1.2025.1.24, features = "nFeature_RNA", split.by = "seurat_clusters")
vplot3 <- VlnPlot(BDL1.2025.1.24, features = "percent.mt", split.by = "seurat_clusters")
cowplot::plot_grid(vplot1, vplot2, vplot3, nrow = 3)
dev.off()

png("BDL1.2025.1.24_featureplot.png", res = 200, unit = "in", height=8, width=11)
FeaturePlot(BDL1.2025.1.24, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), label = TRUE)
dev.off()

#ADJUSTMENT MADE 3/1/2024: CLUSTER 4 GIVES LITTLE INFO, WILL REMOVE

BDL1.2025.1.24<- subset(BDL1.2025.1.24, idents = c("4"), invert = T)

## The object is currently sitting at 32883 features across 6523 samples within 1 assay 

BDL1.2025.1.24<- NormalizeData(BDL1.2025.1.24, verbose=FALSE)
BDL1.2025.1.24<- FindVariableFeatures(BDL1.2025.1.24, verbose=FALSE, nfeatures=2000)
all.genes <- rownames(BDL1.2025.1.24)
BDL1.2025.1.24<- ScaleData(BDL1.2025.1.24, features = VariableFeatures(BDL1.2025.1.24))
BDL1.2025.1.24<- RunPCA(BDL1.2025.1.24, features = VariableFeatures(BDL1.2025.1.24))

setwd("/home/hk738/Satoshi_Data/BASC_Ventilated_Retest/PCA/PCA_4/")

png("BDL1.2025.1.24_PCA1.png", res=200, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 1:9, cells = 500, balanced = TRUE)
dev.off() 

png("BDL1.2025.1.24_PCA2.png", res=200, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA3.png", res=200, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA4.png", res=200, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 28:36, cells = 500, balanced = TRUE)
dev.off()

png("BDL1.2025.1.24_PCA5.png", res=200, unit="in", height=8, width=11)
DimHeatmap(BDL1.2025.1.24, dims = 37:45, cells = 500, balanced = TRUE)
dev.off()

dims.use <- c(1:5, 7:9, 10:12, 14:16, 18:19, 20:26)

BDL1.2025.1.24<- RunUMAP(BDL1.2025.1.24, 
                         dims=dims.use,
                         n.neighbors = 20,
                         min.dist = .15
)

DimPlot(BDL1.2025.1.24, shuffle=TRUE, raster=FALSE, pt.size=1, label=FALSE)

BDL1.2025.1.24 <- FindNeighbors(BDL1.2025.1.24,
                                reduction="pca",
                                k.param = 20,
                                dims=dims.use
                                
)

BDL1.2025.1.24 <- FindClusters(BDL1.2025.1.24,
                               resolution= 2
                               
)

setwd("/home/hk738/Satoshi_Data/BASC_Ventilated_Retest/UMAP/UMAP_4/")

png("BDL1.2025.1.24_UMAP_3.png", res=200, unit="in", height=8, width=11)
DimPlot(BDL1.2025.1.24, shuffle=TRUE, raster=FALSE, pt.size=1, label=TRUE) + NoLegend()
dev.off()

setwd("/home/hk738/Satoshi_Data/BASC_Ventilated_Retest/Plots/Iteration_4/")

png("BDL1.2025.1.24_cluster_VlnPlot.png", res=200, unit="in", height=8, width=11)
vplot1 <- VlnPlot(BDL1.2025.1.24, features = "nCount_RNA", split.by = "seurat_clusters")
vplot2  <- VlnPlot(BDL1.2025.1.24, features = "nFeature_RNA", split.by = "seurat_clusters")
vplot3 <- VlnPlot(BDL1.2025.1.24, features = "percent.mt", split.by = "seurat_clusters")
cowplot::plot_grid(vplot1, vplot2, vplot3, nrow = 3)
dev.off()

png("BDL1.2025.1.24_featureplot.png", res = 200, unit = "in", height=8, width=11)
FeaturePlot(BDL1.2025.1.24, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), label = TRUE)
dev.off()

##Featureplot looks satisfactory. Will find markers to determine cluster number

BDL1.2025.1.24 <- FindClusters(BDL1.2025.1.24,
                               resolution= .2
                               
)
setwd("/home/hk738/Satoshi_Data/BASC_Ventilated_Retest")

#save(BDL1.2025.1.24, file = "BDL1.2025.1.24.2024.2.27.Robj")

save(BDL1.2025.1.24, file = "BDL1.2025.1.24.2024.3.1.Robj")

setwd("/home/hk738/Satoshi_Data/BASC_Ventilated_Retest/Markers/")

BDL1.2025.1.24.markers <- FindAllMarkers(BDL1.2025.1.24, only.pos=TRUE, logfc.threshold=0.25)
BDL1.2025.1.24.markers$ratio <-BDL1.2025.1.24.markers$pct.1/BDL1.2025.1.24.markers$pct.2
BDL1.2025.1.24.markers$power <- BDL1.2025.1.24.markers$avg_log2FC*BDL1.2025.1.24.markers$ratio

write.table(BDL1.2025.1.24.markers, file="BDL1.markers.2024.3.1.txt", sep="\t", row.names=FALSE)

FeaturePlot(BDL1.2025.1.24, features = "Krt16")

###Exploring the cluster markers
setwd("/home/hk738/Satoshi_Data/BASC_Ventilated_Retest/Features/")

###Cluster 2

png("BDL1.2025.1.24_featureplot_cluster2.png", res = 200, unit = "in", height=8, width=11)
k1 <- FeaturePlot(BDL1.2025.1.24, features = "Aldh1a1", label = T)
k2 <- FeaturePlot(BDL1.2025.1.24, features = "Ms4a8", label = T)
k3 <- FeaturePlot(BDL1.2025.1.24, features = "Abca1", label = T)
k4 <- FeaturePlot(BDL1.2025.1.24, features = "MGC105649", label = T)
k5 <- FeaturePlot(BDL1.2025.1.24, features = "Mal", label = T)
k6 <- FeaturePlot(BDL1.2025.1.24, features = "Scgb1a1", label = T)
k7 <- FeaturePlot(BDL1.2025.1.24, features = "Tgfb2", label = T)
k8 <- FeaturePlot(BDL1.2025.1.24, features = "Krt19", label = T)
k9 <- FeaturePlot(BDL1.2025.1.24, features = "Pon3", label = T)
k10 <- FeaturePlot(BDL1.2025.1.24, features = "Ctla2a", label = T)
k11 <- FeaturePlot(BDL1.2025.1.24, features = "Elf3", label = T)
k12 <- FeaturePlot(BDL1.2025.1.24, features = "Cysrt1", label = T)

cowplot::plot_grid(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12)
dev.off()

#Cluster 1
clusterzero <- FindMarkers(BDL1.2025.1.24, ident.1 = "0", ident.2 = "1")
clusterzero$ratio <-clusterzero$pct.1/clusterzero$pct.2
clusterzero$power <- clusterzero$avg_log2FC*clusterzero$ratio

png("BDL1.2025.1.24_featureplot_cluster0.png", res = 200, unit = "in", height=8, width=11)
k1 <- FeaturePlot(BDL1.2025.1.24, features = "Ccne2", label = T)
k2 <- FeaturePlot(BDL1.2025.1.24, features = "Cdc6", label = T)
k3 <- FeaturePlot(BDL1.2025.1.24, features = "Dtl", label = T)
k4 <- FeaturePlot(BDL1.2025.1.24, features = "Ung", label = T)
k5 <- FeaturePlot(BDL1.2025.1.24, features = "Kif20b", label = T)
k6 <- FeaturePlot(BDL1.2025.1.24, features = "Kif11", label = T)
k7 <- FeaturePlot(BDL1.2025.1.24, features = "Kifc1", label = T)
k8 <- FeaturePlot(BDL1.2025.1.24, features = "Spc25", label = T)
k9 <- FeaturePlot(BDL1.2025.1.24, features = "Aurkb", label = T)
k10 <- FeaturePlot(BDL1.2025.1.24, features = "Incenp", label = T)
#k11 <- FeaturePlot(BDL1.2025.1.24, features = "Tpx2", label = T)
#k12 <- FeaturePlot(BDL1.2025.1.24, features = "Cenpf", label = T)

cowplot::plot_grid(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10)
dev.off()

#Cluster 1

png("BDL1.2025.1.24_featureplot_cluster1.png", res = 200, unit = "in", height=8, width=11)
k1 <- FeaturePlot(BDL1.2025.1.24, features = "Top2a", label = T)
k2 <- FeaturePlot(BDL1.2025.1.24, features = "Nusap1", label = T)
k3 <- FeaturePlot(BDL1.2025.1.24, features = "Aurka", label = T)
k4 <- FeaturePlot(BDL1.2025.1.24, features = "Plk1", label = T)
k5 <- FeaturePlot(BDL1.2025.1.24, features = "Kif20b", label = T)
k6 <- FeaturePlot(BDL1.2025.1.24, features = "Kif11", label = T)
k7 <- FeaturePlot(BDL1.2025.1.24, features = "Kifc1", label = T)
k8 <- FeaturePlot(BDL1.2025.1.24, features = "Spc25", label = T)
k9 <- FeaturePlot(BDL1.2025.1.24, features = "Aurkb", label = T)
k10 <- FeaturePlot(BDL1.2025.1.24, features = "Incenp", label = T)
#k11 <- FeaturePlot(BDL1.2025.1.24, features = "Tpx2", label = T)
#k12 <- FeaturePlot(BDL1.2025.1.24, features = "Cenpf", label = T)

cowplot::plot_grid(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10)
dev.off()

#Cluster 3

#Cluster 3 looks extremely strange. Can we separate it out?
cluster3<-FindMarkers(BDL1.2025.1.24, ident.1 = "3", ident.2 = c("0", "2"))
cluster3$ratio <-cluster3$pct.1/cluster3$pct.2
cluster3$power <- cluster3$avg_log2FC*cluster3$ratio

png("BDL1.2025.1.24_featureplot_cluster3.png", res = 200, unit = "in", height=8, width=11)
k1 <- FeaturePlot(BDL1.2025.1.24, features = "Tubb3")
k2 <- FeaturePlot(BDL1.2025.1.24, features = "Mrpl43")
k3 <- FeaturePlot(BDL1.2025.1.24, features = "Crabp2")
k4 <- FeaturePlot(BDL1.2025.1.24, features = "Cd63")
k5 <- FeaturePlot(BDL1.2025.1.24, features = "Dsc3")
k6 <- FeaturePlot(BDL1.2025.1.24, features = "Scgb1a1")
k7 <- FeaturePlot(BDL1.2025.1.24, features = "Tgfb2")
k8 <- FeaturePlot(BDL1.2025.1.24, features = "Krt19")
k9 <- FeaturePlot(BDL1.2025.1.24, features = "Pon3")
k10 <- FeaturePlot(BDL1.2025.1.24, features = "Ctla2a")
k11 <- FeaturePlot(BDL1.2025.1.24, features = "Elf3")
k12 <- FeaturePlot(BDL1.2025.1.24, features = "Cysrt1")

cowplot::plot_grid(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12)
dev.off()
