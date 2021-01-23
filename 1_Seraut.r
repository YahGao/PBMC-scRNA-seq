setwd("/yhgao/pbmc")

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)

CO.data <- Read10X(data.dir = "/yhgao/pbmc/CO/")
T1.data <- Read10X(data.dir = "/yhgao/pbmc/T1/")
T2.data <- Read10X(data.dir = "/yhgao/pbmc/T2/")
T3.data <- Read10X(data.dir = "/yhgao/pbmc/T3/")

CO <- CreateSeuratObject(counts = CO.data, project = "CO", min.cells = 3, min.features = 200)
CO[["percent.mt"]] <- PercentageFeatureSet(CO, pattern = "^MT-")
CO <- subset(CO, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 10000 & percent.mt < 20)
CO <- NormalizeData(CO, normalization.method = "LogNormalize", scale.factor = 10000)
CO <- FindVariableFeatures(CO, selection.method = "vst", nfeatures = 2000)
T1 <- CreateSeuratObject(counts = T1.data, project = "T1", min.cells = 3, min.features = 200)
T1[["percent.mt"]] <- PercentageFeatureSet(T1, pattern = "^MT-")
T1 <- subset(T1, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 10000 & percent.mt < 20)
T1 <- NormalizeData(T1, normalization.method = "LogNormalize", scale.factor = 10000)
T1 <- FindVariableFeatures(T1, selection.method = "vst", nfeatures = 2000)
T2 <- CreateSeuratObject(counts = T2.data, project = "T2", min.cells = 3, min.features = 200)
T2[["percent.mt"]] <- PercentageFeatureSet(T2, pattern = "^MT-")
T2 <- subset(T2, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 10000 & percent.mt < 20)
T2 <- NormalizeData(T2, normalization.method = "LogNormalize", scale.factor = 10000)
T2 <- FindVariableFeatures(T2, selection.method = "vst", nfeatures = 2000)
T3 <- CreateSeuratObject(counts = T3.data, project = "T3", min.cells = 3, min.features = 200)
T3[["percent.mt"]] <- PercentageFeatureSet(T3, pattern = "^MT-")
T3 <- subset(T3, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 10000 & percent.mt < 20)
T3 <- NormalizeData(T3, normalization.method = "LogNormalize", scale.factor = 10000)
T3 <- FindVariableFeatures(T3, selection.method = "vst", nfeatures = 2000)

cattle <- merge(CO, y = c(T1, T2, T3), add.cell.ids = c("CO","T1", "T2", "T3"), project = "pbmc")
cattle.list <- SplitObject(cattle, split.by = "orig.ident")
cattle.list <- lapply(X = cattle.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
cattle.anchors <- FindIntegrationAnchors(object.list = cattle.list, dims = 1:20)
cattle.combined <- IntegrateData(anchorset = cattle.anchors, dims = 1:20)
DefaultAssay(cattle.combined) <- "integrated"
cattle.combined <- ScaleData(cattle.combined, verbose = FALSE)
cattle.combined <- RunPCA(cattle.combined, npcs = 30, verbose = FALSE)
cattle.combined <- JackStraw(cattle.combined, num.replicate = 100)
cattle.combined <- ScoreJackStraw(cattle.combined, dims = 1:20)

JackStrawPlot(cattle.combined, dims = 1:15) +
  theme(legend.text = element_text(size=8),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size=9,face="bold"),
        axis.title.y = element_text(size=9,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=0.9)))
ElbowPlot(cattle.combined) 

cattle.combined <- RunUMAP(cattle.combined, reduction = "pca", dims = 1:12)
cattle.combined <- FindNeighbors(cattle.combined, reduction = "pca", dims = 1:12)
cattle.combined <- FindClusters(cattle.combined, resolution = 0.05)
saveRDS(cattle.combined, "cattle_tutorial_new.rds")

# Fig. S1a
tiff(file="Fig. S1a.tiff", width = 10, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(cattle.combined, reduction = "umap", group.by = "orig.ident", pt.size = 0.5) +
  theme(legend.text = element_text(size=8),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width=unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size=7,face="bold"),
        axis.title.y = element_text(size=7,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=0.9)))  
dev.off()

# Fig. S1b
tiff(file="Fig. S1b.tiff", width = 14, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(cattle.combined, reduction = "umap", split.by = "orig.ident", pt.size = 0.5, ncol=2) +
  theme(legend.text = element_text(size=8),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size=7,face="bold"),
        axis.title.y = element_text(size=7,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=0.9)))
dev.off()

# Fig. 1a
new.cluster.ids <- c("C0: T", "C1: B", "C2: Monocytes", "C3: HSC", "C4: NK", "C5: Macrophages", "C6: DC")
names(new.cluster.ids) <- levels(cattle.combined)
cattle.combined <- RenameIdents(cattle.combined, new.cluster.ids)
tiff(file="Fig. 1a.tiff", width = 8.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(cattle.combined, reduction = "umap", label = T, pt.size = 0.5) +
  theme(legend.text = element_text(size=8),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size=7,face="bold"),
        axis.title.y = element_text(size=7,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=0.9))) + NoLegend()
dev.off()

# Fig. S1b
write.table(cattle.combined$seurat_clusters, "cattle_combined_cluster.txt", quote = F, col.names = F, sep = "\t")
cattle.combined1 <- read.table("cattle_combined_cluster1.txt", header = T, row.names = 1)
colnames(cattle.combined1) <- c("CO", "T1", "T2", "T3")
tiff(file="Fig. 1b.tiff", width = 4, height = 7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
barplot(as.matrix(100*cattle.combined1), legend = rownames(cattle.combined1), border = NA,
        col = c('#f8766d','#c49a00','#53b400','#00c094','#00b6eb','#a58aff','#fb61d7'),
        cex.axis = 1, cex.names = 0.7, ylim = c(0, 100),
        las = 1, width = 0.8, space = 0.25, beside = FALSE, 
        args.legend = list(x = 'right', bty = 'n', inset = -0.73, cex = 1, y.intersp = -0.8, 
                           x.intersp = 0.7, text.width = 1.8, ncol=1))
mtext('Percentage of cells (%)', cex = 1, side = 2, line = 2.3)
dev.off()

# Fig. S1c & # Fig. S2a
c0_marker <- c("CD4", "CD5", "LEF1")
tiff(file="c0_marker_fp.tiff", width = 25.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c0_marker, max.cutoff = 3, ncol = 3)
dev.off()
tiff(file="c0_marker_vp.tiff", width = 25.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c0_marker, pt.size = 0, combine = T, ncol = 3)
dev.off()
c1_marker <- c("MS4A1", "CD74", "CD79B")
tiff(file="c1_marker_fp.tiff", width = 25.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c1_marker, max.cutoff = 3, ncol = 3)
dev.off()
tiff(file="c1_marker_vp.tiff", width = 25.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c1_marker, pt.size = 0, combine = T, ncol = 3)
dev.off()
c2_marker <- c("CD14", "IL1B", "S100A12")
tiff(file="c2_marker_fp.tiff", width = 25.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c2_marker, max.cutoff = 3, ncol = 3)
dev.off()
tiff(file="c2_marker_vp.tiff", width = 25.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c2_marker, pt.size = 0, combine = T, ncol = 3)
dev.off()
c3_marker <- c("GATA3", "CD48")
tiff(file="c3_marker_fp.tiff", width = 17, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c3_marker, max.cutoff = 3, ncol = 2)
dev.off()
tiff(file="c3_marker_vp.tiff", width = 17, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c3_marker, pt.size = 0, combine = T, ncol = 2)
dev.off()
c4_marker <- c("GNLY", "NKG7", "CTSW", "PRF1", "IL2RB")
tiff(file="c4_marker_fp.tiff", width = 42.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c4_marker, max.cutoff = 3, ncol = 5)
dev.off()
tiff(file="c4_marker_vp.tiff", width = 42.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c4_marker, pt.size = 0, combine = T, ncol = 5)
dev.off()
c5_marker <- c("DUSP5", "TNF")
tiff(file="c5_marker_fp.tiff", width = 17, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c5_marker, max.cutoff = 3, ncol = 2)
dev.off()
tiff(file="c5_marker_vp.tiff", width = 17, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c5_marker, pt.size = 0, combine = T, ncol = 2)
dev.off()
c6_marker <- c("FCER1A", "NR4A3")
tiff(file="c6_marker_fp.tiff", width = 17, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c6_marker, max.cutoff = 3, ncol = 2)
dev.off()
tiff(file="c6_marker_vp.tiff", width = 17, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c6_marker, pt.size = 0, combine = T)
dev.off()

# Fig. 1d
new.cluster.ids <- c(paste("Treg\n","TEM\n","TCM\n","CD4+\n","CD8+",sep=""), paste("Naive B\n","Memory B",sep=""), paste("CD14 Mono\n","CD16 Mono\n","Neutrophils",sep=""), "HSC", "NK", "Macrophages", "DC")
names(new.cluster.ids) <- levels(cattle.combined)
cattle.combined <- RenameIdents(cattle.combined, new.cluster.ids)
tiff(file="Fig. 1d.tiff", width = 8.5, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(cattle.combined, reduction = "umap", label = T, pt.size = 0.5) +
  theme(legend.text = element_text(size=8),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size=7,face="bold"),
        axis.title.y = element_text(size=7,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=0.9))) + NoLegend()
dev.off()


# Fig. S2b
cattle.combined <- ScaleData(cattle.combined, verbose = FALSE, assay = "RNA")
tiff(file="c0_marker_fp_detail.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c("SAMHD1", "CD4", "CD27", "TCF7", "TMEM204", "RCAN3"), max.cutoff = 3, ncol = 3)
dev.off()
tiff(file="c0_marker_vp_detail.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c("SAMHD1", "CD4", "CD27", "TCF7", "TMEM204", "RCAN3"), pt.size = 0, combine = T, ncol = 3)
dev.off()
tiff(file="c1_marker_fp_detail.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c("CD19", "CD79A", "CD79B", "MS4A1", "TNFRSF13C", "CD37"), max.cutoff = 3, ncol = 3)
dev.off()
tiff(file="c1_marker_vp_detail.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c("CD19", "CD79A", "CD79B", "MS4A1", "TNFRSF13C", "CD37"), pt.size = 0, combine = T, ncol = 3)
dev.off()
tiff(file="c2_marker_fp_detail.tiff", width = 21, height = 6, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c("CD14", "CD68", "S100A12"), max.cutoff = 3, ncol = 3)
dev.off()
tiff(file="c2_marker_vp_detail.tiff", width = 21, height = 6, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c("CD14", "CD68", "S100A12"), pt.size = 0, combine = T, ncol = 3)
dev.off()
tiff(file="c3_marker_fp_detail1.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle.combined, features = c("CXCL2", "IL1B", "S100A8", "S100A9", "TLR2", "TREM1"), max.cutoff = 3, ncol = 3)
dev.off()
tiff(file="c3_marker_vp_detail1.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle.combined, features = c("CXCL2", "IL1B", "S100A8", "S100A9", "TLR2", "TREM1"), pt.size = 0, combine = T, ncol = 3)
dev.off()


# Identify conserved cell type markers
DefaultAssay(cattle.combined) <- "RNA"
cluster0.markers <- FindConservedMarkers(cattle.combined, ident.1 = 0, grouping.var = "orig.ident", verbose = FALSE)
cluster1.markers <- FindConservedMarkers(cattle.combined, ident.1 = 1, grouping.var = "orig.ident", verbose = FALSE)
cluster2.markers <- FindConservedMarkers(cattle.combined, ident.1 = 2, grouping.var = "orig.ident", verbose = FALSE)
cluster3.markers <- FindConservedMarkers(cattle.combined, ident.1 = 3, grouping.var = "orig.ident", verbose = FALSE)
cluster4.markers <- FindConservedMarkers(cattle.combined, ident.1 = 4, grouping.var = "orig.ident", verbose = FALSE)
cluster5.markers <- FindConservedMarkers(cattle.combined, ident.1 = 5, grouping.var = "orig.ident", verbose = FALSE)
cluster6.markers <- FindConservedMarkers(cattle.combined, ident.1 = 6, grouping.var = "orig.ident", verbose = FALSE)
