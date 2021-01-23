setwd("/yhgao/rumen/")
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)
library(SeuratData)

#=======
# Seurat
#=======
cattle_c.data <- Read10X(data.dir = "c:/YHGAO/work/usda/6-scRNA-sd/analysis/C/")
cattle_c <- CreateSeuratObject(counts = cattle_c.data, project = "cattle_c", min.cells = 3, min.features = 200)
cattle_c[["percent.mt"]] <- PercentageFeatureSet(cattle_c, pattern = "^MT-")
cattle_c <- subset(cattle_c, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 10000 & percent.mt < 15)
cattle_c <- NormalizeData(cattle_c, normalization.method = "LogNormalize", scale.factor = 10000)
cattle_c <- FindVariableFeatures(cattle_c, selection.method = "vst", nfeatures = 2000)

data("ifnb")
ifnb.list <- SplitObject(ifnb, split.by = "stim")
ctrl <- subset(ifnb, idents = "IMMUNE_CTRL")
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

cattle_human <- merge(cattle_c, y = ctrl, add.cell.ids = c("Cattle", "Human"), project = "pbmc")
cattle_human.list <- SplitObject(cattle_human, split.by = "orig.ident")
cattle_human.list <- lapply(X = cattle_human.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

cattle_human.anchors <- FindIntegrationAnchors(object.list = cattle_human.list, dims = 1:20)
cattle_human.combined <- IntegrateData(anchorset = cattle_human.anchors, dims = 1:20)
DefaultAssay(cattle_human.combined) <- "integrated"
cattle_human.combined <- ScaleData(cattle_human.combined, verbose = FALSE)
cattle_human.combined <- RunPCA(cattle_human.combined, npcs = 30, verbose = FALSE)
cattle_human.combined <- JackStraw(cattle_human.combined, num.replicate = 100)
cattle_human.combined <- ScoreJackStraw(cattle_human.combined, dims = 1:20)
ElbowPlot(cattle_human.combined) + geom_vline(aes(xintercept=10), col="black") 
cattle_human.combined <- RunUMAP(cattle_human.combined, reduction = "pca", dims = 1:13)
cattle_human.combined <- FindNeighbors(cattle_human.combined, reduction = "pca", dims = 1:13)
cattle_human.combined <- FindClusters(cattle_human.combined, resolution = 0.18)

tiff(file="Fig. 6a.tiff", width = 10, height = 7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(cattle_human.combined, reduction = "umap", group.by = "orig.ident", pt.size = 0.5,
        cols = c("#f86e65", "#00bfc4")) +
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

tiff(file="Fig. 6b.tiff", width = 13, height = 7.5, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(cattle_human.combined, reduction = "umap", split.by = "orig.ident", pt.size = 0.5, 
        ncol=3, label = T, label.size = 3) +
  theme(legend.text = element_text(size=8),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size=7,face="bold"),
        axis.title.y = element_text(size=7,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=0.9)))
dev.off()

tiff(file="Fig. 6c.tiff", width = 8.5, height = 7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(cattle_human.combined, reduction = "umap", label = T, pt.size = 0.5, label.size = 3) +
  theme(legend.text = element_text(size=8),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size=7,face="bold"),
        axis.title.y = element_text(size=7,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=0.9)))
dev.off()


cattle_human.combined_data <- as.data.frame(cattle_human.combined@assays$integrated@data)
cattle_human.combined_data <- as.data.frame(t(cattle_human.combined_data))
cattle <- cattle_human.combined_data[1:7053,]
cattle$cells <- rownames(cattle)
human <- cattle_human.combined_data [7054:13601,]
human$cells <- rownames(human)
cattle_human.combined_cluster <- as.data.frame(cattle_human.combined$seurat_clusters)
cattle_human.combined_cluster$cells <- rownames(cattle_human.combined_cluster)
colnames(cattle_human.combined_cluster) <- c("cluster", "cells")
cattle_counts <- merge(cattle_human.combined_cluster, cattle)
rownames(cattle_counts) <- cattle_counts$cells
cattle_counts <- cattle_counts[,-1]
cattle_counts1 <- aggregate(x=cattle_counts[,2:2001], by = list(cattle_counts$cluster), FUN=median)
rownames(cattle_counts1) <- cattle_counts1$Group.1
cattle_counts1 <- cattle_counts1[,-1]
cattle_counts1 <- as.data.frame(t(cattle_counts1))
colnames(cattle_counts1) <- c("cc0","cc1","cc2","cc3","cc4","cc5","cc6","cc7","cc8","cc9")
cattle_counts1$genes <- rownames(cattle_counts1)
human_counts <- merge(cattle_human.combined_cluster, human)
rownames(human_counts) <- human_counts$cells
human_counts <- human_counts[,-1]
human_counts1 <- aggregate(x=human_counts[,2:2001], by = list(human_counts$cluster), FUN=median)
rownames(human_counts1) <- human_counts1$Group.1
human_counts1 <- human_counts1[,-1]
human_counts1 <- as.data.frame(t(human_counts1))
colnames(human_counts1) <- c("hc0","hc1","hc2","hc3","hc4","hc5","hc6","hc7","hc8")
human_counts1$genes <- rownames(human_counts1)
c_h_counts <- merge(cattle_counts1, human_counts1)
rownames(c_h_counts) <- c_h_counts$genes
c_h_counts <- c_h_counts[,-1]
c_h_counts_cor <- cor(c_h_counts)
library(pheatmap)
tiff(file="Fig. 6d.tiff", width = 16, height = 15, units = "cm", res = 600, pointsize = 8,compression= "lzw")
pheatmap(c_h_counts_cor, fontsize = 15)
dev.off()


#========
# Azimuth
#========
library(SeuratDisk)
library(ggplot2)
library(patchwork)

reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")
pbmc <- SCTransform(pbmc, verbose = FALSE)
anchors <- FindTransferAnchors(reference = reference, query = pbmc, normalization.method = "SCT", reference.reduction = "spca", dims = 1:50)
pbmc <- MapQuery(anchorset = anchors, query = pbmc, reference = reference, refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2", celltype.l3 = "celltype.l3", predicted_ADT = "ADT"), reference.reduction = "spca", reduction.model = "wnn.umap")
Idents(pbmc) <- 'predicted.celltype.l2'
treg_markers <- FindMarkers(pbmc, ident.1 = "Treg", only.pos = T, logfc.threshold = 0.1)
DefaultAssay(pbmc) <- 'predicted_ADT'
reference$id <- 'reference'
pbmc$id <- 'query'
refquery <- merge(reference, pbmc)
refquery[["spca"]] <- merge(reference[["spca"]], pbmc[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
p4 <- DimPlot(pbmc, group.by = "predicted.celltype.l1", label = T, label.size = 3, repel = T) + NoLegend()
p5 <- DimPlot(pbmc, group.by = "predicted.celltype.l2", label = T, label.size = 3, repel = T) + NoLegend()
p6 <- DimPlot(pbmc, group.by = "predicted.celltype.l3", label = T, label.size = 3, repel = T) + NoLegend()
pdf("Fig. S7.pdf", width=15, height=5)
p4 + p5 + p6
dev.off()



# Fig. S8. marker genes
tiff(file="c046_marker_vp_detail.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle_human.combined, features = c("CD5","LEF1","CD27","RCAN3","GATA3","RHOH"), pt.size = 0, combine = T, ncol = 3)
dev.off()
tiff(file="c046_marker_fp_detail.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle_human.combined, features = c("CD5","LEF1","CD27","RCAN3","GATA3","RHOH"), max.cutoff = 3, ncol = 3)
dev.off()

c1_marker <- read.table("hh.txt", header = F)
tiff(file="c1_marker_vp_detail.tiff", width = 21, height = 6, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle_human.combined, features = c("S100A8","S100A9","VCAN"), pt.size = 0, combine = T, ncol = 3)
dev.off()
tiff(file="c1_marker_fp_detail.tiff", width = 21, height = 6, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle_human.combined, features = c("S100A8","S100A9","VCAN"), max.cutoff = 3, ncol = 3)
dev.off()

c28_marker <- read.table("hh.txt", header = F)
tiff(file="c28_marker_vp_detail.tiff", width = 21, height = 30, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle_human.combined, features = c("CD19","CD79A","CD79B","HHEX","MS4A1","VPREB3","SPIB","LY86","TNFRSF13C","JCHAIN","MZB1","CRELD2","DERL3","SEC61A1"), pt.size = 0, combine = T, ncol = 3)
dev.off()
tiff(file="c28_marker_fp_detail.tiff", width = 21, height = 30, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle_human.combined, features = c("CD19","CD79A","CD79B","HHEX","MS4A1","VPREB3","SPIB","LY86","TNFRSF13C","JCHAIN","MZB1","CRELD2","DERL3","SEC61A1"), max.cutoff = 3, ncol = 3)
dev.off()

c3_marker <- read.table("hh.txt", header = F)
tiff(file="c3_marker_vp_detail.tiff", width = 21, height = 18, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle_human.combined, features = c("GNLY","NKG7","CTSW","PRF1","CCL5","GZMB","HOPX"), pt.size = 0, combine = T, ncol = 3)
dev.off()
tiff(file="c3_marker_fp_detail.tiff", width = 21, height = 18, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle_human.combined, features = c("GNLY","NKG7","CTSW","PRF1","CCL5","GZMB","HOPX"), max.cutoff = 3, ncol = 3)
dev.off()

c5_marker <- read.table("hh.txt", header = F)
tiff(file="c5_marker_vp_detail.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle_human.combined, features = c("FCGR3A","C3AR1","MS4A4A","MS4A7","CYBB","AIF1"), pt.size = 0, combine = T, ncol = 3)
dev.off()
tiff(file="c5_marker_fp_detail.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle_human.combined, features = c("FCGR3A","C3AR1","MS4A4A","MS4A7","CYBB","AIF1"), max.cutoff = 3, ncol = 3)
dev.off()

c7_marker <- read.table("hh.txt", header = F)
tiff(file="c7_marker_vp_detail.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle_human.combined, features = c("CD5","LEF1","CD27","RCAN3","GATA3","RHOH"), pt.size = 0, combine = T, ncol = 3)
dev.off()
tiff(file="c7_marker_fp_detail.tiff", width = 21, height = 12, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle_human.combined, features = c("CD5","LEF1","CD27","RCAN3","GATA3","RHOH"), max.cutoff = 3, ncol = 3)
dev.off()

c9_marker <- read.table("hh.txt", header = F)
tiff(file="c9_marker_vp_detail.tiff", width = 21, height = 6, units = "cm", res = 600, pointsize = 8,compression= "lzw")
VlnPlot(cattle_human.combined, features = c("PTPRE","GATA2","RGS18"), pt.size = 0, combine = T, ncol = 3)
dev.off()
tiff(file="c9_marker_fp_detail.tiff", width = 21, height = 6, units = "cm", res = 600, pointsize = 8,compression= "lzw")
FeaturePlot(cattle_human.combined, features = c("PTPRE","GATA2","RGS18"), max.cutoff = 3, ncol = 3)
dev.off()
