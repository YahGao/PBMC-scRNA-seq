setwd("/yhgao/rumen/")

# Azimuth
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")

pbmc <- readRDS("cattle_tutorial_new.rds")
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

p1 = DimPlot(pbmc, group.by = "predicted.celltype.l1", label = T, label.size = 3, repel = T) + NoLegend()
p2 = DimPlot(pbmc, group.by = "predicted.celltype.l2", label = T, label.size = 3, repel = T) + NoLegend()
p3 = DimPlot(pbmc, group.by = "predicted.celltype.l3", label = T, label.size = 3, repel = T) + NoLegend()
pdf("Fig. S3a.pdf", width=15, height=5)
p1 + p2 + p3
dev.off()


# SingleR
library(SingleR)
library(scRNAseq)
library(scater)
library(cowplot)

blue.se <- BlueprintEncodeData()
cattle_blue <- SingleR(test = as.SingleCellExperiment(cattle.combined), ref = blue.se,labels = blue.se$label.fine)
cattle.combined$labels_blue <- cattle_blue$labels

tiff(file="Fig. S3b1.tiff", width = 27, height = 18, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(cattle.combined, reduction = "umap", label = F, pt.size = 2, group.by = "labels_blue") +
  theme(legend.text = element_text(size=15),
        legend.position = "right",
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 16, angle = 0, hjust = 0.6),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=3), ncol = 1))
dev.off()

tiff(file="Fig. S3b2.tiff", width = 45, height = 30, units = "cm", res = 600, pointsize = 8,compression= "lzw")
DimPlot(cattle.combined, reduction = "umap", label = T, pt.size = 0.5, split.by = "labels_blue", 
        ncol = 6, label.size = 5) +  NoLegend() +
  theme(legend.text = element_text(size=15),
        legend.position = "right",
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        plot.title = element_text(hjust = 0.5, size = 0.3),
        axis.text.x = element_text(size = 10, angle = 0, hjust = 0.6),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=3), ncol = 1))
dev.off()
