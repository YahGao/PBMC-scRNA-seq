setwd("/yhgao/rumen/")
library(Seurat)
library(SingleR)
library(scRNAseq)
library(scater)
library(cowplot)
library(ggplot2)
library(pheatmap)

# Step 1. By cluster
# Cell cycle heatmap by cluster
cell_cycle_g1_s <- read.table("cell_cycle_g1_s.txt", header = F)
cell_cycle_g2_m <- read.table("cell_cycle_g2_m.txt", header = F)
cell_cycle_genes <- rbind(cell_cycle_g1_s, cell_cycle_g2_m)
cell_cycle_genes <- unique(cell_cycle_genes)
cattle.combined_umi <- as.matrix(cattle.combined@assays$RNA@data)
cattle.combined_umi3 <- cattle.combined_umi
cattle.combined_umi3 <- as.data.frame(cattle.combined_umi3)
cattle.combined_umi3$genes <- rownames(cattle.combined_umi3)
colnames(cell_cycle_genes) <- "genes"
cattle.combined_umi3_cycle <- merge(cell_cycle_genes, cattle.combined_umi3)
rownames(cattle.combined_umi3_cycle) <- cattle.combined_umi3_cycle$genes
cattle.combined_umi3_cycle <- cattle.combined_umi3_cycle[,-1]
cattle.combined_umi3_cycle <- t(cattle.combined_umi3_cycle)
cattle.combined_umi3_cycle <- as.data.frame(cattle.combined_umi3_cycle)
cattle.combined_umi3_cycle$cells <- rownames(cattle.combined_umi3_cycle)
cattle.combined_cluster <- as.data.frame(cattle.combined$seurat_clusters)
cattle.combined_cluster$cells <- rownames(cattle.combined_cluster)
colnames(cattle.combined_cluster) <- c("cluster", "cells")
cattle.combined_umi3_cycle_cluster <- merge(cattle.combined_cluster, cattle.combined_umi3_cycle)
cattle.combined_umi3_cycle_cluster <- cattle.combined_umi3_cycle_cluster[order(cattle.combined_umi3_cycle_cluster[,2]),]
rownames(cattle.combined_umi3_cycle_cluster) <- cattle.combined_umi3_cycle_cluster$cells
cattle.combined_umi3_cycle_cluster <- cattle.combined_umi3_cycle_cluster[,-1]
cattle.combined_umi3_cycle_cluster1 <- as.data.frame(cattle.combined_umi3_cycle_cluster$cluster)
cattle.combined_umi3_cycle_cluster2 <- cattle.combined_umi3_cycle_cluster[,-1]
cattle.combined_umi3_cycle_cluster3 <- apply(cattle.combined_umi3_cycle_cluster2, 1, mean)
cattle.combined_umi3_cycle_cluster3 <- as.data.frame(cattle.combined_umi3_cycle_cluster3)
cattle.combined_umi3_cycle_cluster4 <- cbind(cattle.combined_umi3_cycle_cluster1, cattle.combined_umi3_cycle_cluster3)
colnames(cattle.combined_umi3_cycle_cluster4) <- c("cluster", "mean_cc_genes")
sdata <- split(cattle.combined_umi3_cycle_cluster4,cattle.combined_umi3_cycle_cluster4$cluster)
result <- lapply(sdata,function(x) x[order(x[,2], decreasing = T),])
result0 <- result$`0`
result1 <- result$`1`
result2 <- result$`2`
result3 <- result$`3`
result4 <- result$`4`
result5 <- result$`5`
result6 <- result$`6`
result_all <- rbind(result0,result1,result2,result3,result4,result5,result6)
result_all$cells <- rownames(result_all)
cattle.combined_umi3_cycle_cluster$cells <- rownames(cattle.combined_umi3_cycle_cluster)
result_all_umi <- merge(result_all, cattle.combined_umi3_cycle_cluster, sort = F)
result_all_umi_cluster <- data.frame(result_all_umi$cluster)
rownames(result_all_umi_cluster) <- result_all_umi$cells
colnames(result_all_umi_cluster) <- "cluster"
result_all_umi1 <- result_all_umi
rownames(result_all_umi1) <- result_all_umi1$cells
result_all_umi1 <- result_all_umi1[,-c(1:3)]
result_all_umi1 <- t(result_all_umi1)
result_all_umi1 <- log2(result_all_umi1 + 1)
result_all_umi1 <- result_all_umi1 - 1
result_all_umi2 <- as.data.frame(result_all_umi1)
result_all_umi2$genes <- rownames(result_all_umi2)
result_all_umi2_cc <- merge(cell_cycle_genes, result_all_umi2, sort = F)
rownames(result_all_umi2_cc) <- result_all_umi2_cc$genes
result_all_umi2_cc <- result_all_umi2_cc[,-1]
tiff(file="Fig. 2a.tiff", width = 24, height = 15, units = "cm", res = 600, pointsize = 8,compression= "lzw")
pheatmap(result_all_umi2_cc,show_colnames =F,show_rownames = T, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("white", "red"))(100),
         annotation_col=result_all_umi_cluster,annotation_legend = F,annotation_names_col=F,
         fontsize_row = 5, legend_breaks = c(-1:1), legend_labels = c("-1","0","1"))
dev.off()


# Cell cycle index by cluster
cattle.combined_umi <- as.matrix(cattle.combined@assays$RNA@data)
cattle.combined_umi <- as.data.frame(cattle.combined_umi)
cattle.combined_umi$genes <- rownames(cattle.combined_umi)
cc_cattle.combined_umi <- merge(cell_cycle_genes, cattle.combined_umi)
rownames(cc_cattle.combined_umi) <- cc_cattle.combined_umi[,1]
cc_cattle.combined_umi <- cc_cattle.combined_umi[,-1]
cc_cattle.combined_umi <- t(cc_cattle.combined_umi)
cc_cattle.combined_umi <- as.data.frame(cc_cattle.combined_umi)
cc_cattle.combined_umi$cells <- rownames(cc_cattle.combined_umi)
cattle.combined_cluster <- as.data.frame(cattle.combined$seurat_clusters)
colnames(cattle.combined_cluster) <- "cluster"
cattle.combined_cluster$cells <- rownames(cattle.combined_cluster)
cattle.combined_cluster_cells <- merge(cattle.combined_cluster, cc_cattle.combined_umi)
cattle.combined_cluster_cells <- cattle.combined_cluster_cells[order(cattle.combined_cluster_cells[,2]),]
rownames(cattle.combined_cluster_cells) <- cattle.combined_cluster_cells$cells
cattle.combined_cluster_cells <- cattle.combined_cluster_cells[,-1]
cattle.combined_cluster_cells1 <- cattle.combined_cluster_cells
cattle.combined_cluster_cells1$mean <- rowMeans(cattle.combined_cluster_cells1[,2:94])
cattle.combined_cluster_cells1 <- cattle.combined_cluster_cells1[,-c(2:94)]
cattle.combined_cluster_cells2 <- mean(cattle.combined_cluster_cells1$mean)
cattle.combined_cluster_cells2_1 <- rep(cattle.combined_cluster_cells2, times = 26141)
cattle.combined_cluster_cells2_1 <- as.data.frame(cattle.combined_cluster_cells2_1)
cattle.combined_cluster_cells3 <- cbind(cattle.combined_cluster_cells1, cattle.combined_cluster_cells2_1)
colnames(cattle.combined_cluster_cells3) <- c("cluster", "mean", "mean1")
cattle.combined_cluster_cells3$judge <- cattle.combined_cluster_cells3$mean > cattle.combined_cluster_cells3$mean1
cattle.combined_cluster_cells4 <- aggregate(cattle.combined_cluster_cells3[,4], list(cattle.combined_cluster_cells3[,1]), table)
cattle.combined_cluster_cells4
cattle.combined_cluster_cells4_1 <- c(6111,4083,1003,564,650,187,19)  
cattle.combined_cluster_cells5 <- aggregate(cattle.combined_cluster_cells3[,1], list(cattle.combined_cluster_cells3[,1]), length)
cattle.combined_cluster_cells5 <- cbind(cattle.combined_cluster_cells5, cattle.combined_cluster_cells4_1)
colnames(cattle.combined_cluster_cells5) <- c("cluster", "total", "ture")
cattle.combined_cluster_cells5$index <- cattle.combined_cluster_cells5$ture/cattle.combined_cluster_cells5$total
write.table(cattle.combined_cluster_cells5, "cluster_cells.txt", quote = F, row.names = F, sep = "\t")


# Step 2. By ident
# Cell cycle heatmap by ident
cattle.combined_umi3_cycle_cells <- cattle.combined_umi3_cycle
cattle.combined_umi3_cycle_cells <- cattle.combined_umi3_cycle_cells[,-94]
cattle.combined_umi3_cycle_cells1 <- apply(cattle.combined_umi3_cycle_cells, 1, mean)
cattle.combined_umi3_cycle_cells1 <- as.data.frame(cattle.combined_umi3_cycle_cells1)
cattle.combined_umi3_cycle_cells1$cells <- rownames(cattle.combined_umi3_cycle_cells1)
colnames(cattle.combined_umi3_cycle_cells1) <- c("cells_mean", "cells")
cattle.combined_ident <- as.data.frame(cattle.combined$orig.ident)
cattle.combined_ident$cells <- rownames(cattle.combined_ident)
colnames(cattle.combined_ident) <- c("ident", "cells")
cattle.combined_umi3_cycle_cells2 <- merge(cattle.combined_ident, cattle.combined_umi3_cycle_cells1)
cattle.combined_umi3_cycle_cells3 <- merge(cattle.combined_umi3_cycle_cells2, cattle.combined_umi3_cycle)
rownames(cattle.combined_umi3_cycle_cells3) <- cattle.combined_umi3_cycle_cells3$cells
cattle.combined_umi3_cycle_cells3 <- cattle.combined_umi3_cycle_cells3[,-1]
sdata1 <- split(cattle.combined_umi3_cycle_cells3 ,cattle.combined_umi3_cycle_cells3$ident)
result1 <- lapply(sdata1,function(x) x[order(x[,2], decreasing = T),])
result1_c <- result1$cattle_c
result1_t1 <- result1$cattle_t1
result1_t2 <- result1$cattle_t2
result1_t3 <- result1$cattle_t3
result1_combined <- rbind(result1_c, result1_t1, result1_t2, result1_t3)
result1_combined_ident <- as.data.frame(result1_combined$ident)
rownames(result1_combined_ident) <- rownames(result1_combined)
colnames(result1_combined_ident) <- "ident"
result1_combined1 <- result1_combined
result1_combined1 <- result1_combined1[,-c(1:2)]
result1_combined1 <- t(result1_combined1)
result1_combined1 <- as.data.frame(result1_combined1)
result1_combined1$genes <- rownames(result1_combined1)
result1_combined2 <- merge(cell_cycle_genes, result1_combined1, sort = F)
rownames(result1_combined2) <- result1_combined2$genes
result1_combined2 <- result1_combined2[,-1]
result1_combined2 <- log2(result1_combined2 + 1)
result1_combined2 <- result1_combined2 -1
tiff(file="Fig. S4.tiff", width = 24, height = 15, units = "cm", res = 600, pointsize = 8,compression= "lzw")
pheatmap(result1_combined2,show_colnames =F,show_rownames = T, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("white", "red"))(100),
         annotation_col=result1_combined_ident,annotation_legend = F,annotation_names_col=F,
         fontsize_row = 5, legend_breaks = c(-1:1), legend_labels = c("-1","0","1"))
dev.off()


# Cell cycle index by ident
cattle.combined_umi$genes <- rownames(cattle.combined_umi)
cc_cattle.combined_umi <- merge(cell_cycle_genes, cattle.combined_umi)
rownames(cc_cattle.combined_umi) <- cc_cattle.combined_umi[,1]
cc_cattle.combined_umi <- cc_cattle.combined_umi[,-1]
cc_cattle.combined_umi <- t(cc_cattle.combined_umi)
cc_cattle.combined_umi <- as.data.frame(cc_cattle.combined_umi)
cc_cattle.combined_umi$cells <- rownames(cc_cattle.combined_umi)
cattle.combined_ident <- as.data.frame(cattle.combined$orig.ident)
colnames(cattle.combined_ident) <- "ident"
cattle.combined_ident$cells <- rownames(cattle.combined_ident)
cattle.combined_ident_cells <- merge(cattle.combined_ident, cc_cattle.combined_umi)
cattle.combined_ident_cells <- cattle.combined_ident_cells[order(cattle.combined_ident_cells[,2]),]
rownames(cattle.combined_ident_cells) <- cattle.combined_ident_cells$cells
cattle.combined_ident_cells <- cattle.combined_ident_cells[,-1]
cattle.combined_ident_cells1 <- cattle.combined_ident_cells
cattle.combined_ident_cells1$mean <- rowMeans(cattle.combined_ident_cells1[,2:94])
cattle.combined_ident_cells1 <- cattle.combined_ident_cells1[,-c(2:94)]
cattle.combined_ident_cells2 <- mean(cattle.combined_ident_cells1$mean)
cattle.combined_ident_cells2_1 <- rep(cattle.combined_ident_cells2, times = 26141)
cattle.combined_ident_cells2_1 <- as.data.frame(cattle.combined_ident_cells2_1)
cattle.combined_ident_cells3 <- cbind(cattle.combined_ident_cells1, cattle.combined_ident_cells2_1)
colnames(cattle.combined_ident_cells3) <- c("ident", "mean", "mean1")
cattle.combined_ident_cells3$judge <- cattle.combined_ident_cells3$mean > cattle.combined_ident_cells3$mean1
cattle.combined_ident_cells4 <- aggregate(cattle.combined_ident_cells3[,4], list(cattle.combined_ident_cells3[,1]), table)
cattle.combined_ident_cells4
cattle.combined_ident_cells4_1 <- c(3337, 4375, 3266, 1639)  
cattle.combined_ident_cells5 <- aggregate(cattle.combined_ident_cells3[,1], list(cattle.combined_ident_cells3[,1]), length)
cattle.combined_ident_cells5 <- cbind(cattle.combined_ident_cells5, cattle.combined_ident_cells4_1)
colnames(cattle.combined_ident_cells5) <- c("ident", "total", "ture")
cattle.combined_ident_cells5$index <- cattle.combined_ident_cells5$ture/cattle.combined_ident_cells5$total
write.table(cattle.combined_ident_cells5, "ident_cells.txt", quote = F, row.names = F, sep = "\t")
