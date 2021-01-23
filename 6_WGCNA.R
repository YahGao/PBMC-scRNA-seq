setwd("/yhgao/rumen/")

#============================#
# Seurat object construction #
#============================#
library(Seurat)

# Filter cells
datadf <- as.matrix(cattle.combined@assays$RNA@data)
idd1 <- cattle.combined@meta.data
Inter.id1 <- cbind(rownames(idd1),idd1$seurat_clusters)
rownames(Inter.id1) <- rownames(idd1)
colnames(Inter.id1) <- c("CellID","Celltype")
Inter.id1 <- as.data.frame(Inter.id1)
Inter1 <- datadf[,Inter.id1$CellID]
Inter2 <- as.matrix(Inter1)

pseudocell.size = 50 # 100 test, you can change the size as you need
new_ids_list1 = list()
for (i in 1:length(levels(Inter.id1$Celltype))) {
  cluster_id = levels(Inter.id1$Celltype)[i]
  cluster_cells <- rownames(Inter.id1[Inter.id1$Celltype == cluster_id,])
  cluster_size <- length(cluster_cells)     
  pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
  pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
  names(pseudo_ids) <- sample(cluster_cells)    
  new_ids_list1[[i]] <- pseudo_ids      
}

new_ids <- unlist(new_ids_list1)
new_ids <- as.data.frame(new_ids)
new_ids_length <- table(new_ids)
new_colnames <- rownames(new_ids)
all.data <- datadf[, as.character(new_colnames)]
all.data <- t(all.data)
new.data <- aggregate(list(all.data[, 1:length(all.data[1,])]), list(name = new_ids[,1]), FUN=mean)
rownames(new.data) <- new.data$name
new.data <- new.data[, -1]
new_ids_length <- as.matrix(new_ids_length)
short <- which(new_ids_length < 50) # change this number as the same number of pseudocell.size
new_good_ids <- as.matrix(new_ids_length[-short,])
result <- t(new.data)[, rownames(new_good_ids)]

# Filter genes
cattle.combined <- FindVariableFeatures(cattle.combined, nfeatures = 2000) # change this number as you need
Cluster1 <- result[intersect(Seurat::VariableFeatures(cattle.combined), rownames(result)),]


#===========#
#   WGCNA   #
#===========#
library(WGCNA)

type <- "unsigned"  
corType <- "pearson"
corFnc <- ifelse(corType == "pearson", cor, bicor)
maxPOutliers <- ifelse(corType == "pearson", 1, 0.05)
robustY <- ifelse(corType == "pearson", T, F)
dataExpr <- as.matrix(Cluster1)
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),] 
dataExpr <- as.data.frame(t(dataExprVar))

gsg <- goodSamplesGenes(dataExpr, verbose = 3)
nGenes <- ncol(dataExpr)
nSamples <- nrow(dataExpr)
sampleTree <- hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")

powers <- c(c(1:10), seq(from = 12, to =30, by = 2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
cex1 <- 0.9
power <- sft$powerEstimate
softPower <- power

cor <- WGCNA::cor
net <- blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                        TOMType = "unsigned", minModuleSize = 10,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs=TRUE, corType = corType, 
                        maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                        saveTOMFileBase = paste0("dataExpr", ".tom"),
                        verbose = 3)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

pdf("Fig. 3a.pdf",width=10,height=4)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module colors", dendroLabels = F, hang = 0.03, addGuide = T, guideHang = 0.05)
dev.off()


MEs <- net$MEs
MEs_col <- MEs
library(stringr)
colnames(MEs_col) <- paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col <- orderMEs(MEs_col)
MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")
merge <- mergeCloseModules(dataExpr, moduleColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors;
mergedMEs <- merge$newMEs;
which.module <- "turquoise"; 
ME <- mergedMEs[, paste("ME", which.module, sep="")]
load(net$TOMFiles, verbose = T)
TOM <- as.matrix(TOM)
dissTOM <- 1-TOM
plotTOM <- dissTOM^7
diag(plotTOM) <- NA
myrumen <- CreateSeuratObject(Cluster1)
mergedMEs_name <- rownames(mergedMEs)
mergedMEs_name <- substr(mergedMEs_name, 1,1) 

for (i in 1:7) {
  mergedMEs_all <- mergedMEs[mergedMEs_name==i,]
  myrumen_all <- subset(myrumen, idents = i)
  moduleTraitCor_noFP1 <- cor(mergedMEs_all, myrumen_all@meta.data$nFeature_RNA, use = "p")
  moduleTraitPvalue_noFP1 <- corPvalueStudent(moduleTraitCor_noFP1, nSamples)
  textMatrix_noFP1 <- paste(signif(moduleTraitCor_noFP1, 2), "\n(", signif(moduleTraitPvalue_noFP1, 1), ")", sep = "")
  moduleTraitCor_noFP2 <- cor(mergedMEs_all, myrumen_all@meta.data$nCount_RNA, use = "p")
  moduleTraitPvalue_noFP2 <- corPvalueStudent(moduleTraitCor_noFP2, nSamples)
  textMatrix_noFP2 <- paste(signif(moduleTraitCor_noFP2, 2), "\n(", signif(moduleTraitPvalue_noFP2, 1), ")", sep = "")
  moduleTraitCor_noFP <- cbind(moduleTraitCor_noFP1, moduleTraitCor_noFP2)
  moduleTraitCor_noFP <- as.data.frame(moduleTraitCor_noFP)
  colnames(moduleTraitCor_noFP) <- c("nFeature_RNA", "nCount_RNA")
  textMatrix_noFP <- cbind(textMatrix_noFP1, textMatrix_noFP2)
  pdf(file=paste("cattle_wgcna_module",i-1,".pdf", sep = ""), width = 5, height = 7)
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor_noFP, 
                 xLabels = colnames(moduleTraitCor_noFP), 
                 yLabels = names(mergedMEs_all), 
                 ySymbols = names(mergedMEs_all), 
                 colorLabels = FALSE, 
                 colors = blueWhiteRed(50), 
                 textMatrix = textMatrix_noFP,
                 setStdMargins = FALSE, 
                 cex.text = 1, 
                 zlim = c(-1,1), 
                 main = paste("c",i-1, sep="")) 
  dev.off()
}

# Export the genes information
moduleColors_uniq <- unique(moduleColors)
for (color in moduleColors_uniq) {
  module <- color
  probes <- colnames(dataExpr)
  inModule <- (moduleColors == module)
  modProbes <- probes[inModule]
  write.table(modProbes, paste("module_gene_", module, ".txt", sep = ""), quote = F, row.names = F, col.names = F)
}


#=========================#
#   Enrichment analysis   #
#=========================#
library(clusterProfiler)
library(org.Bt.eg.db)

#"black"    "blue"      "brown"      "green"      "grey" "pink"  "red"    "turquoise"       "yellow"

module <- "yellow"
module_gene <- read.table(paste("module_gene_", module, ".txt", sep = ""), header = F)
module_gene$V1 <- as.character(module_gene$V1)
module_gene_id <- bitr(module_gene$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Bt.eg.db")
module_gene_list <- module_gene_id$ENTREZID
module_gene_list[duplicated(module_gene_list)]
module_gene_go <- enrichGO(module_gene_list, OrgDb = org.Bt.eg.db, ont='ALL', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID')
tiff(file = paste("module_gene_", module, "_go1.tiff", sep = ""), width = 20, height = 15, units = "cm", res = 600, pointsize = 8,compression= "lzw")
barplot(module_gene_go, showCategory=20, drop=T, font.size = 20) + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 18))
dev.off()

tiff(file = paste("module_gene_", module, "_go2.tiff", sep = ""), width = 23, height = 11, units = "cm", res = 600, pointsize = 8,compression= "lzw")
dotplot(module_gene_go, showCategory=50, font.size =20) + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 18))
dev.off()

module_gene_kegg <- enrichKEGG(module_gene_list, organism = 'bta', keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
tiff(file = paste("module_gene_", module, "_kegg.tiff", sep = ""), width = 30, height = 21, units = "cm", res = 600, pointsize = 8, compression= "lzw")
dotplot(module_gene_kegg, showCategory = 30, font.size = 20) + theme(legend.text = element_text(size = 15),legend.title = element_text(size = 18))
dev.off()
