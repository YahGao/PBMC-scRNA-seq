setwd("/yhgao/rumen/")

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(SCENIC)

# Step 1. Prepare input file
cells.use <- WhichCells(object = cattle.combined) 
exprMat <- GetAssayData(object = cattle.combined, assay= "RNA", slot = "data")[, cells.use] 
exprMat <- as(Class = 'matrix', object = exprMat)
gene <- cattle.combined$nFeature_RNA
umi <- cattle.combined$nCount_RNA
clus <- cattle.combined$seurat_clusters
cellInfo <- data.frame(gene, umi, clus)

# Step 2. Initialization
cellTypeColumn <- "clus"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")
colVars <- list(CellType=c("0"="forestgreen", "1"="darkorange", "2"="magenta4", "3"="hotpink", "4"="red3", "5"="skyblue", "6"="darkblue"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")
org="hgnc"
dbDir="hg19_mc9nr"
myDatasetTitle="SCENIC example on Human"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle) 
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# Step 3. Construction of gene coexpression network
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions, minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions)


# Step 4. Construct and calculate gene regulation network
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"))
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# Step 5. Convert the network activity to the format of ON/OFF
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
runSCENIC_4_aucell_binarize(scenicOptions) # Set the threshold

nPcs <- c(5)
scenicOptions@settings$seed <- 123 # same seed for all of them
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


# Fig. 2b
library(AUCell)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
thresholds <- loadInt(scenicOptions, "aucell_thresholds")
thresholds <- getThresholdSelected(thresholds)

regulonsCells <- setNames(lapply(names(thresholds), function(x) {
  trh <- thresholds[x]
  names(which(getAUC(regulonAUC)[x, ] > trh))
}), names(thresholds))

regulonActivity <- reshape2::melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity[, 1], regulonActivity[, 2]))
class(binaryRegulonActivity) <- "matrix"
regulonSelection <- loadInt(scenicOptions, "aucell_regulonSelection", ifNotExists = "null", verbose = F)

library(NMF)
cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists = "null")
cellInfo <- data.frame(cellInfo)
colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists = "null")

.openDevHeatmap <- function(fileName, devType)
{
  if(devType!="pdf") 
  {
    if(devType=="png") .openDev(fileName=fileName, devType=devType, width=1200,height=1200)
    if(devType!="png") .openDev(fileName=fileName, devType=devType)
    fileName <- NA
  }else{
    fileName <- paste0(fileName,".pdf")
  }
  return(fileName)
}

.closeDevHeatmap <- function(devType) {if(devType!="pdf") {dev.off()}}

regulon_label1 <- "all"
for (selRegs in regulon_label1) {
  if (length(regulonSelection[[selRegs]]) > 0) {
    regulonSelection[[selRegs]] <- regulonSelection[[selRegs]][which(regulonSelection[[selRegs]] %in%  rownames(binaryRegulonActivity))]
    binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]], , drop = F]
    if (nrow(binaryMat) > 0) {
      fileName <- paste0(getOutName(scenicOptions, "s4_binaryActivityHeatmap"), selRegs)
      fileName <- .openDevHeatmap(fileName = fileName, devType = getSettings(scenicOptions, "devType"))
      rowv <- ifelse(nrow(binaryMat) >= 2, T, NA)
      colv <- ifelse(ncol(binaryMat) >= 2, T, NA)
      NMF::aheatmap(binaryMat, scale = "none", 
                    revC = TRUE, main = selRegs, annCol = cellInfo[colnames(binaryMat),  , drop = F], annColor = colVars, Rowv = rowv, 
                    Colv = colv, color = c("white", "black"), filename = fileName)
      if (getSettings(scenicOptions, "devType") != 
          "pdf") 
        dev.off()
    }}}


cattle_pbmc_cluster <- as.data.frame(cattle.combined$seurat_clusters)
cattle_pbmc_cluster$cells <- rownames(cattle_pbmc_cluster)
colnames(cattle_pbmc_cluster) <- c("clusters", "cells")
cattle_pbmc_cluster <- cattle_pbmc_cluster[order(cattle_pbmc_cluster[,1]),]
binaryMat_col <- as.data.frame(colnames(binaryMat))
colnames(binaryMat_col) <- "cells"
binaryMat_col1 <- merge(cattle_pbmc_cluster, binaryMat_col, sort = F)
binaryMat_col1 <- binaryMat_col1[order(binaryMat_col1[,2]),]

tiff(file="Fig. 2b.tiff", width = 22, height = 18, units = "cm", res = 600, pointsize = 8,compression= "lzw")
aheatmap(binaryMat, scale = "none", revC = T, main = selRegs, labCol = NA,
         annCol = cellInfo[binaryMat_col1$cells, ,drop = T], annColor = colVars, Rowv = rowv, 
         Colv = NA, color = c("white", "black"))
dev.off()
