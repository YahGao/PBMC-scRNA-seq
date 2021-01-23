setwd("/yhgao/rumen/")
library(monocle)
library(Seurat)

cattle_rds <- readRDS('cattle_tutorial.rds')
data <- as(as.matrix(cattle_rds@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = cattle_rds@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
cattle_rds <- monocle_cds
cattle_rds <- estimateSizeFactors(cattle_rds)
cattle_rds <- estimateDispersions(cattle_rds)
cattle_rds <- detectGenes(cattle_rds, min_expr = 3 )
expressed_genes <- row.names(subset(fData(cattle_rds), num_cells_expressed >= 10))
disp_table <- dispersionTable(cattle_rds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cattle_rds <- setOrderingFilter(cattle_rds, unsup_clustering_genes$gene_id)
cattle_rds <- reduceDimension(cattle_rds, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = T)
cattle_rds <- clusterCells(cattle_rds, num_clusters = 2)
diff_test_res <- differentialGeneTest(cattle_rds[expressed_genes,], fullModelFormulaStr = "~percent.mt")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) 
cattle_rds <- setOrderingFilter(cattle_rds, ordering_genes)
cattle_rds <- reduceDimension(cattle_rds, max_components = 2, method = 'DDRTree')
cattle_rds <- orderCells(cattle_rds)


tiff(file="Fig. 2c.tiff", width = 11, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
plot_cell_trajectory(cattle_rds, color_by = "Pseudotime") +
  labs(x = "Componet 1", y="Componet 2") +
  theme(legend.text = element_text(size=15),
        legend.key.size = unit(1.5, 'cm'),
        legend.key.width = unit(0.5,'cm'),
        legend.position = "right",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=7,face="bold"),
        axis.title.y = element_text(size=7,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

tiff(file="Fig. 2d.tiff", width = 12, height = 9, units = "cm", res = 600, pointsize = 8,compression= "lzw")
plot_cell_trajectory(cattle_rds, color_by = "orig.ident") +
  labs(x = "Componet 1", y="Componet 2") +
  theme(legend.text = element_text(size=15),
        legend.key.size = unit(1.5, 'cm'),
        legend.key.width = unit(0.5,'cm'),
        legend.position = "right",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=7,face="bold"),
        axis.title.y = element_text(size=7,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

