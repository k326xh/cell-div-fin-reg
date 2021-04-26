# This script handles outputs from running `cellranger count` on fastq files.
# Load required packages
library(Seurat) # Package version v3.0.2
library(plyr)
library(dplyr) # for mutate, "%>%"
library(Matrix)
library(RColorBrewer) # colorRampPalette(), brewer.pal
library("pheatmap") # for pheatmap()
library(ggplot2) 
library(slingshot) # pseudotime calculation

#### Individual data object constructions ####
# Below is the process of preparing 1dpa1 dataset. Use the same process for all datasets 
data_1dpa1 <- Read10X(data.dir = <raw_cellranger_output_path>) # Replace <> part with actual path
# Set up 1dpa1 object
colnames(x =  data_1dpa1) <- paste('1dpa1', colnames(x =  data_1dpa1), sep = '-')

# Create Seurat object. Set threshold to include only gene features that are detected in at least 5 cells.
obj_1dpa1 <- CreateSeuratObject(data_1dpa1, project = "1dpa1", min.cells = 5)
obj_1dpa1@meta.data$Stage <- "1dpa"

# Estimate mitochondiral contaminations
mito.genes <- grep(pattern = "^mt-", x = rownames(obj_1dpa1), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = obj_1dpa1, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = obj_1dpa1, slot = 'counts'))
## The [[ operator can add columns to object metadata, and is a great place to stash QC stats
obj_1dpa1[['percent.mito']] <- percent.mito

# Filter cells based on number of features captured, number of counts and percentage of mitochondrial reads.
obj_1dpa1 <- subset(x = obj_1dpa1, subset = nFeature_RNA > 500 & nCount_RNA < 40000 & percent.mito < 0.05)

# Normalization and variable gene selection
## sctransform replaces NormalizeData, ScaleData, and FindVariableFeatures.
## Transformed data will be available in the SCT assay, which is set as the default after running sctransform
obj_1dpa1 <- SCTransform(obj_1dpa1, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)

# Regress out cell cycle effect. Cell cycle gene list is converted from regev_lab_cell_cycle_genes.txt
obj_1dpa1 <- CellCycleScoring(obj_1dpa1, s.features = s.genes, g2m.features = g2m.genes)
obj_1dpa1$CC.Difference <- obj_1dpa1$S.Score - obj_1dpa1$G2M.Score
obj_1dpa1 <- ScaleData(obj_1dpa1, vars.to.regress = "CC.Difference", features = rownames(obj_1dpa1))

#### Data integration ####
# Set up list for all data objects for integration
obj.list <- list(obj_1dpa1,
                  obj_1dpa2, 
                  obj_2dpa1,
                  obj_2dpa2, 
                  obj_4dpa1,
                  obj_4dpa2, 
                  obj_samp1,
                  obj_samp2)
# Data integration
obj.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = obj.features, verbose = FALSE)
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = obj.features, verbose = FALSE)
# include all common genes for downstream analysis (integrate for all required genes)
GenesT0Integrate <- Reduce(intersect, lapply(obj.list, rownames))
fin.integrated <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", 
                                features.to.integrate = GenesT0Integrate ,verbose = T)

#### Dimensional reduction and clustering #### 
# Non-linear dim reduction for visualization
## default assay is "integrated" post integration
fin.integrated <- RunPCA(fin.integrated, verbose = T)
fin.integrated <- RunUMAP(fin.integrated, dims = 1:20)

# Linear dim reduction for clustering 
fin.integrated <- FindNeighbors(fin.integrated, reduction = "pca", dims = 1:20)
fin.integrated <- FindClusters(fin.integrated, resolution = 0.1)

# Differential analysis
## Find markers using RNA or SCT data slot
fin.cl.mk.sct <- FindAllMarkers(fin.integrated, assay = "SCT", only.pos = T, logfc.threshold = 0.25, min.pct = 0.25) 
## Take top 30 DEGs for heatmap
fin.cl.mk.sct.top <- fin.cl.mk.sct %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

# Assign cell identity
fin.integrated[["major.cl"]] <- mapvalues(fin.integrated$integrated_snn_res.0.1, 
	from = c("2","5","1","0","4","3"), to = c("Superficial Epithelial","Mucosal-like","Intermediate Epithelial","Basal Epithelial","Hematopoietic","Mesenchymal"))
## Reset levels to match the order from outer to inner layer
fin.integrated[["major.cl"]] <- factor(fin.integrated$major.cl, levels = c("Superficial Epithelial","Mucosal-like","Intermediate Epithelial","Basal Epithelial","Hematopoietic","Mesenchymal"))
## Merge epithelial identities in to Epithelial
fin.integrated[["merged.cl"]] <- mapvalues(fin.integrated$integrated_snn_res.0.1, 
	from = c("0","1","2","3","4","5"), to = c("Epithelial","Epithelial","Epithelial","Mesenchymal","Hematopoietic","Epithelial"))
## Reset levels to match the order from outer to inner layer 
fin.integrated[["merged.cl"]] <- factor(fin.integrated$merged.cl, levels = c("Epithelial","Hematopoietic","Mesenchymal"))
## Merge regenerating stages as regenerating
fin.integrated[["preReg"]] <- mapvalues(fin.integrated$Stage, 
                                         from = c("Preinjury","1dpa","2dpa","4dpa"),
                                         to = c("Preinjury","Regenerating","Regenerating","Regenerating"))
## Manual set stringent G1/G2M/S phase separations
Idents(fin.integrated, cells = colnames(fin.integrated)[fin.integrated$CC.Diff >= 0.1 & fin.integrated$S.Score > 0 & fin.integrated$G2M.Score > 0]) <- "S"
Idents(fin.integrated, cells = colnames(fin.integrated)[fin.integrated$CC.Diff <= -0.1 & fin.integrated$S.Score > 0 & fin.integrated$G2M.Score > 0]) <- "G2M"
## stringent cut off for G1. Use both score below zero (as what used in seurat by default)
Idents(fin.integrated, cells = colnames(fin.integrated)[fin.integrated$S.Score <= 0 & fin.integrated$G2M.Score <= 0]) <- "G1"
## Assign the rest as unselected
Idents(fin.integrated, cells = colnames(fin.integrated)[Idents(fin.integrated) != "S" & Idents(fin.integrated) != "G2M" & Idents(fin.integrated) != "G1"]) <- "Unselected"
fin.integrated[["manualPhase"]] <- Idents(object = fin.integrated)
fin.integrated[["manualPhase"]] <- factor(fin.integrated$manualPhase, levels = c("G1","S","G2M","Unselected"))

# Plot major cell cluster on UMAP axes and split by group. Other metadata options could be used for grouping/splitting.
DimPlot(fin.integrated, 
	group.by = "major.cl", 
	split.by = "Stage",
	cols = brewer.pal(7, "Set1")[c(1:5,7)])
# Plot gene expression distribution
FeaturePlot(object = fin.integrated, 
              features = c("epcam","cdh1","mpeg1.1","cxcr3.2","msx1b","msx3"), 
              min.cutoff = "q5", max.cutoff = "q95", 
              cols = c("lightgrey", "blue"), 
              pt.size = 0.1, 
              ncol = 2)
DoHeatmap(fin.integrated, features = fin.cl.mk.sct.top$gene)
DotPlot(fin.integrated, assay = "SCT", 
        features = rev(c("epcam","cdh1","krt4","agr2","foxp1b","foxp4","tp63","krtt1c19e","mpeg1.1","cxcr3.2","msx1b","twist1a")), 
        cols = c("lightgrey","black"),
        dot.scale = 8) + RotatedAxis()

#### Cell cycle analysis ####
Idents(fin.integrated) <- "manualPhase"
# Subset S phase
S.subset <- subset(x = fin.integrated, idents = "S")

S.subset <- RunPCA(S.subset, npcs = 50, verbose = F)
S.subset <- FindNeighbors(S.subset, reduction = "pca", dims = 1:10) # PC number could change
S.subset <- FindClusters(S.subset, resolution = 0.3)
S.subset <- RunUMAP(S.subset, reduction = "pca", dims = 1:10, n.neighbors = 50, min.dist = 0.5)
DimPlot(S.subset, reduction = "umap") + NoAxes()

#### Major cell population sub-clustering ####
# Hematopoietic 
## Subset H cells
Idents(fin.integrated) <- "merged.cl"
H.subset <- subset(fin.integrated, idents = "Hematopoietic")
DefaultAssay(H.subset) <- "integrated"
H.subset <- RunPCA(H.subset, npcs = 50, verbose = F)
H.subset <- FindNeighbors(H.subset, reduction = "pca", dims = 1:10) 
H.subset <- FindClusters(H.subset, resolution = 0.3)
H.subset <- RunUMAP(H.subset, reduction = "pca", dims = 1:10, n.neighbors = 50, min.dist = 0.2)
DimPlot(H.subset, reduction = "umap", cols = brewer.pal(4,"Dark2")) + NoAxes()
## FindAllMarkers
H.subset.markers <- FindAllMarkers(H.subset, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
H.subset.markers.top <- H.subset.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
H.subset <- ScaleData(H.subset,assay = "SCT", vars.to.regress = c("percent.mito", "nCount_RNA"))
DoHeatmap(H.subset, features = H.subset.markers.top$gene, assay = "SCT") + FontSize(y.text = 0, x.text = 0, main = 0) + NoLegend()

# Epithelial
Epi.subset <- subset(fin.integrated, idents = "Epithelial")
Epi.subset <- SCTransform(Epi.subset, vars.to.regress = c("percent.mito", "nCount_RNA","CC.Diff"), verbose = F, return.only.var.genes = F)
Epi.subset <- RunPCA(Epi.subset, npcs = 50, verbose = F)
Epi.subset <- FindNeighbors(Epi.subset, reduction = "pca", dims = 1:10) 
Epi.subset <- FindClusters(Epi.subset, resolution = 0.2)
Epi.subset <- RunUMAP(Epi.subset, reduction = "pca", dims = 1:10, n.neighbors = 50, min.dist = 0.2)
DimPlot(Epi.subset, reduction = "umap") + NoAxes()
## Plot expression distribution of krt and cldn family genes, and known epithelial marker 
DotPlot(Epi.subset, assay = "RNA",split.by = "preReg", 
        features = rev(c("cldne", "cldnf", "krt1-19d","krt17","cldnh","cldna", "krt93", "krt94", "cldn1", "cldni","krt4","fn1b","tp63","krtt1c19e")), 
        cols = c("red","blue"),
        dot.scale = 8) + RotatedAxis()
## Subset regenerating basal epithelial cells
WE.subset <- subset(Epi.subset, cells = rownames(Epi.subset@meta.data)[Epi.subset$major.cl == "Basal Epithelial" & Epi.subset$preReg == "Regenerating" ])
DefaultAssay(WE.subset) <- "integrated"
WE.subset <- RunPCA(WE.subset, npcs = 50, verbose = F)
WE.subset <- FindNeighbors(WE.subset, reduction = "pca", dims = 1:10) 
WE.subset <- FindClusters(WE.subset, resolution = 0.1)
WE.subset <- RunUMAP(WE.subset, reduction = "pca", dims = 1:10, n.neighbors = 50, min.dist = 0.2)
DimPlot(WE.subset, reduction = "umap") + NoAxes()

# Mesenchymal
Mes.subset <- subset(fin.integrated, idents = "Mesenchymal")
Mes.subset <- SCTransform(Mes.subset, vars.to.regress = c("percent.mito", "nCount_RNA","CC.Diff"), verbose = F, return.only.var.genes = F)
Mes.subset <- RunPCA(Mes.subset, npcs = 50, verbose = F)
Mes.subset <- FindNeighbors(Mes.subset, reduction = "pca", dims = 1:30) # PC number could change
Mes.subset <- FindClusters(Mes.subset, resolution = 0.4)
Mes.subset <- RunUMAP(Mes.subset, reduction = "pca", dims = 1:30, n.neighbors = 50, min.dist = 0.2)
## PC30 + res 0.4 revealed cluster 4 and 8 as epithelial populations. Remove and redo clustering and UMAP calc
Mes.epiRm.subset <- subset(Mes.subset, idents = c("0","1","2","3","5","6","7"))
## Redo clustering and non-linear dimensional reductions.
Mes.epiRm.subset <- FindNeighbors(Mes.epiRm.subset, reduction = "pca", dims = 1:25) # PC number could change
Mes.epiRm.subset <- FindClusters(Mes.epiRm.subset, resolution = 0.6)
Mes.epiRm.subset <- RunUMAP(Mes.epiRm.subset, reduction = "pca", dims = 1:25, n.neighbors = 20, min.dist = 0.2)
DimPlot(Mes.epiRm.subset, group.by = "Stage") + NoAxes()
## Slingshot pseudotime reconstructions
rd <- Mes.epiRm.subset@reductions$umap@cell.embeddings
cl <- Mes.epiRm.subset@active.ident
## Number of clusters
numCl <- length(unique(cl))
sds <- slingshot(rd, cl)
## Number of curves
numCurves <- length(sds@curves)
png("Mesenchymal_withSlingShotPT.png",    
    width = 3*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 5)
plot(rd, col = brewer.pal(numCl,"Paired")[as.factor(cl)], asp = 1, pch = 16)
# For fixed 4 colorful lines
lines(slingCurves(sds)$curve1, lwd = 2, col = "black")
lines(slingCurves(sds)$curve2, lwd = 2, col = "blue")
lines(slingCurves(sds)$curve3, lwd = 2, col = "red")
lines(slingCurves(sds)$curve4, lwd = 2, col = "brown")
dev.off()

