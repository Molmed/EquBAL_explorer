# EquineAsthmaCellAtlas code can ONLY run on Seuratv4 (preferably v.4.3)

library(Seurat)

obj <- readRDS("subset_tcells_raw.rds")

# Renormalizing
DefaultAssay(obj) <-"RNA"
sample_list <- SplitObject(obj, split.by = "orig.ident")
sample_list <- lapply(sample_list, function(x) {
  SCTransform(x, verbose = TRUE, vars.to.regress= "nFeature_RNA")
})
obj <- Reduce(function(x, y) merge(x, y), sample_list)

# Selecting variable features for the merged SCT object (seurat github #4145)
obj_features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 2000)
VariableFeatures(obj[["SCT"]]) <- obj_features

# Filtering
obj <- subset(obj, nFeature_RNA > 300)

# Running PCA
obj <- RunPCA(obj)

# Calculating significant PCs
stdv <- obj[["pca"]]@stdev
sum.stdv <- sum(obj[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

# Running UMAP, clustering
obj <- RunUMAP(obj, reduction = "pca", dims = 1:pcs)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:pcs)
obj <- FindClusters(obj, resolution = seq(0.1, 1.0, by=0.1))

# Filtering cells identified as double-positive
obj <- SetIdent(obj, value = obj$SCT_snn_res.0.1)
obj <- subset(obj, idents = 1, invert = TRUE)

# Rerunning PCA
obj <- RunPCA(obj)

# Calculating significant PCs
stdv <- obj[["pca"]]@stdev
sum.stdv <- sum(obj[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

# Rerunning UMAP, clustering
obj <- RunUMAP(obj, reduction = "pca", dims = 1:pcs)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:pcs)
obj <- FindClusters(obj, resolution = seq(0.1, 0.8, by=0.1))

# Running FindAllMarkers on selected resolutions
DefaultAssay(object = obj) <- "SCT"
Idents(object = obj) <- "SCT_snn_res.0.5"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.5 <- FindAllMarkers(obj, group.by = "SCT_snn_res.0.5", recorrect_umi=FALSE, only.pos = TRUE)
Idents(object = obj) <- "SCT_snn_res.0.6"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.6 <- FindAllMarkers(obj, group.by = "SCT_snn_res.0.6", recorrect_umi=FALSE, only.pos = TRUE)
Idents(object = obj) <- "SCT_snn_res.0.7"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.7 <- FindAllMarkers(obj, group.by = "SCT_snn_res.0.7", recorrect_umi=FALSE, only.pos = TRUE)
Idents(object = obj) <- "SCT_snn_res.0.8"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.8 <- FindAllMarkers(obj, group.by = "SCT_snn_res.0.8", recorrect_umi=FALSE, only.pos = TRUE)

# Annotating celltypes on lvl2
obj <- SetIdent(obj, value = obj$SCT_snn_res.0.7)
obj <- RenameIdents(object=obj, "0"="CD4+ T-cells", 
                    "1"="GZMA+ cytotoxic T-cells", 
                    "2"="γδ T-cells", 
                    "3"="NK cells / NKT cells", 
                    "4"="GZMA+ cytotoxic T-cells", 
                    "5"="GZMK+ cytotoxic T-cells", 
                    "6"="CD4+ ribosomal+ T-cells", 
                    "7"="KLRB1+- TRDC- cells", 
                    "8"="ISG+ T-cells", 
                    "9"="GZMK+ cytotoxic T-cells")
obj$CellType_lvl2 <- Idents(object = obj)

saveRDS(obj, "analyzed_tcells.rds")
