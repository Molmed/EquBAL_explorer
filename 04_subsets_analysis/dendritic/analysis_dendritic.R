# EquineAsthmaCellAtlas code can ONLY run on Seuratv4 (preferably v.4.3)

library(Seurat)

obj <- readRDS("subset_dendritic_raw.rds")

# Setting variable features of a subsetted object (Seurat github #2814)
VariableFeatures(obj[["SCT"]]) <- rownames(obj[["SCT"]]@scale.data)

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
obj <- FindClusters(obj, resolution = seq(0.2, 0.8, by=0.1))

# Running FindAllMarkers on selected resolutions
DefaultAssay(object = obj) <- "SCT"
Idents(object = obj) <- "integrated_snn_res.0.2"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.2 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.2", recorrect_umi=FALSE, only.pos = TRUE)
Idents(object = obj) <- "integrated_snn_res.0.3"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.3 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.3", recorrect_umi=FALSE, only.pos = TRUE)

# Annotating celltypes on lvl2
obj <- SetIdent(obj, value = obj$integrated_snn_res.0.2)
obj <- RenameIdents(object=obj, "0"="conventional DC", 
                    "1"="monocyte derived DC", 
                    "2"="migrating DC", 
                    "3"="migrating DC", 
                    "4"="migrating DC", 
                    "5"="SPP1+ DC", 
                    "6"="monocytes", 
                    "7"="unknown")
obj$CellType_lvl2 <- Idents(object = obj)

saveRDS(obj, "analyzed_dendritic.rds")