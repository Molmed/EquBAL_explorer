# EquineAsthmaCellAtlas code can ONLY run on Seuratv4 (preferably v.4.3)

library(Seurat)

obj <- readRDS("../02_integration/integrated_Dropseq_HIVE.rds")

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

# Running FindAllMarkers on every resolution
DefaultAssay(object = obj) <- "SCT"
Idents(object = obj) <- "integrated_snn_res.0.1"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.1 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.1")
Idents(object = obj) <- "integrated_snn_res.0.2"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.2 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.2")
Idents(object = obj) <- "integrated_snn_res.0.3"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.3 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.3")
Idents(object = obj) <- "integrated_snn_res.0.4"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.4 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.4")
Idents(object = obj) <- "integrated_snn_res.0.5"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.5 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.5")
Idents(object = obj) <- "integrated_snn_res.0.6"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.6 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.6")
Idents(object = obj) <- "integrated_snn_res.0.7"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.7 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.7")
Idents(object = obj) <- "integrated_snn_res.0.8"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.8 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.8")
Idents(object = obj) <- "integrated_snn_res.0.9"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.9 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.9")
Idents(object = obj) <- "integrated_snn_res.1"
obj <- PrepSCTFindMarkers(obj)
markers_res.1.0 <- FindAllMarkers(obj, group.by = "integrated_snn_res.1")

# Finding and removing double positives
clust0 <- subset(obj, subset = integrated_snn_res.0.6 == 0)
clust0 <- RunUMAP(clust0, dims = 1:20)
clust0 <- FindNeighbors(object = clust0, reduction="pca", dims = 1:20)
clust0 <- FindClusters(object = clust0, resolution = seq(0.1, 0.3, by=0.1))
clust0 <- PrepSCTFindMarkers(clust0)
markers_clust0 <- FindAllMarkers(clust0, group.by = "integrated_snn_res.0.3")
clust0_0 <- subset(clust0, subset = integrated_snn_res.0.3 == 0)
cells_to_delete <- colnames(clust0_0)
obj_noDP <- subset(obj, cells = setdiff(Cells(obj), cells_to_delete))

# Reclustering cleaned dataset (res 0.1-1.0)

# Running PCA
obj <- RunPCA(obj_noDP)

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

# Running FindAllMarkers on every resolution
DefaultAssay(object = obj) <- "SCT"
Idents(object = obj) <- "integrated_snn_res.0.1"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.1 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.1")
Idents(object = obj) <- "integrated_snn_res.0.2"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.2 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.2")
Idents(object = obj) <- "integrated_snn_res.0.3"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.3 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.3")
Idents(object = obj) <- "integrated_snn_res.0.4"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.4 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.4")
Idents(object = obj) <- "integrated_snn_res.0.5"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.5 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.5")
Idents(object = obj) <- "integrated_snn_res.0.6"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.6 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.6")
Idents(object = obj) <- "integrated_snn_res.0.7"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.7 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.7")
Idents(object = obj) <- "integrated_snn_res.0.8"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.8 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.8")
Idents(object = obj) <- "integrated_snn_res.0.9"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.9 <- FindAllMarkers(obj, group.by = "integrated_snn_res.0.9")
Idents(object = obj) <- "integrated_snn_res.1"
obj <- PrepSCTFindMarkers(obj)
markers_res.1.0 <- FindAllMarkers(obj, group.by = "integrated_snn_res.1")

# Annotating major celltypes
obj <- SetIdent(obj, value = obj$integrated_snn_res.0.1)
obj <- RenameIdents(object=obj, "0"="Macrophages", "1"="T cells", 
                    "2"="Proliferating Macrophages", "3"="Putative doublets", 
                    "4"="Dendritic cells", "5"="Neutrophils", 
                    "6"="Mast cells", "7"="Dendritic cells", 
                    "8"="Eosinophils")
obj$CellType_lvl1 <- Idents(object = obj)

# Subsetting major celltypes
subset_tcells <- subset(obj, subset = CellType_lvl1  == "T cells")
subset_macrophages <- subset(obj, subset = CellType_lvl1  == "Macrophages")
subset_neutrophils <- subset(obj, subset = CellType_lvl1  == "Neutrophils")
subset_mastcells <- subset(obj, subset = CellType_lvl1  == "Mast cells")
subset_dendritic <- subset(obj, subset = CellType_lvl1  == "Dendritic cells")

saveRDS(obj, "analyzed_Dropseq_HIVE.rds")
saveRDS(subset_tcells, "../04_subsets_analysis/tcells/subset_tcells_raw.rds")
saveRDS(subset_macrophages, "../04_subsets_analysis/macrophages/subset_macrophages_raw.rds")
saveRDS(subset_neutrophils, "../04_subsets_analysis/neutrophils/subset_neutrophils_raw.rds")
saveRDS(subset_mastcells, "../04_subsets_analysis/mastcells/subset_mastcells_raw.rds")
saveRDS(subset_dendritic, "../04_subsets_analysis/dendritic/subset_dendritic_raw.rds")


