# EquineAsthmaCellAtlas code can ONLY run on Seuratv4 (preferably v.4.3)

library(Seurat)

# Reading in objects
dropseq <- readRDS("../01_preprocessing/Dropseq/preprocessed_dropseq.rds")
hive <- readRDS("../01_preprocessing/HIVE/preprocessed_hive.rds")

# Merging objects, setting variable features (default var.features slot empties during merging)
mergedobj <- merge(dropseq, y=hive, project="Equine asthma Cell Atlas")
VariableFeatures(mergedobj[["SCT"]]) <- rownames(mergedobj[["SCT"]]@scale.data)

# Splitting by horse ID: integration is done between individual horses
obj_list <- SplitObject(mergedobj, split.by = "orig.ident")

# Fixing SCTmodel bug in Seuratv4 (SCTransformed objects after merge+split having several models instead of one)
for (i in 1:33) {
  objsmall <- obj_list[[i]]
  empty_models <- sapply(objsmall[['SCT']]@SCTModel.list, function(m) nrow(m@cell.attributes))
  objsmall[['SCT']]@SCTModel.list <- objsmall[['SCT']]@SCTModel.list[empty_models > 0]
  obj_list[[i]] <- objsmall
}

# Integration
var_features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 2000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = var_features)
obj_list <- lapply(X = obj_list, FUN = RunPCA, features = var_features)
anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = var_features, 
                                  reduction = "rpca", normalization.method = "SCT", 
                                  k.anchor = 15)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

saveRDS(integrated, "integrated_Dropseq_HIVE.rds")