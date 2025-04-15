# EquineAsthmaCellAtlas code can ONLY run on Seuratv4 (preferably v.4.3)

library(Seurat)
library(DoubletFinder)

# Reading in count matrices
AS <- Read10X_h5("../../raw_data/Dropseq/AS.h5")
BO <- Read10X_h5("../../raw_data/Dropseq/BO.h5")
CA <- Read10X_h5("../../raw_data/Dropseq/CA.h5")
DO <- Read10X_h5("../../raw_data/Dropseq/DO.h5")
FAN <- Read10X_h5("../../raw_data/Dropseq/FAN.h5")
FS <- Read10X_h5("../../raw_data/Dropseq/FS.h5")
FI <- Read10X_h5("../../raw_data/Dropseq/FI.h5")
FLO <- Read10X_h5("../../raw_data/Dropseq/FLO.h5")
GU <- Read10X_h5("../../raw_data/Dropseq/GU.h5")
HO <- Read10X_h5("../../raw_data/Dropseq/HO.h5")
ME <- Read10X_h5("../../raw_data/Dropseq/ME.h5")
MY <- Read10X_h5("../../raw_data/Dropseq/MY.h5")
NE <- Read10X_h5("../../raw_data/Dropseq/NE.h5")
OD <- Read10X_h5("../../raw_data/Dropseq/OD.h5")
PR <- Read10X_h5("../../raw_data/Dropseq/PR.h5")
QU <- Read10X_h5("../../raw_data/Dropseq/QU.h5")
TI <- Read10X_h5("../../raw_data/Dropseq/TI.h5")
VA <- Read10X_h5("../../raw_data/Dropseq/VA.h5")
VAL <- Read10X_h5("../../raw_data/Dropseq/VAL.h5")
ZI <- Read10X_h5("../../raw_data/Dropseq/ZI.h5")

# Creating seurat objects
Obj_AS<-CreateSeuratObject(counts=AS, project= "AS") 
Obj_BO<-CreateSeuratObject(counts=BO, project= "BO") 
Obj_CA<-CreateSeuratObject(counts=CA, project= "CA") 
Obj_DO<-CreateSeuratObject(counts=DO, project= "DO") 
Obj_FAN<-CreateSeuratObject(counts=FAN, project= "FAN") 
Obj_FS<-CreateSeuratObject(counts=FS, project= "FS") 
Obj_FI<-CreateSeuratObject(counts=FI, project= "FI") 
Obj_FLO<-CreateSeuratObject(counts=FLO, project= "FLO") 
Obj_GU<-CreateSeuratObject(counts=GU, project= "GU") 
Obj_HO<-CreateSeuratObject(counts=HO, project= "HO") 
Obj_ME<-CreateSeuratObject(counts=ME, project= "ME") 
Obj_MY<-CreateSeuratObject(counts=MY, project= "MY") 
Obj_NE<-CreateSeuratObject(counts=NE, project= "NE") 
Obj_OD<-CreateSeuratObject(counts=OD, project= "OD") 
Obj_PR<-CreateSeuratObject(counts=PR, project= "PR") 
Obj_QU<-CreateSeuratObject(counts=QU, project= "QU") 
Obj_TI<-CreateSeuratObject(counts=TI, project= "TI") 
Obj_VA<-CreateSeuratObject(counts=VA, project= "VA") 
Obj_VAL<-CreateSeuratObject(counts=VAL, project= "VAL") 
Obj_ZI<-CreateSeuratObject(counts=ZI, project= "ZI") 

# Merging seurat objects
mergedobj <- merge(Obj_AS, c(Obj_BO, Obj_CA, Obj_DO, Obj_FAN, Obj_FS, 
                             Obj_FI, Obj_FLO, Obj_GU, Obj_HO, Obj_ME, 
                             Obj_MY, Obj_NE, Obj_OD, Obj_PR, Obj_QU,  
                             Obj_TI, Obj_VA, Obj_VAL, Obj_ZI), 
                   add.cell.ids = c("AS", "BO","CA","DO","FAN",
                                    "FS","FI","FLO","GU","HO",
                                    "ME","MY","NE","OD","PR",
                                    "QU","TI","VA","VAL","ZI"))

# Adding mitochondrial and ribosomal proportion data
total_counts_per_cell <- colSums(mergedobj@assays$RNA@counts)
mito_genes <- rownames(mergedobj)[grep("^MT-", rownames(mergedobj))]
mergedobj$percent_mito <- colSums(mergedobj@assays$RNA@counts[mito_genes, ])/total_counts_per_cell
ribo_genes <- rownames(mergedobj)[grep("^RP[SL]", rownames(mergedobj))]
mergedobj$percent_ribo <- colSums(mergedobj@assays$RNA@counts[ribo_genes, ])/total_counts_per_cell
obj <- mergedobj

# Filtering out cells with >8000 UMI, >10% mitochondrial genome, >15% ribosomal genome
obj_filtered <- subset(obj, nCount_RNA <= 8000) 
obj_filtered <- subset(obj_filtered, percent_mito <= 0.10) 
obj_filtered <- subset(obj_filtered, percent_ribo <= 0.15)

# Filtering out genes found in <0.1% of cells (for each horse separately, to preserve biological variability)
subsets <- SplitObject(obj_filtered, split.by = "orig.ident") 
subsets_filtered <- c()
for (i in 1:20) {
  horse_subset <- subsets[[i]]
  total_cells <- dim(horse_subset)[2]
  rarity_threshold <- as.integer(total_cells/1000)
  non_rare_genes <- rownames(horse_subset)[Matrix::rowSums(horse_subset)>rarity_threshold]
  subset_filtered <- subset(horse_subset, features = non_rare_genes)
  subsets_filtered <- c(subsets_filtered, subset_filtered)
}
obj_filtered <- merge(subsets_filtered[[1]], c(subsets_filtered[[2]], 
                                               subsets_filtered[[3]], 
                                               subsets_filtered[[4]],
                                               subsets_filtered[[5]],
                                               subsets_filtered[[6]],
                                               subsets_filtered[[7]], 
                                               subsets_filtered[[8]],
                                               subsets_filtered[[9]],
                                               subsets_filtered[[10]], 
                                               subsets_filtered[[11]],
                                               subsets_filtered[[12]],
                                               subsets_filtered[[13]], 
                                               subsets_filtered[[14]],
                                               subsets_filtered[[15]],
                                               subsets_filtered[[16]], 
                                               subsets_filtered[[17]],
                                               subsets_filtered[[18]],
                                               subsets_filtered[[19]], 
                                               subsets_filtered[[20]]),
                      add.cell.ids = c("AS", "BO","CA","DO", "FAN","FS","FI",
                                       "FLO","GU","HO","ME", "MY","NE","OD",
                                       "PR", "QU","TI", "VA", "VAL","ZI"))

# Filtering out cells with <200 genes
obj_filtered <- subset(obj_filtered, nFeature_RNA >= 200)
obj <- obj_filtered

# Loading local files for two DoubletFinder functions due to incompatibility with Seuratv4
# (see DoubletFinder github #180)
source("../custom_code/paramSweep.R")
source("../custom_code/doubletFinder.R") 

# Running normalization and detecting doublets for each horse
subsets <- SplitObject(obj, split.by = "orig.ident") 
for (i in 1:20) {
  
  # Running SCTransform, PCA
  horse_obj <- SCTransform(subsets[[i]], vars.to.regress= "nFeature_RNA")  
  horse_obj <- RunPCA(horse_obj, nfeatures.print = 10)
  
  # Calculating significant PCs
  stdv <- horse_obj[["pca"]]@stdev
  sum.stdv <- sum(horse_obj[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  
  # Running UMAP, FindNeighbors, FindClusters
  horse_obj <- RunUMAP(horse_obj, dims = 1:pcs)
  horse_obj <- FindNeighbors(object = horse_obj, dims = 1:pcs)              
  horse_obj <- FindClusters(object = horse_obj, resolution = 0.1)
  
  # Running PK identification (no ground-truth)
  sweep.list <- paramSweep(horse_obj, PCs = 1:pcs, num.cores = 1, sct=T) 
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Calculating optimal pk
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  # Estimating homotypic doublet proportion
  annotations <- horse_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(horse_obj@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # Running DoubletFinder
  horse_obj <- doubletFinder(seu = horse_obj, PCs = 1:pcs, pK = optimal.pk, nExp = nExp.poi.adj, sct=T)
  metadata <- horse_obj@meta.data
  colnames(metadata)[11] <- "doublet_finder2"
  horse_obj@meta.data <- metadata
  
  # Subsetting singlets
  horse_obj_singlets <- subset(horse_obj, doublet_finder2 == "Singlet")
  subsets[[i]] <- horse_obj_singlets
  remove(horse_obj_singlets)
}

# Merging singlet horse-sets back into one dataset
obj_singlets <- merge(x = subsets[[1]],
                      y = c(subsets[[2]], subsets[[3]], subsets[[4]],
                            subsets[[5]], subsets[[6]], subsets[[7]],
                            subsets[[8]], subsets[[9]], subsets[[10]],
                            subsets[[11]], subsets[[12]], subsets[[13]],
                            subsets[[14]], subsets[[15]], subsets[[16]],
                            subsets[[17]], subsets[[18]], subsets[[19]],
                            subsets[[20]]), project = "DropSeq Equine Asthma")

# Selecting variable features for the merged SCT object (seurat github #4145)
obj_features <- SelectIntegrationFeatures(object.list = subsets, nfeatures = 2000)
VariableFeatures(obj_singlets[["SCT"]]) <- obj_features

# Cell-Cycle Scoring
obj <- obj_singlets
obj <- SetIdent(obj, value = obj$orig.ident)
obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
obj <- SetIdent(obj, value = obj$orig.ident)

# Adding metadata

# Method
cells <- (dim(obj)[2])
obj$method <- rep("DropSeq", times = cells)

# Case/Control
obj@meta.data[["Phenotype1"]] <- NA
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "AS"] <- "control"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "BO"] <- "control"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "CA"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "DO"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "FAN"] <- "control"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "FS"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "FI"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "FLO"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "GU"] <- "control"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "HO"] <- "control"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "ME"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "MY"] <- "control"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "NE"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "OD"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "PR"] <- "control"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "QU"] <- "control"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "TI"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "VA"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "VAL"] <- "case"  
obj@meta.data[["Phenotype1"]][obj@meta.data[["orig.ident"]] == "ZI"] <- "case"  

# BAL phenotype
obj@meta.data[["BAL_phenotype"]] <- NA
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "AS"] <- "control"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "BO"] <- "neutrophilic_control"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "CA"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "DO"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "FAN"] <- "control"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "FS"] <- "mastocytic_eosinophilic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "FI"] <- "normal_BAL_case"   
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "FLO"] <- "normal_BAL_case"   
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "GU"] <- "control"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "HO"] <- "control"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "ME"] <- "neutrophilic_mixed"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "MY"] <- "control"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "NE"] <- "mastocytic_eosinophilic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "OD"] <- "mastocytic_eosinophilic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "PR"] <- "control"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "QU"] <- "control"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "TI"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "VA"] <- "mastocytic_eosinophilic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "VAL"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "ZI"] <- "neutrophilic_mixed"  

# Age
obj@meta.data[["Age"]] <- NA
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "AS"] <- 13  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "BO"] <- 10   
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "CA"] <- 12  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "DO"] <- 15  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "FAN"] <- 12  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "FS"] <- 18  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "FI"] <- NA   
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "FLO"] <- 9   
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "GU"] <- 17  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "HO"] <- 16  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "ME"] <- 13  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "MY"] <- 14  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "NE"] <- 4  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "OD"] <- 9  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "PR"] <- 18  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "QU"] <- 8  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "TI"] <- 21  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "VA"] <- 6  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "VAL"] <- 25  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "ZI"] <- 14  

# Sex
obj@meta.data[["Sex"]] <- NA
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "AS"] <- "G"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "BO"] <- "M"   
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "CA"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "DO"] <- "G"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "FAN"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "FS"] <- "G"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "FI"] <- "G"   
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "FLO"] <- "G"   
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "GU"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "HO"] <- "G"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "ME"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "MY"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "NE"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "OD"] <- "G"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "PR"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "QU"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "TI"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "VA"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "VAL"] <- "G"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "ZI"] <- "G"  

# Breed
obj@meta.data[["Breed"]] <- NA
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "AS"] <- "Standardbred trotter"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "BO"] <- "Standardbred trotter"   
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "CA"] <- "Swedish Warmblood"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "DO"] <- "Swedish Warmblood"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "FAN"] <- "Standardbred trotter"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "FS"] <- "Swedish Riding Pony"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "FI"] <- "Swedish Warmblood"   
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "FLO"] <- "Swedish Warmblood"   
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "GU"] <- "Standardbred trotter"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "HO"] <- "Standardbred trotter"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "ME"] <- "Connemara"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "MY"] <- "Standardbred trotter"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "NE"] <- "Standardbred trotter"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "OD"] <- "Coldblooded trotter"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "PR"] <- "Standardbred trotter"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "QU"] <- "Standardbred trotter"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "TI"] <- "Icelandic horse"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "VA"] <- "Tinker"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "VAL"] <- "Welsh Mountain"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "ZI"] <- "NS Draught Horse"  

saveRDS(obj, "preprocessed_dropseq.rds")
