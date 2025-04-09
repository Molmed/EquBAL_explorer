# EquineAsthmaCellAtlas code can ONLY run on Seuratv4 (preferably v.4.3)

library(Seurat)
library(DoubletFinder)

# Reading in count matrices
AM <- read.table("../../raw_data/HIVE/FU190-AM_20230505120209_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
BE <- read.table("../../raw_data/HIVE/FU190-BE_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
DE <- read.table("../../raw_data/HIVE/FU190-DE_20230505120209_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
FA <- read.table("../../raw_data/HIVE/FU190-FA_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
FL <- read.table("../../raw_data/HIVE/FU190-FL_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
GN <- read.table("../../raw_data/HIVE/FU190-GN_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
KA20 <- read.table("../../raw_data/HIVE/FU190-KA20_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
KA30 <- read.table("../../raw_data/HIVE/FU190-KA30_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
LE <- read.table("../../raw_data/HIVE/FU190-LE_20230505120209_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
NAT <- read.table("../../raw_data/HIVE/FU190-NAT_20230505120209_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
QU15 <- read.table("../../raw_data/HIVE/FU190-QU15_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
QU30 <- read.table("../../raw_data/HIVE/FU190-QU30_20230505120209_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
ST <- read.table("../../raw_data/HIVE/FU190-ST_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)

# Creating seurat objects
Obj_AM<-CreateSeuratObject(counts=AM, project= "AM") 
Obj_BE <-CreateSeuratObject(counts=BE, project="BE") 
Obj_DE <-CreateSeuratObject(counts=DE, project = "DE") 
Obj_FA <-CreateSeuratObject(counts=FA, project="FA") 
Obj_FL <-CreateSeuratObject(counts=FL, project = "FL") 
Obj_GN <-CreateSeuratObject(counts=GN, project="GN") 
Obj_KA20 <-CreateSeuratObject(counts=KA20, project="KA20") 
Obj_KA30 <-CreateSeuratObject(counts=KA30, project ="KA30") 
Obj_LE <-CreateSeuratObject(counts=LE, project ="LE") 
Obj_NAT <-CreateSeuratObject(counts=NAT, project ="NAT")
Obj_QU15 <-CreateSeuratObject(counts=QU15, project ="QU15")
Obj_QU30 <-CreateSeuratObject(counts=QU30, project ="QU30")
Obj_ST <-CreateSeuratObject(counts=ST, project ="ST")


# Merging seurat objects
mergedobj <- merge(Obj_AM, c(Obj_BE, Obj_DE, Obj_FA, Obj_FL, Obj_GN, 
                             Obj_KA20, Obj_KA30, Obj_LE, Obj_NAT, Obj_QU15, 
                             Obj_QU30, Obj_ST), 
                 add.cell.ids = c("AM", "BE", "DE", "FA", "FL", 
                                  "GN", "KA20", "KA30", "LE", "NAT", 
                                  "QU15", "QU30", "ST"))

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
for (i in 1:13) {
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
                                               subsets_filtered[[13]]),
                      add.cell.ids = c("AM", "BE", "DE", "FA", "FL", 
                                       "GN", "KA20", "KA30", "LE", "NAT", 
                                       "QU15", "QU30", "ST"))

# Filtering out cells with <400 genes
obj_filtered <- subset(obj_filtered, nFeature_RNA >= 400)
obj <- obj_filtered

# Loading local files for two DoubletFinder functions due to incompatibility with Seuratv4
# (see DoubletFinder github #180)
source("../custom_code/paramSweep.R")
source("../custom_code/doubletFinder.R") 

# Running normalization and detecting doublets for each horse
subsets <- SplitObject(obj, split.by = "orig.ident") 
for (i in 1:13) {
  
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
                            subsets[[11]], subsets[[12]], subsets[[13]]), 
                      project = "HIVE Equine Asthma")

# Selecting variable features for the merged SCT object (seurat github #4145)
obj_features <- SelectIntegrationFeatures(object.list = subsets, nfeatures = 2000)
VariableFeatures(obj_singlets[["SCT"]]) <- obj_features
obj <- obj_singlets

# Cell-Cycle Scoring
obj <- SetIdent(obj, value = obj$orig.ident)
obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
obj <- SetIdent(obj, value = obj$orig.ident)

# Adding metadata

# Method
cells <- (dim(obj)[2])
obj$method <- rep("HIVE", times = cells)

# Case/Control
obj$Phenotype1 <- rep("case", times = cells)

# BAL phenotype
obj@meta.data[["BAL_phenotype"]] <- NA
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "AM"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "BE"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "DE"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "FA"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "FL"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "GN"] <- "neutrophilic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "KA20"] <- "mastocytic_eosinophilic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "KA30"] <- "mastocytic_eosinophilic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "LE"] <- "mastocytic_eosinophilic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "NAT"] <- "normal_BAL_case"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "QU15"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "QU30"] <- "mastocytic"  
obj@meta.data[["BAL_phenotype"]][obj@meta.data[["orig.ident"]] == "ST"] <- "mastocytic"  

# Age
obj@meta.data[["Age"]] <- NA
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "AM"] <- 18  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "BE"] <- 12  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "DE"] <- 13  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "FA"] <- 8  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "FL"] <- 15  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "GN"] <- 20  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "KA20"] <- 16  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "KA30"] <- 16  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "LE"] <- 13  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "NAT"] <- 5   
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "QU15"] <- 8  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "QU30"] <- 8  
obj@meta.data[["Age"]][obj@meta.data[["orig.ident"]] == "ST"] <- NA  

# Sex
obj@meta.data[["Sex"]] <- NA
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "AM"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "BE"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "DE"] <- "M" 
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "FA"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "FL"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "GN"] <- "G"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "KA20"] <- "G"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "KA30"] <- "G"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "LE"] <- "G"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "NAT"] <- "M"   
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "QU15"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "QU30"] <- "M"  
obj@meta.data[["Sex"]][obj@meta.data[["orig.ident"]] == "ST"] <- "G"  

# Breed
obj@meta.data[["Breed"]] <- NA
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "AM"] <- "Swedish Warmblood"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "BE"] <- "Swedish Warmblood"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "DE"] <- "Danish Warmblood"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "FA"] <- "Shetland pony"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "FL"] <- "Swedish Warmblood"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "GN"] <- "Icelandic horse"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "KA20"] <- "Icelandic horse"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "KA30"] <- "Icelandic horse"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "LE"] <- "Swedish Warmblood"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "NAT"] <- "Icelandic horse"   
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "QU15"] <- "Swedish Warmblood"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "QU30"] <- "Swedish Warmblood"  
obj@meta.data[["Breed"]][obj@meta.data[["orig.ident"]] == "ST"] <- "Welsh pony"  

saveRDS(obj, "preprocessed_hive.rds")
